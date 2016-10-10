"""
Cython wrapper for the C++ translation of the Voronoi diagram part of Boost's
Polygon library: http://www.boost.org/doc/libs/1_53_0_beta1/libs/polygon/doc/voronoi_diagram.htm
"""
from __future__ import division

import sys as _sys
import struct
import copy as _copy
import unicodedata as _unicodedata
import time as _time
import math

from cython.operator cimport dereference as deref

SILENT = False
def log_action(description):
    if not SILENT:
        print description

log_action("Python binding clipper library")


cdef extern from "Python.h":
    Py_INCREF(object o)
    object Py_BuildValue(char * format, ...)
    object PyBuffer_FromMemory(void * ptr, int size)
    #int PyArg_ParseTuple(object struct,void* ptr)
    char * PyString_AsString(object string)
    int PyArg_VaParse(object args, char *format, ...)
    int PyArg_Parse(object args, char *format, ...)
    int PyObject_AsReadBuffer(object obj, void*buffer, int * buffer_len)
    object PyBuffer_FromObject(object base, int offset, int size)
    object PyBuffer_FromReadWriteObject(object base, int offset, int size)
    PyBuffer_New(object o)

cdef extern from "stdio.h":
    cdef void printf(char * , ...)

cdef extern from "stdlib.h":
    cdef void*malloc(unsigned int size)
    cdef void*free(void * p)
    char * strdup(char * str)
    int strcpy(void * str, void * src)
    int memcpy(void * str, void * src, int size)

cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()
        int size()


cdef extern from "voronoi.hpp":
    cdef struct Point:
        int X
        int Y

    cdef struct Segment:
        Point p0
        Point p1
        
    cdef struct c_Vertex:
        double X
        double Y
        
    cdef struct c_Edge:
        long long start
        long long end
        int isPrimary
        size_t site1
        size_t site2
        int isLinear
        long cell;
        long twin;			

    cdef struct c_Cell:
        size_t cell_identifier;
        size_t site
        int contains_point
        int contains_segment
        int is_open		
        vector[long long] vertices
        vector[long long] edges
        int source_category
				
    cdef cppclass VoronoiDiagram:
        VoronoiDiagram()
        void AddPoint(Point p)
        void AddSegment(Segment s)
        void Construct()
        void GetEdges(vector[c_Vertex] & c_vertices, vector[c_Edge] & c_edges)
        void GetCells(vector[c_Vertex] & c_vertices, vector[c_Edge] & c_edges, vector[c_Cell] & c_cell_edges)
        vector[Point] GetPoints()
        vector[Segment] GetSegments()

####################################
##VORONOY UTILS
####################################
class Vertex:
    X = 0.0
    Y = 0.0
    def __init__(self,x,y):
        self.X = x
        self.Y = y
		
class Edge:
    start = -1
    end = -1
    is_primary = False
    site1 = -1
    site2 = -1
    is_linear = False
    cell = -1
    twin = -1
			
class Cell:
    cell_identifier = -1
    site = -1
    contains_point = False
    contains_segment = False
    is_open = False
	
    vertices = None
    edges = None
	
    source_category = None

    def __init__(self, cell_identifier, site, vertices, edges, source_category):
        self.cell_identifier = cell_identifier
        self.site = site
        self.vertices = vertices
        self.edges = edges
        self.vertices.append(self.vertices[0])
        self.source_category = source_category

class VoronoiException(Exception):
    pass
	
class FocusOnDirectixException(Exception):
    pass
	
class UnsolvableParabolaEquation(Exception):
    pass
	
####################################
##ROTATION
####################################		
def Rotate(point, theta):
    t = -1 * theta
    cos = math.cos(t)
    sin = math.sin(t)
    return [(point[0] * cos) - (point[1] * sin),
	(point[0] * sin) + (point[1] * cos)]
	
def RotateWithShift(point, theta, shift_x, shift_y):
    return Rotate([point[0] - shift_x, point[1] - shift_y], theta)
	
def Unrotate(point, theta, shift_x, shift_y):
    cos = math.cos(theta)
    sin = math.sin(theta)
    return [(point[0] * cos) - (point[1] * sin) + shift_x, (point[0] * sin) + (point[1] * cos) + shift_y]
	
def GetLineAngleInRadians(start_point_x, start_point_y,end_point_x, end_point_y):
    return math.atan2(end_point_y - start_point_y, end_point_x - start_point_x)	
	
	
####################################
##DISTANCE
####################################	
def DistanceSquared(point_start, point_end):
    """Returns the squared length of the line.
    :return: a float representing the squared length of the line.
    """
    return pow(point_end[0] - point_start[0], 2) + pow(point_end[1] - point_start[1], 2)
	
def Distance(point_start, point_end):
    """Returns the length of the line"""
    return math.sqrt(DistanceSquared(point_start, point_end))
	
			
####################################
##PYVORONOI OPERATIONS
####################################		

cdef class Pyvoronoi:
    cdef VoronoiDiagram *thisptr
    cdef int constructed

    inputPoints = []
    inputSegments = []
    outputEdges = []
    outputVertices = []
    outputCells = []

    cdef public int SCALING_FACTOR

    def __cinit__(self, scaling_factor = None):
        """ Creates an instance of the Pyvoronoi class.
        """
        log_action("Creating an VoronoiDiagram instance")
        self.thisptr = new VoronoiDiagram()
        self.constructed = 0

        if scaling_factor != None:
            self.SCALING_FACTOR = scaling_factor
        else:
            self.SCALING_FACTOR = 1
			
        del self.inputPoints[:]
        del self.inputSegments[:]			

    def __dealloc__(self):
        log_action("Deleting the VoronoiDiagram instance")
        del self.thisptr
                        
    def AddPoint(self, point):
        """ Add a point
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Point c_point = self._to_voronoi_point(point)
        self.thisptr.AddPoint(c_point)
        self.inputPoints.append(point)

    def AddSegment(self, segment):
        """ Add a segment
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Segment c_segment = self._to_voronoi_segment(segment)
        self.thisptr.AddSegment(c_segment)
        self.inputSegments.append(segment)
		
    def Construct(self):
        """ Generates the voronoi diagram for the added points and segments. Voronoi cell structure will be generated.
        """
        
        if self.constructed == 1:
            raise VoronoiException('Construct() has already been called')

        self.constructed = 1
		
        self.thisptr.Construct()
        
        cdef vector[c_Edge] c_edges
        cdef vector[c_Vertex] c_vertices
        cdef vector[c_Cell] c_cells
        self.thisptr.GetCells(c_vertices, c_edges, c_cells)		
		
        del self.outputEdges[:] 
        del self.outputVertices[:]
        del self.outputCells[:]
		
		
        #--------------------------------------------------------
        #LOGIC TO PRINT CELL OBJECTS
        #--------------------------------------------------------
        cdef size_t count = c_edges.size()
        for i in range(count):
            edge = Edge()
            edge.start = c_edges[i].start
            edge.end = c_edges[i].end
            edge.is_primary = c_edges[i].isPrimary != False
            edge.is_linear = c_edges[i].isLinear != False
        
            edge.site1 = c_edges[i].site1
            edge.site2 = c_edges[i].site2
			
            edge.cell = c_edges[i].cell
            edge.twin = c_edges[i].twin
        
            self.outputEdges.append(edge)

        count = c_vertices.size()      
        for i in range(count):
            vertex = Vertex(self._from_voronoi_value(c_vertices[i].X),  self._from_voronoi_value(c_vertices[i].Y))
            self.outputVertices.append(vertex)

        count = c_cells.size()

        for i in range(count):
            c_cell = c_cells[i]
            outputCell = Cell(c_cell.cell_identifier, c_cell.site, c_cell.vertices, c_cell.edges, c_cell.source_category)
            outputCell.contains_point = c_cell.contains_point != False
            outputCell.contains_segment = c_cell.contains_segment != False
            outputCell.is_open = c_cell.is_open != False
            self.outputCells.append(outputCell)
				
        cdef vector[Point] points = self.thisptr.GetPoints()
        count = points.size()	
		
				
    def GetVertices(self):
        """ Returns the edges of the voronoi diagram.             
        """

        if self.constructed == 0:
            raise VoronoiException('Construct() has not been called')

        return self.outputVertices

    def GetEdges(self):
        """ Returns the edges of the voronoi diagram.             
        """

        if self.constructed == 0:
            raise VoronoiException('Construct() has not been called')
        
        return self.outputEdges
		
    def GetCells(self):
        """ Return pointer to the edges that make the cellss            
        """	
        if self.constructed == 0:
            raise VoronoiException('Construct() has not been called')
                
        return self.outputCells

    def GetPoints(self):
        """ Returns the points added to the voronoi diagram
        """
        return self.inputPoints

    def GetSegments(self):
        """ Returns the segments added to the voronoi diagram
        """
        return self.inputSegments
		

    def ReturnCurvedSiteInformation(self, edge):
        """Return the index of the point side and the segment site associated  with a segment index
        """		
        twinEdge = self.outputEdges[edge.twin]	
		
        cell = self.outputCells[edge.cell]
        twinCell = self.outputCells[twinEdge.cell]

        pointSite = self.RetrievePoint(cell) if cell.contains_point == True else self.RetrievePoint(twinCell)
        segmentSite = self.RetrieveSegment(twinCell) if cell.contains_point == True else self.RetrieveSegment(cell)
			
        return [pointSite, segmentSite]
        
    def DiscretizeCurvedEdge(self, index, max_dist, parabola_equation_tolerance = 1):
        if(max_dist <= 0):
            raise Exception("Max distance must be greater than 0. Value passed: {0}".format(max_dist))
        edge = self.outputEdges[index]
        sites = self.ReturnCurvedSiteInformation(edge)
        pointSite = sites[0]
        segmentSite = sites[1]
		
        edgeStartVertex = self.outputVertices[edge.start]
        edgeEndVertex = self.outputVertices[edge.end]
        return self.Discretize(pointSite,segmentSite, [edgeStartVertex.X,edgeStartVertex.Y], [edgeEndVertex.X, edgeEndVertex.Y], max_dist, parabola_equation_tolerance)


    def RetrievePoint(self, cell):
        """Retrive the input point associated with a cell.
		:param cell: the cell that contains a point. The point can be either a input point or the end point of an input segment.
        """
        if(cell.source_category == 0):
            return self.inputPoints[cell.site]

        input_segment = self.RetrieveSegment(cell)
        if(cell.source_category == 1):
            return input_segment[0]
        else:			
            return input_segment[1]
		
    def RetrieveSegment(self, cell):
        """Retrive the input segment associated with a cell.
        """
        return self.inputSegments[cell.site - len(self.inputPoints)]
		
    def GetParabolaY(self, x, focus, directrix_y):
        """
		Solve the parabola equation for a given value on the x-axis and return the associated value on the y-axis. 
		This equation assumes that the directix is parallel to the x-axis. 
        Parabola equation are different if the directix is parallel to the y-axis.
		:param x: the x-value used to solve the equation.
        :param focus: the focus point used for solving the equation of the parabola.
        :param directix: the directix value used for solving the equation of the parabola.	
        :return: the associated value on the y-axis.
        """
        return (pow(x - focus[0], 2) + pow(focus[1], 2) - pow(directrix_y, 2)) / (2 * (focus[1] - directrix_y));

    def CheckUnsolvableParabolaEquation(self, boost_x, boost_y, focus, directix, tolerance):
        """
        Compare the y-coordinate of a point on the parabola returned by Boost with the computed value. 
		The function will return an exception if the difference between the computed y-value and the y-value returned by Boost.
		The computed point will be returned otherwise.
        :param boost_x: the x-value of the point parabola returned by boost.
        :param boost_y: the y-value of the point parabola returned by boost.
        :param focus: the focus point used for solving the equation of the parabola.
        :param directix: the directix value used for solving the equation of the parabola.
		:param tolerance: the distance allowed between the point computed by boost and the point computed by the equation.
		:return: the point on the parabola computed using the value of boost_x.
		"""	
        computed_point_y = self.GetParabolaY(boost_x, focus, directix)
        delta = computed_point_y - boost_y if computed_point_y > boost_y else boost_y - computed_point_y
        if delta > tolerance:
            raise UnsolvableParabolaEquation("The computed Y on the parabola for the starting / ending point is different from the rotated point returned by Boost. Difference: {0}. Maximum tolerance: {1}".format(delta, tolerance))
        return [boost_x, computed_point_y]
	
    def Discretize(self, point, segment, parabola_start, parabola_end, max_dist, parabola_equation_tolerance):
        """
        Interpolate points on a parabola. The points are garanteed to be closer than the value of the parameter max_dist.
        :param point: The input point associated with the cell or the neighbour cell. The point is used as the focus in the equation of the parabola.
        :param segment: The input segment associated with the cell or the neighbour cell. The point is used as the directix in the equation of the parabola.
        :param max dist: The maximum distance between 2 vertices on the discretized geometry.
        :param parabola_equation_tolerance: The maximum difference allowed between the y coordinate returned by Boost, and the equation of the parabola.
		:return: the list of points on the parabola.
		"""
	
		#Test if the input point is on the input line. If yes, it is impossible to compute a parabola.
        if(point[0] == segment[0][0] and point[1] == segment[0][1]) or (point[0] == segment[1][0] and point[1] == segment[1][1]):
            raise FocusOnDirectixException()
		
		###############################
        #Rotate the input information
		#Rotation is done so that the 
		#Input segment is aligned with the x-axis.
		#Input point becomes the focus of the parabola equation.
		#Input segment becomes the directix of the parabola equation.
		###############################
        shift_x = min(segment[0][0], segment[1][0])
        shift_y = min(segment[0][1], segment[1][1])	
        angle = GetLineAngleInRadians(
            segment[0][0],
            segment[0][1],
            segment[1][0],
            segment[1][1],			
        )
		
        focus_rotated = RotateWithShift(point, angle, shift_x, shift_y)
        directix_start_rotated = RotateWithShift(segment[0], angle, shift_x, shift_y)
        directix_end_rotated = RotateWithShift(segment[1], angle, shift_x, shift_y)
        parabola_start_rotated = RotateWithShift(parabola_start, angle, shift_x, shift_y)
        parabola_end_rotated = RotateWithShift(parabola_end, angle, shift_x, shift_y)
				
		###############################
        #Validate the equation for the 
		#first and last point given by 
		#Boost
		###############################
        directix = directix_start_rotated[1]
        parabola_start_rotated_check = self.CheckUnsolvableParabolaEquation(parabola_start_rotated[0], parabola_start_rotated[1], focus_rotated, directix, parabola_equation_tolerance)
        parabola_end_rotated_check = self.CheckUnsolvableParabolaEquation(parabola_end_rotated[0], parabola_end_rotated[1], focus_rotated, directix, parabola_equation_tolerance)

		
		###############################
        #Compute the intermediate points
		###############################
        densified_rotated  = [parabola_start_rotated_check]
        previous = densified_rotated[0]
        next  = [parabola_end_rotated_check]
		
        while (len(next) > 0):	
            current = next[-1]
            distance = Distance(current, previous)
			
            if distance > max_dist:
                mid_point_x = (previous[0] + current[0]) / 2
                mid_point = [mid_point_x, self.GetParabolaY(mid_point_x, focus_rotated, directix)]			
                next.append(mid_point)
            else:
                densified_rotated.append(current)
                next.pop()
                previous = current			
                

        densified_rotated[0] = parabola_start_rotated
        densified_rotated[-1] = parabola_end_rotated		
		###############################
        #Unrotate the computed intermediate points
		###############################
        densified = map(lambda x: Unrotate(x, angle, shift_x, shift_y), densified_rotated)
        return densified
				

    cdef Segment _to_voronoi_segment(self, object py_segment):
        return Segment(self._to_voronoi_point(py_segment[0]), self._to_voronoi_point(py_segment[1]))

    cdef Point _to_voronoi_point(self, object py_point):
        return Point(self._to_voronoi_int(py_point[0]), self._to_voronoi_int(py_point[1]))

    cdef int _to_voronoi_int(self, val):
        return val * self.SCALING_FACTOR

    cdef double _to_voronoi_double(self, val):
        return val * <double>self.SCALING_FACTOR

    cdef double _from_voronoi_value(self, val):
        return val / <double>self.SCALING_FACTOR