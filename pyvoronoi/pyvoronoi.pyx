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

SILENT = True
def log_action(description):
    if not SILENT:
        print(description)


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
        int isLinear
        long long cell;
        long long twin;

    cdef struct c_Cell:
        long long cell_identifier;
        long long site
        int contains_point
        int contains_segment
        int is_open
        int is_degenerate
        vector[long long] vertices
        vector[long long] edges
        int source_category

    cdef cppclass VoronoiDiagram:
        VoronoiDiagram()
        void AddPoint(Point p)
        void AddSegment(Segment s)
        void Construct() nogil
        vector[Point] GetPoints()
        vector[Segment] GetSegments()
        vector[int] GetIntersectingSegments()
        vector[int] GetDegenerateSegments()
        vector[int] GetPointsOnSegments()
        void MapVertexIndexes()
        void MapEdgeIndexes()
        void MapCellIndexes()

        long long CountVertices()
        long long CountEdges()
        long long CountCells()

        Point GetPoint(int index)
        Segment GetSegment(int index)
        c_Vertex GetVertex(long long index)
        c_Edge GetEdge(long long index)
        c_Cell GetCell(long long index)



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
    is_linear = False
    cell = -1
    twin = -1

    def __init__(self, start, end, cell, twin):
        self.start = start
        self.end = end
        self.cell = cell
        self.twin = twin

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
        self.source_category = source_category
        self.vertices = vertices
        self.edges = edges
        if len(self.vertices) > 0:
            self.vertices.append(self.vertices[0])


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

    cdef readonly list inputPoints
    cdef readonly list inputSegments

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

        self.inputPoints = []
        self.inputSegments = []

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
        self.inputPoints.append([c_point.X, c_point.Y])

    def AddSegment(self, segment):
        """ Add a segment
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Segment c_segment = self._to_voronoi_segment(segment)
        self.thisptr.AddSegment(c_segment)
        self.inputSegments.append([[c_segment.p0.X, c_segment.p0.Y], [c_segment.p1.X, c_segment.p1.Y]])

    def Construct(self):
        """ Generates the voronoi diagram for the added points and segments. Voronoi cell structure will be generated.
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has already been called')

        self.constructed = 1
        with nogil:
            self.thisptr.Construct()

        self.thisptr.MapVertexIndexes()
        self.thisptr.MapEdgeIndexes()
        self.thisptr.MapCellIndexes()

    def GetPoint(self, index):
        return self.thisptr.GetPoint(index)

    def GetSegment(self, index):
        return self.thisptr.GetSegment(index)

    def GetVertex(self, index):
        """
        """
        if index < 0 or index >= self.CountVertices():
            raise IndexError(index)
        c_vertex = self.thisptr.GetVertex(index)
        return Vertex(c_vertex.X / self.SCALING_FACTOR, c_vertex.Y / self.SCALING_FACTOR)

    def GetEdge(self, index):
        if index < 0 or index >= self.CountEdges():
            raise IndexError(index)
        c_edge =  self.thisptr.GetEdge(index)
        edge = Edge(c_edge.start, c_edge.end, c_edge.cell, c_edge.twin)
        edge.is_primary = c_edge.isPrimary != False
        edge.is_linear = c_edge.isLinear != False
        return edge

    def GetCell(self, index):
        if index < 0 or index >= self.CountCells():
            raise IndexError(index)
        c_cell = self.thisptr.GetCell(index)
        cell = Cell(c_cell.cell_identifier, c_cell.site, c_cell.vertices, c_cell.edges, c_cell.source_category)
        cell.contains_point = c_cell.contains_point != False
        cell.contains_segment = c_cell.contains_segment != False
        cell.is_degenerate = c_cell.is_degenerate != False
        cell.is_open = c_cell.is_open != False
        return cell

    def CountVertices(self):
        return self.thisptr.CountVertices()

    def CountEdges(self):
        return self.thisptr.CountEdges()

    def CountCells(self):
        return self.thisptr.CountCells()

    def GetPoints(self):
        """ Returns the points added to the voronoi diagram
        """
        return self.inputPoints

    def GetSegments(self):
        """ Returns the segments added to the voronoi diagram
        """
        return self.inputSegments

    def GetVertices(self):
        count = self.CountVertices()
        output = []
        for index in  range(count):
            output.append(self.GetVertex(index))
        return output

    def GetIntersectingSegments(self):
        """
        Returns the indexes of segments that intersect another segment. The indexes are returned as a list.
        """
        return self.thisptr.GetIntersectingSegments()

    def GetDegenerateSegments(self):
        """
        Return the indexes of segments which has identical coordinates for its start point and end point. The indexes are returned as a list.
        """
        return self.thisptr.GetDegenerateSegments()

    def GetPointsOnSegments(self):
        """
        Return the indexes of points located on a segments. Connection at any of the end points is disregarded. The indexes are returned as a list.
        """
        return self.thisptr.GetPointsOnSegments()

    def GetEdges(self):
        count = self.CountEdges()
        output = []
        for index in range(count):
            output.append(self.GetEdge(index))
        return output

    def GetCells(self):
        count = self.CountCells()
        output = []
        for index in range(count):
            output.append(self.GetCell(index))
        return output

    def ReturnCurvedSiteInformation(self, edge):
        """Return the index of the point side and the segment site associated  with a segment index
        """
        twinEdge = self.GetEdge(edge.twin)

        cell = self.GetCell(edge.cell)
        twinCell = self.GetCell(twinEdge.cell)

        pointSite = self.RetrieveScaledPoint(cell) if cell.contains_point == True else self.RetrieveScaledPoint(twinCell)
        segmentSite = self.RetriveScaledSegment(twinCell) if cell.contains_point == True else self.RetriveScaledSegment(cell)

        return [pointSite, segmentSite]

    def DiscretizeCurvedEdge(self, index, max_dist, parabola_equation_tolerance = 0.0001):
        if(max_dist <= 0):
            raise ValueError("Max distance must be greater than 0. Value passed: {0}".format(max_dist))

        if(parabola_equation_tolerance < 0):
            raise ValueError("Parabola equation tolerance must be greater than 0 or equal to 0. Value passed: {0}".format(parabola_equation_tolerance))

        edge = self.GetEdge(index)
        sites = self.ReturnCurvedSiteInformation(edge)
        pointSite = sites[0]
        segmentSite = sites[1]

        edgeStartVertex = self.GetVertex(edge.start)
        edgeEndVertex = self.GetVertex(edge.end)
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


    def RetrieveScaledPoint(self, cell):
        non_scaled_point = self.RetrievePoint(cell)
        return [
			non_scaled_point[0] / self.SCALING_FACTOR,
			non_scaled_point[1] / self.SCALING_FACTOR
		]

    def RetriveScaledSegment(self, cell):
        non_scaled_segment = self.RetrieveSegment(cell)
        return [
			[non_scaled_segment[0][0] / self.SCALING_FACTOR, non_scaled_segment[0][1] / self.SCALING_FACTOR],
			[non_scaled_segment[1][0] / self.SCALING_FACTOR,non_scaled_segment[1][1] / self.SCALING_FACTOR]
		]

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
        return round(val * self.SCALING_FACTOR)

    cdef double _to_voronoi_double(self, val):
        return val * <double>self.SCALING_FACTOR

    cdef double _from_voronoi_value(self, val):
        return val / <double>self.SCALING_FACTOR
