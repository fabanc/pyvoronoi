"""
Cython wrapper for the C++ translation of the Voronoi diagram part of Boost's
Polygon library: http://www.boost.org/doc/libs/1_53_0_beta1/libs/polygon/doc/voronoi_diagram.htm
"""
SILENT = False

def log_action(description):
    if not SILENT:
        print description

log_action("Python binding clipper library")


import sys as _sys
import struct
import copy as _copy
import unicodedata as _unicodedata
import time as _time

from cython.operator cimport dereference as deref

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
				
    cdef cppclass VoronoiDiagram:
        VoronoiDiagram()
        void AddPoint(Point p)
        void AddSegment(Segment s)
        void Construct()
        void GetEdges(vector[c_Vertex] & c_vertices, vector[c_Edge] & c_edges)
        void GetCells(vector[c_Vertex] & c_vertices, vector[c_Edge] & c_edges, vector[c_Cell] & c_cell_edges)
        vector[Point] GetPoints()
        vector[Segment] GetSegments()

        
	
class Vertex:
    X = 0.0
    Y = 0.0

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

    def AddSegment(self, segment):
        """ Add a segment
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Segment c_segment = self._to_voronoi_segment(segment)
        self.thisptr.AddSegment(c_segment)

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
        print "Get cells information..."
        self.thisptr.GetCells(c_vertices, c_edges, c_cells)		
        print "Done!"
		
        del self.inputPoints[:]
        del self.inputSegments[:]
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
            vertex = Vertex()
            vertex.X = self._from_voronoi_value(c_vertices[i].X)
            vertex.Y = self._from_voronoi_value(c_vertices[i].Y)

            self.outputVertices.append(vertex)

        count = c_cells.size()

        for i in range(count):
            c_cell = c_cells[i]
            outputCell = Cell(c_cell.cell_identifier, c_cell.site, c_cell.vertices, c_cell.edges)
            outputCell.contains_point = c_cell.contains_point != False
            outputCell.contains_segment = c_cell.contains_segment != False
            outputCell.is_open = c_cell.is_open != False
            outputCell.source_category = c_cell.source_category
            self.outputCells.append(outputCell)
			
        self.inputPoints = self.GetPoints()
        self.inputSegments = self.GetSegments()
				
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
        cdef vector[Point] points = self.thisptr.GetPoints()
        cdef size_t count = points.size()
        outputPoints = []

        for i in range(count):
            point = [self._from_voronoi_value(points[i].X), self._from_voronoi_value(points[i].Y)]
            outputPoints.append(point)

        return outputPoints

    def GetSegments(self):
        """ Returns the segments added to the voronoi diagram
        """
        cdef vector[Segment] segments = self.thisptr.GetSegments()
        cdef size_t count = segments.size()
        outputSegments = []

        for i in range(count):
            segment = []
            startPoint = [self._from_voronoi_value(segments[i].p0.X), self._from_voronoi_value(segments[i].p0.Y)]
            segment.append(startPoint)

            endPoint = [self._from_voronoi_value(segments[i].p1.X), self._from_voronoi_value(segments[i].p1.Y)]
            segment.append(endPoint)
            outputSegments.append(segment)

        return outputSegments
		


    def RetrievePoint(self, cell):
        #input point
        if(cell.source_category == 0):
            return self.inputPoints[cell.site]
        elif(cell.source_category == 1):
            return self.inputSegments.startPoint
        else:			
            return self.inputSegments.endPoint
		
    def RetrieveSegment(self, cell):
        """Retrive the segment associated with a cell.
        """
        return cell.site - len(self.inputPoints) if cell.contains_segment == 1 else None
		
    def GetParabolaY(self, x,a,b):
        return ((x - a) * (x - a) + b * b) / (b + b)		
	
    def Discretize(self, point, segment, max_dist, discretization):
        """
        Discretize
            :param point: The input point associated with the cell or the neighbour cell
            :param segment: The input segment associated with the cell or the neighbour cell
            :param max dist: The maximum distance between 2 vertices on the discretized geometry
            :param discretization: The curved output edge		
		"""
        low_segment_x = segment.x1 if segment.x1 < segment.x2 else segment.x2
        low_segment_y = segment.y1 if segment.y1 < segment.y2 else segment.y2
    
        max_segment_x = segment.x2 if segment.x1 < segment.x2 else segment.x1
        max_segment_y = segment.y2 if segment.y1 < segment.y2 else segment.y1
    
    
        # Apply the linear transformation to move start point of the segment to
        # the point with coordinates (0, 0) and the direction of the segment to
        # coincide the positive direction of the x-axis.
        segm_vec_x = max_segment_x - low_segment_x
        segm_vec_y = max_segment_y - low_segment_y
        sqr_segment_length = segm_vec_x * segm_vec_x + segm_vec_y * segm_vec_y;
    
        #Compute x-coordinates of the endpoints of the edge
        ##in the transformed space.
        projection_start = sqr_segment_length * self.GetPointProjection(discretization.start, segment);
        projection_end = sqr_segment_length * self.GetPointProjection(discretization.end, segment);
    
        #Compute parabola parameters in the transformed space.
        #Parabola has next representation:
        ##f(x) = ((x-rot_x)^2 + rot_y^2) / (2.0*rot_y).
        point_vec_x = point.x - segment.x1 if segment.x1 < segment.x2 else point.x - segment.x2;
        point_vec_y = point.y - segment.y1 if segment.y1 < segment.y2 else point.y - segment.y2;
        rot_x = segm_vec_x * point_vec_x + segm_vec_y * point_vec_y;
        rot_y = segm_vec_x * point_vec_y - segm_vec_y * point_vec_x;
    
        #Save the last point.
        last_point = discretization.end;
        discretization_array = [discretization.start]
    
        #Use stack to avoid recursion.
        point_stack = [projection_end]
        cur_x = projection_start
        cur_y = self.GetParabolaY(cur_x, rot_x, rot_y)
    
        #Adjust max_dist parameter in the transformed space.
        max_dist_transformed = max_dist * max_dist * sqr_segment_length;
    
        while (len(point_stack) != 0):
            #Get last element on the stack
            new_x = point_stack[-1]
            new_y = self.GetParabolaY(new_x, rot_x, rot_y)
    
            #Compute coordinates of the point of the parabola that is
            #furthest from the current line segment.
            mid_x = (new_y - cur_y) / (new_x - cur_x) * rot_y + rot_x
            mid_y = self.GetParabolaY(mid_x, rot_x, rot_y)
    
            dist = (new_y - cur_y) * (mid_x - cur_x) - (new_x - cur_x) * (mid_y - cur_y)
            dist = dist * dist / ((new_y - cur_y) * (new_y - cur_y) + (new_x - cur_x) * (new_x - cur_x))
    
            if (dist <= max_dist_transformed):
                #Distance between parabola and line segment is less than max_dist.
                point_stack.pop()
                inter_x = (segm_vec_x * new_x - segm_vec_y * new_y) / sqr_segment_length + low_segment_x
                inter_y = (segm_vec_x * new_y + segm_vec_y * new_x) / sqr_segment_length + low_segment_y
                vertex = Vertex()
                vertex.X = inter_x
                vertex.Y = inter_y
                discretization_array.append(vertex)
                cur_x = new_x
                cur_y = new_y
            else:
                point_stack.append
    
        discretization[-1] = last_point
		
		
    def GetPointProjection(self, point, segment):
        """
        Get normalized length of the distance between:
          1) point projection onto the segment
          2) start point of the segment
        Return this length divided by the segment length. This is made to avoid
        sqrt computation during transformation from the initial space to the
        transformed one and vice versa. The assumption is made that projection of
        the point lies between the start-point and endpoint of the segment.	
        """
        segment_vec_x = segment.x2 - segment.x1 if segment.x1 < segment.x2 else segment.x1 - segment.x2
        segment_vec_y = segment.y2 - segment.y1 if segment.y1 < segment.y2 else segment.y1 - segment.y2
        point_vec_x = point.x - segment.x1 if segment.x1 < segment.x2 else point.x - segment.x2
        point_vec_y = point.y - segment.y1 if segment.y1 < segment.y2 else point.y - segment.y2
        sqr_segment_length = segment_vec_x * segment_vec_x + segment_vec_y * segment_vec_y
        vec_dot = segment_vec_x * point_vec_x + segment_vec_y * point_vec_y
        return vec_dot / sqr_segment_length
    
    

		
		

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