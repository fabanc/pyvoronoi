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

    cdef struct c_Edge2:
        double x1
        double y1
        double x2
        double y2		
        int isPrimary
        size_t site
        int isLinear
		
    cdef struct c_CellEdge:
        size_t cellId
        size_t source_index
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
        void GetCellEdges(vector[c_Vertex] & c_vertices, vector[c_Edge] & c_edges, vector[c_CellEdge] & c_cell_edges)
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
	
class Edge2:
    x1 = -1
    y1 = -1
    x2 = -1
    y2 = -1	
    is_primary = False
    site = -1
    is_linear = False

    def __init__(self, x1,y1,x2,y2, is_primary, site, is_linear):
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.is_primary = is_primary
        self.site = site
        self.is_linear = is_linear
		
class VoronoiCell:
    cellId = -1
    source_index = -1
    contains_point = 0
    contains_segment = 0
    is_open = 0
	
    vertices = None
    segments = None

    def __init__(self,cellId, source_index, vertices, segments):
        self.cellId = cellId
        self.source_index = source_index
        self.vertices = vertices
        self.segments = segments
        self.vertices.append(self.vertices[0])

class VoronoiException(Exception):
    pass

cdef class Pyvoronoi:
    cdef VoronoiDiagram *thisptr
    cdef int constructed

    outputEdges = []
    outputVertices = []
	
    outputCellEdges = []
    outputCellEdges2 = []
    outputCellVertices = []	
    outputCells = []
    outputSortedCells = []	

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
        """ Generates the voronoi diagram for the added points and segments.
        """
        
        if self.constructed == 1:
            raise VoronoiException('Construct() has already been called')

        self.constructed = 1
        print "Construct..."
        self.thisptr.Construct()
        
        cdef vector[c_Edge] c_edges
        cdef vector[c_Vertex] c_vertices
		
        cdef vector[c_Vertex] c_vertices2
        cdef vector[c_Edge] c_edges2
        cdef vector[c_CellEdge] c_cell_edges		
        print "Get cells information..."
        self.thisptr.GetCellEdges(c_vertices2, c_edges2, c_cell_edges)		
        print "Done!"
		
        #self.thisptr.GetEdges(c_vertices, c_edges)
		

        #cdef size_t count = c_edges.size()
        
        del self.outputEdges[:] 
        del self.outputVertices[:]
		
        del self.outputCellEdges[:] 
        del self.outputCellEdges2[:] 
        del self.outputCellVertices[:]		
        del self.outputCells[:]
        del self.outputSortedCells[:]
		
		
        #for i in range(count):
        #    edge = Edge()
        #    edge.start = c_edges[i].start
        #    edge.end = c_edges[i].end
        #    edge.is_primary = c_edges[i].isPrimary != False
        #    edge.is_linear = c_edges[i].isLinear != False
        #
        #    edge.site1 = c_edges[i].site1
        #    edge.site2 = c_edges[i].site2
        #
        #    self.outputEdges.append(edge)
        #
        #count = c_vertices.size()
        #
        #for i in range(count):
        #    vertex = Vertex()
        #    vertex.X = self._from_voronoi_value(c_vertices[i].X)
        #    vertex.Y = self._from_voronoi_value(c_vertices[i].Y)
        #
        #    self.outputVertices.append(vertex)
			
		

        #--------------------------------------------------------
        #LOGIC TO PRINT CELL OBJECTS
        #--------------------------------------------------------
        count = c_edges2.size()
        print "Cell edges size = " + str(count)
        for i in range(count):
            edge = Edge()
            edge.start = c_edges2[i].start
            edge.end = c_edges2[i].end
            edge.is_primary = c_edges2[i].isPrimary != False
            edge.is_linear = c_edges2[i].isLinear != False
        
            edge.site1 = c_edges2[i].site1
            edge.site2 = c_edges2[i].site2
        
            self.outputCellEdges.append(edge)
			
					
        count = c_vertices2.size()
        print "Cell vertex size = " + str(count)        
        for i in range(count):
            vertex = Vertex()
            vertex.X = self._from_voronoi_value(c_vertices2[i].X)
            vertex.Y = self._from_voronoi_value(c_vertices2[i].Y)

            self.outputCellVertices.append(vertex)
				
        count = c_cell_edges.size()
        print "Number of cells: {0}".format(count)
		
		
        for i in range(count):
            c_cell = c_cell_edges[i]
            self.outputCells.append(VoronoiCell(c_cell.cellId, c_cell.source_index, c_cell.vertices, c_cell.edges))
            self.outputCells[-1].contains_point = c_cell.contains_point != False
            self.outputCells[-1].contains_segment = c_cell.contains_segment != False
            self.outputCells[-1].is_open = c_cell.is_open != False		
			
        #print "Building Cells"
        #for i in range(count):
        #    cedge = c_cell_edges[i]
        #    #print "Edge index: {0}".format(cedge.edgeId)
        #    edge = self.outputCellEdges[cedge.edgeId]
        #    #print "Site: {0}, Edge Index: {1}, Start: {2}, End: {3}, Is Linear: {4}".format(cedge.cellId,i,edge.start, edge.end, edge.is_linear)
        #    if (i==0):
        #        cell = VoronoiCell(cellId)
        #        cell.segments.append(i)
        #        self.outputCells.append(cell)
        #        cellId = cellId + 1				
        #    else:
        #        if(self.outputCells[-1].cellId != cedge.cellId):
        #            cell = VoronoiCell(cedge.cellId)
        #            cell.segments.append(i)
        #            self.outputCells.append(cell)
        #            cellId = cellId + 1
        #        else:
        #            self.outputCells[-1].segments.append(i)
        #
		#			
        #print "Sorting Cells"
        #sortedCells = []
        #for c in self.outputCells:
        #    print "Cell: {0}".format(c.cellId)
        #    sortedSegments = self.SortCell(c.segments)
        #    if(sortedSegments != None):
        #        c.segments = sortedSegments
        #        sortedCells.append(c)
		#		
        #self.outputCells = sortedCells
				
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
		
    def GetCellVertices(self):
        """ Returns the edges of the voronoi diagram.             
        """

        if self.constructed == 0:
            raise VoronoiException('Construct() has not been called')

        return self.outputCellVertices

    def GetCellEdges(self):
        """ Returns the edges of the voronoi diagram.             
        """

        if self.constructed == 0:
            raise VoronoiException('Construct() has not been called')
        
        return self.outputCellEdges

    def GetCells(self):
        """ Return pointer to the edges that make the cellss            
        """	
        return self.outputCells

    def SortCell(self, cellSegments):
	
        #Default checks
        if len(cellSegments) == 0:
            return None
        if len(cellSegments) == 1:
            if(self.outputCellEdges[cellSegments[0]].start == -1 or self.outputCellEdges[cellSegments[0]].end == -1):
                return None
            else:
                return cellSegments
				
        #Sorting logic		
        sortedList = [cellSegments.pop()]
			
        while len(cellSegments) > 0:
            count = len(sortedList)
            lastSortedEdge = self.outputCellEdges[sortedList[-1]]
            #print "  Last Sorted Edge. Start: {0}, End: {1}".format(lastSortedEdge.start, lastSortedEdge.end)
            for i in range(len(cellSegments)):
                edge = self.outputCellEdges[cellSegments[i]]
                #print "     Edge. Start: {0}, End: {1}".format(edge.start, edge.end)
                if edge.start == -1 or edge.end == -1:
                     return None
                if lastSortedEdge.end == edge.start or lastSortedEdge.end == edge.end:
                    sortedList.append(cellSegments.pop(i))
                    break
            if len(sortedList) == count:
                return None#raise Exception ("Error while sorting segments")	

        #for e in sortedList:
        #    tEdge = self.outputCellEdges[e]
        #    print "Sorted List Element. Start element: {0}, End element: {1}".format(tEdge.start, tEdge.end)

        return sortedList			

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