"""
Cython wrapper for the C++ translation of the Voronoi diagram part of Boost's
Polygon library: http://www.boost.org/doc/libs/1_53_0_beta1/libs/polygon/doc/voronoi_diagram.htm
"""
SILENT = True
"""
Voronoi library operates with integer coordinates. To preserve the degree of
floating point precision use the SCALING_FACTOR with which all the coordinates and
relevant properties will be multiplied before used with the Voronoi library.
More info: http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Classes/ClipperOffset/Properties/ArcTolerance.htm"""
SCALING_FACTOR = 1

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
        int hasStart
        int hasEnd
        c_Vertex start
        c_Vertex end
        int isPrimary
        int site1
        int site2
        
    cdef cppclass VoronoiDiagram:
        VoronoiDiagram()
        void AddPoint(Point p)
        void AddSegment(Segment s)
        void Construct()
        vector[c_Edge] GetEdges()
        vector[Point] GetPoints()
        vector[Segment] GetSegments()
        

class Vertex:
    X = 0.0
    Y = 0.0

class Edge:
    has_start = False
    has_end = False
    start = Vertex()
    end = Vertex()
    is_primary = False
    site1 = -1
    site2 = -1

class VoronoiException(Exception):
    pass

cdef class Pyvoronoi:
    cdef VoronoiDiagram *thisptr
    cdef int constructed

    def __cinit__(self):
        """ Creates an instance of the Pyvoronoi class.
        """
        log_action("Creating an VoronoiDiagram instance")
        self.thisptr = new VoronoiDiagram()
        self.constructed = 0

    def __dealloc__(self):
        log_action("Deleting the VoronoiDiagram instance")
        del self.thisptr

    def AddPoint(self, point):
        """ Add a point
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Point c_point = _to_voronoi_point(point)
        self.thisptr.AddPoint(c_point)

    def AddSegment(self, segment):
        """ Add a segment
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Segment c_segment = _to_voronoi_segment(segment)
        self.thisptr.AddSegment(c_segment)

    def Construct(self):
        """ Generates the voronoi diagram for the added points and segments.
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has already been called')

        self.constructed = 1
        self.thisptr.Construct()

    def GetEdges(self):
        """ Returns the edges of the voronoi diagram.             
        """

        if self.constructed == 0:
            raise VoronoiException('Construct() has not been called')

        cdef vector[c_Edge] edges = self.thisptr.GetEdges()
        cdef size_t count = edges.size()
        outputEdges = []

        for i in range(count):
            edge = Edge()
            edge.has_start = edges[i].hasStart != False
            edge.has_end = edges[i].hasEnd != False

            start = Vertex()
            start.X = _from_voronoi_value(edges[i].start.X)
            start.Y = _from_voronoi_value(edges[i].start.Y)

            end = Vertex()
            end.X = _from_voronoi_value(edges[i].end.X)
            end.Y = _from_voronoi_value(edges[i].end.Y)

            edge.start = start
            edge.end = end
            edge.is_primary = edges[i].isPrimary

            edge.site1 = edges[i].site1
            edge.site2 = edges[i].site2

            outputEdges.append(edge)

        return outputEdges

    def GetPoints(self):
        """ Returns the points added to the voronoi diagram
        """
        cdef vector[Point] points = self.thisptr.GetPoints()
        cdef size_t count = points.size()
        outputPoints = []

        for i in range(count):
            point = [_from_voronoi_value(points[i].X), _from_voronoi_value(points[i].Y)]
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
            startPoint = [_from_voronoi_value(segments[i].p0.X), _from_voronoi_value(segments[i].p0.Y)]
            segment.append(startPoint)

            endPoint = [_from_voronoi_value(segments[i].p1.X), _from_voronoi_value(segments[i].p1.Y)]
            segment.append(endPoint)
            outputSegments.append(segment)

        return outputSegments


cdef Segment _to_voronoi_segment(object py_segment):
    return Segment(_to_voronoi_point(py_segment[0]), _to_voronoi_point(py_segment[1]))

cdef Point _to_voronoi_point(object py_point):
    return Point(_to_voronoi_int(py_point[0]), _to_voronoi_int(py_point[1]))

cdef int _to_voronoi_int(val):
    return val * SCALING_FACTOR

cdef double _to_voronoi_double(val):
    return val * <double>SCALING_FACTOR

cdef double _from_voronoi_value(val):
    return val / <double>SCALING_FACTOR