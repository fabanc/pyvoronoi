"""
Cython wrapper for the C++ translation of the Voronoi diagram part of Boost's
Polygon library: https://www.boost.org/doc/libs/1_80_0/libs/polygon/doc/voronoi_diagram.htm
"""
import struct
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
        vector[unsigned long long] GetIntersectingSegments()
        vector[unsigned long long] GetDegenerateSegments()
        vector[unsigned long long] GetPointsOnSegments()
        void MapVertexIndexes()
        void MapEdgeIndexes()
        void MapCellIndexes()

        long long CountPoints()
        long long CountSegments()
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
    """
    This class represents any vertex generated when constructing Boost Voronoi output. Those vertices are referenced by edges
    and cells.
    """
    def __init__(self,x: float,y: float) -> None:
        """
        Constructor
        :param x: The coordinate of the vertex on the x-axis.
        :type x: float
        :param y: The coordinate of the vertex on the y-axis.
        :type y: float
        """
        self.X = x
        self.Y = y

class Edge:
    """
    This class represents any edge generated when constructing Boost Voronoi output. Those edges are referenced by cells. Think
    of an edge as a linear component. It can be either straight lines or parabolas when returned by Boost Voronoi. An edge can
    only belong to one cell. That means that if two cells share a border, two edges will exist. Those two edges will be twins.
    """
    def __init__(self, start: int, end: int, cell: int, twin: int) -> None:
        """
        Constructor
        :param start: The index of the vertex at the start of the edge.
        :type start: int
        :param end: The index of the vertex at the start of the edge.
        :type end: int
        :param cell: The index of the cell this edge belong to.
        :type cell: int
        :param twin: The index of the twin edge.
        :type twin: int
        """
        self.start = start
        self.end = end
        self.cell = cell
        self.twin = twin

class Cell:
    """
    This class represents any cell generated when constructing Boost Voronoi output. A Voronoi cell represents a region of the Voronoi diagram bounded by the Voronoi edges.
    Voronoi cells are built around an input site that can be either site or segment. When a cells is on the border of the Voronoi solution, it might have
    infinite edges, meaning edges that do not have coordinates at one of their end point.

    The list of source categories for each cell can be found on the boost website.
     - SOURCE_CATEGORY_SINGLE_POINT: 0
     - SOURCE_CATEGORY_SEGMENT_START_POINT: 1
     - SOURCE_CATEGORY_SEGMENT_END_POINT: 2
     - SOURCE_CATEGORY_INITIAL_SEGMENT: 3
     - SOURCE_CATEGORY_REVERSE_SEGMENT: 4
     - SOURCE_CATEGORY_GEOMETRY_SHIFT: 5
     - SOURCE_CATEGORY_BITMASK: 6
    """
    def __init__(self, cell_identifier: int, site: int, vertices: list[int], edges: list[int], source_category: int) -> None:
        """
        :param cell_identifier:The identifier of the cell.
        :type cell_identifier: int
        :param site: The index of the input site that was used. Segments are added in the list of sites after points.
        :type site: int
        :param vertices: The list of indexes that make up the cell.
        :type vertices: list[int]
        :param edges: The list of edges that make up the cell.
        :type edges: list[int]
        :param source_category: The source category associated with the cell and its input site.
        :type source_category: int
        """
        self.cell_identifier = cell_identifier
        self.site = site
        self.source_category = source_category
        self.vertices = vertices
        self.edges = edges
        if len(self.vertices) > 0:
            self.vertices.append(self.vertices[0])

def _point_dict_to_point_array(p):
    """
    Used to read the input point from Boost Voronoi.
    :param p: The point as returned by Boost Voronoi
    :type p: Dictionary
    :return: A list representing the input point. The first element represents the X-coordinate. The second element represents the Y coordinate.
    :rtype: list[int, int]
    """
    return [p['X'], p['Y']]

def _segment_dict_to_point_array(s):
    """
    Used to read an input segment from Boost Voronoi.
    :param s: The segment as returned by Boost Voronoi
    :type s: Dictionary
    :return: A list representing the input segment. The first element represents the start point, the second element represents the end point..
    :rtype: list[int, int]
    """
    return[
        [s['p0']['X'], s['p0']['Y']],
        [s['p1']['X'], s['p1']['Y']]
    ]

class VoronoiException(Exception):
    """Indicate a problem generating a output with Boost Voronoi. This might indicate problem with your input data."""
    pass

class FocusOnDirectixException(Exception):
    """When computing a parabola between a point and a line, the input point must not be on the input line. If this scenario occurs, this exception is raised."""
    pass

class UnsolvableParabolaEquation(Exception):
    """Indicates a scenario where a parabola cannot be computed."""
    pass

####################################
##ROTATION
####################################
# TODO This function should be moved in a separate module. Conceptually, they have noting to do with Boost Voronoi and could be pure python
def rotate(point: list[float, float], theta: float) ->  list[float, float]:
    """Rotate a point around the origin using an angle

    :param point: the input point
    :type point: list[float, float]
    :param theta: the angle of rotation
    :type theta: float
    :return: the rotated point.
    :rtype: list[float, float]
    """
    t = -1 * theta
    cos = math.cos(t)
    sin = math.sin(t)
    return [(point[0] * cos) - (point[1] * sin),
	(point[0] * sin) + (point[1] * cos)]

# TODO This function should be moved in a separate module. Conceptually, they have noting to do with Boost Voronoi and could be pure python
def rotate_with_shift(point: list[float, float], theta: float, shift_x: float, shift_y: float) -> list[float, float]:
    """Rotate a point around the origin and shift it along the x-axis and the y-axis

    :param point: the input point
    :type point: list[float, float]
    :param theta: the angle of rotation around the origin
    :type theta: float
    :param shift_x: the translation to apply on the x-value.
    :type shift_x: float
    :param shift_y: the translation to apply on the y-value.
    :type shift_y: float
    :return: the rotated point.
    :rtype: list[float, float]
    """
    return rotate([point[0] - shift_x, point[1] - shift_y], theta)

# TODO This function should be moved in a separate module. Conceptually, they have noting to do with Boost Voronoi and could be pure python
def unrotate_with_shift(point: list[float, float], theta: float, shift_x: float, shift_y: float) -> list[float, float]:
    """ Undo the rotation and shift done by the function rotate_with_shift

    :param point: the input point
    :type point: list[float, float]
    :param theta: the angle of rotation around the origin used during the initial rotation
    :type theta: float
    :param shift_x: the translation to apply on the x-value to undo.
    :type shift_x: float
    :param shift_y: the translation to apply on the y-value to undo.
    :type shift_y: float
    :return: the unrotated point.
    :rtype: list[float, float]
    """
    cos = math.cos(theta)
    sin = math.sin(theta)
    return [(point[0] * cos) - (point[1] * sin) + shift_x, (point[0] * sin) + (point[1] * cos) + shift_y]

# TODO This function should be moved in a separate module. Conceptually, they have noting to do with Boost Voronoi and could be pure python
def get_line_angle_in_radians(start_point_x: float, start_point_y: float, end_point_x: float, end_point_y: float) -> float:
    """Get the orientation of a line in radians

    :param start_point_x: The x coordinate of the start point.
    :type start_point_x: float
    :param start_point_y: The y coordinate of the start point.
    :type start_point_y: float
    :param end_point_x: The x coordinate of the end point.
    :type end_point_x: float
    :param end_point_y: The y coordinate of the end point.
    :type end_point_y: float
    :return: The orientation of the line from the start point to the end point in radians
    :rtype: float
    """
    return math.atan2(end_point_y - start_point_y, end_point_x - start_point_x)


####################################
##DISTANCE
####################################
# TODO This function should be moved in a separate module. Conceptually, they have noting to do with Boost Voronoi and could be pure python
def get_distance_squared(point_start: list[float, float], point_end: list[float, float]) -> float:
    """Returns the squared distance between two points.

    :param point_start: the start point of the line as a list of coordinates
    :type point_start: list[float, float]
    :param point_end: the end point of the line as a list of coordinates
    :type point_end: list[float, float]
    :return: the squared distance between the two points
    :rtype: float
    """
    return pow(point_end[0] - point_start[0], 2) + pow(point_end[1] - point_start[1], 2)

# TODO This function should be moved in a separate module. Conceptually, they have noting to do with Boost Voronoi and could be pure python
def get_distance(point_start: list[float, float], point_end: list[float, float]) -> float:
    """Returns the distance between two points.

    :param point_start: the start point of the line.
    :type point_start: list[float, float]
    :param point_end: the end point of the line.
    :type point_end: list[float, float]
    :return: the distance between the two points
    :rtype: float
    """
    return math.sqrt(get_distance_squared(point_start, point_end))


####################################
##PYVORONOI OPERATIONS
####################################

cdef class Pyvoronoi:
    """
    The wrapper class around Boost Voronoi. Add input point or segments, and then call Construct to generate Voronoi
    cells.
    """
    cdef VoronoiDiagram *thisptr
    cdef int constructed

    cdef public int SCALING_FACTOR

    def __cinit__(self, scaling_factor:int = None) -> None:
        """
        Creates an instance of the Pyvoronoi class.
        :param scaling_factor: The scaling factor that can be used to multiply the coordinate of input points or segment. If null, the value will be forced to 1.
        :type scaling_factor: int
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

    def AddPoint(self, point:list[float, float])->None:
        """Add a point to the Voronoi builder. The coordinate will be multiplied by the factor, then rounded as an integer.

        :param point: A list with 2 elements. The first element represents the X-coordinate of the point.  The second element represents the Y-coordinate of the point.
        :type point: list[float, float]
        :return: None
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Point c_point = self._to_voronoi_point(point)
        self.thisptr.AddPoint(c_point)

    def AddSegment(self, segment: list[list[float, float], list[float, float]]) -> None:
        """Add a segment to the Voronoi builder. The segment is made of an array of coordinates. The coordinates will be multiplied by the factor, then rounded as an integer.

        :param segment: A list that contains two elements. Each element is an array containing the [X, Y] coordinate for each of the segment end points.
        :type segment: list[list[float, float], list[float, float]]
        :return: None
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has been called, can\'t add more elements')

        cdef Segment c_segment = self._to_voronoi_segment(segment)
        self.thisptr.AddSegment(c_segment)

    def Construct(self):
        """Generates the voronoi diagram for the added points and segments. Voronoi cell structure will be generated.
        Calling this method will prevent adding any new input point or segment.
        """

        if self.constructed == 1:
            raise VoronoiException('Construct() has already been called')

        self.constructed = 1
        with nogil:
            self.thisptr.Construct()

        self.thisptr.MapVertexIndexes()
        self.thisptr.MapEdgeIndexes()
        self.thisptr.MapCellIndexes()

    def GetPoint(self, index: int) -> list[int, int]:
        """Return an input point used to generate the voronoi diagram.

        :param index: The index of the point to retrieve. The function CountPoints can be used to retrieve the number of input points stored in memory.
        :type index: int
        :return: A list with two elements. The first element is the X coordinate of the input point. The second element is the Y coordinate.
        :rtype: list[int, int]
        """
        return _point_dict_to_point_array(self.thisptr.GetPoint(index))

    def GetSegment(self, index: int) -> list[list[int, int], list[int, int]]:
        """Return an input point segment to generate the voronoi diagram.

        :param index: The index of the segment to retrieve. The function CountSegments can be used to retrieve the number of input segments stored in memory.
        :type index: int
        :return: A list with two elements. The first element is a list representing the start point of the segment. The second element is a list representing the end point of the segment.
        :rtype: list[int, int]
        """
        return _segment_dict_to_point_array(self.thisptr.GetSegment(index))

    def GetVertex(self, index: int) -> Vertex:
        """Return  the output vertex at a given index. The list of vertex is generated upon calling Construct.

        :param index: The index of the vertex to retrieve.
        :type index: int
        :return: The matching vertex.
        :rtype: Vertex
        """
        if index < 0 or index >= self.CountVertices():
            raise IndexError(index)
        c_vertex = self.thisptr.GetVertex(index)
        return Vertex(c_vertex.X / self.SCALING_FACTOR, c_vertex.Y / self.SCALING_FACTOR)

    def GetEdge(self, index: int) -> Edge:
        """Return  the edge at a given index. The list of edge is generated upon calling Construct.

        :param index: The index of the edge to retrieve.
        :type index: int
        :return: The matching edge.
        :rtype: Edge
        """
        if index < 0 or index >= self.CountEdges():
            raise IndexError(index)
        c_edge =  self.thisptr.GetEdge(index)
        edge = Edge(c_edge.start, c_edge.end, c_edge.cell, c_edge.twin)
        edge.is_primary = c_edge.isPrimary != False
        edge.is_linear = c_edge.isLinear != False
        return edge

    def GetCell(self, index: int) -> Cell:
        """Return  the cell at a given index. The list of cells is generated upon calling Construct.

        :param index: The index of the cell to retrieve.
        :type index: int
        :return: The matching cell.
        :rtype: Cell
        """
        if index < 0 or index >= self.CountCells():
            raise IndexError(index)
        c_cell = self.thisptr.GetCell(index)
        cell = Cell(c_cell.cell_identifier, c_cell.site, c_cell.vertices, c_cell.edges, c_cell.source_category)
        cell.contains_point = c_cell.contains_point != False
        cell.contains_segment = c_cell.contains_segment != False
        cell.is_degenerate = c_cell.is_degenerate != False
        cell.is_open = c_cell.is_open != False
        return cell

    def CountPoints(self):
        """Returns the number of input points stored in memory to solve the Voronoi problem.

        :return: The number of input points.
        :rtype: int
        """
        return self.thisptr.CountPoints()

    def CountSegments(self):
        """Returns the number of input segments stored in memory to solve the Voronoi problem.

        :return: The number of input segments.
        :rtype: int
        """
        return self.thisptr.CountSegments()

    def CountVertices(self):
        """Returns the number of output Vertices generated as part of the Voronoi solution after having called Construct.

        :return: The number of vertices.
        :rtype: int
        """
        return self.thisptr.CountVertices()

    def CountEdges(self):
        """Returns the number of output Edges generated as part of the Voronoi solution after having called Construct.

        :return: The number of edges.
        :rtype: int
        """
        return self.thisptr.CountEdges()

    def CountCells(self):
        """Returns the number of output Cells generated as part of the Voronoi solution after having called Construct.

        :return: The number of cells.
        :rtype: int
        """
        return self.thisptr.CountCells()

    def GetPoints(self):
        """Iterate through the input points added to the voronoi builder

        :return: A generator that iterates through the list of input points.
        :rtype: Generator[list[int, int]]
        """
        for p in self.thisptr.GetPoints():
            yield _point_dict_to_point_array(p)


    def GetSegments(self):
        """Iterate through the input segments added to the voronoi builder

        :return: A generator that iterates through the list of input segments.
        :rtype: Generator[list[int, int]]
        """
        for s in self.thisptr.GetSegments():
            yield _segment_dict_to_point_array(s)

    def GetIntersectingSegments(self):
        """Input Data Validation
            Returns the indexes of segments that intersect another segment beyond sharing an end point.
            Those segments can prevent the voronoi algorithm from solving, or generate a corrupted output.
            The indexes are returned as a list.

        :return: A list of indexes.
        :rtype: list[int]
        """
        return self.thisptr.GetIntersectingSegments()

    def GetDegenerateSegments(self):
        """Input Data Validation
            Return the indexes of segments which has identical coordinates for its start point and end point.
            Those segments can prevent the voronoi algorithm from solving, or generate a corrupted output.
            The indexes are returned as a list.

        :return: A list of indexes.
        :rtype: list[int]
        """
        return self.thisptr.GetDegenerateSegments()

    def GetPointsOnSegments(self):
        """Input Data Validation
            Return the indexes of points located on a segments. Connection at any of the end points is disregarded.
            Those situations can prevent the voronoi algorithm from solving, or generate a corrupted output.
            The indexes are returned as a list.

        :return: A list of indexes.
        :rtype: list[int]
        """
        return self.thisptr.GetPointsOnSegments()

    def GetVertices(self):
        """Get the list of vertices generated after calling construct. This returns a duplicated list of the output Vertices.
        Consider using EnumerateVertices instead.

        :return: A copy of the output vertices.
        :rtype: list[Vertex]
        """
        count = self.CountVertices()
        output = []
        for index in  range(count):
            output.append(self.GetVertex(index))
        return output

    def EnumerateVertices(self):
        """Iterate through the list of output vertices generated after calling Construct.

        :return: A generator iterating through the output vertices. Each object is tuple with two elements.
            The first one is the index. The second element is the vertex.
        :rtype: Generator[(int, Vertex)]
        """
        for index in range(self.CountVertices()):
            yield index, self.GetVertex(index)

    def GetEdges(self):
        """Get the list of edges generated after calling construct. This returns a duplicated list of the output edges.
        Consider using EnumerateEdges instead.

        :return: A copy of the output edges.
        :rtype: list[Edge]
        """
        count = self.CountEdges()
        output = []
        for index in range(count):
            output.append(self.GetEdge(index))
        return output

    def EnumerateEdges(self):
        """Iterate through the list of output edges generated after calling Construct.

        :return: A generator iterating through the output edges. Each object is tuple with two elements.
            The first one is the index. The second element is the edge.
        :rtype: Generator[(int, Edge)]
        """
        for index in range(self.CountEdges()):
            yield index, self.GetEdge(index)

    def GetCells(self):
        """Get the list of cells generated after calling construct. This returns a duplicated list of the output cells.
        Consider using EnumerateCells instead.

        :return: A copy of the output cells.
        :rtype: list[Cell]
        """
        count = self.CountCells()
        output = []
        for index in range(count):
            output.append(self.GetCell(index))
        return output

    def EnumerateCells(self):
        """Iterate through the list of output cells generated after calling Construct.

        :return: A generator iterating through the output edges. Each object is tuple with two elements.
            The first one is the index. The second element is the cell.
        :rtype: Generator[(int, Cell)]
        """
        for index in range(self.CountCells()):
            yield index, self.GetCell(index)

    def ReturnCurvedSiteInformation(self, edge: Edge) -> list[int, int]:
        """Returns the index of the input point site and the segment site associated  with a segment index.

        :param edge: The edge for which to retrieve information.
        :type edge: Edge
        :return: A list with two elements.  The first element is the input point. The second element is the input segment.
        :rtype: list[list[float, float], list[list[float, float], list[float, float]]]
        """
        twinEdge = self.GetEdge(edge.twin)

        cell = self.GetCell(edge.cell)
        twinCell = self.GetCell(twinEdge.cell)

        pointSite = self.RetrieveScaledPoint(cell) if cell.contains_point == True else self.RetrieveScaledPoint(twinCell)
        segmentSite = self.RetrieveScaledSegment(twinCell) if cell.contains_point == True else self.RetrieveScaledSegment(cell)

        return [pointSite, segmentSite]

    def DiscretizeCurvedEdge(self, index: int, max_dist: float, parabola_equation_tolerance = 0.0001) -> map[list[float, float]]:
        """Returns a list of point that represent the parabola given the index of an edge. This is a convenience wrapper
        around the function Discretize

        :param index: The index of the edge discretize.
        :type index: int
        :param max_dist: The maximum distance between points.
        :type max_dist: float
        :param parabola_equation_tolerance: The maximum difference allowed between the y coordinate returned by Boost, and the equation of the parabola.
        :type parabola_equation_tolerance: float
        :return: The list of points on the parabola returned as a map.
        :rtype: map[list[float, float]]
        """
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


    def RetrievePoint(self, cell: Cell) -> list[int, int]:
        """Retrieve the input point associated with a cell. The point coordinates are as used by Boost Voronoi: they have been multiplied by the factor.

		:param cell: the cell that contains a point. The point can either an input point or the end point of an input segment.
		:type cell: Cell
		:return: An input point
		:rtype: list[int, int]
        """
        if(cell.source_category == 0):
            return self.GetPoint(cell.site)

        input_segment = self.RetrieveSegment(cell)
        if(cell.source_category == 1):
            return input_segment[0]
        else:
            return input_segment[1]

    def RetrieveSegment(self, cell: Cell) -> list[list[int, int], list[int, int]]:
        """Retrieve the input segment associated with a cell. The segment coordinates are as used by Boost Voronoi: they have been multiplied by the factor.

		:param cell: the cell that contains a segment.
		:type cell: Cell
		:return: An input segment
		:rtype: list[list[int, int], list[int, int]]
        """
        i = cell.site - self.thisptr.CountPoints()
        s = self.GetSegment(i)
        return s

    def RetrieveScaledPoint(self, cell: Cell) -> list[float, float]:
        """Retrieve the input point associated with a cell. The point coordinates are returned as passed to Boost Voronoi, before applying the factor.

		:param cell: the cell that contains a point. The point can either an input point or the end point of an input segment.
		:type cell: Cell
		:return: An input point
		:rtype: list[int, int]
        """
        non_scaled_point = self.RetrievePoint(cell)
        return [
			non_scaled_point[0] / self.SCALING_FACTOR,
			non_scaled_point[1] / self.SCALING_FACTOR
		]

    def RetrieveScaledSegment(self, cell: Cell) -> list[list[float, float], list[float, float]]:
        """Retrieve the input segment associated with a cell. The segment coordinates are returned as passed to Boost Voronoi, before applying the factor.

		:param cell: the cell that contains a segment.
		:type cell: Cell
		:return: An input segment
		:rtype: list[list[int, int], list[int, int]]
        """
        non_scaled_segment = self.RetrieveSegment(cell)
        return [
			[non_scaled_segment[0][0] / self.SCALING_FACTOR, non_scaled_segment[0][1] / self.SCALING_FACTOR],
			[non_scaled_segment[1][0] / self.SCALING_FACTOR,non_scaled_segment[1][1] / self.SCALING_FACTOR]
		]

    def GetParabolaY(self, x: float, focus: list[float, float], directrix_y: float) -> float:
        """Used for interpolating parabolas
            Solve the parabola equation for a given value on the x-axis and return the associated value on the y-axis. This equation assumes that the directix is parallel to the x-axis. Parabola equation are different if the directix is parallel to the y-axis.

        :param x: The position of a point on the X-axis. The Y value is derived based on this coordinate.
        :type x: float
        :param focus: The focus point used to solve the equation.
        :type focus: list[float, float]
        :param directrix_y: The value on the Y-axis of the directix.
        :type directrix_y: float
        :return: the value on the y-axis derivated for x
        :rtype: float
        """
        return (pow(x - focus[0], 2) + pow(focus[1], 2) - pow(directrix_y, 2)) / (2 * (focus[1] - directrix_y))

    def CheckUnsolvableParabolaEquation(self, boost_x: float, boost_y: float, focus: list[float, float], directix: float, tolerance: float) -> list[float, float]:
        """Used for interpolating parabolas
            Compare the y-coordinate of a point on the parabola returned by Boost with the computed value.
		    The function will return an exception if the difference between the computed y-value and the y-value returned by Boost.
		    The computed point will be returned otherwise.

        :param boost_x: the x-value of the point parabola returned by boost.
        :type boost_x: float
        :param boost_y: the y-value of the point parabola returned by boost.
        :type boost_y: float
        :param focus: the focus point used for solving the equation of the parabola.
        :type focus: list[float, float]
        :param directix: the directix value used for solving the equation of the parabola.
        :type directix: float
		:param tolerance: the distance allowed between the point computed by boost and the point computed by the equation.
		:type tolerance: float
		:raises UnsolvableParabolaEquation: Indicates a scenario where a parabola cannot be computed.
		:return: the point on the parabola computed using the value of boost_x.
		:rtype: list[float, float]
		"""

        computed_point_y = self.GetParabolaY(boost_x, focus, directix)
        delta = computed_point_y - boost_y if computed_point_y > boost_y else boost_y - computed_point_y
        if delta > tolerance:
            raise UnsolvableParabolaEquation("The computed Y on the parabola for the starting / ending point is different from the rotated point returned by Boost. Difference: {0}. Maximum tolerance: {1}".format(delta, tolerance))
        return [boost_x, computed_point_y]

    def Discretize(self, point: list[float, float], segment: list[list[float, float], list[float, float]], parabola_start: list[float, float], parabola_end: list[float, float], max_dist: float, parabola_equation_tolerance: float) -> map[list[float, float]]:
        """Interpolate points on a parabola. The points are garanteed to be closer than the value of the parameter max_dist.

        :param point: The input point associated with the cell or the neighbour cell. The point is used as the focus in the equation of the parabola.
        :type point: list[float, float]
        :param segment: The input segment associated with the cell or the neighbour cell. The point is used as the directix in the equation of the parabola.
        :type segment: list[float, float]
        :param parabola_start: The starting point of the parabola.
        :type parabola_start: list[float, float]
        :param parabola_end: The end point of the parabola
        :type parabola_end: list[float, float]
        :param max_dist: The maximum distance between 2 vertices on the discretized geometry.
        :type max_dist: float
        :param parabola_equation_tolerance: The maximum difference allowed between the y coordinate returned by Boost, and the equation of the parabola.
        :type parabola_equation_tolerance: float
        :return: The list of points on the parabola.
        :rtype: map[list[float, float]]
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
        angle = get_line_angle_in_radians(
            segment[0][0],
            segment[0][1],
            segment[1][0],
            segment[1][1],
        )

        focus_rotated = rotate_with_shift(point, angle, shift_x, shift_y)
        directix_start_rotated = rotate_with_shift(segment[0], angle, shift_x, shift_y)
        directix_end_rotated = rotate_with_shift(segment[1], angle, shift_x, shift_y)
        parabola_start_rotated = rotate_with_shift(parabola_start, angle, shift_x, shift_y)
        parabola_end_rotated = rotate_with_shift(parabola_end, angle, shift_x, shift_y)

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
            distance = get_distance(current, previous)

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
        densified = map(lambda x: unrotate_with_shift(x, angle, shift_x, shift_y), densified_rotated)
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
