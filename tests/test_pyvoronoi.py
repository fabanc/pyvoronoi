#!/usr/bin/python
"""
Tests for Pyvoronoi wrapper library.
"""

from __future__ import print_function
from unittest import TestCase, main

import pyvoronoi

class TestPyvoronoiModule(TestCase):
    def test_has_classes(self):
        self.assertTrue(hasattr(pyvoronoi, 'Pyvoronoi'))

    def test_scaling_factor(self):
        pv = pyvoronoi.Pyvoronoi()
        self.assertTrue(pv.SCALING_FACTOR == 1)

        pv = pyvoronoi.Pyvoronoi(10)
        self.assertTrue(pv.SCALING_FACTOR == 10)

        pv = pyvoronoi.Pyvoronoi(1)
        self.assertTrue(pv.SCALING_FACTOR == 1)

class TestPyvoronoiAdd(TestCase):
    def test_add_point(self):
        factor = 10
        inputPoint = [0.5, 1]
        pv = pyvoronoi.Pyvoronoi(factor)
        pv.AddPoint(inputPoint)
        points = pv.GetPoints()
        self.assertTrue(len(points) == 1)
        self.assertTrue(points[0][0] == inputPoint[0] * factor)
        self.assertTrue(points[0][1] == inputPoint[1] * factor)

    def test_add_segment(self):
        factor = 10
        segment = [[0.5, 1], [0, 2]]
        pv = pyvoronoi.Pyvoronoi(factor)
        pv.AddSegment(segment)
        segments = pv.GetSegments()
        self.assertTrue(len(segments) == 1)
        self.assertTrue(segments[0] == [
            [segment[0][0] * factor, segment[0][1] * factor],
            [segment[1][0] * factor, segment[1][1] * factor],
        ])

    def test_add_point_after_construct(self):
        pv = pyvoronoi.Pyvoronoi()
        pv.Construct()
        self.assertRaises(pyvoronoi.VoronoiException, pv.AddPoint, [0, 0])

    def test_add_segment_rounding(self):
        factor = 10
        segment = [[0.59, 1.09], [0, 2]]
        pv = pyvoronoi.Pyvoronoi(factor)
        pv.AddSegment(segment)
        segments = pv.GetSegments()
        self.assertTrue(len(segments) == 1)
        self.assertTrue(segments[0] == [
            [round(segment[0][0] * factor), round(segment[0][1] * factor)],
            [segment[1][0] * factor, segment[1][1] * factor],
        ])

class TestPyvoronoiConstruct(TestCase):
    def test_square(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddSegment([[0, 0], [0, 1]])
        pv.AddSegment([[0, 1], [1, 1]])
        pv.AddSegment([[1, 1], [1, 0]])
        pv.AddSegment([[1, 0], [0, 0]])
        pv.Construct()
        edges = pv.GetEdges()
        vertices = pv.GetVertices()
        cells = pv.GetCells()
        self.assertTrue(len(cells) == 8)
        self.assertTrue(len([i for i in edges if i.is_primary == True]) == 8)
        self.assertTrue(len(vertices) == 5)
        self.assertTrue(len(cells[0].edges) == 2)

    def test_rectangle(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddSegment([[0, 0], [0, 2]])
        pv.AddSegment([[0, 2], [1, 2]])
        pv.AddSegment([[1, 2], [1, 0]])
        pv.AddSegment([[1, 0], [0, 0]])
        pv.Construct()
        edges = pv.GetEdges()
        vertices = pv.GetVertices()
        cells = pv.GetCells()
        self.assertTrue(len(cells) == 8)
        self.assertTrue(len([i for i in edges if i.is_primary == True]) == 10)
        self.assertTrue(len(vertices) == 6)
        self.assertTrue(len(list(filter(lambda e: edges[e].is_primary, cells[1].edges))) == 3)
        self.assertTrue(len(list(filter(lambda e: edges[e].is_primary, cells[3].edges))) == 2)

    def test_twins(self):
        """
        Validate that the twin attribute is consistent.
        :return: 
        """
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.AddSegment([[0,0],[0,10]])
        pv.AddSegment([[0,0],[10,0]])
        pv.AddSegment([[0,10],[10,10]])
        pv.AddSegment([[10,0],[10,10]])
        pv.Construct()
        edges = pv.GetEdges()
        for i in range(len(edges)):
            edge = edges[i]
            self.assertTrue(edges[edge.twin].twin == i)
        cells = pv.GetCells()

    def test_cells_vertices_duplication(self):
        """
        Validate that the first and last vertex are the same on cells.
        :return:
        """
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.AddSegment([[0,0],[0,10]])
        pv.AddSegment([[0,0],[10,0]])
        pv.AddSegment([[0,10],[10,10]])
        pv.AddSegment([[10,0],[10,10]])
        pv.Construct()
        cells = pv.GetCells()
        vertices = pv.GetVertices()
        cell = cells[5]
        self.assertNotEqual(cell.vertices[-2], cell.vertices[-1])
        self.assertEqual(cell.vertices[0], cell.vertices[-1])

    def test_vertex_reference_for_edges(self):
        """
        Test the node edge have both ends not referencing a vertex.
        :return: 
        """
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.AddSegment([[0,0],[0,10]])
        pv.AddSegment([[0,0],[10,0]])
        pv.AddSegment([[0,10],[10,10]])
        pv.AddSegment([[10,0],[10,10]])
        pv.Construct()
        edges = pv.GetEdges()
        cells = pv.GetCells()
        for i in range(len(edges)):
            edge = edges[i]
            if cells[edge.cell].is_open == False:
                self.assertTrue(edge.start != -1 and edge.end != -1)

    def test_edge_vertices_indexes(self):
        """
        Test the node edge have both ends not referencing a vertex.
        :return: 
        """
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddSegment([[0,0],[0,10]])
        pv.AddSegment([[0,0],[10,0]])
        pv.AddSegment([[0,10],[10,10]])
        pv.AddSegment([[10,0],[10,10]])
        pv.Construct()
        edges = pv.GetEdges()

        vertices_count = 0

        for i in range(len(edges)):
            edge = edges[i]
            vertices = [edge.start, edge.end]
            for v in vertices:
                if v == -1:
                    vertices_count += 1

        self.assertTrue(vertices_count == 16)

    def test_vertex_reference_for_cells(self):
        """
        Test that cells reference at least one vertex.
        :return: 
        """
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.AddSegment([[0,0],[0,10]])
        pv.AddSegment([[0,0],[10,0]])
        pv.AddSegment([[0,10],[10,10]])
        pv.AddSegment([[10,0],[10,10]])
        pv.Construct()
        cells = pv.GetCells()
        for i in range(len(cells)):
            cell = cells[i]
            if not cell.is_degenerate:
                valid_vertices = [v for v in cell.vertices if v != -1]
                self.assertTrue(len(valid_vertices) > 0)

    def test_input_point(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.Construct()
        self.assertTrue(1 == len(pv.GetPoints()))

    def test_discretize(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.AddSegment([[0,0],[0,10]])
        pv.AddSegment([[0,0],[10,0]])
        pv.AddSegment([[0,10],[10,10]])
        pv.AddSegment([[10,0],[10,10]])
        pv.Construct()

        vertices = pv.GetVertices()
        edges = pv.GetEdges()
        cells = pv.GetCells()

        testEdgeIndex = -1
        for i in range(len(edges)):
            if cells[edges[edges[i].twin].cell].source_category == 0 and cells[edges[i].cell].site == 1:
                testEdgeIndex = i

        #Test edge index. This should be the edge going from  {2.92893218813452, 2.92893218813452} To {2.92893218813452, 7.07106781186548}
        testEdge = edges[testEdgeIndex]
        sites = pv.ReturnCurvedSiteInformation(testEdge)
        self.assertTrue(sites[0] == [5,5])
        self.assertTrue(sites[1] == [[0,0],[0,10]])
        startVertex = vertices[testEdge.start]
        endVertex = vertices[testEdge.end]
        points = list(pv.DiscretizeCurvedEdge(testEdgeIndex, 3))


        #Validate the  start point
        self.assertAlmostEqual(points[0][0],startVertex.X)
        self.assertAlmostEqual(points[0][1], startVertex.Y)
        self.assertAlmostEqual(points[1][0], 2.5)
        self.assertAlmostEqual(points[1][1], 5)
        self.assertAlmostEqual(points[-1][0], endVertex.X)
        self.assertAlmostEqual(points[-1][1], endVertex.Y)

    def test_discretize_with_invalid_distance(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5, 5])
        pv.AddSegment([[0, 0], [0, 10]])
        pv.AddSegment([[0, 0], [10, 0]])
        pv.AddSegment([[0, 10], [10, 10]])
        pv.AddSegment([[10, 0], [10, 10]])
        pv.Construct()

        with self.assertRaises(ValueError):
            pv.DiscretizeCurvedEdge(2, -1)

    def test_get_raises_indexerror(self):
        pv = pyvoronoi.Pyvoronoi(1)
        with self.assertRaises(IndexError):
            pv.GetEdge(0) # shouldn't crash
        with self.assertRaises(IndexError):
            pv.GetVertex(0) # shouldn't crash
        with self.assertRaises(IndexError):
            pv.GetCell(0) # shouldn't crash
            
    def test_objects_dont_share_data(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5, 5])
        pv.AddSegment([[0, 0], [0, 10]])

        pv2 = pyvoronoi.Pyvoronoi(1)
        pv2.AddPoint([9, 9])
        pv2.AddSegment([[1, 1], [1, 9]])

        pv.Construct()
        pv2.Construct()

        self.assertEqual([[5, 5]], pv.GetPoints())
        self.assertEqual([[[0, 0], [0, 10]]], pv.GetSegments())
        self.assertEqual([[9, 9]], pv2.GetPoints())
        self.assertEqual([[[1, 1], [1, 9]]], pv2.GetSegments())


class TestInputSegmentIntersects(TestCase):
    def test_true_intersection_1(self):
        """
         Test that our 2 intersecting segments are detected among two other segment that do not intersect anything.
         :return:
        """

        pv = pyvoronoi.Pyvoronoi(1)

        # Those first two segments do not intersect
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [10, 10]])

        # Those two segments intersect but do not intersect the first two segments
        pv.AddSegment([[0, 0], [10, 0]])
        pv.AddSegment([[5, -5], [5, 10]])

        # The output should be our two intersecting segments indexed, without duplicating any idenfiers
        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(2, len(intersecting_segments))
        self.assertEqual(2, intersecting_segments[0])
        self.assertEqual(3, intersecting_segments[1])

    def test_true_intersection_2(self):
        """
         Test that our 2 intersecting segments are detected
         :return:
        """

        pv = pyvoronoi.Pyvoronoi(1)

        pv.AddSegment([[10, 10], [-10, -10]])
        pv.AddSegment([[-10, 10], [10, -10]])

        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(2, len(intersecting_segments))
        self.assertEqual(0, intersecting_segments[0])
        self.assertEqual(1, intersecting_segments[1])

    def test_end_intersection_3(self):
        """
         Test that our 2 intersecting segments are detected
         :return:
        """

        pv = pyvoronoi.Pyvoronoi(1)

        pv.AddSegment([[0, 0], [10, 0]])
        pv.AddSegment([[5, 0], [5, 10]])

        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(2, len(intersecting_segments))


    def test_true_intersection_at_ends_is_disregarded_horizontal(self):
        """
        Test that the intersection test returns false when line intersects at endpoints.
        In that case, they touch, but do not intersect.
        :return:
        """
        pv = pyvoronoi.Pyvoronoi(1)

        # Those first two segments not intersect
        pv.AddSegment([[0, 0], [10, 0]])
        pv.AddSegment([[-10, 0], [0, 0]])

        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(0, len(intersecting_segments))

    def test_true_intersection_at_ends_is_disregarded_vertical(self):
        """
        Test that the intersection test returns false when line intersects at endpoints.
        In that case, they touch, but do not intersect.
        :return:
        """
        pv = pyvoronoi.Pyvoronoi(1)

        # Those first two segments not intersect
        pv.AddSegment([[0, 0], [0, 10]])
        pv.AddSegment([[0, -10], [0, 0]])

        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(0, len(intersecting_segments))


    def test_true_intersection_at_ends_is_disregarded(self):
        """
        Test that the intersection test returns false when line intersects at endpoints.
        In that case, they touch, but do not intersect.
        :return:
        """
        pv = pyvoronoi.Pyvoronoi(1)

        # Those first two segments not intersect
        pv.AddSegment([[0, 0], [5, 5]])
        pv.AddSegment([[5, 5], [10, 10]])

        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(0, len(intersecting_segments))


    def test_colinearity_intersects_horizontal(self):
        pv = pyvoronoi.Pyvoronoi(1)

        # Those two segments overlap on 0,0 --> 5,0
        pv.AddSegment([[0, 0], [10, 0]])
        pv.AddSegment([[-10, 0], [5, 0]])

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [10, 10]])

        # Validate that the two segments that intersects on 0 --> 5 intersect
        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(2, len(intersecting_segments))
        self.assertEqual(0, intersecting_segments[0])
        self.assertEqual(1, intersecting_segments[1])

    def test_colinearity_intersects_vertical(self):
        pv = pyvoronoi.Pyvoronoi(1)

        # Those two segments overlap on 0,0 --> 5,0
        pv.AddSegment([[0, 0], [0, 10]])
        pv.AddSegment([[0, -10], [0, 5]])

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [10, 10]])

        # Validate that the two segments that intersects on 0 --> 5 intersect
        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(2, len(intersecting_segments))
        self.assertEqual(0, intersecting_segments[0])
        self.assertEqual(1, intersecting_segments[1])


    def test_colinearity_intersects(self):
        pv = pyvoronoi.Pyvoronoi(1)

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [10, 10]])

        # Those two segments overlap
        pv.AddSegment([[0, 0], [2, 2]])
        pv.AddSegment([[1, 1], [3, 3]])


        # Validate that the two segments that intersects on 0 --> 5 intersect
        intersecting_segments = pv.GetIntersectingSegments()
        self.assertEqual(2, len(intersecting_segments))
        self.assertEqual(2, intersecting_segments[0])
        self.assertEqual(3, intersecting_segments[1])

class TestDegeneratedInputSegment(TestCase):
    def test_no_degenerate_segment(self):
        pv = pyvoronoi.Pyvoronoi(1)

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [10, 10]])

        invalid_segments = pv.GetDegenerateSegments()
        self.assertEqual(0, len(invalid_segments))

    def test_degenerate_segment(self):
        pv = pyvoronoi.Pyvoronoi(1)

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [6, 6]])

        invalid_segments = pv.GetDegenerateSegments()
        # s = pv.GetSegment(invalid_segments[0])
        self.assertEqual(1, len(invalid_segments))

class TestInputPointOnInputSegment(TestCase):
    def test_no_point_on_segment(self):
        pv = pyvoronoi.Pyvoronoi(1)

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [10, 10]])
        pv.AddPoint([0,0])

        invalid_points = pv.GetPointsOnSegments()
        self.assertEqual(0, len(invalid_points))

    def test_point_on_segment_end_point(self):
        """
        Pyvoronoi does not consider the point on the line if it is equal to an end point.
        :return:
        """
        pv = pyvoronoi.Pyvoronoi(1)

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6, 6], [10, 10]])
        pv.AddPoint([10, 10])

        invalid_points = pv.GetPointsOnSegments()
        self.assertEqual(0, len(invalid_points))



    def test_point_on_segment_factor10(self):
        pv = pyvoronoi.Pyvoronoi(10)

        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-6, -6], [-10, -10]])
        pv.AddSegment([[6.6, 6.6], [10.1, 10.1]])
        pv.AddPoint([0, 0])
        pv.AddPoint([7.7, 7.7])

        invalid_points = pv.GetPointsOnSegments()
        self.assertEqual(1, len(invalid_points))
        self.assertEqual(1, invalid_points[0])

    def test_scenario(self):
        pv = pyvoronoi.Pyvoronoi(1)
        # Those two segments not intersect or overlap anything
        pv.AddSegment([[-8433001, 5672399], [-8418599, 5672399]])
        pv.AddSegment([[-8433001, 5672399], [-8433001, 5687401]])
        # pv.AddSegment([[-8418599, 5687401], [-8418599, 5672399]])



        invalid_segments = pv.GetIntersectingSegments()
        self.assertEqual(0, len(invalid_segments))

def run_tests():
    main()

if __name__ == '__main__':
    run_tests()
