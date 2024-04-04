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
        print(points)
        self.assertTrue(len(points) == 1)
        self.assertTrue(points[0][0] == inputPoint[0] * factor)
        self.assertTrue(points[0][1] == inputPoint[1] * factor)

    def test_add_segment(self):
        factor = 10
        segment = [[0.5, 1], [0, 2]]
        pv = pyvoronoi.Pyvoronoi(factor)
        pv.AddSegment(segment)
        segments = pv.GetSegments()
        print(segments)
        self.assertTrue(len(segments) == 1)
        self.assertEqual(
            segments[0],
            [
                [segment[0][0] * factor, segment[0][1] * factor],
                [segment[1][0] * factor, segment[1][1] * factor],
            ]
        )

    def test_add_point_after_construct(self):
        pv = pyvoronoi.Pyvoronoi()
        pv.Construct()
        self.assertRaises(pyvoronoi.VoronoiException, pv.AddPoint, [0, 0])

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
        #self.assertTrue(len([i for i in edges if i.is_primary == True]) == 4)Should have been changes when fabanc added the cell concept
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
        #self.assertTrue(len([i for i in edges if i.is_primary == True]) == 5) Should have been changes when fabanc added the cell concept
        self.assertTrue(len([i for i in edges if i.is_primary == True]) == 10)
        self.assertTrue(len(vertices) == 6)
        self.assertTrue(len(list(filter(lambda e: edges[e].is_primary, cells[1].edges))) == 3)
        self.assertTrue(len(list(filter(lambda e: edges[e].is_primary, cells[3].edges))) == 2)

    def test_square(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddSegment([[0, 0], [0, 2]])
        pv.AddSegment([[0, 2], [2, 2]])
        pv.AddSegment([[2, 2], [2, 0]])
        pv.AddSegment([[2, 0], [0, 0]])
        pv.Construct()
        edges = pv.GetEdges()
        vertices = pv.GetVertices()
        cells = pv.GetCells()

        self.assertTrue(len(cells) == 8)
        self.assertEqual(len([i for i in edges if i.is_primary == True]),8)
        self.assertTrue(len(vertices) == 5)
        self.assertTrue(len(list(filter(lambda e: edges[e].is_primary, cells[1].edges))) == 2)


    def test_square_middle_point_first(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([1,1])
        pv.AddSegment([[0, 0], [0, 2]])
        pv.AddSegment([[0, 2], [2, 2]])
        pv.AddSegment([[2, 2], [2, 0]])
        pv.AddSegment([[2, 0], [0, 0]])
        pv.Construct()
        edges = pv.GetEdges()
        vertices = pv.GetVertices()
        cells = pv.GetCells()

        r = pv.RetrieveSegment(cells[1])
        self.assertEqual(r, [[0, 0], [0, 2]])

    def test_square_middle_point_last(self):
        pv = pyvoronoi.Pyvoronoi(1)

        pv.AddSegment([[0, 0], [0, 2]])
        pv.AddSegment([[0, 2], [2, 2]])
        pv.AddSegment([[2, 2], [2, 0]])
        pv.AddSegment([[2, 0], [0, 0]])
        pv.AddPoint([1, 1])
        pv.Construct()
        edges = pv.GetEdges()
        vertices = pv.GetVertices()
        cells = pv.GetCells()


        r = pv.RetrieveSegment(cells[1])
        self.assertEqual(r, [[0, 0], [0, 2]])


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

    def test_iterators(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.AddSegment([[0,0],[0,10]])
        pv.AddSegment([[0,0],[10,0]])
        pv.AddSegment([[0,10],[10,10]])
        pv.AddSegment([[10,0],[10,10]])
        pv.Construct()

        vertices = pv.GetVertices()
        vertex_generator = pv.IterateVertices()
        for v in vertices:
            v_generator = next(vertex_generator)
            self.assertEqual(v, v_generator)

        edges = pv.GetEdges()
        edge_generator = pv.IterateEdges()
        for e in edges:
            e_generator = next(edge_generator)
            self.assertEqual(e, e_generator)

        cells = pv.GetCells()
        cell_generator = pv.IterateCells()
        for c in cells:
            c_generator = next(cell_generator)
            self.assertEqual(c, c_generator)

    def test_cells_vertices_duplication(self):
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

def run_tests():
    main()

if __name__ == '__main__':
    run_tests()
