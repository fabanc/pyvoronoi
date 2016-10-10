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
        pv = pyvoronoi.Pyvoronoi(10)
        pv.AddPoint([0.5, 1])
        points = pv.GetPoints()
        print(points)
        self.assertTrue(len(points) == 1)
        self.assertTrue(points[0][0] == 0.5)
        self.assertTrue(points[0][1] == 1)

    def test_add_segment(self):
        pv = pyvoronoi.Pyvoronoi(10)
        segment = [[0.5, 1], [0, 2]]
        pv.AddSegment(segment)
        segments = pv.GetSegments()
        print(segments)
        self.assertTrue(len(segments) == 1)
        self.assertTrue(segments[0] == segment)

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
        self.assertTrue(len(filter(lambda e: edges[e].is_primary, cells[1].edges)) == 3)
        self.assertTrue(len(filter(lambda e: edges[e].is_primary, cells[3].edges)) == 2)

    def test_twins(self):
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

    def test_input_point(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddPoint([5,5])
        pv.Construct()
        self.assertTrue(1 == len(pv.inputPoints))

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
        points = pv.DiscretizeCurvedEdge(testEdgeIndex, 3)


        #Validate the  start point
        self.assertAlmostEquals(points[0][0],startVertex.X)
        self.assertAlmostEquals(points[0][1], startVertex.Y)
        self.assertAlmostEquals(points[1][0], 2.5)
        self.assertAlmostEquals(points[1][1], 5)
        self.assertAlmostEquals(points[-1][0], endVertex.X)
        self.assertAlmostEquals(points[-1][1], endVertex.Y)




def run_tests():
    main()

if __name__ == '__main__':
    run_tests()