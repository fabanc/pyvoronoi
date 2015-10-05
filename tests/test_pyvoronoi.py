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
        self.assertTrue(len([i for i in edges if i.is_primary == True]) == 4)
        self.assertTrue(len(vertices) == 5)

    def test_rectangle(self):
        pv = pyvoronoi.Pyvoronoi(1)
        pv.AddSegment([[0, 0], [0, 2]])
        pv.AddSegment([[0, 2], [1, 2]])
        pv.AddSegment([[1, 2], [1, 0]])
        pv.AddSegment([[1, 0], [0, 0]])
        pv.Construct()
        edges = pv.GetEdges()
        vertices = pv.GetVertices()
        self.assertTrue(len([i for i in edges if i.is_primary == True]) == 5)
        self.assertTrue(len(vertices) == 6)

def run_tests():
    main()

if __name__ == '__main__':
    run_tests()