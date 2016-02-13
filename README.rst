==========
 pyvoronoi
==========

A wrapper for Boost's Voronoi diagram library

Install
=======

Dependencies
------------

Cython dependency is optional. Cpp sources generated with Cython are available in releases.

Note on using the setup.py:

setup.py operates in 2 modes that are based on the presence of the dev file in the root of the project.

* When dev is **present**, Cython will be used to compile the .pyx sources. This is the development mode (as you get it in the git repository).

* When dev is **absent**, C/C++ compiler will be used to compile the .cpp sources (that were prepared in in the development mode). This is the distribution mode (as you get it on PyPI).

This way the package can be used without or with an incompatible version of Cython.

The idea comes from Matt Shannon's bandmat library.

From PyPI
---------

Cython not required.

``pip install pyvoronoi``

From source
-----------

Cython required.

Clone the repository:

``git clone git@github.com:Voxel8/pyvoronoi.git``

Install:

``python setup.py install``

After every modification of .pyx files compile with Cython:

``python setup.py build_ext --inplace``

Using
=====

Create a new instance, passing the scaling factor into the constructor:
``` 
import pyvoronoi
pv = pyvoronoi.Pyvoronoi(10)
```

Since the voronoi library uses integer representation for points, the scaling factor chosen must be high enough
to avoid roundoff error when converting from point coordinates to integers.

Add points and segments:

.. code:: python

	pv.AddPoint([0, 0])
	pv.AddSegment([[1,5],[2,2]])

Call ``Construct()`` and get the edges and vertices:

.. code:: python

	pv.Construct()
	edges = pv.GetEdges()
	vertices = pv.GetVertices()

Edges have the following properties:

* ``start, end`` contain the indices of the start and end vertices or ``-1`` if the edge is infinite at that end.
* ``is_primary`` is true if the edge is not coincident with any of the source inputs.
* ``site1, site2`` are the indices of the sites which generated this edge. Sites are indexed as points, then segments, so if there are 5 points and 3 segments, a site index of 7 means the last segment, a site index of 2 means the third point.

The function construct does not return any cell. In order to get cells, use the following function:

.. code:: python
	pv = pyvoronoi.Pyvoronoi(100)
	pv.AddSegment([[0.1,0.8],[0.3,0.6]])
	pv.AddSegment([[0.3,0.6],[0.4,0.6]])
	pv.AddSegment([[0.4,0.6],[0.4,0.5]])
	pv.AddSegment([[0.4,0.6],[0.4,0.7]])
	pv.AddSegment([[0.4,0.7],[0.5,0.8]])
	pv.AddSegment([[0.4,0.7],[0.5,0.6]])
	pv.AddSegment([[0.5,0.6],[0.7,0.7]])

	for i in range(len(pfroms)):
		pv.AddSegment([pfroms[i],ptos[i]])

	pv.ConstructWithCells()
	edges = pv.GetCellEdges()
	vertices = pv.GetCellVertices()
		
	print "Cell Edges"
	cells = pv.GetCells()
	print "Cell Count: " + str(len(cells))
	for c in cells:
		print "Cell ID: {0}. Contains point: {1}. Contains segment: {2}. Is open: {3}, Site Index: {4}".format(c.cellId, c.contains_point, c.contains_segment, c.is_open, c.source_index)#Works fine
		print ",".join(map(str,c.vertices))
		for sIndex in c.segments:
			print "Start Index: {0}, End Index = {1}".format(edges[sIndex].start, edges[sIndex].end)#Fail with error AttributeError: 'dict' object has no attribute 'x1'
		print "\n"

Note that when using the method ConstructWithCells instead of Construct , the object are retrieved using different methods:

* GetCells() --> GetCellVertices()
* GetEdges() --> GetCellEdges()

You can also retrieve object that belong to the class VoronoiCell using the method GetCells()

.. code:: python
class VoronoiCell:
    cellId = -1
    source_index = -1
    contains_point = 0
    contains_segment = 0
    is_open = 0
	
    vertices = None
    segments = None	
		

License
=======

-  Pyvoronoi is available under `MIT
   license <http://opensource.org/licenses/MIT>`__.
-  The core Voronoi library is available under `Boost Software
   License <http://www.boost.org/LICENSE_1_0.txt>`__. Freeware for both
   open source and commercial applications.