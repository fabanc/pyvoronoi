==========
 pyvoronoi
==========

A wrapper for Boost's Voronoi diagram library

=======
Install
=======

Dependencies
============

Cython dependency is optional. Cpp sources generated with Cython are available in releases.

Note on using the setup.py:

setup.py operates in 2 modes that are based on the presence of the dev file in the root of the project.

* When dev is **present**, Cython will be used to compile the .pyx sources. This is the development mode (as you get it in the git repository).

* When dev is **absent**, C/C++ compiler will be used to compile the .cpp sources (that were prepared in in the development mode). This is the distribution mode (as you get it on PyPI).

This way the package can be used without or with an incompatible version of Cython.

The idea comes from Matt Shannon's bandmat library.

From source
===========

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


License
=======

-  Pyvoronoi is available under `MIT
   license <http://opensource.org/licenses/MIT>`__.
-  The core Voronoi library is available under `Boost Software
   License <http://www.boost.org/LICENSE_1_0.txt>`__. Freeware for both
   open source and commercial applications.