# pyvoronoi
A wrapper for Boost's Voronoi diagram library

# Install
## Dependencies

Cython dependency is optional. Cpp sources generated with Cython are available in releases.

Note on using the setup.py:

setup.py operates in 2 modes that are based on the presence of the dev file in the root of the project.

* When dev is **present**, Cython will be used to compile the .pyx sources. This is the development mode (as you get it in the git repository).

* When dev is **absent**, C/C++ compiler will be used to compile the .cpp sources (that were prepared in in the development mode). This is the distribution mode (as you get it on PyPI).

This way the package can be used without or with an incompatible version of Cython.

The idea comes from Matt Shannon's bandmat library.

## From source

Cython required.

Clone the repository:

git clone git@github.com:Voxel8/pyvoronoi.git

Install:

```python setup.py install```

After every modification of .pyx files compile with Cython:

```python setup.py build_ext --inplace```

## Using

Create a new instance:
``` 
import pyvoronoi
pv = pyvoronoi.Pyvoronoi()
```

Add points and segments:

```
pv.AddPoint([0, 0])
pv.AddSegment([[1,5],[2,2]])
```

Call ```Construct()``` and get the edges:
``` 
pv.Construct()
edges = pv.GetEdges()
```

Edges have the following properties:

* ```has_start, has_end``` indicate whether the edge has a start and end point, correspondingly. 
Edges can be infinite or one-ended. In this case the corresponding ```start/end``` property is invalid.
* ```start, end``` contain the X, Y coordinates of the start and end points if ```has_start/has_end``` is set.
* ```is_primary``` is true if the edge is not coincident with one of the source inputs.
* ```site1, site2``` are the indices of the sites which generated this edge. Sites are indexed as points, then segments,
so if there are 5 points and 3 segments, a site index of 7 means the last segment, a site index of 2 means the third 
point.

## License
Copyright 2014 Voxel8, Inc.

This program can not be copied and/or distributed without the permission of Voxel8, Inc.
