
# pyvoronoi

A wrapper for Boost's Voronoi diagram library. The full documentation of the Boost Voronoi API is available [here](https://www.boost.org/doc/libs/1_75_0/libs/polygon/doc/voronoi_main.htm).

## Documentation

The documentation for Pyvoronoi is available [here](https://pyvoronoi.readthedocs.io/en/latest/#).

The documentation is built with Sphinx. See the `docs` folder and the file `requirements.txt` to check the requirements if you want to build it locally .

## Change log
 * 1.2.4: Removed C++ warning during compile time by addressing possible data loss when switching between integer types.
 * 1.2.3: (skipped)
 * 1.2.2: Added API documentation using Sphynx and Read The Doc theme.
 * 1.2.1: Introduction of type hints for all function and their documentation. This enables intellisense and code completion.
 * 1.2.0: Improved memory management when returning Boost Voronoi output.
 * 1.1.9: Improved memory management when reading Boost Voronoi input.
 * 1.1.8: (skipped)
 * 1.1.7: Added method to validate the input data served to Boost Voronoi.

## Install

The installation have been tested on Windows and Linux Ubuntu. If you notice any issue on Mac, reach out to us, we are interested in making sure it works for you.

Windows users will need Microsoft Visual C++ installed on their machine. You can find information about the version needed on [this](https://wiki.python.org/moin/WindowsCompilers) link. Python version from 3.5 to 3.12 rely on Visual C++ 14.x.

### Dependencies

Cython dependency is optional. Cpp sources generated with Cython are available in releases.

Note on using the setup.py:

setup.py operates in 2 modes that are based on the presence of the dev file in the root of the project.

* When dev is **present**, Cython will be used to compile the .pyx sources. This is the development mode (as you get it in the git repository).

* When dev is **absent**, C/C++ compiler will be used to compile the .cpp sources (that were prepared in in the development mode). This is the distribution mode (as you get it on PyPI).

This way the package can be used without or with an incompatible version of Cython.

The idea comes from Matt Shannon's bandmat library.

### From PyPI


Cython not required.

``pip install pyvoronoi``

### From source


Cython required.

Clone the repository:

``git clone https://github.com/fabanc/pyvoronoi.git``

Install:

``python setup.py install``

After every modification of .pyx files compile with Cython:

``python setup.py build_ext --inplace``

Note in order to build the wheels, you will need to also install ``wheel``

``pip install wheel``

### Using


Create a new instance, passing the scaling factor into the constructor:
```
import pyvoronoi
pv = pyvoronoi.Pyvoronoi(10)
```

Since the voronoi library uses integer representation for points, the scaling factor chosen must be high enough
to avoid roundoff error when converting from point coordinates to integers.

Add points and segments:

```python
pv.AddPoint([0, 0])
pv.AddSegment([[1,5],[2,2]])
```

Call ``Construct()`` and get the edges and vertices:


```python
pv.Construct()
edges = pv.GetEdges()
vertices = pv.GetVertices()
cells = pv.GetCells()
```

Note that vertices, edges, and cells, can be accessed individually. The methods above are just convenience wrappers around
the following functions:

* GetVertex

* GetEdge

* Get Cell

```python
def GetVertices(self):
    count = self.CountVertices()
    output = []
    for index in  range(count):
        output.append(self.GetVertex(index))
    return output
```

```python
def GetEdges(self):
    count = self.CountEdges()
    output = []
    for index in range(count):
        output.append(self.GetEdge(index))
    return output
```

```python
def GetCells(self):
    count = self.CountCells()
    output = []
    for index in range(count):
        output.append(self.GetCell(index))
    return output
```

Note code above duplicates output data, which might be a problem for big datasets. As of 1.2.0, there is a more memory-friendly way to do that.

```python
    def GetVertices(self):
        count = self.CountVertices()
        output = []
        for index in  range(count):
            output.append(self.GetVertex(index))
        return output

    def EnumerateVertices(self):
        for index in range(self.CountVertices()):
            yield index, self.GetVertex(index)
            
    def EnumerateEdges(self):
        for index in range(self.CountEdges()):
            yield index, self.GetEdge(index)

    def EnumerateCells(self):
        for index in range(self.CountCells()):
            yield index, self.GetCell(index)
```

They return the same information as the list, but fetch each output object from the Boost C++ API. This is much more memory-friendly.

```python
   import pyvoronoi
   pv = pyvoronoi.Pyvoronoi(1)
   pv.AddPoint([5,5])
   pv.AddSegment([[0,0],[0,10]])
   pv.AddSegment([[0,0],[10,0]])
   pv.AddSegment([[0,10],[10,10]])
   pv.AddSegment([[10,0],[10,10]])
   pv.Construct()
   
   for index, vertex in pv.EnumerateVertices():
      # ... do things with vertex.X or vertex.Y
   
   for index, edge in pv.EnumerateEdges():
      # ... do things with the current edge
   
   for index, cell in pv.EnumerateCells():
      # ... do things with the current cell
```

Vertices have the following properties:

* ``X``: the position on the X-axis of the vertex.
* ``Y``: the position on the Y-axis of the vertex.

Edges have the following properties:

* ``start, end`` contain the indices of the start and end vertices or ``-1`` if the edge is infinite at that end.
* ``is_primary`` is true if the edge is not coincident with any of the source inputs.
* ``is_linear`` is true if the edge is linear (not curved).
* ``cell`` is the identifier of the cell this segment is part of.
* ``twin`` is the identifier of the twin segment as defined in the boost voronoi API.

Cells have the following properties:

* ``cell_identifier`` is the index of the cell.
* ``site`` is the index of the site which generated this cell (same as site1, site2 on the edges).
* ``contains_point`` is true if the site was generated by a point.
* ``contains_segment`` is true if the site was generated by a segment.
* ``is_open`` is true if any of the cell's edges is infinite.
* ``is_degenerate`` is true if the cell doesn't have an incident edge. Can happen if a few input segments share a common endpoint.
* ``vertices`` contains indices into the vertex array.
* ``edges`` contains indices into the edge array.

They have also a few instance methods. All those are instance methods of the class `Pyvoronoi`. Those methods takes a cell object as a parameter:
* ``RetrieveScaledPoint`` retrives information about the input point at the origin of a Voronoi Cell, when the center of the cell is a point. This method removes the scaling factor.
* ``RetrieveScaledSegment`` retrives information about the input segment at the origin of a Voronoi Cell, when the center of the cell is a segment.This method removes the scaling factor.
* ``RetrievePoint`` retrives information about the input point at the origin of a Voronoi Cell, when the center of the cell is a point. This method uses the scaling factor and show the coordinates as used by the voronoi builder in Boost.
* ``RetrieveSegment`` retrives information about the input segment at the origin of a Voronoi Cell, when the center of the cell is a segment.This method uses the scaling factor and show the coordinates as used by the voronoi builder in Boost.


```python
pv = pyvoronoi.Pyvoronoi(100)
pv.AddSegment([[0.1,0.8],[0.3,0.6]])
pv.AddSegment([[0.3,0.6],[0.4,0.6]])
pv.AddSegment([[0.4,0.6],[0.4,0.5]])
pv.AddSegment([[0.4,0.6],[0.4,0.7]])
pv.AddSegment([[0.4,0.7],[0.5,0.8]])
pv.AddSegment([[0.4,0.7],[0.5,0.6]])
pv.AddSegment([[0.5,0.6],[0.7,0.7]])

pv.Construct()
edges = pv.GetEdges()
vertices = pv.GetVertices()
cells = pv.GetCells()
print("Cell Count: {0}".format(len(cells)))
for c in cells:
    print("Cell contains point: {0}. Contains segment: {1}. Is open: {2}, Site Index: {3}".format(c.contains_point, c.contains_segment, c.is_open, c.site))
    print(",".join(map(str,c.vertices)))
    for sIndex in c.edges:
        print("Start Index: {0}, End Index = {1}".format(edges[sIndex].start, edges[sIndex].end))
```

Some output edges returned by the boost voronoi API are suposed to be curved. In the C++ API, it is up to you to code it. Luckily, you can do it in python using the following the function DiscretizeCurvedEdge.
The sample below shows you how to do that:

```python
for cIndex in range(len(cells)):
    cell = cells[cIndex]
    if cell.is_open == False:
        for i in range(len(cell.edges)):
            e = edges[cell.edges[i]]
            startVertex = vertices[e.start]
            endVertex = vertices[e.end]

            max_distance  = distance([startVertex.X, startVertex.Y], [endVertex.X, endVertex.Y]) / 10
            if startVertex != -1 and endVertex != -1:
                if(e.is_linear == True):
                    array = [[startVertex.X, startVertex.Y],[endVertex.X, endVertex.Y]]
                else:
                    points = pv.DiscretizeCurvedEdge(i, max_distance)
                    for p in points:
                        print "{0},{1}".format(p[0], p[1])
```

The curve interpolation code can return 2 exceptions.

* FocusOnDirectixException: this happens when the input point is on the segment side. In that cases, it makes no sense to interpolate a parabola between those two geometries since a parabola equation is supposed to find an equidistant point between the two geometries.

* UnsolvableParabolaEquation: there are cases where the point returned by boost does not fit with the parabola equation (for a same position on the x-axis, we get 2 different points, both equidistant). Understanding this issue is still under investigation. It is possible to mitigate this issue by setting an optional 3rd parameter of the function DiscretizeCurvedEdge). A higher value means more tolerance to this exception. The recommended value would be 1 / Scaling Factor.

### Data validation

According to the Boost Voronoi Documentation [here](https://www.boost.org/doc/libs/1_84_0/libs/polygon/doc/voronoi_main.htm)

> Input points and segments should not overlap except their endpoints. This means that input point should not lie inside the input segment and input segments should not intersect except their endpoints.

As of version 1.1.7 Pyvoronoi gives you 3 method to validate your input points and segments.

* GetPointsOnSegments: this function returns the list of indexes of all the input points located anywhere on a segment. Segments end points are disregarded.
* GetDegenerateSegments: this function returns the list of indexes of all degenerate segments. Degenerate segments use the same coordinates for their first and last point.
* GetIntersectingSegments: this function returns the list of indexes of all the segments that intersect another segment. Intersections between segments at endpoints only are disregarded.

Those function are can be handy if you are using a factor greater than 1 since the code validates the data after the factor has been applied. In other words, the coordinates tested are the coordinates used to solve the Voronoi problem. 

#### Example 1

```python
     pv = pyvoronoi.Pyvoronoi(1)

     # Those two segments do not intersect or overlap anything
     pv.AddSegment([[-6, -6], [-10, -10]])
     pv.AddSegment([[6, 6], [10, 10]])
     
     # The second point is located on the second segment
     pv.AddPoint([0, 0])
     pv.AddPoint([7, 7])
        
     # Will return [1] as the second point is on the second segment
     invalid_points = pv.GetPointsOnSegments()
```

#### Example 2

```python
     pv = pyvoronoi.Pyvoronoi(1)

     # Those two segments overlap on 0,0 --> 5,0
     pv.AddSegment([[0, 0], [10, 0]])
     pv.AddSegment([[-10, 0], [5, 0]])

     # Those two segments not intersect or overlap anything
     pv.AddSegment([[-6, -6], [-10, -10]])
     pv.AddSegment([[6, 6], [10, 10]])

    # Will return [0, 1] since the first two segments overlap
     intersecting_segments = pv.GetIntersectingSegments()
```

### Retrieving input geometries

Along with the validation data, you can inspect the data passed to pyvoronoi using a few convenience methods. Note that the coordinates returned are the coordinates after pyvoronoi applies the factor.
The coordinates you see are the coordinates used solve the Voronoi problem. 

As of version 1.1.9, you can access the input geometries used by pyvoronoi using

* GetPoint(index): returns the coordinates of the input point as pair of coordinate [x, y]
* GetSegment(index): returns the coordinates of the input segment pair of pair of coordinate [[x1, y1], [x2, y2]]
* CountPoints(): returns the number of input points passed to pyvoronoi
* CountSegments(): returns the number of input segments passed to pyvoronoi

Pyvoronoi also provides two generator to iterate through input points and segments:
* GetPoints()
* GetSegments()

A good example on how to use this code can be found in this unit test:

```python
    def test_retrieve_input(self):
        pv = pyvoronoi.Pyvoronoi(1)

        p1 = [5, 5]
        p2 = [0, 0]
        s1 = [[10, 10], [20, 20]]
        s2 = [[10, 10], [20, 20]]
        s3 = [[-10, -10], [-20, -20]]

        pv.AddPoint(p1)
        pv.AddPoint(p2)
        pv.AddSegment(s1)
        pv.AddSegment(s2)
        pv.AddSegment(s3)


        pv.Construct()

        self.assertTrue(2 == len(list(pv.GetPoints())))
        self.assertEqual(2, pv.CountPoints())
        self.assertTrue(3 == len(list(pv.GetSegments())))
        self.assertEqual(3, pv.CountSegments())
        self.assertEqual(p1, pv.GetPoint(0))
        self.assertEqual(p2, pv.GetPoint(1))
        self.assertEqual(s1, pv.GetSegment(0))
        self.assertEqual(s2, pv.GetSegment(1))
```

#### Example:

```python

```


# License

-  Pyvoronoi is available under `MIT
   license <http://opensource.org/licenses/MIT>`__.
-  The core Voronoi library is available under `Boost Software
   License <http://www.boost.org/LICENSE_1_0.txt>`__. Freeware for both
   open source and commercial applications.

# Development

## Build tools

This project uses [cibuildwheel](https://github.com/pypa/cibuildwheel) to build wheels on multiple platforms.

### Stubfile generation

I used CythonPEG.
https://github.com/RaubCamaioni/CythonPEG


### Documentation

I use Sphinx. No particular reason except I like how it works and the output is user-friendly, especially with the RTD theme.

 * Documentation: https://www.sphinx-doc.org/en/master/tutorial/index.html
 * Read the doc theme: https://sphinx-rtd-tutorial.readthedocs.io/en/latest/index.html

To generate the documentation:

```commandline
python.exe setup.py build_ext --inplace
D:\arcgis-pro-envs\pyvoronoi\Scripts\sphinx-build -M html docs/source/ docs/build/ -E -a
```

$ (sudo) pip install sphinx
$ (sudo) pip install sphinx-rtd-theme
D:\arcgis-pro-envs\pyvoronoi\Scripts\sphinx-build -M html docs/source/ docs/build/ -E -a
