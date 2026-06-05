# Quickstart

This page walks through the core workflow: build a primitive, run an algorithm,
and visualize the result. Every snippet is self-contained.

## Import

The package is imported as `cgeom`. Primitives live at the top level, while
algorithms and plotting helpers live in submodules:

```python
import cgeom
from cgeom import Point, Segment, Circle, Polygon, Line
from cgeom.algorithms import (
    ConvexHull,
    DelaunayTriangulation,
    VoronoiDiagram,
    MinimumCircle,
    SegmentIntersection,
    PolygonTriangulation,
)
from cgeom.visualization import plot_convex_hull, plot_voronoi
```

## Primitives

Primitives are immutable Pydantic models and accept several construction forms:

```python
p1 = Point(x=1.0, y=2.0)      # keyword
p2 = Point([3, 4])            # from a list
p3 = Point((5, 6))            # from a tuple

p1.distance_to(p2)           # Euclidean distance

seg = Segment([[0, 0], [3, 4]])
seg.length                   # 5.0
seg.midpoint                 # Point(x=1.5, y=2.0)

circle = Circle(center=Point(0, 0), radius=5.0)
circle.area
circle.contains_point(Point(3, 4))

triangle = Polygon([[0, 0], [4, 0], [0, 3]])
triangle.area
triangle.perimeter
```

See the [Elements guide](guide/elements) for the full primitive API.

## Convex hull

```python
points = [(326, 237), (373, 209), (378, 265), (443, 241),
          (396, 231), (416, 270), (361, 335), (324, 297)]

hull = ConvexHull(points)
plot_convex_hull(hull)

print("Hull vertex indices:", hull.get_indexes())
```

## Voronoi diagram

```python
import numpy as np

points = np.loadtxt("examples/points1.txt")

voronoi = VoronoiDiagram(points)
cells = voronoi.build_voronoi_diagram()
plot_voronoi(voronoi, cells)
```

## Delaunay triangulation

```python
points = [
    (326, 237), (373, 209), (378, 265), (443, 241), (396, 231),
    (416, 270), (361, 335), (324, 297), (400, 306), (454, 315),
]

dt = DelaunayTriangulation(points)
triangles = dt.triangulate()
edges = dt.get_edges()

from cgeom.visualization import plot_delaunay
plot_delaunay(dt, title="Delaunay with Circumcircles", show_circumcircles=True)
```

## Segment intersection

```python
segments = [
    [[0, 0], [4, 4]],   # diagonal /
    [[0, 4], [4, 0]],   # diagonal \
    [[0, 2], [4, 2]],   # horizontal
    [[2, 0], [2, 4]],   # vertical
]

si = SegmentIntersection(segments)
intersections = si.find_intersections()      # Bentley–Ottmann sweep line
pairs = si.get_intersection_pairs()

from cgeom.visualization import plot_intersections
plot_intersections(si)
```

Continue to the [User guide](guide/algorithms) for the algorithm-by-algorithm
walkthrough, or jump to the [API reference](api/index).
