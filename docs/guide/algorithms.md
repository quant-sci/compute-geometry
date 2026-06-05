# Algorithms

The `cgeom.algorithms` module bundles classical computational-geometry
algorithms. Each is a class constructed from a set of points (or segments) and
validated through the matching Pydantic input model.

```python
from cgeom.algorithms import (
    ConvexHull,
    DelaunayTriangulation,
    VoronoiDiagram,
    MinimumCircle,
    SegmentIntersection,
    PolygonTriangulation,
)
```

## ConvexHull

Computes the convex hull of a point set with the Gift Wrapping (Jarvis march)
algorithm.

```python
points = [(326, 237), (373, 209), (378, 265), (443, 241),
          (396, 231), (416, 270), (361, 335), (324, 297)]

hull = ConvexHull(points)
hull.convex_hull()        # ordered hull vertices [[x, y], ...]
hull.get_indexes()        # indices of the input points on the hull
```

Construction requires at least 3 non-collinear points.

## DelaunayTriangulation

Builds the Delaunay triangulation of a point set.

```python
dt = DelaunayTriangulation(points)
dt.triangulate()          # list of index triples
dt.get_edges()            # unique edges
dt.get_circumcircles()    # circumcircle of each triangle
```

## VoronoiDiagram

Builds the Voronoi diagram — the dual of the Delaunay triangulation.

```python
import numpy as np

points = np.loadtxt("examples/points1.txt")
voronoi = VoronoiDiagram(points)
cells = voronoi.build_voronoi_diagram()
```

## MinimumCircle

Finds the smallest enclosing circle of a point set.

```python
mc = MinimumCircle()
mc.minimum_circle([(0, 0), (1, 0), (0, 1), (1, 1)])   # [[cx, cy], radius]
```

## SegmentIntersection

Reports all pairwise intersections of a set of segments. The primary method uses
the Bentley–Ottmann sweep-line algorithm; a brute-force method is provided for
verification.

```python
segments = [
    [[0, 0], [4, 4]],
    [[0, 4], [4, 0]],
    [[0, 2], [4, 2]],
    [[2, 0], [2, 4]],
]

si = SegmentIntersection(segments)
si.find_intersections()              # sweep line
si.find_intersections_brute_force()  # verification
si.get_intersection_pairs()          # (i, j, point) tuples
```

## PolygonTriangulation

Triangulates a simple polygon.

```python
pt = PolygonTriangulation([[0, 0], [4, 0], [4, 4], [0, 4]])
pt.triangulate()
```

## Input validation

Every algorithm validates its input through a Pydantic model in
`cgeom.elements.models` (for example, `ConvexHullInput`,
`DelaunayTriangulationInput`). Invalid input — wrong shape, non-numeric values,
too few points, or degenerate configurations — raises
`pydantic.ValidationError` at construction time.

See the [API reference](../api/algorithms) for full signatures.
