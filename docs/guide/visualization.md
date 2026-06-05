# Visualization

The `cgeom.visualization` module provides Matplotlib-based plotting helpers, one
per algorithm. They share a polished, minimalist navy/blue palette and sensible
defaults, so a single call produces a publication-ready figure.

```python
from cgeom.visualization import (
    plot_convex_hull,
    plot_delaunay,
    plot_intersections,
    plot_min_circle,
    plot_min_circle_random,
    plot_triangulation,
    plot_voronoi,
)
```

## Plotting helpers

Each helper takes the corresponding algorithm object (already constructed) and
draws its result.

```python
from cgeom.algorithms import ConvexHull, DelaunayTriangulation, SegmentIntersection

# Convex hull
plot_convex_hull(ConvexHull(points))

# Delaunay triangulation (optionally with circumcircles)
dt = DelaunayTriangulation(points)
dt.triangulate()
plot_delaunay(dt, title="Delaunay with Circumcircles", show_circumcircles=True)

# Segment intersections
plot_intersections(SegmentIntersection(segments))
```

## Voronoi diagrams

`plot_voronoi` takes the diagram object together with the cells returned by
`build_voronoi_diagram`:

```python
voronoi = VoronoiDiagram(points)
cells = voronoi.build_voronoi_diagram()
plot_voronoi(voronoi, cells)
```

## Minimum enclosing circle

```python
plot_min_circle(min_circle)        # plot a computed minimum circle
plot_min_circle_random(n=20)       # generate random points and plot
```

See the [API reference](../api/visualization) for full signatures.
