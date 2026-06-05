# compute-geometry

A research-focused computational geometry library for Python, providing a set of
composable primitives and algorithms for solving geometric problems.

`compute-geometry` ships geometric **primitives** (points, segments, lines, circles
and polygons) validated with Pydantic, a collection of classical **algorithms**
(convex hull, Voronoi diagrams, Delaunay and polygon triangulation, minimum
enclosing circle and segment intersection), and a polished **visualization** layer
built on Matplotlib.

```python
import numpy as np
from cgeom.algorithms import VoronoiDiagram
from cgeom.visualization import plot_voronoi

points = np.loadtxt("examples/points1.txt")

voronoi = VoronoiDiagram(points)
cells = voronoi.build_voronoi_diagram()
plot_voronoi(voronoi, cells)
```

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
quickstart
```

```{toctree}
:maxdepth: 2
:caption: User guide

guide/elements
guide/algorithms
guide/visualization
```

