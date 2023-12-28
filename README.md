<h1 align="center">
<img src="logo.png" width="400">
</h1>

Welcome to the **compute-geometry** library, a comprehensive computational geometry library for Python. This library is designed to provide a set of tools and algorithms for solving geometric problems, making it a valuable resource for developers and researchers working in areas such as computer graphics, robotics, and geographic information systems.

## Features

- **Geometry Primitives:** Efficient implementations of fundamental geometric primitives, including points, lines, and polygons.
  
- **Algorithms:** A variety of computational geometry algorithms, such as convex hull computation, line segment intersection, and Voronoi diagram generation.

- **Geometry Operations:** Functions for performing common geometric operations, such as distance calculations, area computations, and geometric transformations.

- **Visualization:** Functions for visualizing geometric objects and algorithms using Matplotlib.

- **Documentation:** Comprehensive documentation with usage examples to help you quickly integrate the library into your projects.

## Installation

You can install the **compute-geometry** library using `pip`:

```bash
pip install compute-geometry
```

## Getting Started

To get started with the _compute-geometry_ library, import it into your Python script or project:

```python
import compute-geometry as cgeom
```

Explore the library's documentation to understand the available functionality and how to use it effectively.

## Examples

Here are some examples to demonstrate how to use the Geometry library:

```python
# Example 1: Compute the convex hull of a set of points

from cgeom.algorithms import ConvexHull

# create a list of points
points = [(326, 237),(373, 209), (378, 265), (443, 241), (396, 231), (416, 270), (361, 335), (324, 297)]

# create a convex hull object with the list of points
convex_hull = ConvexHull(points)

# plot the convex hull
convex_hull.plot()

# print the indexes of the points that form the convex hull
print('Convex Hull: ', convex_hull.get_indexes())


# Example 2: Compute Voronoi diagram of a set of points

from cgeom.voronoi import VoronoiDiagram

# load a set of points
points = np.loadtxt("./points1.txt")

# create a voronoi diagram object
voronoi = VoronoiDiagram(points)

# build the voronoi diagram
cells = voronoi.build_voronoi_diagram()

# plot the voronoi diagram
voronoi.plot_voronoi(cells)

```

## Contributing

Contributions to the _compute-geometry_ library are welcome! If you have additional algorithms, improvements, or bug fixes, please submit a pull request. Be sure to follow the [contribution guidelines](CONTRIBUTING.md) outlined in the repository.

## License

This library is licensed under the [MIT License](LICENSE), allowing you to use, modify, and distribute it for both commercial and non-commercial purposes.

Start exploring the world of computational geometry with the _compute-geometry_ library in Python!