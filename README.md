<h1 align="center">
<img src="logo.png" width="400">
</h1>

Welcome to the **scikit-geometry** library, a comprehensive computational geometry library for Python. This library is designed to provide a set of tools and algorithms for solving geometric problems, making it a valuable resource for developers and researchers working in areas such as computer graphics, robotics, and geographic information systems.

## Features

- **Geometry Primitives:** Efficient implementations of fundamental geometric primitives, including points, lines, and polygons.
  
- **Algorithms:** A variety of computational geometry algorithms, such as convex hull computation, line segment intersection, and Voronoi diagram generation.

- **Geometry Operations:** Functions for performing common geometric operations, such as distance calculations, area computations, and geometric transformations.

- **Documentation:** Comprehensive documentation with usage examples to help you quickly integrate the library into your projects.

## Installation

You can install the Geometry library using `pip`:

```bash
pip install scikit-geometry
```

## Getting Started

To get started with the _scikit-geometry_ library, import it into your Python script or project:

```python
import scikit-geometry as geometry
```

Explore the library's documentation to understand the available functionality and how to use it effectively.

## Examples

Here are some examples to demonstrate how to use the Geometry library:

```python
# Example 1: Compute the convex hull of a set of points
points = [(0, 0), (1, 1), (2, 2), (0, 2)]
convex_hull = geometry.convex_hull(points)
print("Convex Hull:", convex_hull.points)

convex_hull.plot()

# Example 2: Check if two line segments intersect
segment1 = [(1, 1), (4, 4)]
segment2 = [(3, 1), (1, 3)]
print(geometry.segments_intersect(segment1, segment2))
```

## Contributing

Contributions to the _scikit-geometry_ library are welcome! If you have additional algorithms, improvements, or bug fixes, please submit a pull request. Be sure to follow the [contribution guidelines](CONTRIBUTING.md) outlined in the repository.

## License

This library is licensed under the [MIT License](LICENSE), allowing you to use, modify, and distribute it for both commercial and non-commercial purposes.

Start exploring the world of computational geometry with the _scikit-geometry_ library in Python!