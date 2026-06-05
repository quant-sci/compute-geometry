# Elements

Geometric primitives are the building blocks of the library. They are immutable
[Pydantic](https://docs.pydantic.dev) models, so they validate their inputs on
construction and are hashable and safe to share.

All primitives are re-exported from the top-level package:

```python
from cgeom import Point, Line, Segment, Circle, Polygon
```

```{image} ../../public/elements.png
:alt: The geometric primitives
:align: center
```

## Point

A 2D point with `x` and `y` coordinates. It accepts keyword arguments, a list, a
tuple, or a NumPy array, and supports destructuring and indexing.

```python
import numpy as np

p1 = Point(x=1.0, y=2.0)
p2 = Point([3, 4])
p3 = Point(np.array([7, 8]))

x, y = p1                 # destructuring
p1[0], p1[1]              # indexing

p1.distance_to(p2)        # Euclidean distance
p3.to_numpy()             # -> np.ndarray([7., 8.])
p3.to_list()              # -> [7.0, 8.0]
```

Coordinates must be finite — `NaN` or `inf` raises a validation error.

## Line

An infinite line defined by two points. Exposes `slope`, `y_intercept`, the
implicit `coefficients` `(a, b, c)`, and a `contains_point` test. Vertical lines
report `slope` and `y_intercept` as `None`.

```python
line = Line([[0, 0], [2, 4]])
line.slope                # 2.0
line.y_intercept          # 0.0
line.coefficients         # (a, b, c)
line.contains_point(Point(1, 2))   # True

vertical = Line([[3, 0], [3, 5]])
vertical.slope            # None
```

## Segment

A finite line segment between two endpoints. Exposes `length` and `midpoint`.

```python
seg = Segment([[0, 0], [3, 4]])
seg.length                # 5.0
seg.midpoint              # Point(x=1.5, y=2.0)
seg.to_list()             # [[0, 0], [3, 4]]
```

## Circle

A circle defined by a center `Point` and a `radius`. Exposes `area`,
`circumference`, and `contains_point`.

```python
c1 = Circle(center=Point(0, 0), radius=5.0)
c1.area                   # 78.54...
c1.circumference
c1.contains_point(Point(3, 4))    # True (on the boundary)

c2 = Circle([[1, 2], 3])  # from [[cx, cy], r]
```

## Polygon

A polygon defined by an ordered list of vertices. Exposes `area`, `perimeter`,
and `num_vertices`, and is iterable.

```python
triangle = Polygon([[0, 0], [4, 0], [0, 3]])
triangle.area             # 6.0
triangle.perimeter        # 12.0
triangle.num_vertices     # 3

for vertex in triangle:
    print(vertex)

triangle.to_numpy()
```

## Interoperability with algorithms

Algorithm outputs map cleanly back onto primitives:

```python
from cgeom.algorithms import MinimumCircle, ConvexHull

# MinimumCircle returns [[cx, cy], radius] -> Circle
mc = MinimumCircle()
circle = Circle(mc.minimum_circle([(0, 0), (1, 0), (0, 1), (1, 1)]))

# ConvexHull returns [[x, y], ...] -> Polygon
ch = ConvexHull([(326, 237), (373, 209), (378, 265), (443, 241)])
hull_polygon = Polygon(ch.convex_hull())
```

See the [API reference](../api/elements) for full signatures.
