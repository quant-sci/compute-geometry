"""Immutable geometric primitive classes built on Pydantic v2."""

import math
from typing import Iterator

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator


class Point(BaseModel):
    """A 2D point with coordinates (x, y)."""

    model_config = ConfigDict(frozen=True)
    x: float
    y: float

    def __init__(self, *args, **kwargs):
        if len(args) == 2 and not kwargs:
            super().__init__(x=args[0], y=args[1])
        elif len(args) == 1 and not kwargs:
            obj = type(self).model_validate(args[0])
            super().__init__(**obj.model_dump())
        else:
            super().__init__(**kwargs)

    @model_validator(mode="before")
    @classmethod
    def _coerce_input(cls, data):
        if isinstance(data, dict):
            return data
        if isinstance(data, np.ndarray):
            data = data.tolist()
        if isinstance(data, (list, tuple)) and len(data) == 2:
            return {"x": data[0], "y": data[1]}
        return data

    @model_validator(mode="after")
    def _validate_finite(self):
        if not (math.isfinite(self.x) and math.isfinite(self.y)):
            raise ValueError("Point coordinates must be finite (no NaN or inf)")
        return self

    def distance_to(self, other: "Point") -> float:
        """Euclidean distance to another point."""
        return math.hypot(self.x - other.x, self.y - other.y)

    def to_numpy(self) -> np.ndarray:
        """Return as numpy array [x, y]."""
        return np.array([self.x, self.y])

    def to_list(self) -> list[float]:
        """Return as [x, y]."""
        return [self.x, self.y]

    def __getitem__(self, index: int) -> float:
        if index == 0:
            return self.x
        if index == 1:
            return self.y
        raise IndexError(f"Point index {index} out of range (0..1)")

    def __len__(self) -> int:
        return 2

    def __iter__(self) -> Iterator[float]:
        yield self.x
        yield self.y

    def __repr__(self) -> str:
        return f"Point(x={self.x}, y={self.y})"


class Line(BaseModel):
    """An infinite line defined by two distinct points."""

    model_config = ConfigDict(frozen=True)
    point1: Point
    point2: Point

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and not kwargs:
            obj = type(self).model_validate(args[0])
            super().__init__(**obj.model_dump())
        else:
            super().__init__(**kwargs)

    @model_validator(mode="before")
    @classmethod
    def _coerce_input(cls, data):
        if isinstance(data, dict):
            return data
        if isinstance(data, (list, tuple)) and len(data) == 2:
            return {"point1": data[0], "point2": data[1]}
        return data

    @model_validator(mode="after")
    def _validate_distinct(self):
        if self.point1.x == self.point2.x and self.point1.y == self.point2.y:
            raise ValueError("Line requires two distinct points")
        return self

    @property
    def coefficients(self) -> tuple[float, float, float]:
        """General form (a, b, c) normalized so a² + b² = 1."""
        a = self.point1.y - self.point2.y
        b = self.point2.x - self.point1.x
        c = self.point1.x * self.point2.y - self.point2.x * self.point1.y
        norm = math.hypot(a, b)
        return (a / norm, b / norm, c / norm)

    @property
    def slope(self) -> float | None:
        """Slope of the line, or None if vertical."""
        dx = self.point2.x - self.point1.x
        if dx == 0:
            return None
        return (self.point2.y - self.point1.y) / dx

    @property
    def y_intercept(self) -> float | None:
        """Y-intercept, or None if vertical."""
        s = self.slope
        if s is None:
            return None
        return self.point1.y - s * self.point1.x

    def contains_point(self, point: "Point", tol: float = 1e-10) -> bool:
        """Check whether a point lies on this line."""
        a, b, c = self.coefficients
        return abs(a * point.x + b * point.y + c) <= tol

    def __repr__(self) -> str:
        return f"Line(point1={self.point1!r}, point2={self.point2!r})"


class Segment(BaseModel):
    """A line segment between two distinct points."""

    model_config = ConfigDict(frozen=True)
    start: Point
    end: Point

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and not kwargs:
            obj = type(self).model_validate(args[0])
            super().__init__(**obj.model_dump())
        else:
            super().__init__(**kwargs)

    @model_validator(mode="before")
    @classmethod
    def _coerce_input(cls, data):
        if isinstance(data, dict):
            return data
        if isinstance(data, (list, tuple)) and len(data) == 2:
            return {"start": data[0], "end": data[1]}
        return data

    @model_validator(mode="after")
    def _validate_distinct(self):
        if self.start.x == self.end.x and self.start.y == self.end.y:
            raise ValueError("Segment requires two distinct endpoints")
        return self

    @property
    def length(self) -> float:
        """Euclidean length of the segment."""
        return self.start.distance_to(self.end)

    @property
    def midpoint(self) -> Point:
        """Midpoint of the segment."""
        return Point(x=(self.start.x + self.end.x) / 2,
                     y=(self.start.y + self.end.y) / 2)

    def to_list(self) -> list[list[float]]:
        """Return as [[x1, y1], [x2, y2]]."""
        return [self.start.to_list(), self.end.to_list()]

    def __repr__(self) -> str:
        return f"Segment(start={self.start!r}, end={self.end!r})"


class Circle(BaseModel):
    """A circle defined by a center point and radius."""

    model_config = ConfigDict(frozen=True)
    center: Point
    radius: float

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and not kwargs:
            obj = type(self).model_validate(args[0])
            super().__init__(**obj.model_dump())
        else:
            super().__init__(**kwargs)

    @model_validator(mode="before")
    @classmethod
    def _coerce_input(cls, data):
        if isinstance(data, dict):
            return data
        if isinstance(data, (list, tuple)) and len(data) == 2:
            first, second = data
            is_point_like = isinstance(first, (list, tuple, np.ndarray, Point))
            is_scalar = isinstance(second, (int, float, np.floating))
            if is_point_like and is_scalar:
                return {"center": first, "radius": second}
        return data

    @model_validator(mode="after")
    def _validate_radius(self):
        if not (math.isfinite(self.radius) and self.radius > 0):
            raise ValueError("Circle radius must be positive and finite")
        return self

    @property
    def area(self) -> float:
        """Area of the circle."""
        return math.pi * self.radius ** 2

    @property
    def circumference(self) -> float:
        """Circumference of the circle."""
        return 2 * math.pi * self.radius

    def contains_point(self, point: "Point", tol: float = 1e-10) -> bool:
        """Check whether a point lies inside or on the circle."""
        return self.center.distance_to(point) <= self.radius + tol

    def to_list(self) -> list:
        """Return as [[cx, cy], r]."""
        return [self.center.to_list(), self.radius]

    def __repr__(self) -> str:
        return f"Circle(center={self.center!r}, radius={self.radius})"


class Polygon(BaseModel):
    """A polygon defined by an ordered sequence of vertices."""

    model_config = ConfigDict(frozen=True)
    vertices: tuple[Point, ...]

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and not kwargs:
            obj = type(self).model_validate(args[0])
            super().__init__(**obj.model_dump())
        else:
            super().__init__(**kwargs)

    @model_validator(mode="before")
    @classmethod
    def _coerce_input(cls, data):
        if isinstance(data, dict):
            return data
        if isinstance(data, np.ndarray):
            return {"vertices": data.tolist()}
        if isinstance(data, (list, tuple)):
            return {"vertices": data}
        return data

    @model_validator(mode="after")
    def _validate_min_vertices(self):
        if len(self.vertices) < 3:
            raise ValueError("Polygon requires at least 3 vertices")
        return self

    @property
    def num_vertices(self) -> int:
        """Number of vertices."""
        return len(self.vertices)

    @property
    def area(self) -> float:
        """Signed area using the shoelace formula (positive for CCW)."""
        n = len(self.vertices)
        s = 0.0
        for i in range(n):
            j = (i + 1) % n
            s += self.vertices[i].x * self.vertices[j].y
            s -= self.vertices[j].x * self.vertices[i].y
        return s / 2.0

    @property
    def perimeter(self) -> float:
        """Perimeter of the polygon."""
        n = len(self.vertices)
        return sum(
            self.vertices[i].distance_to(self.vertices[(i + 1) % n])
            for i in range(n)
        )

    def to_numpy(self) -> np.ndarray:
        """Return as numpy array of shape (n, 2)."""
        return np.array([v.to_list() for v in self.vertices])

    def to_list(self) -> list[list[float]]:
        """Return as [[x, y], ...]."""
        return [v.to_list() for v in self.vertices]

    def __iter__(self) -> Iterator["Point"]:
        return iter(self.vertices)

    def __len__(self) -> int:
        return len(self.vertices)

    def __getitem__(self, index: int) -> "Point":
        return self.vertices[index]

    def __repr__(self) -> str:
        return f"Polygon(vertices={list(self.vertices)!r})"
