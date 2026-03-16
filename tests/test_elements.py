"""Tests for geometric primitive classes."""

import math

import numpy as np
import pytest
from pydantic import ValidationError

from cgeom.elements.elements import Circle, Line, Point, Polygon, Segment


class TestPoint:
    def test_keyword_construction(self):
        p = Point(x=1.0, y=2.0)
        assert p.x == 1.0
        assert p.y == 2.0

    def test_from_list(self):
        p = Point([3, 4])
        assert p.x == 3.0
        assert p.y == 4.0

    def test_from_tuple(self):
        p = Point((5, 6))
        assert p.x == 5.0
        assert p.y == 6.0

    def test_from_ndarray(self):
        p = Point(np.array([7.0, 8.0]))
        assert p.x == 7.0
        assert p.y == 8.0

    def test_getitem(self):
        p = Point(x=1, y=2)
        assert p[0] == 1.0
        assert p[1] == 2.0
        with pytest.raises(IndexError):
            p[2]

    def test_iter(self):
        p = Point(x=3, y=4)
        x, y = p
        assert x == 3.0
        assert y == 4.0

    def test_len(self):
        assert len(Point(x=0, y=0)) == 2

    def test_distance_to(self):
        p1 = Point(x=0, y=0)
        p2 = Point(x=3, y=4)
        assert p1.distance_to(p2) == pytest.approx(5.0)

    def test_to_numpy(self):
        arr = Point(x=1, y=2).to_numpy()
        assert isinstance(arr, np.ndarray)
        np.testing.assert_array_equal(arr, [1.0, 2.0])

    def test_to_list(self):
        assert Point(x=1, y=2).to_list() == [1.0, 2.0]

    def test_immutability(self):
        p = Point(x=1, y=2)
        with pytest.raises(ValidationError):
            p.x = 5

    def test_hashable(self):
        p1 = Point(x=1, y=2)
        p2 = Point(x=1, y=2)
        assert hash(p1) == hash(p2)
        assert {p1, p2} == {p1}

    def test_rejects_nan(self):
        with pytest.raises(ValidationError):
            Point(x=float("nan"), y=1)

    def test_rejects_inf(self):
        with pytest.raises(ValidationError):
            Point(x=1, y=float("inf"))

    def test_repr(self):
        assert repr(Point(x=1, y=2)) == "Point(x=1.0, y=2.0)"


class TestLine:
    def test_keyword_construction(self):
        line = Line(point1=Point(x=0, y=0), point2=Point(x=1, y=1))
        assert line.point1 == Point(x=0, y=0)

    def test_from_raw_lists(self):
        line = Line([[0, 0], [1, 1]])
        assert line.point1 == Point(x=0, y=0)
        assert line.point2 == Point(x=1, y=1)

    def test_coefficients_normalized(self):
        line = Line([[0, 0], [1, 0]])  # horizontal
        a, b, c = line.coefficients
        assert a**2 + b**2 == pytest.approx(1.0)
        assert c == pytest.approx(0.0)

    def test_slope_normal(self):
        line = Line([[0, 0], [2, 4]])
        assert line.slope == pytest.approx(2.0)

    def test_slope_vertical(self):
        line = Line([[1, 0], [1, 5]])
        assert line.slope is None

    def test_y_intercept(self):
        line = Line([[0, 3], [1, 5]])
        assert line.y_intercept == pytest.approx(3.0)

    def test_y_intercept_vertical(self):
        line = Line([[1, 0], [1, 5]])
        assert line.y_intercept is None

    def test_contains_point(self):
        line = Line([[0, 0], [1, 1]])
        assert line.contains_point(Point(x=5, y=5)) is True
        assert line.contains_point(Point(x=1, y=0)) is False

    def test_rejects_identical_points(self):
        with pytest.raises(ValidationError):
            Line([[1, 1], [1, 1]])

    def test_repr(self):
        line = Line([[0, 0], [1, 1]])
        assert "Line" in repr(line)


class TestSegment:
    def test_keyword_construction(self):
        seg = Segment(start=Point(x=0, y=0), end=Point(x=3, y=4))
        assert seg.start == Point(x=0, y=0)

    def test_from_raw_lists(self):
        seg = Segment([[0, 0], [3, 4]])
        assert seg.end == Point(x=3, y=4)

    def test_length(self):
        seg = Segment([[0, 0], [3, 4]])
        assert seg.length == pytest.approx(5.0)

    def test_midpoint(self):
        seg = Segment([[0, 0], [4, 6]])
        mid = seg.midpoint
        assert mid.x == pytest.approx(2.0)
        assert mid.y == pytest.approx(3.0)

    def test_to_list(self):
        assert Segment([[1, 2], [3, 4]]).to_list() == [[1.0, 2.0], [3.0, 4.0]]

    def test_rejects_zero_length(self):
        with pytest.raises(ValidationError):
            Segment([[1, 1], [1, 1]])

    def test_immutability(self):
        seg = Segment([[0, 0], [1, 1]])
        with pytest.raises(ValidationError):
            seg.start = Point(x=5, y=5)


class TestCircle:
    def test_keyword_construction(self):
        c = Circle(center=Point(x=0, y=0), radius=5.0)
        assert c.center == Point(x=0, y=0)
        assert c.radius == 5.0

    def test_from_list_format(self):
        c = Circle([[1, 2], 3.0])
        assert c.center == Point(x=1, y=2)
        assert c.radius == 3.0

    def test_area(self):
        c = Circle(center=Point(x=0, y=0), radius=1.0)
        assert c.area == pytest.approx(math.pi)

    def test_circumference(self):
        c = Circle(center=Point(x=0, y=0), radius=1.0)
        assert c.circumference == pytest.approx(2 * math.pi)

    def test_contains_point_inside(self):
        c = Circle(center=Point(x=0, y=0), radius=5.0)
        assert c.contains_point(Point(x=1, y=1)) is True

    def test_contains_point_on_boundary(self):
        c = Circle(center=Point(x=0, y=0), radius=5.0)
        assert c.contains_point(Point(x=5, y=0)) is True

    def test_contains_point_outside(self):
        c = Circle(center=Point(x=0, y=0), radius=5.0)
        assert c.contains_point(Point(x=6, y=0)) is False

    def test_to_list(self):
        c = Circle([[1, 2], 3.0])
        assert c.to_list() == [[1.0, 2.0], 3.0]

    def test_rejects_non_positive_radius(self):
        with pytest.raises(ValidationError):
            Circle(center=Point(x=0, y=0), radius=0)
        with pytest.raises(ValidationError):
            Circle(center=Point(x=0, y=0), radius=-1)

    def test_rejects_inf_radius(self):
        with pytest.raises(ValidationError):
            Circle(center=Point(x=0, y=0), radius=float("inf"))

    def test_interop_with_minimum_circle_output(self):
        """Circle accepts the exact format MinimumCircle returns: [[cx,cy], r]."""
        raw = [[2.5, 3.0], 4.123]
        c = Circle(raw)
        assert c.center.x == pytest.approx(2.5)
        assert c.radius == pytest.approx(4.123)


class TestPolygon:
    def test_from_list_of_lists(self):
        poly = Polygon([[0, 0], [4, 0], [4, 3], [0, 3]])
        assert poly.num_vertices == 4

    def test_from_ndarray(self):
        arr = np.array([[0, 0], [1, 0], [0, 1]])
        poly = Polygon(arr)
        assert poly.num_vertices == 3

    def test_area_ccw(self):
        poly = Polygon([[0, 0], [4, 0], [4, 3], [0, 3]])
        assert poly.area == pytest.approx(12.0)

    def test_area_cw(self):
        poly = Polygon([[0, 0], [0, 3], [4, 3], [4, 0]])
        assert poly.area == pytest.approx(-12.0)

    def test_perimeter(self):
        poly = Polygon([[0, 0], [4, 0], [4, 3], [0, 3]])
        assert poly.perimeter == pytest.approx(14.0)

    def test_iteration(self):
        pts = [[0, 0], [1, 0], [0, 1]]
        poly = Polygon(pts)
        for i, v in enumerate(poly):
            assert v == Point(pts[i])

    def test_indexing(self):
        poly = Polygon([[0, 0], [1, 0], [0, 1]])
        assert poly[0] == Point(x=0, y=0)
        assert poly[2] == Point(x=0, y=1)

    def test_len(self):
        poly = Polygon([[0, 0], [1, 0], [0, 1]])
        assert len(poly) == 3

    def test_to_numpy(self):
        poly = Polygon([[0, 0], [1, 0], [0, 1]])
        arr = poly.to_numpy()
        assert arr.shape == (3, 2)
        np.testing.assert_array_equal(arr[0], [0.0, 0.0])

    def test_to_list(self):
        poly = Polygon([[1, 2], [3, 4], [5, 6]])
        assert poly.to_list() == [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]

    def test_rejects_fewer_than_3(self):
        with pytest.raises(ValidationError):
            Polygon([[0, 0], [1, 1]])

    def test_hashable(self):
        p1 = Polygon([[0, 0], [1, 0], [0, 1]])
        p2 = Polygon([[0, 0], [1, 0], [0, 1]])
        assert hash(p1) == hash(p2)

    def test_keyword_construction(self):
        poly = Polygon(vertices=[(0, 0), (1, 0), (0, 1)])
        assert poly.num_vertices == 3
