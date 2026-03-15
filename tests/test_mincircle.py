import matplotlib
matplotlib.use("Agg")

import math
import random
import pytest
from cgeom.algorithms import MinimumCircle


def _all_enclosed(points, center, radius, tol=1e-6):
    """Check every point lies within (or on) the circle."""
    for p in points:
        dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
        if dist > radius + tol:
            return False
    return True


class TestMinimumCircleTwoPoints:
    """Two points — circle diameter equals the distance between them."""

    def setup_method(self):
        self.mc = MinimumCircle()
        self.points = [[0, 0], [4, 0]]

    def test_center(self):
        circle = self.mc.minimum_circle(list(self.points))
        assert circle[0][0] == pytest.approx(2.0, abs=1e-6)
        assert circle[0][1] == pytest.approx(0.0, abs=1e-6)

    def test_radius(self):
        circle = self.mc.minimum_circle(list(self.points))
        assert circle[1] == pytest.approx(2.0, abs=1e-6)


class TestMinimumCircleThreePoints:
    """Three points on a known circumscribed circle."""

    def setup_method(self):
        self.mc = MinimumCircle()
        # equilateral triangle centered at origin, circumradius = 2
        self.points = [
            [2.0, 0.0],
            [-1.0, math.sqrt(3)],
            [-1.0, -math.sqrt(3)],
        ]

    def test_all_points_enclosed(self):
        pts = [p[:] for p in self.points]
        circle = self.mc.minimum_circle(pts)
        assert _all_enclosed(self.points, circle[0], circle[1])

    def test_radius_is_circumradius(self):
        pts = [p[:] for p in self.points]
        circle = self.mc.minimum_circle(pts)
        assert circle[1] == pytest.approx(2.0, abs=1e-4)


class TestMinimumCircleWithInterior:
    """Points with some strictly inside — hull points determine the circle."""

    def setup_method(self):
        self.mc = MinimumCircle()
        self.points = [
            [0, 0], [10, 0], [10, 10], [0, 10],  # square
            [5, 5], [3, 7],                         # interior
        ]

    def test_all_enclosed(self):
        pts = [p[:] for p in self.points]
        circle = self.mc.minimum_circle(pts)
        assert _all_enclosed(self.points, circle[0], circle[1])

    def test_radius_upper_bound(self):
        """Radius should not exceed the diagonal / 2 of the bounding box."""
        pts = [p[:] for p in self.points]
        circle = self.mc.minimum_circle(pts)
        half_diag = math.sqrt(10**2 + 10**2) / 2
        assert circle[1] <= half_diag + 1e-4


class TestMinimumCircleHeuristic:
    """Heuristic should also enclose all points."""

    def setup_method(self):
        self.mc = MinimumCircle()
        random.seed(42)
        self.points = [[random.uniform(-10, 10), random.uniform(-10, 10)] for _ in range(20)]

    def test_heuristic_encloses_all(self):
        circle = MinimumCircle.minimum_circle_heuristic(self.points)
        assert _all_enclosed(self.points, circle[0], circle[1])

    def test_heuristic_positive_radius(self):
        circle = MinimumCircle.minimum_circle_heuristic(self.points)
        assert circle[1] > 0


class TestMinimumCircleCollinear:
    """Collinear points — circle should span the full segment."""

    def setup_method(self):
        self.mc = MinimumCircle()
        self.points = [[0, 0], [2, 0], [4, 0], [6, 0]]

    def test_encloses_all(self):
        pts = [p[:] for p in self.points]
        circle = self.mc.minimum_circle(pts)
        assert _all_enclosed(self.points, circle[0], circle[1])

    def test_radius_is_half_span(self):
        pts = [p[:] for p in self.points]
        circle = self.mc.minimum_circle(pts)
        assert circle[1] == pytest.approx(3.0, abs=1e-4)
