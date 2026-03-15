import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from cgeom.algorithms import ConvexHull


class TestConvexHullTriangle:
    """Three points already forming a convex hull."""

    def setup_method(self):
        self.points = [[0, 0], [4, 0], [2, 3]]
        self.ch = ConvexHull(self.points)

    def test_hull_contains_all_points(self):
        hull = self.ch.convex_hull()
        assert len(hull) == 3

    def test_find_low_right(self):
        low = self.ch.find_low_right()
        assert low[1] == 0  # lowest y

    def test_get_indexes_length(self):
        indexes = self.ch.get_indexes()
        assert len(indexes) == 3


class TestConvexHullSquare:
    """Four points forming a square."""

    def setup_method(self):
        self.points = [[0, 0], [4, 0], [4, 4], [0, 4]]
        self.ch = ConvexHull(self.points)

    def test_hull_has_four_vertices(self):
        hull = self.ch.convex_hull()
        assert len(hull) == 4

    def test_hull_vertices_are_original_points(self):
        hull = self.ch.convex_hull()
        original = {(p[0], p[1]) for p in self.points}
        hull_set = {(p[0], p[1]) for p in hull}
        assert hull_set == original


class TestConvexHullWithInteriorPoints:
    """Points with some strictly inside the hull."""

    def setup_method(self):
        self.points = [
            [0, 0], [10, 0], [10, 10], [0, 10],  # corners
            [5, 5], [3, 3], [7, 2],                # interior
        ]
        self.ch = ConvexHull(self.points)

    def test_hull_excludes_interior(self):
        hull = self.ch.convex_hull()
        assert len(hull) == 4

    def test_interior_points_not_on_hull(self):
        hull = self.ch.convex_hull()
        hull_set = {(p[0], p[1]) for p in hull}
        for interior in [[5, 5], [3, 3], [7, 2]]:
            assert (interior[0], interior[1]) not in hull_set

    def test_get_indexes_match_hull(self):
        hull = self.ch.convex_hull()
        indexes = self.ch.get_indexes()
        for i, idx in enumerate(indexes):
            assert list(self.ch.points[idx]) == hull[i]


class TestConvexHullLargerSet:
    """Larger random-ish point set with known hull."""

    def setup_method(self):
        # circle points form the hull, interior ones don't
        angles = np.linspace(0, 2 * np.pi, 8, endpoint=False)
        circle = [[float(10 * np.cos(a)), float(10 * np.sin(a))] for a in angles]
        interior = [[1, 1], [-1, 2], [0, 0], [3, -2]]
        self.points = circle + interior
        self.ch = ConvexHull(self.points)

    def test_hull_size(self):
        hull = self.ch.convex_hull()
        assert len(hull) == 8

    def test_find_low_right_is_minimum_y(self):
        low = self.ch.find_low_right()
        min_y = min(p[1] for p in self.points)
        assert low[1] == pytest.approx(min_y)
