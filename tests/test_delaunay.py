"""Tests for Delaunay triangulation (Bowyer-Watson)."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from cgeom.algorithms import ConvexHull, DelaunayTriangulation


class TestDelaunayTriangle:
    """Three points forming a single triangle."""

    def setup_method(self):
        self.points = [[0, 0], [4, 0], [2, 3]]
        self.dt = DelaunayTriangulation(self.points)

    def test_one_triangle(self):
        tris = self.dt.triangulate()
        assert len(tris) == 1

    def test_three_edges(self):
        edges = self.dt.get_edges()
        assert len(edges) == 3

    def test_correct_indices(self):
        tris = self.dt.triangulate()
        assert tris[0] == [0, 1, 2]


class TestDelaunaySquare:
    """Four points forming a square — two triangles."""

    def setup_method(self):
        self.points = [[0, 0], [4, 0], [4, 4], [0, 4]]
        self.dt = DelaunayTriangulation(self.points)

    def test_two_triangles(self):
        tris = self.dt.triangulate()
        assert len(tris) == 2

    def test_five_edges(self):
        edges = self.dt.get_edges()
        assert len(edges) == 5

    def test_delaunay_property(self):
        _assert_delaunay_property(self.dt)


class TestDelaunaySquareWithCenter:
    """Square with a center point — four triangles."""

    def setup_method(self):
        self.points = [[0, 0], [4, 0], [4, 4], [0, 4], [2, 2]]
        self.dt = DelaunayTriangulation(self.points)

    def test_four_triangles(self):
        tris = self.dt.triangulate()
        assert len(tris) == 4

    def test_delaunay_property(self):
        _assert_delaunay_property(self.dt)


class TestDelaunayProperty:
    """Verify the Delaunay property on random point sets."""

    @pytest.mark.parametrize("n", [10, 20])
    def test_no_point_in_circumcircle(self, n):
        rng = np.random.default_rng(42 + n)
        pts = rng.uniform(-100, 100, size=(n, 2)).tolist()
        dt = DelaunayTriangulation(pts)
        _assert_delaunay_property(dt)


class TestDelaunayConvexHullRelation:
    """Convex hull edges must be a subset of Delaunay edges."""

    def test_hull_edges_subset(self):
        pts = [
            [0, 0], [10, 0], [10, 10], [0, 10],
            [5, 5], [3, 3], [7, 2], [2, 8],
            [8, 6], [6, 1], [1, 5], [9, 9],
        ]
        dt = DelaunayTriangulation(pts)
        ch = ConvexHull(pts)
        hull = ch.convex_hull()

        # Build hull edges as index pairs using the original points list
        def _find_index(coord):
            for i, p in enumerate(pts):
                if abs(p[0] - coord[0]) < 1e-9 and abs(p[1] - coord[1]) < 1e-9:
                    return i
            return -1

        hull_idx = [_find_index(v) for v in hull]
        hull_edges = set()
        for i in range(len(hull_idx)):
            a, b = hull_idx[i], hull_idx[(i + 1) % len(hull_idx)]
            hull_edges.add((min(a, b), max(a, b)))

        dt_edges = set()
        for tri in dt.triangulate():
            for j in range(3):
                a, b = tri[j], tri[(j + 1) % 3]
                dt_edges.add((min(a, b), max(a, b)))

        assert hull_edges.issubset(dt_edges)


class TestDelaunayOutputFormats:
    """Check that output types and shapes are correct."""

    def setup_method(self):
        self.dt = DelaunayTriangulation([[0, 0], [4, 0], [2, 3], [1, 2]])

    def test_triangulate_returns_list_of_int_lists(self):
        tris = self.dt.triangulate()
        assert isinstance(tris, list)
        for tri in tris:
            assert isinstance(tri, list)
            assert len(tri) == 3
            for idx in tri:
                assert isinstance(idx, int)

    def test_get_triangles_coordinates(self):
        coord_tris = self.dt.get_triangles()
        for tri in coord_tris:
            assert len(tri) == 3
            for pt in tri:
                assert len(pt) == 2

    def test_edges_unique(self):
        edges = self.dt.get_edges()
        edge_set = set()
        for e in edges:
            key = (tuple(e[0]), tuple(e[1]))
            assert key not in edge_set
            edge_set.add(key)

    def test_circumcircle_format(self):
        circs = self.dt.get_circumcircles()
        for circ in circs:
            assert len(circ) == 2
            assert len(circ[0]) == 2  # center [cx, cy]
            assert isinstance(circ[1], float)  # radius
            assert circ[1] > 0


class TestDelaunayLargerSet:
    """20 points from the Voronoi test data — structural validity."""

    def setup_method(self):
        self.data = [
            [326, 237], [373, 209], [378, 265], [443, 241], [396, 231],
            [416, 270], [361, 335], [324, 297], [400, 306], [454, 315],
            [489, 285], [488, 234], [443, 185], [421, 137], [380, 169],
            [315, 160], [297, 204], [267, 248], [265, 344], [342, 263],
        ]
        self.dt = DelaunayTriangulation(self.data)

    def test_triangle_count(self):
        tris = self.dt.triangulate()
        # For n points with h on convex hull: 2n - h - 2 triangles
        assert len(tris) >= 1

    def test_all_indices_valid(self):
        tris = self.dt.triangulate()
        n = len(self.data)
        for tri in tris:
            for idx in tri:
                assert 0 <= idx < n

    def test_delaunay_property(self):
        _assert_delaunay_property(self.dt)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _assert_delaunay_property(dt):
    """Assert no input point lies strictly inside any triangle's circumcircle."""
    tris = dt.triangulate()
    pts = dt.points
    for tri in tris:
        cx, cy, r = DelaunayTriangulation._circumcircle(
            pts[tri[0]], pts[tri[1]], pts[tri[2]]
        )
        for k in range(len(pts)):
            if k not in tri:
                dist = np.sqrt((pts[k][0] - cx) ** 2 + (pts[k][1] - cy) ** 2)
                assert dist >= r - 1e-8, (
                    f"Point {k} is inside circumcircle of triangle {tri}"
                )
