import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from cgeom.algorithms import PolygonTriangulation


class TestTriangulationTriangle:
    """A triangle needs zero diagonals."""

    def setup_method(self):
        poly = np.array([[0, 0], [4, 0], [2, 3]], dtype=float)
        self.tri = PolygonTriangulation(poly)

    def test_no_diagonals(self):
        assert len(self.tri.diagonals) == 0


class TestTriangulationSquare:
    """A convex quadrilateral needs exactly 1 diagonal."""

    def setup_method(self):
        # CCW order
        poly = np.array([[0, 0], [4, 0], [4, 4], [0, 4]], dtype=float)
        self.tri = PolygonTriangulation(poly)

    def test_one_diagonal(self):
        assert len(self.tri.diagonals) == 1

    def test_diagonal_connects_polygon_vertices(self):
        diag = self.tri.diagonals[0]
        vertices = {(0, 0), (4, 0), (4, 4), (0, 4)}
        assert (diag[0][0], diag[0][1]) in vertices
        assert (diag[1][0], diag[1][1]) in vertices


class TestTriangulationPentagon:
    """A convex pentagon needs n-3 = 2 diagonals."""

    def setup_method(self):
        # regular pentagon (CCW)
        angles = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        poly = np.column_stack([np.cos(angles), np.sin(angles)])
        self.tri = PolygonTriangulation(poly)

    def test_two_diagonals(self):
        assert len(self.tri.diagonals) == 2


class TestTriangulationHexagon:
    """A convex hexagon needs n-3 = 3 diagonals."""

    def setup_method(self):
        angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
        poly = np.column_stack([np.cos(angles), np.sin(angles)])
        self.tri = PolygonTriangulation(poly)

    def test_three_diagonals(self):
        assert len(self.tri.diagonals) == 3


class TestTriangulationDiagonalCount:
    """For any simple polygon with n vertices, triangulation produces n-3 diagonals."""

    @pytest.mark.parametrize("n", [4, 5, 6, 7, 8])
    def test_n_minus_3_diagonals(self, n):
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        poly = np.column_stack([np.cos(angles), np.sin(angles)])
        tri = PolygonTriangulation(poly)
        assert len(tri.diagonals) == n - 3


class TestTriangulationLShaped:
    """Concave (L-shaped) polygon."""

    def setup_method(self):
        # L-shape in CCW order
        poly = np.array([
            [0, 0], [4, 0], [4, 2], [2, 2], [2, 4], [0, 4]
        ], dtype=float)
        self.tri = PolygonTriangulation(poly)

    def test_diagonal_count(self):
        # 6 vertices -> 3 diagonals
        assert len(self.tri.diagonals) == 3

    def test_diagonals_use_polygon_vertices(self):
        vertices = {(0, 0), (4, 0), (4, 2), (2, 2), (2, 4), (0, 4)}
        for diag in self.tri.diagonals:
            assert (diag[0][0], diag[0][1]) in vertices
            assert (diag[1][0], diag[1][1]) in vertices
