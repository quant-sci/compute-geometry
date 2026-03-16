"""Tests for pydantic input validation on all algorithm models."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
import warnings
from pydantic import ValidationError

from cgeom.elements.models import (
    ConvexHullInput,
    MinimumCircleInput,
    PolygonTriangulationInput,
    VoronoiDiagramInput,
)
from cgeom.algorithms import (
    ConvexHull,
    MinimumCircle,
    PolygonTriangulation,
    VoronoiDiagram,
)


# ---------------------------------------------------------------------------
# ConvexHullInput
# ---------------------------------------------------------------------------

class TestConvexHullValidation:

    def test_too_few_points(self):
        with pytest.raises(ValidationError, match="at least 3 points"):
            ConvexHullInput(points=[[0, 0], [1, 1]])

    def test_collinear_points(self):
        with pytest.raises(ValidationError, match="collinear"):
            ConvexHullInput(points=[[0, 0], [1, 1], [2, 2]])

    def test_non_numeric(self):
        with pytest.raises(ValidationError, match="numeric"):
            ConvexHullInput(points=[["a", "b"], ["c", "d"], ["e", "f"]])

    def test_wrong_shape_3d(self):
        with pytest.raises(ValidationError, match="shape"):
            ConvexHullInput(points=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    def test_duplicate_warning(self):
        with pytest.warns(UserWarning, match="Duplicate"):
            ConvexHullInput(points=[[0, 0], [1, 0], [0, 1], [0, 0]])

    def test_valid_list(self):
        result = ConvexHullInput(points=[[0, 0], [4, 0], [2, 3]])
        assert len(result.points) == 3

    def test_valid_ndarray(self):
        arr = np.array([[0, 0], [4, 0], [2, 3]], dtype=float)
        result = ConvexHullInput(points=arr)
        assert len(result.points) == 3

    def test_valid_tuples(self):
        result = ConvexHullInput(points=[(0, 0), (4, 0), (2, 3)])
        assert len(result.points) == 3

    def test_integration_rejects_bad_input(self):
        with pytest.raises(ValidationError):
            ConvexHull([[0, 0], [1, 1]])


# ---------------------------------------------------------------------------
# MinimumCircleInput
# ---------------------------------------------------------------------------

class TestMinimumCircleValidation:

    def test_too_few_points(self):
        with pytest.raises(ValidationError, match="at least 2 points"):
            MinimumCircleInput(points=[[0, 0]])

    def test_non_numeric(self):
        with pytest.raises(ValidationError, match="numeric"):
            MinimumCircleInput(points=[["x", "y"], ["a", "b"]])

    def test_valid_two_points(self):
        result = MinimumCircleInput(points=[[0, 0], [4, 0]])
        assert len(result.points) == 2

    def test_integration_minimum_circle(self):
        mc = MinimumCircle()
        with pytest.raises(ValidationError):
            mc.minimum_circle([[0, 0]])

    def test_integration_heuristic(self):
        with pytest.raises(ValidationError):
            MinimumCircle.minimum_circle_heuristic([[0, 0]])


# ---------------------------------------------------------------------------
# PolygonTriangulationInput
# ---------------------------------------------------------------------------

class TestPolygonTriangulationValidation:

    def test_too_few_vertices(self):
        with pytest.raises(ValidationError, match="at least 3 vertices"):
            PolygonTriangulationInput(poly=[[0, 0], [1, 1]])

    def test_collinear_vertices(self):
        with pytest.raises(ValidationError, match="collinear"):
            PolygonTriangulationInput(poly=[[0, 0], [1, 0], [2, 0]])

    def test_non_numeric(self):
        with pytest.raises(ValidationError, match="numeric"):
            PolygonTriangulationInput(poly=[["a", "b"], ["c", "d"], ["e", "f"]])

    def test_duplicate_warning(self):
        with pytest.warns(UserWarning, match="Duplicate"):
            PolygonTriangulationInput(
                poly=[[0, 0], [4, 0], [2, 3], [0, 0]]
            )

    def test_valid_ndarray(self):
        arr = np.array([[0, 0], [4, 0], [2, 3]], dtype=float)
        result = PolygonTriangulationInput(poly=arr)
        assert len(result.poly) == 3

    def test_integration_rejects_bad_input(self):
        with pytest.raises(ValidationError):
            PolygonTriangulation([[0, 0], [1, 1]])


# ---------------------------------------------------------------------------
# VoronoiDiagramInput
# ---------------------------------------------------------------------------

class TestVoronoiDiagramValidation:

    def test_too_few_points(self):
        with pytest.raises(ValidationError, match="at least 2 points"):
            VoronoiDiagramInput(points=[[0, 0]])

    def test_same_y_first_two(self):
        with pytest.raises(ValidationError, match="same y-coordinate"):
            VoronoiDiagramInput(points=[[0, 5], [3, 5]])

    def test_non_numeric(self):
        with pytest.raises(ValidationError, match="numeric"):
            VoronoiDiagramInput(points=[["a", "b"], ["c", "d"]])

    def test_valid_input(self):
        result = VoronoiDiagramInput(points=[[0, 0], [4, 2]])
        assert len(result.points) == 2

    def test_integration_rejects_same_y(self):
        with pytest.raises(ValidationError):
            VoronoiDiagram([[0, 5], [3, 5]])
