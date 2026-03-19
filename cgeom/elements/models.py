"""Pydantic validation models for computational geometry algorithm inputs."""

import warnings
from typing import List

import numpy as np
from pydantic import BaseModel, model_validator


def _to_point_list(raw) -> List[List[float]]:
    """Accept list/tuple/ndarray, validate shape (n, 2), return List[List[float]]."""
    if raw is None:
        raise ValueError("Points data is required")

    if isinstance(raw, np.ndarray):
        arr = raw.astype(float)
    else:
        try:
            arr = np.array(raw, dtype=float)
        except (ValueError, TypeError):
            raise ValueError("Points must contain only numeric values")

    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"Points must have shape (n, 2), got shape {arr.shape}")

    return arr.tolist()


def _has_duplicates(points: List[List[float]]) -> bool:
    """Set-based duplicate check on point coordinates."""
    seen = set()
    for p in points:
        key = (p[0], p[1])
        if key in seen:
            return True
        seen.add(key)
    return False


def _all_collinear(points: List[List[float]]) -> bool:
    """Cross-product test: returns True if all points lie on the same line."""
    if len(points) < 3:
        return True
    x0, y0 = points[0]
    x1, y1 = points[1]
    for i in range(2, len(points)):
        x2, y2 = points[i]
        cross = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
        if abs(cross) > 1e-10:
            return False
    return True


class ConvexHullInput(BaseModel):
    """Validation model for ConvexHull inputs."""

    points: List[List[float]]

    @model_validator(mode="before")
    @classmethod
    def validate_input(cls, data):
        if isinstance(data, dict):
            raw = data.get("points")
        else:
            raw = data

        points = _to_point_list(raw)

        if len(points) < 3:
            raise ValueError("ConvexHull requires at least 3 points")

        if _all_collinear(points):
            raise ValueError(
                "All points are collinear; convex hull is undefined"
            )

        if _has_duplicates(points):
            warnings.warn(
                "Duplicate points detected; they will be included but may "
                "not affect the hull",
                UserWarning,
            )

        return {"points": points}


class MinimumCircleInput(BaseModel):
    """Validation model for MinimumCircle inputs."""

    points: List[List[float]]

    @model_validator(mode="before")
    @classmethod
    def validate_input(cls, data):
        if isinstance(data, dict):
            raw = data.get("points")
        else:
            raw = data

        points = _to_point_list(raw)

        if len(points) < 2:
            raise ValueError("MinimumCircle requires at least 2 points")

        if _has_duplicates(points):
            warnings.warn("Duplicate points detected", UserWarning)

        return {"points": points}


class PolygonTriangulationInput(BaseModel):
    """Validation model for PolygonTriangulation inputs."""

    poly: List[List[float]]
    poly_name: str = "Polygon"

    @model_validator(mode="before")
    @classmethod
    def validate_input(cls, data):
        if not isinstance(data, dict):
            raise ValueError(
                "PolygonTriangulationInput expects keyword arguments "
                "(poly=..., poly_name=...)"
            )

        raw = data.get("poly")
        poly_name = data.get("poly_name", "Polygon")

        points = _to_point_list(raw)

        if len(points) < 3:
            raise ValueError(
                "Polygon triangulation requires at least 3 vertices"
            )

        if _all_collinear(points):
            raise ValueError(
                "All vertices are collinear; polygon is degenerate"
            )

        if _has_duplicates(points):
            warnings.warn(
                "Duplicate vertices detected in polygon", UserWarning
            )

        return {"poly": points, "poly_name": poly_name}


class VoronoiDiagramInput(BaseModel):
    """Validation model for VoronoiDiagram inputs."""

    points: List[List[float]]

    @model_validator(mode="before")
    @classmethod
    def validate_input(cls, data):
        if isinstance(data, dict):
            raw = data.get("points")
        else:
            raw = data

        points = _to_point_list(raw)

        if len(points) < 2:
            raise ValueError("Voronoi diagram requires at least 2 points")

        if abs(points[0][1] - points[1][1]) < 1e-10:
            raise ValueError(
                "First two points have the same y-coordinate, which causes "
                "division by zero in the Voronoi construction. "
                "Reorder points so the first two have different y-coordinates."
            )

        if _has_duplicates(points):
            warnings.warn("Duplicate points detected", UserWarning)

        return {"points": points}


def _to_segment_list(raw) -> List[List[List[float]]]:
    """Accept list/tuple/ndarray of segments, validate shape (n, 2, 2)."""
    if raw is None:
        raise ValueError("Segments data is required")

    # Handle list of Segment objects
    from cgeom.elements.elements import Segment
    if isinstance(raw, (list, tuple)) and len(raw) > 0 and isinstance(raw[0], Segment):
        raw = [s.to_list() for s in raw]

    if isinstance(raw, np.ndarray):
        arr = raw.astype(float)
    else:
        try:
            arr = np.array(raw, dtype=float)
        except (ValueError, TypeError):
            raise ValueError("Segments must contain only numeric values")

    if arr.ndim != 3 or arr.shape[1] != 2 or arr.shape[2] != 2:
        raise ValueError(
            f"Segments must have shape (n, 2, 2), got shape {arr.shape}"
        )

    return arr.tolist()


class SegmentIntersectionInput(BaseModel):
    """Validation model for SegmentIntersection inputs."""

    segments: List[List[List[float]]]

    @model_validator(mode="before")
    @classmethod
    def validate_input(cls, data):
        if isinstance(data, dict):
            raw = data.get("segments")
        else:
            raw = data

        segments = _to_segment_list(raw)

        if len(segments) < 2:
            raise ValueError(
                "SegmentIntersection requires at least 2 segments"
            )

        for i, seg in enumerate(segments):
            x1, y1 = seg[0]
            x2, y2 = seg[1]
            if abs(x1 - x2) < 1e-12 and abs(y1 - y2) < 1e-12:
                raise ValueError(
                    f"Segment {i} has zero length (identical endpoints)"
                )

        seen = set()
        for seg in segments:
            key = (tuple(seg[0]), tuple(seg[1]))
            rkey = (tuple(seg[1]), tuple(seg[0]))
            if key in seen or rkey in seen:
                warnings.warn(
                    "Duplicate segments detected", UserWarning
                )
                break
            seen.add(key)

        return {"segments": segments}


class DelaunayTriangulationInput(BaseModel):
    """Validation model for DelaunayTriangulation inputs."""

    points: List[List[float]]

    @model_validator(mode="before")
    @classmethod
    def validate_input(cls, data):
        if isinstance(data, dict):
            raw = data.get("points")
        else:
            raw = data

        points = _to_point_list(raw)

        if len(points) < 3:
            raise ValueError(
                "Delaunay triangulation requires at least 3 points"
            )

        if _all_collinear(points):
            raise ValueError(
                "All points are collinear; triangulation is undefined"
            )

        if _has_duplicates(points):
            warnings.warn(
                "Duplicate points detected; they will be included but may "
                "not affect the triangulation",
                UserWarning,
            )

        return {"points": points}
