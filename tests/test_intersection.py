"""Tests for line segment intersection (Bentley-Ottmann and brute force)."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from cgeom.algorithms import SegmentIntersection


def _round_pts(pts):
    """Round intersection points for stable comparison."""
    return sorted((round(p[0], 6), round(p[1], 6)) for p in pts)


class TestIntersectionSimpleCross:
    """Two segments forming an X — one intersection."""

    def setup_method(self):
        self.segs = [[[0, 0], [4, 4]], [[0, 4], [4, 0]]]
        self.si = SegmentIntersection(self.segs)

    def test_one_intersection(self):
        pts = self.si.find_intersections()
        assert len(pts) == 1

    def test_intersection_at_center(self):
        pts = self.si.find_intersections()
        assert abs(pts[0][0] - 2.0) < 1e-6
        assert abs(pts[0][1] - 2.0) < 1e-6

    def test_brute_force_agrees(self):
        pts_bo = _round_pts(self.si.find_intersections())
        pts_bf = _round_pts(self.si.find_intersections_brute_force())
        assert pts_bo == pts_bf


class TestIntersectionParallel:
    """Two parallel segments — no intersection."""

    def setup_method(self):
        self.segs = [[[0, 0], [4, 0]], [[0, 2], [4, 2]]]
        self.si = SegmentIntersection(self.segs)

    def test_no_intersection(self):
        pts = self.si.find_intersections()
        assert len(pts) == 0


class TestIntersectionSharedEndpoint:
    """Two segments sharing one endpoint."""

    def setup_method(self):
        self.segs = [[[0, 0], [2, 2]], [[2, 2], [4, 0]]]
        self.si = SegmentIntersection(self.segs)

    def test_one_intersection(self):
        pts = self.si.find_intersections()
        assert len(pts) == 1

    def test_at_shared_point(self):
        pts = self.si.find_intersections()
        assert abs(pts[0][0] - 2.0) < 1e-6
        assert abs(pts[0][1] - 2.0) < 1e-6


class TestIntersectionTOnSegment:
    """T-junction: endpoint of one segment on interior of another."""

    def setup_method(self):
        self.segs = [[[0, 0], [4, 0]], [[2, -2], [2, 0]]]
        self.si = SegmentIntersection(self.segs)

    def test_one_intersection(self):
        pts = self.si.find_intersections()
        assert len(pts) == 1

    def test_at_t_junction(self):
        pts = self.si.find_intersections()
        assert abs(pts[0][0] - 2.0) < 1e-6
        assert abs(pts[0][1] - 0.0) < 1e-6


class TestIntersectionVertical:
    """Vertical segment crossing horizontal — one intersection."""

    def setup_method(self):
        self.segs = [[[0, 2], [4, 2]], [[2, 0], [2, 4]]]
        self.si = SegmentIntersection(self.segs)

    def test_one_intersection(self):
        pts = self.si.find_intersections()
        assert len(pts) == 1

    def test_at_crossing(self):
        pts = self.si.find_intersections()
        assert abs(pts[0][0] - 2.0) < 1e-6
        assert abs(pts[0][1] - 2.0) < 1e-6


class TestIntersectionMultiple:
    """Multiple segments with known intersection count."""

    def setup_method(self):
        # Star pattern — 4 segments with multiple crossings
        self.segs = [
            [[0, 0], [4, 4]],
            [[0, 4], [4, 0]],
            [[0, 2], [4, 2]],
            [[2, 0], [2, 4]],
        ]
        self.si = SegmentIntersection(self.segs)

    def test_correct_count(self):
        pts = self.si.find_intersections()
        # All pairs intersect except parallel ones; 4 segments = 6 pairs, all cross
        # (0,1) at (2,2), (0,2) at (1,1)? Let's just check brute force agrees
        pts_bf = self.si.find_intersections_brute_force()
        assert len(pts) == len(pts_bf)

    def test_brute_force_agreement(self):
        pts_bo = _round_pts(self.si.find_intersections())
        pts_bf = _round_pts(self.si.find_intersections_brute_force())
        assert pts_bo == pts_bf


class TestIntersectionNoIntersection:
    """Non-crossing segments — 0 intersections."""

    def setup_method(self):
        self.segs = [
            [[0, 0], [1, 0]],
            [[2, 2], [3, 2]],
            [[5, 5], [6, 5]],
        ]
        self.si = SegmentIntersection(self.segs)

    def test_no_intersections(self):
        pts = self.si.find_intersections()
        assert len(pts) == 0

    def test_brute_force_agrees(self):
        pts_bf = self.si.find_intersections_brute_force()
        assert len(pts_bf) == 0


class TestIntersectionBruteForceAgreement:
    """Parametrized random sets — Bentley-Ottmann matches brute force."""

    @pytest.mark.parametrize("seed", [42, 123, 999])
    def test_agreement(self, seed):
        rng = np.random.default_rng(seed)
        n = 10
        segs = rng.uniform(-50, 50, size=(n, 2, 2)).tolist()
        si = SegmentIntersection(segs)
        pts_bo = _round_pts(si.find_intersections())
        pts_bf = _round_pts(si.find_intersections_brute_force())
        assert pts_bo == pts_bf


class TestIntersectionOutputFormats:
    """Check that output types and shapes are correct."""

    def setup_method(self):
        self.segs = [[[0, 0], [4, 4]], [[0, 4], [4, 0]]]
        self.si = SegmentIntersection(self.segs)

    def test_find_intersections_returns_list_of_lists(self):
        pts = self.si.find_intersections()
        assert isinstance(pts, list)
        for pt in pts:
            assert isinstance(pt, list)
            assert len(pt) == 2
            for v in pt:
                assert isinstance(v, float)

    def test_pairs_format(self):
        pairs = self.si.get_intersection_pairs()
        assert isinstance(pairs, list)
        for item in pairs:
            assert len(item) == 3
            seg_i, seg_j, pt = item
            assert isinstance(seg_i, int)
            assert isinstance(seg_j, int)
            assert isinstance(pt, list)
            assert len(pt) == 2

    def test_segments_format(self):
        segs = self.si.get_segments()
        assert isinstance(segs, list)
        for seg in segs:
            assert len(seg) == 2
            for pt in seg:
                assert len(pt) == 2


class TestIntersectionLargerSet:
    """15+ segments — structural validity and brute-force agreement."""

    def setup_method(self):
        rng = np.random.default_rng(7)
        self.segs = rng.uniform(-100, 100, size=(18, 2, 2)).tolist()
        self.si = SegmentIntersection(self.segs)

    def test_structural_validity(self):
        pts = self.si.find_intersections()
        assert isinstance(pts, list)
        for pt in pts:
            assert len(pt) == 2

    def test_brute_force_agreement(self):
        pts_bo = _round_pts(self.si.find_intersections())
        pts_bf = _round_pts(self.si.find_intersections_brute_force())
        assert pts_bo == pts_bf
