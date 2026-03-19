"""Line segment intersection using Bentley-Ottmann sweep line and brute force."""

import heapq
from enum import IntEnum

import numpy as np


# ---------------------------------------------------------------------------
# Event types — ordering matters for the heap (left < intersection < right)
# ---------------------------------------------------------------------------

class _EventType(IntEnum):
    LEFT = 0
    INTERSECTION = 1
    RIGHT = 2


# ---------------------------------------------------------------------------
# Geometric helpers
# ---------------------------------------------------------------------------

_EPS = 1e-9


def _segment_intersection(s1, s2):
    """Parametric intersection of two segments. Returns [x, y] or None."""
    x1, y1 = s1[0]
    x2, y2 = s1[1]
    x3, y3 = s2[0]
    x4, y4 = s2[1]

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    if abs(denom) < _EPS:
        # Parallel or collinear — check for overlap
        return _collinear_intersection(s1, s2)

    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom

    if -_EPS <= t <= 1 + _EPS and -_EPS <= u <= 1 + _EPS:
        ix = x1 + t * (x2 - x1)
        iy = y1 + t * (y2 - y1)
        return [ix, iy]
    return None


def _collinear_intersection(s1, s2):
    """Check if two collinear segments overlap, return one endpoint of overlap or None."""
    x1, y1 = s1[0]
    x2, y2 = s1[1]
    x3, y3 = s2[0]
    x4, y4 = s2[1]

    # Check if segments are on the same line
    cross = (x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)
    if abs(cross) > _EPS:
        return None

    # Project onto the longer axis
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    if dx >= dy:
        a_lo, a_hi = min(x1, x2), max(x1, x2)
        b_lo, b_hi = min(x3, x4), max(x3, x4)
    else:
        a_lo, a_hi = min(y1, y2), max(y1, y2)
        b_lo, b_hi = min(y3, y4), max(y3, y4)

    lo = max(a_lo, b_lo)
    hi = min(a_hi, b_hi)

    if lo > hi + _EPS:
        return None

    # Return the start of overlap as the intersection point
    if dx >= dy:
        t = (lo - x1) / (x2 - x1) if abs(x2 - x1) > _EPS else 0
        return [x1 + t * (x2 - x1), y1 + t * (y2 - y1)]
    else:
        t = (lo - y1) / (y2 - y1) if abs(y2 - y1) > _EPS else 0
        return [x1 + t * (x2 - x1), y1 + t * (y2 - y1)]


def _y_at_x(seg, x):
    """Y-coordinate where segment crosses the sweep line at x."""
    x1, y1 = seg[0]
    x2, y2 = seg[1]
    if abs(x2 - x1) < _EPS:
        return min(y1, y2)
    t = (x - x1) / (x2 - x1)
    return y1 + t * (y2 - y1)


# ---------------------------------------------------------------------------
# Sweep-line status structure
# ---------------------------------------------------------------------------

class _SweepLineStatus:
    """Maintains sorted list of active segments ordered by y at current sweep x."""

    def __init__(self):
        self._segments = []  # list of segment indices
        self._sweep_x = 0.0
        self._seg_data = None  # reference to segments array

    def set_context(self, segments):
        self._seg_data = segments

    def set_sweep_x(self, x):
        self._sweep_x = x

    def _sort_key(self, seg_idx):
        return _y_at_x(self._seg_data[seg_idx], self._sweep_x)

    def insert(self, seg_idx):
        key = self._sort_key(seg_idx)
        # Binary search for insertion point
        lo, hi = 0, len(self._segments)
        while lo < hi:
            mid = (lo + hi) // 2
            if self._sort_key(self._segments[mid]) < key - _EPS:
                lo = mid + 1
            else:
                hi = mid
        self._segments.insert(lo, seg_idx)

    def remove(self, seg_idx):
        try:
            self._segments.remove(seg_idx)
        except ValueError:
            pass

    def _find_pos(self, seg_idx):
        for i, s in enumerate(self._segments):
            if s == seg_idx:
                return i
        return -1

    def swap(self, seg_a, seg_b):
        pos_a = self._find_pos(seg_a)
        pos_b = self._find_pos(seg_b)
        if pos_a >= 0 and pos_b >= 0:
            self._segments[pos_a] = seg_b
            self._segments[pos_b] = seg_a

    def neighbors(self, seg_idx):
        """Return (above, below) neighbor indices or None."""
        pos = self._find_pos(seg_idx)
        if pos < 0:
            return None, None
        above = self._segments[pos - 1] if pos > 0 else None
        below = self._segments[pos + 1] if pos < len(self._segments) - 1 else None
        return above, below


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class SegmentIntersection:
    """Find intersection points among a set of 2D line segments.

    Provides both a Bentley-Ottmann sweep line method and a brute-force
    O(n^2) pairwise check.

    Attributes:
        segments (np.ndarray): Array of segments with shape (n, 2, 2).
    """

    def __init__(self, segments):
        """Initialize with a set of 2D line segments.

        Args:
            segments: Collection of segments as a list, tuple, or numpy array.
                Each segment is [[x1, y1], [x2, y2]].

        Raises:
            pydantic.ValidationError: If fewer than 2 segments, zero-length
                segment, non-numeric, or wrong shape.
        """
        from cgeom.elements.models import SegmentIntersectionInput
        validated = SegmentIntersectionInput(segments=segments)
        self.segments = np.array(validated.segments)

    def find_intersections(self):
        """Find all intersection points using Bentley-Ottmann sweep line.

        Returns:
            list[list[float]]: Unique intersection points as [[x, y], ...].
        """
        segs = self._normalize_segments()
        n = len(segs)

        if n < 2:
            return []

        # Build event queue
        events = []
        for i in range(n):
            lx, ly = segs[i][0]
            rx, ry = segs[i][1]
            heapq.heappush(events, (lx, _EventType.LEFT, ly, i, -1))
            heapq.heappush(events, (rx, _EventType.RIGHT, ry, i, -1))

        status = _SweepLineStatus()
        status.set_context(segs)

        found_points = {}  # (round_x, round_y) -> [x, y]
        found_pairs = set()  # frozenset(i, j) to avoid duplicate events

        def _check_and_add(seg_a, seg_b):
            if seg_a is None or seg_b is None:
                return
            pair = frozenset((seg_a, seg_b))
            if pair in found_pairs:
                return
            pt = _segment_intersection(segs[seg_a], segs[seg_b])
            if pt is not None:
                rkey = (round(pt[0], 9), round(pt[1], 9))
                if rkey not in found_points:
                    found_points[rkey] = pt
                found_pairs.add(pair)
                heapq.heappush(events, (
                    pt[0], _EventType.INTERSECTION, pt[1], seg_a, seg_b
                ))

        while events:
            x, etype, y, s1, s2 = heapq.heappop(events)
            status.set_sweep_x(x)

            if etype == _EventType.LEFT:
                status.insert(s1)
                above, below = status.neighbors(s1)
                _check_and_add(s1, above)
                _check_and_add(s1, below)

            elif etype == _EventType.RIGHT:
                above, below = status.neighbors(s1)
                _check_and_add(above, below)
                status.remove(s1)

            elif etype == _EventType.INTERSECTION:
                status.swap(s1, s2)
                # After swap, check new neighbors
                above_s1, below_s1 = status.neighbors(s1)
                above_s2, below_s2 = status.neighbors(s2)
                _check_and_add(s1, above_s1)
                _check_and_add(s1, below_s1)
                _check_and_add(s2, above_s2)
                _check_and_add(s2, below_s2)

        return list(found_points.values())

    def find_intersections_brute_force(self):
        """Find all intersection points using O(n^2) pairwise check.

        Returns:
            list[list[float]]: Unique intersection points as [[x, y], ...].
        """
        segs = self.segments.tolist()
        n = len(segs)
        found = {}

        for i in range(n):
            for j in range(i + 1, n):
                pt = _segment_intersection(segs[i], segs[j])
                if pt is not None:
                    rkey = (round(pt[0], 9), round(pt[1], 9))
                    if rkey not in found:
                        found[rkey] = pt

        return list(found.values())

    def get_intersection_pairs(self):
        """Return intersection details with segment indices.

        Returns:
            list[tuple[int, int, list[float]]]: Each entry is
                (seg_i, seg_j, [x, y]).
        """
        segs = self.segments.tolist()
        n = len(segs)
        pairs = []

        for i in range(n):
            for j in range(i + 1, n):
                pt = _segment_intersection(segs[i], segs[j])
                if pt is not None:
                    pairs.append((i, j, pt))

        return pairs

    def get_segments(self):
        """Return segments as a plain list.

        Returns:
            list[list[list[float]]]: Segments as [[[x1,y1],[x2,y2]], ...].
        """
        return self.segments.tolist()

    def plot(self, title="Segment Intersections"):
        """Deprecated: use cgeom.visualization.plot_intersections() instead."""
        import warnings
        warnings.warn(
            "SegmentIntersection.plot() is deprecated. "
            "Use cgeom.visualization.plot_intersections(si_obj, title) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        from cgeom.visualization import plot_intersections
        plot_intersections(self, title)

    def _normalize_segments(self):
        """Return segments normalized left-to-right (or bottom-to-top for vertical)."""
        segs = self.segments.tolist()
        normalized = []
        for seg in segs:
            (x1, y1), (x2, y2) = seg
            if x1 > x2 + _EPS or (abs(x1 - x2) < _EPS and y1 > y2):
                normalized.append([[x2, y2], [x1, y1]])
            else:
                normalized.append([[x1, y1], [x2, y2]])
        return normalized
