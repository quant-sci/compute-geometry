import numpy as np


class DelaunayTriangulation:
    """Delaunay triangulation of a 2D point set using the Bowyer-Watson algorithm.

    Attributes:
        points (np.ndarray): Array of input points with shape (n, 2).
    """

    def __init__(self, points):
        """Initialize with a set of 2D points.

        Args:
            points: Collection of 2D points as a list, tuple, or numpy array
                with shape (n, 2).

        Raises:
            pydantic.ValidationError: If fewer than 3 points, all collinear,
                non-numeric, or wrong shape.
        """
        from cgeom.elements.models import DelaunayTriangulationInput
        validated = DelaunayTriangulationInput(points=points)
        self.points = np.array(validated.points)
        self._triangles = None

    def triangulate(self):
        """Compute the Delaunay triangulation using Bowyer-Watson.

        Returns:
            list[list[int]]: List of triangles, each a sorted list of three
                point indices into ``self.points``.
        """
        if self._triangles is not None:
            return self._triangles

        pts = self.points
        n = len(pts)

        # Build a super-triangle that encloses all points
        min_x = pts[:, 0].min()
        max_x = pts[:, 0].max()
        min_y = pts[:, 1].min()
        max_y = pts[:, 1].max()
        dx = max_x - min_x
        dy = max_y - min_y
        delta = max(dx, dy, 1.0)
        mid_x = (min_x + max_x) / 2
        mid_y = (min_y + max_y) / 2
        margin = 20.0

        # Super-triangle vertices stored at indices n, n+1, n+2
        super_pts = np.array([
            [mid_x - margin * delta, mid_y - delta],
            [mid_x + margin * delta, mid_y - delta],
            [mid_x, mid_y + margin * delta],
        ])
        all_pts = np.vstack([pts, super_pts])

        # Each triangle is a frozenset of 3 indices
        triangles = {frozenset([n, n + 1, n + 2])}

        for i in range(n):
            # Find bad triangles whose circumcircle contains point i
            bad = set()
            for tri in triangles:
                a, b, c = tri
                if self._in_circumcircle(all_pts[i], all_pts[a], all_pts[b], all_pts[c]):
                    bad.add(tri)

            # Find boundary polygon: edges belonging to exactly one bad triangle
            edge_count = {}
            for tri in bad:
                verts = list(tri)
                for j in range(3):
                    edge = frozenset([verts[j], verts[(j + 1) % 3]])
                    edge_count[edge] = edge_count.get(edge, 0) + 1

            boundary = [edge for edge, cnt in edge_count.items() if cnt == 1]

            # Remove bad triangles
            triangles -= bad

            # Create new triangles from point i to each boundary edge
            for edge in boundary:
                v1, v2 = edge
                triangles.add(frozenset([i, v1, v2]))

        # Remove triangles that reference super-triangle vertices
        super_indices = {n, n + 1, n + 2}
        triangles = {tri for tri in triangles if not tri & super_indices}

        self._triangles = [sorted(tri) for tri in triangles]
        self._triangles.sort()
        return self._triangles

    def get_triangles(self):
        """Return triangles as coordinate triples.

        Returns:
            list[list[list[float]]]: Each triangle is [[x,y], [x,y], [x,y]].
        """
        tris = self.triangulate()
        return [[self.points[i].tolist() for i in tri] for tri in tris]

    def get_edges(self):
        """Return unique edges of the triangulation.

        Returns:
            list[list[list[float]]]: Each edge is [[x1,y1], [x2,y2]].
        """
        tris = self.triangulate()
        seen = set()
        edges = []
        for tri in tris:
            for j in range(3):
                edge = (tri[j], tri[(j + 1) % 3])
                key = (min(edge), max(edge))
                if key not in seen:
                    seen.add(key)
                    edges.append([self.points[key[0]].tolist(),
                                  self.points[key[1]].tolist()])
        return edges

    def get_circumcircles(self):
        """Return the circumcircle for each triangle.

        Returns:
            list[list]: Each entry is [[cx, cy], radius].
        """
        tris = self.triangulate()
        result = []
        for tri in tris:
            cx, cy, r = self._circumcircle(
                self.points[tri[0]], self.points[tri[1]], self.points[tri[2]]
            )
            result.append([[cx, cy], r])
        return result

    def plot(self, title="Delaunay Triangulation"):
        """Deprecated: use cgeom.visualization.plot_delaunay() instead."""
        import warnings
        warnings.warn(
            "DelaunayTriangulation.plot() is deprecated. "
            "Use cgeom.visualization.plot_delaunay(dt_obj, title) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        from cgeom.visualization import plot_delaunay
        plot_delaunay(self, title)

    @staticmethod
    def _circumcircle(a, b, c):
        """Compute the circumcircle of triangle (a, b, c).

        Returns:
            tuple: (cx, cy, radius).
        """
        ax, ay = float(a[0]), float(a[1])
        bx, by = float(b[0]), float(b[1])
        cx, cy = float(c[0]), float(c[1])

        D = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        ux = ((ax * ax + ay * ay) * (by - cy)
              + (bx * bx + by * by) * (cy - ay)
              + (cx * cx + cy * cy) * (ay - by)) / D
        uy = ((ax * ax + ay * ay) * (cx - bx)
              + (bx * bx + by * by) * (ax - cx)
              + (cx * cx + cy * cy) * (bx - ax)) / D
        r = np.sqrt((ax - ux) ** 2 + (ay - uy) ** 2)
        return ux, uy, r

    @staticmethod
    def _in_circumcircle(p, a, b, c):
        """Determinant-based test: is point p inside circumcircle of (a, b, c)?

        Orientation-independent: works regardless of vertex winding order.
        Returns True if p is strictly inside (with small epsilon tolerance).
        """
        ax, ay = float(a[0]) - float(p[0]), float(a[1]) - float(p[1])
        bx, by = float(b[0]) - float(p[0]), float(b[1]) - float(p[1])
        cx, cy = float(c[0]) - float(p[0]), float(c[1]) - float(p[1])

        det = (ax * ax + ay * ay) * (bx * cy - cx * by) \
            - (bx * bx + by * by) * (ax * cy - cx * ay) \
            + (cx * cx + cy * cy) * (ax * by - bx * ay)

        # The sign of det depends on triangle orientation (CCW vs CW).
        # Compute orientation and flip if clockwise.
        orient = (float(b[0]) - float(a[0])) * (float(c[1]) - float(a[1])) \
               - (float(b[1]) - float(a[1])) * (float(c[0]) - float(a[0]))
        if orient < 0:
            det = -det

        return det > 1e-10
