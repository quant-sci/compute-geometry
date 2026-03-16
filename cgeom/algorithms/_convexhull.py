import numpy as np
import math

class ConvexHull:
    """Class to find the convex hull of a set of points using the Gift Wrapping algorithm""

    Reference:
        R.A. Jarvis, On the identification of the convex hull of a finite set of points in the plane, Information Processing Letters,
        Volume 2, Issue 1, 1973, Pages 18-21, ISSN 0020-0190, https://doi.org/10.1016/0020-0190(73)90020-3.

    Attributes:
        points (np.array): Array of points
    """
    def __init__(self, points):
        """Initialize the convex hull with a set of 2D points.

        Args:
            points: Collection of 2D points as a list, tuple, or numpy array
                with shape (n, 2).

        Raises:
            pydantic.ValidationError: If fewer than 3 points, all collinear,
                non-numeric, or wrong shape.
        """
        from cgeom.elements.models import ConvexHullInput
        validated = ConvexHullInput(points=points)
        self.points = np.array(validated.points)

    def find_low_right(self):
        """Find the point with the lowest y-coordinate (rightmost if tied).

        Returns:
            list: The [x, y] coordinates of the lowest-rightmost point.
        """
        lower_y = float("inf")
        for point in self.points:
            if point[1] < lower_y:
                lower_y = point[1]
                lower = point
        return list(lower)

    def low_angle(self, b, c):
        """Find the point that makes the smallest angle from vector c->b.

        Used by the Gift Wrapping algorithm to select the next hull vertex.

        Args:
            b: Current point on the hull boundary.
            c: Previous point (defines the reference direction).

        Returns:
            list: The [x, y] coordinates of the point with the lowest angle.
        """
        low_angle = float("inf")
        for a in self.points:
            if a[0] != b[0] or a[1] != b[1]:
                if a[0] != c[0] or a[1] != c[1]:
                    u = [c[0] - b[0], c[1] - b[1]]
                    v = [a[0] - b[0], a[1] - b[1]]
                    mag_u = math.sqrt((u[0]) ** 2 + (u[1]) ** 2)
                    mag_v = math.sqrt((v[0]) ** 2 + (v[1]) ** 2)
                    u = [u[0] / mag_u, u[1] / mag_u]
                    v = [v[0] / mag_v, v[1] / mag_v]
                    angle = math.acos(u[0] * v[0] + u[1] * v[1])
                    orient = b[0] * c[1] + b[1] * a[0] + c[0] * a[1] - a[0] * c[1] - a[1] * b[0] - c[0] * b[1]
                    if orient < 0:
                        angle = 2 * math.pi - angle
                    if angle < low_angle:
                        low_angle = angle
                        low_point = a
        return list(low_point)

    def convex_hull(self):
        """Compute the convex hull using the Gift Wrapping (Jarvis march) algorithm.

        Returns:
            list[list]: Ordered list of [x, y] vertices forming the convex hull.
        """
        a = self.find_low_right()
        b = self.low_angle(a, [1, a[1]])
        init = a
        ch = [a]
        while b[0] != init[0] or b[1] != init[1]:
            ch.append(b)
            a = b
            b = self.low_angle(a, ch[-2])
        return ch

    def plot(self, title='Convex Hull for a set of points'):
        """Deprecated: use cgeom.visualization.plot_convex_hull() instead."""
        import warnings
        warnings.warn(
            "ConvexHull.plot() is deprecated. "
            "Use cgeom.visualization.plot_convex_hull(hull_obj, title) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        from cgeom.visualization import plot_convex_hull
        plot_convex_hull(self, title)

    def get_indexes(self):
        """Get the indices of the convex hull vertices in the original points array.

        Returns:
            list[int]: Indices into ``self.points`` for each hull vertex.
        """
        indexes = []
        for element in self.convex_hull():
            j = 0
            while element[0] != self.points[j][0] or element[1] != self.points[j][1]:
                j += 1
            indexes.append(j)
        return(indexes)
