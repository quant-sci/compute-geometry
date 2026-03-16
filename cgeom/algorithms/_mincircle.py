import random
import numpy as np
from math import sqrt

class MinimumCircle:
    """Find the minimum enclosing circle for a set of 2D points.

    Provides both an exact randomized incremental algorithm and a
    faster heuristic approximation.
    """

    def __init__(self):
        pass

    @staticmethod
    def minimum_circle_heuristic(num_points):
        """
        Heuristic algorithm to find the minimum circle that contains all the points
        :param num_points: list of points
        :return: circle
        """
        from cgeom.elements.models import MinimumCircleInput
        validated = MinimumCircleInput(points=num_points)
        num_points = [p[:] for p in validated.points]
        X_max = max([abs(x) for x, y in num_points])
        X_max = [[x, y] for x, y in num_points if abs(x) == X_max][0]
        X_min = min([abs(x) for x, y in num_points])
        X_min = [[x, y] for x, y in num_points if abs(x) == X_min][0]
        Y_max = max([abs(y) for x, y in num_points])
        Y_max = [[x, y] for x, y in num_points if abs(y) == Y_max][0]
        Y_min = min([abs(y) for x, y in num_points])
        Y_min = [[x, y] for x, y in num_points if abs(y) == Y_min][0]
        points = [X_max, X_min, Y_max, Y_min]

        d_max = 0
        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                d = sqrt((points[i][0] - points[j][0]) ** 2 + (points[i][1] - points[j][1]) ** 2)
                if d > d_max:
                    d_max = d
                    Pi = points[i]
                    Pj = points[j]

        circle = [[0, 0], 0]
        circle[0] = [(Pi[0] + Pj[0]) / 2, (Pi[1] + Pj[1]) / 2]
        circle[1] = d_max / 2
        for point in num_points:
            d_mod = sqrt((point[0] - circle[0][0]) ** 2 + (point[1] - circle[0][1]) ** 2)
            if d_mod > circle[1]:
                d_norm = [(point[0] - circle[0][0]) / d_mod, (point[1] - circle[0][1]) / d_mod]
                circle[0][0] = circle[0][0] + ((d_mod - circle[1]) / 2) * d_norm[0]
                circle[0][1] = circle[0][1] + ((d_mod - circle[1]) / 2) * d_norm[1]
                circle[1] = (d_mod + circle[1]) / 2

        return circle

    @staticmethod
    def randomized_permutation(num_points):
        """
        Randomized permutation of the points
        :param num_points: list of points
        :return: list of points
        """
        for k in range(len(num_points) - 1, 0, -1):
            r = random.randrange(0, k)
            temp = num_points[k]
            num_points[k] = num_points[r]
            num_points[r] = temp
        return num_points

    @staticmethod
    def minimum_circle_points(num_points, q1, q2):
        """
        Compute the minimum circle that contains the points in num_points and the points q1 and q2
        :param num_points: list of points
        :param q1: point
        :param q2: point
        :return: circle
        """
        circle = [[0.0, 0.0], 0.0]
        circle[0][0] = (q1[0] + q2[0]) / 2
        circle[0][1] = (q1[1] + q2[1]) / 2
        circle[1] = (sqrt((q2[0] - q1[0]) ** 2 + (q2[1] - q1[1]) ** 2)) / 2
        for i in range(0, len(num_points)):
            d = sqrt((circle[0][0] - num_points[i][0]) ** 2 + (circle[0][1] - num_points[i][1]) ** 2)
            if d > circle[1]:
                D = 2 * ((q2[0] - q1[0]) * (num_points[i][1] - q1[1]) - (num_points[i][0] - q1[0]) * (q2[1] - q1[1]))
                circleX = ((num_points[i][1] - q1[1]) * (q2[0] ** 2 + q2[1] ** 2 - q1[0] ** 2 - q1[1] ** 2) - (q2[1] - q1[1]) * (
                            num_points[i][0] ** 2 + num_points[i][1] ** 2 - q1[0] ** 2 - q1[1] ** 2)) / D
                circleY = ((q2[0] - q1[0]) * (num_points[i][0] ** 2 + num_points[i][1] ** 2 - q1[0] ** 2 - q1[1] ** 2) - (num_points[i][0] - q1[0]) * (
                            q2[0] ** 2 + q2[1] ** 2 - q1[0] ** 2 - q1[1] ** 2)) / D
                circle[0] = [circleX, circleY]
                circle[1] = sqrt((q1[0] - circle[0][0]) ** 2 + (q1[1] - circle[0][1]) ** 2)
        return circle

    @staticmethod
    def min_circle_with_point(num_points, q):
        """
        Compute the minimum circle that contains the points in num_points and the point q
        :param num_points: list of points
        :param q: point
        :return: circle
        """
        circle = [[0.0, 0.0], 0.0]
        circle[0][0] = (num_points[0][0] + q[0]) / 2
        circle[0][1] = (num_points[0][1] + q[1]) / 2
        circle[1] = (sqrt((num_points[0][0] - q[0]) ** 2 + (num_points[0][1] - q[1]) ** 2)) / 2
        for i in range(1, len(num_points)):
            d = sqrt((circle[0][0] - num_points[i][0]) ** 2 + (circle[0][1] - num_points[i][1]) ** 2)
            if d > circle[1]:
                circle = MinimumCircle.minimum_circle_points(num_points[0:i], num_points[i], q)
        return circle

    def minimum_circle(self, num_points):
        """
        Compute the minimum circle that contains all the points in num_points
        :param num_points: list of points
        :return: circle
        """
        from cgeom.elements.models import MinimumCircleInput
        validated = MinimumCircleInput(points=num_points)
        num_points = [p[:] for p in validated.points]
        random_perm = self.randomized_permutation(num_points)
        circle = [[0.0, 0.0], 0.0]
        circle[0][0] = (random_perm[0][0] + random_perm[1][0]) / 2
        circle[0][1] = (random_perm[0][1] + random_perm[1][1]) / 2
        circle[1] = (sqrt((random_perm[1][0] - random_perm[0][0]) ** 2 + (random_perm[1][1] - random_perm[0][1]) ** 2)) / 2
        for i in range(2, len(random_perm)):
            d = sqrt((circle[0][0] - random_perm[i][0]) ** 2 + (circle[0][1] - random_perm[i][1]) ** 2)
            if d > circle[1]:
                circle = self.min_circle_with_point(random_perm[0:i], random_perm[i])
        return circle

    def plot_min_circle_random(self, sizes, path=None, show=False):
        """Deprecated: use cgeom.visualization.plot_min_circle_random() instead."""
        import warnings
        warnings.warn(
            "MinimumCircle.plot_min_circle_random() is deprecated. "
            "Use cgeom.visualization.plot_min_circle_random(mc_obj, sizes, path, show) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        from cgeom.visualization import plot_min_circle_random
        plot_min_circle_random(self, sizes, path, show)

    def plot_min_circle(self, data, path=None, show=False):
        """Deprecated: use cgeom.visualization.plot_min_circle() instead."""
        import warnings
        warnings.warn(
            "MinimumCircle.plot_min_circle() is deprecated. "
            "Use cgeom.visualization.plot_min_circle(mc_obj, data, path, show) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        from cgeom.visualization import plot_min_circle
        plot_min_circle(self, data, path, show)
