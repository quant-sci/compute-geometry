import matplotlib.pyplot as plt
import seaborn as sns
plt.rc('text', usetex=True)
plt.rc('font', size=12)
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.2, rc={"lines.linewidth": 2.5})

hull_color = "lightslategray"
points_color = "navy"

import numpy as np
import math
import time

np.random.seed(42)

class JarvisMarchGift:
    """Class to find the convex hull of a set of points using the Gift Wrapping algorithm""
    
    Reference: 
        R.A. Jarvis, On the identification of the convex hull of a finite set of points in the plane, Information Processing Letters,
        Volume 2, Issue 1, 1973, Pages 18-21, ISSN 0020-0190, https://doi.org/10.1016/0020-0190(73)90020-3.

    Attributes:
        points (np.array): Array of points
    """
    def __init__(self, points):
        self.points = np.array(points)

    def find_low_right(self):
        lower_y = float("inf")
        for point in self.points:
            if point[1] < lower_y:
                lower_y = point[1]
                lower = point
        return list(lower)

    def low_angle(self, b, c):
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
        a = self.find_low_right()
        b = self.low_angle(a, [1, a[1]])
        init = a
        ch = [a]
        while b[0] != init[0] or b[1] != init[1]:
            ch.append(b)
            a = b
            b = self.low_angle(a, ch[-2])
        return ch

    def plot_convex_hull(self, title):
        plt.set_title(title)
        plt.grid(True)
        plt.scatter(self.points[:, 0], self.points[:, 1])
        ch = self.convex_hull()
        for i in range(len(ch)):
            if i != len(ch) - 1:
                plt.plot([ch[i][0], ch[i + 1][0]], [ch[i][1], ch[i + 1][1]], color=hull_color)
            else:
                plt.plot([ch[i][0], ch[0][0]], [ch[i][1], ch[0][1]], color=hull_color)

    def get_indexes(self):
        indexes = []
        for element in self.convex_hull():
            j = 0
            while element[0] != self.points[j][0] or element[1] != self.points[j][1]:
                j += 1
            indexes.append(j)
        return(indexes)