import random
import numpy as np
import pandas as pd
import time
from math import sqrt

import matplotlib.pyplot as plt
import seaborn as sns

plt.rc('text', usetex=True)
plt.rc('font', size=12)
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5})

circle_color = "lightslategray"
points_color = "navy"

class MinimumCircle:
    def __init__(self):
        pass
    
    @staticmethod
    def minimum_circle_heuristic(num_points):
        """
        Heuristic algorithm to find the minimum circle that contains all the points
        :param num_points: list of points
        :return: circle
        """
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

    def plot_min_circle_random(self, sizes, path=None, show = False):

        time_heuristic = []
        time_min_circle = []
        for num_points in sizes:
            x = np.random.standard_normal(num_points)
            y = np.random.standard_normal(num_points)

            points = []
            for i in range(num_points):
                points.append([x[i], y[i]])

            start_time = time.time()
            min_circle = self.minimum_circle(points)
            end_time = time.time()
            time_min_circle.append(end_time - start_time)

            start_time = time.time()
            min_circle_heuristic = self.minimum_circle_heuristic(points)
            end_time = time.time()
            time_heuristic.append(end_time - start_time)

            print("MinCircle Center/Radius: " + str(min_circle[0]) + " / " + str(min_circle[1]))
            print("Heuristic Center/Radius: " + str(min_circle_heuristic[0]) + " / " + str(min_circle_heuristic[1]))

            min_circle_plt = plt.Circle(min_circle[0], min_circle[1], color=circle_color, clip_on=False, fill=False)
            min_circle_heuristic_plt = plt.Circle(min_circle_heuristic[0], min_circle_heuristic[1], color=circle_color, clip_on=False, fill=False)
            
            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
            ax[0].set_title('MinCircle, size = ' + str(num_points))
            ax[0].set_aspect(1)
            ax[0].scatter(x, y, color=points_color)
            ax[0].set_xlim((-10, 10))
            ax[0].set_ylim((-10, 10))
            ax[0].add_artist(min_circle_plt)
            ax[0].plot(min_circle[0][0], min_circle[0][1], '+', color='orangered')

            ax[1].set_title('Minimum Circle (Heuristic), size = ' + str(num_points))
            ax[1].set_aspect(1)
            ax[1].scatter(x, y, color=points_color)
            ax[1].set_xlim((-10, 10))
            ax[1].set_ylim((-10, 10))
            ax[1].add_artist(min_circle_heuristic_plt)
            ax[1].plot(min_circle[0][0], min_circle[0][1], '+', color="orangered")
            
            if path is not None:
                plt.savefig(path + f'plot_min_circle_{str(num_points)}_points' + ".pdf")
            if show is True:
                plt.show()
        fig = plt.figure(figsize=(5, 5))
        plt.plot(sizes, time_min_circle, "indianred")
        plt.plot(sizes, time_heuristic, "seagreen")
        plt.xlabel('Input size')
        plt.ylabel('time (s)')
        plt.legend(["Minimum Circle", "Minimum Circle (Heuristic)"])
        plt.grid(True)
        if path is not None:
            plt.savefig(path + "time.pdf")
        if show is True:
            plt.show()

    def plot_min_circle(self, data, path=None, show=False):
        time_heuristic = []
        time_min_circle = []

        start_time = time.time()
        min_circle = self.minimum_circle(data)
        end_time = time.time()
        time_min_circle.append(end_time - start_time)

        start_time = time.time()
        min_circle_heuristic = self.minimum_circle_heuristic(data)
        end_time = time.time()
        time_heuristic.append(end_time - start_time)
        
        print("MinCircle Center/Radius: " + str(min_circle[0]) + " / " + str(min_circle[1]))
        print("Heuristic Center/Radius: " + str(min_circle_heuristic[0]) + " / " + str(min_circle_heuristic[1]))

        min_circle_plt = plt.Circle(min_circle[0], min_circle[1], color=circle_color, clip_on=False, fill=False)
        min_circle_heuristic_plt = plt.Circle(min_circle_heuristic[0], min_circle_heuristic[1], color=circle_color, clip_on=False, fill=False)
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.set_title('MinCircle, points.txt')
        ax1.set_aspect(1)
        ax1.scatter(data[:, 0], data[:, 1], color=points_color)
        ax1.set_xlim((0, 800))
        ax1.set_ylim((0, 800))
        ax1.add_artist(min_circle_plt)
        ax1.plot(min_circle[0][0], min_circle[0][1], '+', color="orangered")

        ax2.set_title('Heuristc, points.txt')
        ax2.set_aspect(1)
        ax2.scatter(data[:, 0], data[:, 1], color=points_color)
        ax2.set_xlim((0, 800))
        ax2.set_ylim((0, 800))
        ax2.add_artist(min_circle_heuristic_plt)
        ax2.plot(min_circle[0][0], min_circle[0][1], '+', color="orangered")
        
        if path is not None:
            plt.savefig(path + f'plot_min_circle' + ".pdf")
        if show is True:
            plt.show()

        fig = plt.figure(figsize=(5, 5))
        
        times = pd.DataFrame({"time_min_circle": time_min_circle, "time_heuristic": time_heuristic})
        times.plot(kind = 'bar', color = ["indianred", "seagreen"])
        plt.legend(["Minimum Circle", "Minimum Circle (Heuristic)"])
        plt.grid(True)
        
        if path is not None:
            plt.savefig(path + "time.pdf")
        if show is True:
            plt.show()