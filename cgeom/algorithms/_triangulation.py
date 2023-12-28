import numpy as np
import matplotlib.pyplot as plt
import math
import random

class PolygonTriangulation:
    def __init__(self, poly, poly_name="Polygon"):
        self.poly = poly
        self.poly_name = poly_name
        self.diagonals = self.Triangulation()

    def is_ear(self, poly, vertex):  # Function to define if a given vertex is an ear or not

        # First, we start testing if the line lays inside the polygon by analysing the triangle made by vertex, vertex+1, vertex-1

        p1 = poly[vertex]
        if vertex == len(poly) - 1: # If the vertex is in the last position, vertex + 1 is the first vertex
            p2 = poly[0]
        else:
            p2 = poly[vertex + 1]
        p3 = poly[vertex - 1]

        v1 = [p2[0] - p1[0], p2[1] - p1[1]]
        v2 = [p3[0] - p2[0], p3[1] - p2[1]]
        v3 = [p1[0] - p3[0], p1[1] - p3[1]]

        # Now we need to find the triangle's orientation

        orient = v1[0] * v2[1] + v3[0] * v1[1] + v2[0] * v3[1] - v3[0] * v2[1] - v3[1] * v1[0] - v2[0] * v1[1]

        if orient <= 0:
            return False
        else:

            # Finally, we need to check if any other vertex is inside the triangle made by vertex-1, vertex and vertex+1

            for i, p in enumerate(poly):
                D = p1[0] * p2[1] + p2[0] * p3[1] + p3[0] * p1[1] - p2[1] * p3[0] - p1[0] * p3[1] - p1[1] * p2[0]
                D1 = p[0] * p2[1] + p2[0] * p3[1] + p3[0] * p[1] - p2[1] * p3[0] - p3[1] * p[0] - p[1] * p2[0]
                D2 = p1[0] * p[1] + p[0] * p3[1] + p3[0] * p1[1] - p[1] * p3[0] - p3[1] * p1[0] - p1[1] * p[0]
                D3 = p1[0] * p2[1] + p2[0] * p[1] + p[0] * p1[1] - p2[1] * p[0] - p[1] * p1[0] - p1[1] * p2[0]
                lambda1 = D1 / D
                lambda2 = D2 / D
                lambda3 = D3 / D
                if lambda1 > 0 and lambda2 > 0 and lambda3 > 0: return False
            return True

    def Triangulation(self):    # Main function to triangulate the polygon
        poly = self.poly
        diag = []   # List of diagonals
        vertex = 0
        while len(poly) > 3:    # The loop runs until only one triangle is left
            if vertex >= len(poly):     # When the loop reaches the last vertex
                vertex = 0
            if self.is_ear(poly, vertex):   # Test if the vertex is an ear
                if vertex == len(poly) - 1:     # If the vertex is in the last position, vertex + 1 is the first vertex
                    diag.append([list(poly[vertex - 1]), list(poly[0])])
                else:
                    diag.append([list(poly[vertex - 1]), list(poly[vertex + 1])])
                poly = np.delete(poly, vertex, 0)   # Remove the vertex
            vertex += 1
        return diag

    def get_diag_vertexes(self):
        diag_vertexes = []
        for i, diagonal in enumerate(self.diagonals):
            for vertex in diagonal:
                index = 0
                while vertex[0] != self.poly[index][0] or vertex[1] != self.poly[index][1]:
                    index += 1
                diag_vertexes.append(index)
            diag_vertexes[i] = [diag_vertexes[-2],diag_vertexes[-1]]
            del diag_vertexes[-1]
        return(diag_vertexes)

    def plot_triangulation(self):
        fig, ax = plt.subplots()
        ax.set_aspect(1)
        ax.set_title(self.poly_name)
        x, y = zip(*np.append(self.poly, [self.poly[0]], axis=0))  # Adding the first vertex at the end, for ploting purposes
        ax.plot(x, y, color = "blue")      # Plotting Polygon
        for diag in self.diagonals:     # Plotting triangulation
            ax.plot([diag[0][0], diag[1][0]], [diag[0][1], diag[1][1]], color = "blue")
        plt.show()

#   Function to generate random polygons

def generatePolygon( ctrX, ctrY, aveRadius, irregularity, spikeyness, numVerts ):

    irregularity = clip( irregularity, 0,1 ) * 2*math.pi / numVerts
    spikeyness = clip( spikeyness, 0,1 ) * aveRadius

    # generate n angle steps
    angleSteps = []
    lower = (2*math.pi / numVerts) - irregularity
    upper = (2*math.pi / numVerts) + irregularity
    sum = 0
    for i in range(numVerts) :
        tmp = random.uniform(lower, upper)
        angleSteps.append( tmp )
        sum = sum + tmp

    # normalize the steps so that point 0 and point n+1 are the same
    k = sum / (2*math.pi)
    for i in range(numVerts) :
        angleSteps[i] = angleSteps[i] / k

    # now generate the points
    points = []
    angle = random.uniform(0, 2*math.pi)
    for i in range(numVerts) :
        r_i = clip( random.gauss(aveRadius, spikeyness), 0, 2*aveRadius )
        x = ctrX + r_i*math.cos(angle)
        y = ctrY + r_i*math.sin(angle)
        points.append( (int(x),int(y)) )

        angle = angle + angleSteps[i]

    return points

def clip(x, min, max) :
    if( min > max ) :  return x
    elif( x < min ) :  return min
    elif( x > max ) :  return max
    else :             return x