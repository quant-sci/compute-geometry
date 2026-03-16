import numpy as np
import matplotlib.pyplot as plt
import math
import random

class PolygonTriangulation:
    """Triangulate a simple polygon using the ear-clipping algorithm.

    Attributes:
        poly (np.ndarray): Array of polygon vertices with shape (n, 2).
        poly_name (str): Display name for the polygon.
        diagonals (list): List of diagonals produced by the triangulation.
    """

    def __init__(self, poly, poly_name="Polygon"):
        """Initialize and triangulate a polygon.

        Args:
            poly: Polygon vertices as a list, tuple, or numpy array with shape (n, 2).
                Vertices must be in counter-clockwise order.
            poly_name: Display name used in plot titles.

        Raises:
            pydantic.ValidationError: If fewer than 3 vertices, all collinear,
                non-numeric, or wrong shape.
        """
        from cgeom.elements.models import PolygonTriangulationInput
        validated = PolygonTriangulationInput(poly=poly, poly_name=poly_name)
        self.poly = np.array(validated.poly)
        self.poly_name = validated.poly_name
        self.diagonals = self.Triangulation()

    def is_ear(self, poly, vertex):
        """Determine whether the given vertex is an ear of the polygon.

        A vertex is an ear if the triangle formed by it and its two neighbours
        is oriented counter-clockwise and contains no other polygon vertices.

        Args:
            poly: Current polygon vertex array.
            vertex: Index of the vertex to test.

        Returns:
            bool: True if the vertex is an ear.
        """
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

    def Triangulation(self):
        """Triangulate the polygon by iteratively clipping ears.

        Returns:
            list: List of diagonals, where each diagonal is a pair of [x, y] points.
        """
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
        """Get the vertex indices of each triangulation diagonal.

        Returns:
            list[list[int]]: Pairs of indices into ``self.poly`` for each diagonal.
        """
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
        """Plot the polygon and its triangulation diagonals."""
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
    """Generate a random simple polygon.

    Args:
        ctrX: X-coordinate of the polygon centre.
        ctrY: Y-coordinate of the polygon centre.
        aveRadius: Average radius from centre to vertices.
        irregularity: Parameter in [0, 1] controlling angular variance.
        spikeyness: Parameter in [0, 1] controlling radial variance.
        numVerts: Number of vertices.

    Returns:
        list[tuple[int, int]]: Vertices of the generated polygon.
    """
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
    """Clamp *x* to the range [min, max].

    Args:
        x: Value to clamp.
        min: Lower bound.
        max: Upper bound.

    Returns:
        The clamped value.
    """
    if( min > max ) :  return x
    elif( x < min ) :  return min
    elif( x > max ) :  return max
    else :             return x