# polygon with scipy
from scipy.spatial import ConvexHull
import numpy as np

# Define the vertices of the polygon as a list of (x, y) coordinates
#vertices = [(1, 2), (3, 1), (1, 6), (7, 1)]

# create 1000 random vertices
vertices = np.random.rand(1000, 2)


# Create a convex hull using the vertices
hull = ConvexHull(vertices)

# Get the list of hull vertices in the correct order
hull_vertices = hull.vertices.tolist()
hull_vertices.sort()

# plot the polygon
import matplotlib.pyplot as plt
x = [vertices[i][0] for i in hull_vertices]
x.append(vertices[hull_vertices[0]][0])
y = [vertices[i][1] for i in hull_vertices]
y.append(vertices[hull_vertices[0]][1])
plt.plot(x, y)
plt.show()

for i in hull_vertices:
    print(vertices[i][0], vertices[i][1])