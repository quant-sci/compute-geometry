from cgeom.algorithms import DelaunayTriangulation
from cgeom.visualization import plot_delaunay

# create a list of points
points = [
    (326, 237), (373, 209), (378, 265), (443, 241), (396, 231),
    (416, 270), (361, 335), (324, 297), (400, 306), (454, 315),
]

# create a Delaunay triangulation object
dt = DelaunayTriangulation(points)

# compute the triangulation
triangles = dt.triangulate()
print("Triangles (index triples):", triangles)

# get the edges
edges = dt.get_edges()
print("Number of edges:", len(edges))

# get circumcircles
circumcircles = dt.get_circumcircles()
print("Number of circumcircles:", len(circumcircles))

# plot the triangulation
plot_delaunay(dt)

# plot with circumcircles
plot_delaunay(dt, title="Delaunay with Circumcircles", show_circumcircles=True)
