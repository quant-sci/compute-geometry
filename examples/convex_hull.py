from cgeom.algorithms import ConvexHull
from cgeom.visualization import plot_convex_hull

# create a list of points
points = [(326, 237),(373, 209), (378, 265), (443, 241), (396, 231), (416, 270), (361, 335), (324, 297)]

# create a convex hull object with the list of points
convex_hull = ConvexHull(points)

# plot the convex hull
plot_convex_hull(convex_hull)

# print the indexes of the points that form the convex hull
print("Convex Hull vertex indices:", convex_hull.get_indexes())