import numpy as np
from cgeom.algorithms import VoronoiDiagram
from cgeom.visualization import plot_voronoi

# load a set of points
points = np.loadtxt("examples/points1.txt")

# create a voronoi diagram object
voronoi = VoronoiDiagram(points)

# build the voronoi diagram
cells = voronoi.build_voronoi_diagram()

# plot the voronoi diagram
plot_voronoi(voronoi, cells)