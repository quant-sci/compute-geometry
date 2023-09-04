import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Define the vertices of the polygon as a list of (x, y) coordinates
vertices = [(1, 2), (3, 4), (5, 6), (7, 8)]

# Create a figure and axis
fig, ax = plt.subplots()

# Create a polygon patch and add it to the axis
polygon = patches.Polygon(vertices, closed=True, edgecolor='r', facecolor='none')
ax.add_patch(polygon)

# Set axis limits
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

# Show the plot
plt.gca().set_aspect('equal', adjustable='box')  # Equal aspect ratio
plt.show()