from cgeom.algorithms import SegmentIntersection
from cgeom.visualization import plot_intersections

# create a list of segments: X cross + horizontal + vertical
segments = [
    [[0, 0], [4, 4]],   # diagonal /
    [[0, 4], [4, 0]],   # diagonal \
    [[0, 2], [4, 2]],   # horizontal
    [[2, 0], [2, 4]],   # vertical
]

# create a SegmentIntersection object
si = SegmentIntersection(segments)

# find intersections using Bentley-Ottmann sweep line
intersections = si.find_intersections()
print("Intersection points:", intersections)
print("Number of intersections:", len(intersections))

# find intersections using brute force (for verification)
bf_intersections = si.find_intersections_brute_force()
print("Brute force intersections:", bf_intersections)

# get detailed intersection pairs (which segments intersect)
pairs = si.get_intersection_pairs()
for seg_i, seg_j, pt in pairs:
    print(f"  Segment {seg_i} x Segment {seg_j} at ({pt[0]:.2f}, {pt[1]:.2f})")

# plot the result
plot_intersections(si)
