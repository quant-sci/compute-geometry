import numpy as np
import matplotlib.pyplot as plt
from cgeom import Point, Line, Segment, Circle, Polygon
from cgeom.algorithms import MinimumCircle, ConvexHull

# --- Point ---
# Multiple construction forms
p1 = Point(x=1.0, y=2.0)       # keyword
p2 = Point([3, 4])              # from list
p3 = Point((5, 6))              # from tuple
p4 = Point(np.array([7, 8]))    # from numpy array

# Destructuring and indexing
x, y = p1
print(f"p1 destructured: x={x}, y={y}")
print(f"p1[0]={p1[0]}, p1[1]={p1[1]}")

# Distance and conversions
print(f"Distance p1->p2: {p1.distance_to(p2):.4f}")
print(f"p3 as numpy: {p3.to_numpy()}")
print(f"p4 as list:  {p4.to_list()}")

# --- Line ---
# Construction from raw lists
line = Line([[0, 0], [2, 4]])
print(f"\nLine slope: {line.slope}")
print(f"Line y-intercept: {line.y_intercept}")
print(f"Line coefficients (a,b,c): {line.coefficients}")
print(f"Line contains (1, 2)? {line.contains_point(Point(1, 2))}")

# Vertical line
vertical = Line([[3, 0], [3, 5]])
print(f"Vertical line slope: {vertical.slope}")       # None
print(f"Vertical y-intercept: {vertical.y_intercept}") # None

# --- Segment ---
seg = Segment([[0, 0], [3, 4]])
print(f"\nSegment length: {seg.length}")
print(f"Segment midpoint: {seg.midpoint}")
print(f"Segment as list: {seg.to_list()}")

# --- Circle ---
# Keyword construction
c1 = Circle(center=Point(0, 0), radius=5.0)
print(f"\nCircle area: {c1.area:.4f}")
print(f"Circle circumference: {c1.circumference:.4f}")
print(f"Circle contains (3, 4)? {c1.contains_point(Point(3, 4))}")
print(f"Circle contains (10, 0)? {c1.contains_point(Point(10, 0))}")

# Construction from [[cx, cy], r]
c2 = Circle([[1, 2], 3])
print(f"c2 center: {c2.center}, radius: {c2.radius}")

# --- Polygon ---
# From list of lists
triangle = Polygon([[0, 0], [4, 0], [0, 3]])
print(f"\nTriangle area: {triangle.area}")
print(f"Triangle perimeter: {triangle.perimeter:.4f}")
print(f"Triangle vertices: {triangle.num_vertices}")

# Iteration
for i, vertex in enumerate(triangle):
    print(f"  vertex {i}: {vertex}")

# Conversion
print(f"Triangle as numpy:\n{triangle.to_numpy()}")

# --- Algorithm interop ---
# Circle from MinimumCircle output
points = [(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
mc = MinimumCircle()
raw_circle = mc.minimum_circle(points)  # returns [[cx, cy], radius]
circle = Circle(raw_circle)
print(f"\nMinimumCircle -> Circle: {circle}")

# Polygon from ConvexHull output
hull_points = [(326, 237), (373, 209), (378, 265), (443, 241),
               (396, 231), (416, 270), (361, 335), (324, 297)]
ch = ConvexHull(hull_points)
hull_vertices = ch.convex_hull()  # returns [[x, y], ...]
hull_polygon = Polygon(hull_vertices)
print(f"ConvexHull -> Polygon: {hull_polygon.num_vertices} vertices, area={hull_polygon.area:.1f}")

# --- Grid plot of all elements ---
# Match project palette (cgeom.visualization._plotting)
_CMAP = plt.cm.Greys
_C0 = _CMAP(0.85)
_C1 = _CMAP(0.60)
_C2 = _CMAP(0.40)
_C3 = _CMAP(0.20)
_GRID = "#e8e8e8"

def _style_ax(ax, title=None):
    if title:
        ax.set_title(title, fontsize=12, fontweight="500", pad=10)
    ax.set_facecolor("white")
    ax.grid(False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_color(_GRID)
        ax.spines[sp].set_linewidth(0.7)
    ax.tick_params(colors="#555555", labelsize=8, length=3, width=0.6)

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
fig.patch.set_facecolor("white")

# 1) Points
ax = axes[0, 0]
_style_ax(ax, "Point")
for p, label in [(p1, "keywords"), (p2, "list"), (p3, "tuple"), (p4, "numpy")]:
    ax.scatter(p.x, p.y, c=[_C0], s=50, zorder=3, linewidths=0)
    ax.annotate(label, (p.x, p.y), textcoords="offset points",
                xytext=(6, 6), fontsize=8, color="#555555")
# show distance between p1 and p2
ax.plot([p1.x, p2.x], [p1.y, p2.y], color=_C2, linewidth=1, linestyle="--", zorder=2)
ax.annotate(f"d={p1.distance_to(p2):.2f}", xy=((p1.x+p2.x)/2, (p1.y+p2.y)/2),
            textcoords="offset points", xytext=(6, -10), fontsize=8, color="#555555")

# 2) Lines
ax = axes[0, 1]
_style_ax(ax, "Line")
ax.set_xlim(-1, 5)
ax.set_ylim(-2, 8)
# regular line y=2x through (0,0) and (2,4)
lx = np.array([-0.5, 4.0])
ax.plot(lx, line.slope * lx + line.y_intercept, color=_C0, linewidth=1.6, zorder=2)
ax.scatter([line.point1.x, line.point2.x], [line.point1.y, line.point2.y],
           c=[_C1], s=44, zorder=4, edgecolors="white", linewidths=0.5)
ax.annotate(f"slope={line.slope}", xy=(2.5, 5), fontsize=8, color="#555555")
# vertical line at x=3
ax.axvline(x=3, color=_C2, linewidth=1.6, linestyle="--", zorder=2)
ax.annotate("slope=None", xy=(3.1, 6.5), fontsize=8, color="#555555")

# 3) Segment
ax = axes[0, 2]
_style_ax(ax, "Segment")
ax.set_aspect(1)
ax.plot([seg.start.x, seg.end.x], [seg.start.y, seg.end.y],
        color=_C0, linewidth=1.8, zorder=2, solid_capstyle="round")
ax.scatter([seg.start.x, seg.end.x], [seg.start.y, seg.end.y],
           c=[_C1], s=44, zorder=4, edgecolors="white", linewidths=0.5)
mid = seg.midpoint
ax.scatter(mid.x, mid.y, c=[_C2], s=36, zorder=4, marker="D", edgecolors="white", linewidths=0.5)
ax.annotate(f"len={seg.length:.1f}", xy=((seg.start.x+seg.end.x)/2, (seg.start.y+seg.end.y)/2),
            textcoords="offset points", xytext=(8, -8), fontsize=8, color="#555555")
ax.annotate("midpoint", xy=(mid.x, mid.y), textcoords="offset points",
            xytext=(8, 6), fontsize=8, color="#555555")

# 4) Circle
ax = axes[1, 0]
_style_ax(ax, "Circle")
ax.set_aspect(1)
circle_patch = plt.Circle((c1.center.x, c1.center.y), c1.radius,
                           facecolor=(*_C2[:3], 0.10), edgecolor=_C1,
                           linewidth=1.4, zorder=2)
ax.add_patch(circle_patch)
ax.plot(c1.center.x, c1.center.y, "+", color=_C3,
        markersize=8, markeredgewidth=1.5, zorder=4)
# point inside
pin = Point(3, 4)
pout = Point(10, 0)
ax.scatter(pin.x, pin.y, c=[_C0], s=36, zorder=4, linewidths=0)
ax.annotate("inside", (pin.x, pin.y), textcoords="offset points",
            xytext=(6, 6), fontsize=8, color="#555555")
margin = c1.radius * 0.2
ax.set_xlim(c1.center.x - c1.radius - margin, c1.center.x + c1.radius + margin)
ax.set_ylim(c1.center.y - c1.radius - margin, c1.center.y + c1.radius + margin)

# 5) Polygon
ax = axes[1, 1]
_style_ax(ax, "Polygon")
ax.set_aspect(1)
verts = triangle.to_list()
xs = [v[0] for v in verts] + [verts[0][0]]
ys = [v[1] for v in verts] + [verts[0][1]]
ax.fill(xs, ys, facecolor=(*_C2[:3], 0.12), edgecolor=_C1, linewidth=1.6, zorder=2)
ax.scatter([v[0] for v in verts], [v[1] for v in verts],
           c=[_C0], s=44, zorder=4, edgecolors="white", linewidths=0.5)
ax.annotate(f"area={triangle.area:.1f}", xy=(1.0, 1.0), fontsize=8, color="#555555")

# 6) Algorithm interop — ConvexHull polygon + MinimumCircle
ax = axes[1, 2]
_style_ax(ax, "Algorithm interop")
ax.set_aspect(1)
# convex hull polygon
hv = hull_polygon.to_list()
hxs = [v[0] for v in hv] + [hv[0][0]]
hys = [v[1] for v in hv] + [hv[0][1]]
ax.fill(hxs, hys, facecolor=(*_C2[:3], 0.12), edgecolor=_C1, linewidth=1.4, zorder=2)
# all points
all_pts = np.array(hull_points)
ax.scatter(all_pts[:, 0], all_pts[:, 1], c=[_C0], s=30, zorder=3, linewidths=0)
# hull vertices
ax.scatter([v[0] for v in hv], [v[1] for v in hv],
           c=[_C1], s=44, zorder=4, edgecolors="white", linewidths=0.5)
# minimum enclosing circle from the interop section
mc2 = MinimumCircle()
raw2 = mc2.minimum_circle(hull_points)
circ2 = Circle(raw2)
cp = plt.Circle((circ2.center.x, circ2.center.y), circ2.radius,
                facecolor="none", edgecolor=_C2, linewidth=1.2,
                linestyle="--", zorder=1)
ax.add_patch(cp)
ax.plot(circ2.center.x, circ2.center.y, "+", color=_C3,
        markersize=8, markeredgewidth=1.5, zorder=4)
pad = circ2.radius * 0.15
ax.set_xlim(circ2.center.x - circ2.radius - pad, circ2.center.x + circ2.radius + pad)
ax.set_ylim(circ2.center.y - circ2.radius - pad, circ2.center.y + circ2.radius + pad)

fig.suptitle("Geometric Primitives", fontsize=14, fontweight="600", y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
