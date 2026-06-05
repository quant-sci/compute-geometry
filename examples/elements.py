import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from cgeom import Circle, Line, Point, Polygon, Segment
from cgeom.algorithms import ConvexHull, MinimumCircle

# --- Point ---
# Multiple construction forms
p1 = Point(x=1.0, y=2.0)  # keyword
p2 = Point([3, 4])  # from list
p3 = Point((5, 6))  # from tuple
p4 = Point(np.array([7, 8]))  # from numpy array

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
print(f"Vertical line slope: {vertical.slope}")  # None
print(f"Vertical y-intercept: {vertical.y_intercept}")  # None

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
hull_points = [
    (326, 237),
    (373, 209),
    (378, 265),
    (443, 241),
    (396, 231),
    (416, 270),
    (361, 335),
    (324, 297),
]
ch = ConvexHull(hull_points)
hull_vertices = ch.convex_hull()  # returns [[x, y], ...]
hull_polygon = Polygon(hull_vertices)
print(
    f"ConvexHull -> Polygon: {hull_polygon.num_vertices} vertices, area={hull_polygon.area:.1f}"
)

# ---------------------------------------------------------------------------
# Grid plot of all elements
# Mirrors the project palette (cgeom.visualization._plotting)
# ---------------------------------------------------------------------------
_INK = "#002855"  # Prussian Blue — primary structure, edges
_ACCENT = "#0466c8"  # Smart Blue — focal data points, emphasis
_CHARCOAL = "#023e7d"  # Regal Navy — secondary edges, lines
_STEEL = "#5c677d"  # Blue Slate — labels
_SILVER = "#7d8597"  # Slate Grey — annotations, dashes
_ASH = "#979dac"  # Lavender Grey — pill borders, ticks
_MIST = "#e6e9f1"  # faint lavender — grid lines
_FILL = "#eef1f6"  # panel/shape fill
_BG = "#ffffff"  # figure background

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [
            "Inter",
            "Helvetica Neue",
            "Helvetica",
            "Arial",
            "DejaVu Sans",
        ],
        "axes.unicode_minus": False,
        "figure.dpi": 150,
        "savefig.dpi": 220,
    }
)


def _style_ax(ax, index, title):
    """Apply a clean, modern panel look to an axes object."""
    # small index "eyebrow" above a left-aligned lightweight title
    ax.text(
        0.0,
        1.20,
        f"{index:02d}",
        transform=ax.transAxes,
        fontsize=8,
        fontweight="700",
        color=_SILVER,
        ha="left",
        va="bottom",
    )
    ax.set_title(
        f"{title}", loc="left", fontsize=12.5, fontweight="500", color=_INK, pad=10
    )

    ax.set_facecolor(_BG)
    ax.set_axisbelow(True)
    ax.grid(True, color=_MIST, linewidth=0.7, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(colors=_STEEL, labelsize=7, length=0, width=0, pad=4)


def _pill(
    ax, x, y, text, dx=8, dy=8, ha="left", va="center", color=_STEEL, weight="500"
):
    """Draw a small rounded label badge anchored to a data point."""
    ax.annotate(
        text,
        (x, y),
        textcoords="offset points",
        xytext=(dx, dy),
        fontsize=7.5,
        color=color,
        ha=ha,
        va=va,
        fontweight=weight,
        zorder=6,
        bbox=dict(boxstyle="round,pad=0.32", fc="white", ec=_ASH, lw=0.7),
    )


fig, axes = plt.subplots(2, 3, figsize=(15, 9.2))
fig.patch.set_facecolor(_BG)

# 1) Points -----------------------------------------------------------------
ax = axes[0, 0]
_style_ax(ax, 1, "Point")
# distance connector first (behind points)
ax.plot(
    [p1.x, p2.x],
    [p1.y, p2.y],
    color=_SILVER,
    linewidth=1.0,
    linestyle=(0, (3, 3)),
    zorder=2,
)
for p, label in [(p1, "keyword"), (p2, "list"), (p3, "tuple"), (p4, "numpy")]:
    ax.scatter(p.x, p.y, c=_ACCENT, s=46, zorder=4, edgecolors="white", linewidths=0.9)
    _pill(ax, p.x, p.y, label, dx=9, dy=9, ha="left")
_pill(
    ax,
    (p1.x + p2.x) / 2,
    (p1.y + p2.y) / 2,
    f"d = {p1.distance_to(p2):.2f}",
    dx=10,
    dy=-12,
    color=_CHARCOAL,
)
ax.margins(0.18)

# 2) Line -------------------------------------------------------------------
ax = axes[0, 1]
_style_ax(ax, 2, "Line")
ax.set_xlim(-1, 5)
ax.set_ylim(-2, 8)
lx = np.array([-0.6, 4.2])
ax.plot(
    lx,
    line.slope * lx + line.y_intercept,
    color=_INK,
    linewidth=1.8,
    zorder=3,
    solid_capstyle="round",
)
ax.scatter(
    [line.point1.x, line.point2.x],
    [line.point1.y, line.point2.y],
    c=_ACCENT,
    s=42,
    zorder=5,
    edgecolors="white",
    linewidths=0.9,
)
ax.axvline(x=3, color=_SILVER, linewidth=1.6, linestyle=(0, (4, 4)), zorder=2)
_pill(ax, 2.0, 4.0, f"slope = {line.slope:.0f}", dx=10, dy=-2, color=_CHARCOAL)
_pill(ax, 3.0, 6.6, "slope = None", dx=10, dy=0, color=_STEEL)

# 3) Segment ----------------------------------------------------------------
ax = axes[0, 2]
_style_ax(ax, 3, "Segment")
ax.set_aspect("equal")
ax.plot(
    [seg.start.x, seg.end.x],
    [seg.start.y, seg.end.y],
    color=_INK,
    linewidth=2.0,
    zorder=3,
    solid_capstyle="round",
)
ax.scatter(
    [seg.start.x, seg.end.x],
    [seg.start.y, seg.end.y],
    c=_ACCENT,
    s=42,
    zorder=5,
    edgecolors="white",
    linewidths=0.9,
)
mid = seg.midpoint
ax.scatter(
    mid.x, mid.y, c="white", s=46, zorder=5, marker="o", edgecolors=_INK, linewidths=1.4
)
_pill(ax, mid.x, mid.y, "midpoint", dx=11, dy=10, color=_CHARCOAL)
_pill(
    ax,
    (seg.start.x + seg.end.x) / 2,
    (seg.start.y + seg.end.y) / 2,
    f"len = {seg.length:.1f}",
    dx=12,
    dy=-14,
    color=_STEEL,
)
ax.margins(0.18)

# 4) Circle -----------------------------------------------------------------
ax = axes[1, 0]
_style_ax(ax, 4, "Circle")
ax.set_aspect("equal")
ax.add_patch(
    plt.Circle(
        (c1.center.x, c1.center.y),
        c1.radius,
        facecolor=_FILL,
        edgecolor=_CHARCOAL,
        linewidth=1.6,
        zorder=2,
    )
)
# radius indicator
edge = Point(
    c1.center.x + c1.radius * np.cos(np.pi / 5),
    c1.center.y + c1.radius * np.sin(np.pi / 5),
)
ax.plot(
    [c1.center.x, edge.x],
    [c1.center.y, edge.y],
    color=_SILVER,
    linewidth=1.1,
    linestyle=(0, (3, 3)),
    zorder=3,
)
ax.plot(
    c1.center.x,
    c1.center.y,
    "+",
    color=_CHARCOAL,
    markersize=9,
    markeredgewidth=1.6,
    zorder=4,
)
pin = Point(3, 4)
ax.scatter(pin.x, pin.y, c=_ACCENT, s=46, zorder=5, edgecolors="white", linewidths=0.9)
_pill(ax, pin.x, pin.y, "contains", dx=9, dy=9, color=_CHARCOAL)
_pill(
    ax,
    (c1.center.x + edge.x) / 2,
    (c1.center.y + edge.y) / 2,
    f"r = {c1.radius:.0f}",
    dx=6,
    dy=8,
    color=_STEEL,
)
margin = c1.radius * 0.22
ax.set_xlim(c1.center.x - c1.radius - margin, c1.center.x + c1.radius + margin)
ax.set_ylim(c1.center.y - c1.radius - margin, c1.center.y + c1.radius + margin)

# 5) Polygon ----------------------------------------------------------------
ax = axes[1, 1]
_style_ax(ax, 5, "Polygon")
ax.set_aspect("equal")
verts = triangle.to_list()
xs = [v[0] for v in verts] + [verts[0][0]]
ys = [v[1] for v in verts] + [verts[0][1]]
ax.fill(xs, ys, facecolor=_FILL, edgecolor=_INK, linewidth=1.8, zorder=2)
ax.scatter(
    [v[0] for v in verts],
    [v[1] for v in verts],
    c=_ACCENT,
    s=42,
    zorder=5,
    edgecolors="white",
    linewidths=0.9,
)
cx = sum(v[0] for v in verts) / len(verts)
cy = sum(v[1] for v in verts) / len(verts)
_pill(
    ax, cx, cy, f"area = {triangle.area:.1f}", dx=0, dy=0, ha="center", color=_CHARCOAL
)
ax.margins(0.18)

# 6) Algorithm interop — ConvexHull polygon + MinimumCircle -----------------
ax = axes[1, 2]
_style_ax(ax, 6, "Algorithm interop")
ax.set_aspect("equal")
hv = hull_polygon.to_list()
hxs = [v[0] for v in hv] + [hv[0][0]]
hys = [v[1] for v in hv] + [hv[0][1]]
mc2 = MinimumCircle()
circ2 = Circle(mc2.minimum_circle(hull_points))
ax.add_patch(
    plt.Circle(
        (circ2.center.x, circ2.center.y),
        circ2.radius,
        facecolor="none",
        edgecolor=_SILVER,
        linewidth=1.3,
        linestyle=(0, (4, 4)),
        zorder=1,
    )
)
ax.fill(hxs, hys, facecolor=_FILL, edgecolor=_CHARCOAL, linewidth=1.6, zorder=2)
all_pts = np.array(hull_points)
ax.scatter(
    all_pts[:, 0],
    all_pts[:, 1],
    c=_STEEL,
    s=26,
    zorder=3,
    edgecolors="white",
    linewidths=0.7,
)
ax.scatter(
    [v[0] for v in hv],
    [v[1] for v in hv],
    c=_ACCENT,
    s=42,
    zorder=5,
    edgecolors="white",
    linewidths=0.9,
)
ax.plot(
    circ2.center.x,
    circ2.center.y,
    "+",
    color=_CHARCOAL,
    markersize=9,
    markeredgewidth=1.6,
    zorder=4,
)
_pill(
    ax,
    circ2.center.x,
    circ2.center.y + circ2.radius,
    "min circle",
    dx=0,
    dy=10,
    ha="center",
    color=_STEEL,
)
_pill(ax, hv[0][0], hv[0][1], "hull", dx=-6, dy=-12, ha="right", color=_CHARCOAL)
pad = circ2.radius * 0.18
ax.set_xlim(circ2.center.x - circ2.radius - pad, circ2.center.x + circ2.radius + pad)
ax.set_ylim(circ2.center.y - circ2.radius - pad, circ2.center.y + circ2.radius + pad)

# Header --------------------------------------------------------------------
fig.text(
    0.5,
    0.975,
    "Geometric Primitives",
    ha="center",
    va="top",
    fontsize=17,
    fontweight="600",
    color=_INK,
)
fig.text(
    0.5,
    0.938,
    "Composable building blocks of the compute-geometry library",
    ha="center",
    va="top",
    fontsize=9.5,
    color=_STEEL,
)

plt.tight_layout(rect=[0.015, 0.02, 0.985, 0.90], h_pad=3.2, w_pad=2.6)
fig.savefig(
    "public/elements.png", dpi=220, facecolor=_BG, bbox_inches="tight", pad_inches=0.25
)
print("\nSaved public/elements.png")
plt.show()
