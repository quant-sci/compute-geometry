import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ---------------------------------------------------------------------------
# Global rcParams — polished minimalist defaults
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Inter", "Helvetica Neue", "Helvetica",
                           "Arial", "DejaVu Sans"],
    "font.size":          9,
    "axes.unicode_minus": False,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.15,
})

# ---------------------------------------------------------------------------
# Grayscale palette — hand-picked hex grays
# ---------------------------------------------------------------------------
_INK      = "#1a1a1a"   # near-black — primary data (points, main edges)
_CHARCOAL = "#4a4a4a"   # dark gray — secondary edges, lines
_STEEL    = "#7a7a7a"   # medium gray — tertiary, tick labels, axis labels
_SILVER   = "#a8a8a8"   # light-medium — annotations, dashes, subtle markers
_ASH      = "#d4d4d4"   # light gray — spines, ticks (if ever needed)
_MIST     = "#ededed"   # very faint — fill base color
_BG       = "#fafafa"   # off-white — figure & axes background

# ---------------------------------------------------------------------------
# Visual sizing constants
# ---------------------------------------------------------------------------
_PT_SIZE      = 28      # primary scatter point size
_PT_ACCENT    = 40      # emphasized scatter point size
_PT_EDGE_W    = 0.8     # white edge width on scatter points
_LINE_W       = 1.0     # primary edge/line width
_LINE_W_THIN  = 0.7     # secondary/dashed line width
_FILL_ALPHA   = 0.06    # polygon/triangle fill opacity
_CIRCLE_ALPHA = 0.08    # circle fill opacity

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _hex_to_rgb(hex_color):
    """Convert ``'#rrggbb'`` to an ``(r, g, b)`` float tuple."""
    h = hex_color.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _style_ax(ax, title=None):
    """Apply a minimal look to an axes object."""
    if title:
        ax.set_title(title, fontsize=12, fontweight="300",
                     color=_INK, pad=14)
    ax.set_facecolor(_BG)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(colors=_STEEL, labelsize=7.5, length=0, width=0)


def _new_fig(figsize=(5.5, 5.5)):
    """Create a figure with an off-white background."""
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(_BG)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.90, bottom=0.10)
    return fig, ax


# ---------------------------------------------------------------------------
# Convex Hull
# ---------------------------------------------------------------------------

def plot_convex_hull(hull_obj, title="Convex Hull"):
    """Plot the points and their convex hull.

    Args:
        hull_obj: A ``ConvexHull`` instance.
        title: Title for the matplotlib figure.
    """
    fig, ax = _new_fig()
    _style_ax(ax, title)
    ax.set_aspect("equal")

    ch = hull_obj.convex_hull()
    xs = [p[0] for p in ch] + [ch[0][0]]
    ys = [p[1] for p in ch] + [ch[0][1]]
    ax.fill(xs, ys,
            facecolor=(*_hex_to_rgb(_MIST), _FILL_ALPHA),
            edgecolor=_CHARCOAL, linewidth=_LINE_W, zorder=2)

    ax.scatter(
        hull_obj.points[:, 0], hull_obj.points[:, 1],
        c=_STEEL, s=_PT_SIZE, zorder=3,
        edgecolors="white", linewidths=_PT_EDGE_W,
    )
    ax.scatter(
        [p[0] for p in ch], [p[1] for p in ch],
        c=_INK, s=_PT_ACCENT, zorder=4,
        edgecolors="white", linewidths=_PT_EDGE_W,
    )

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Minimum Circle (random sizes)
# ---------------------------------------------------------------------------

def plot_min_circle_random(mc_obj, sizes, path=None, show=False):
    """Plot minimum circles for randomly generated point sets of varying sizes.

    Args:
        mc_obj: A ``MinimumCircle`` instance.
        sizes: List of point-set sizes to generate and evaluate.
        path: Directory path to save PDF figures. None to skip saving.
        show: If True, display each figure interactively.
    """
    time_heuristic = []
    time_min_circle = []

    for num_points in sizes:
        x = np.random.standard_normal(num_points)
        y = np.random.standard_normal(num_points)

        points = [[x[i], y[i]] for i in range(num_points)]

        start_time = time.time()
        min_circle = mc_obj.minimum_circle(points)
        time_min_circle.append(time.time() - start_time)

        start_time = time.time()
        min_circle_heuristic = mc_obj.minimum_circle_heuristic(points)
        time_heuristic.append(time.time() - start_time)

        print(f"MinCircle  Center/Radius: {min_circle[0]} / {min_circle[1]}")
        print(f"Heuristic  Center/Radius: {min_circle_heuristic[0]} / {min_circle_heuristic[1]}")

        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        fig.patch.set_facecolor(_BG)

        for idx, (circ, label) in enumerate([
            (min_circle, f"Exact (n={num_points})"),
            (min_circle_heuristic, f"Heuristic (n={num_points})"),
        ]):
            ax = axes[idx]
            _style_ax(ax, label)
            ax.set_aspect("equal")
            ax.scatter(x, y, c=_INK, s=16, zorder=3,
                       edgecolors="white", linewidths=_PT_EDGE_W)
            circle_patch = plt.Circle(
                circ[0], circ[1],
                facecolor=(*_hex_to_rgb(_MIST), _CIRCLE_ALPHA),
                edgecolor=_CHARCOAL, linewidth=_LINE_W, zorder=1,
            )
            ax.add_patch(circle_patch)
            ax.plot(circ[0][0], circ[0][1], "+", color=_SILVER,
                    markersize=7, markeredgewidth=1.0, zorder=4)
            ax.set_xlim(-10, 10)
            ax.set_ylim(-10, 10)

        plt.tight_layout()
        if path is not None:
            fig.savefig(path + f"plot_min_circle_{num_points}_points.pdf",
                        bbox_inches="tight")
        if show:
            plt.show()

    # Runtime comparison
    fig, ax = _new_fig(figsize=(5.5, 4))
    _style_ax(ax, "Runtime")
    ax.plot(sizes, time_min_circle, color=_INK, linewidth=_LINE_W,
            marker="o", markersize=4, label="Exact")
    ax.plot(sizes, time_heuristic, color=_SILVER, linewidth=_LINE_W,
            marker="o", markersize=4, label="Heuristic")
    ax.set_xlabel("Input size", fontsize=8.5, color=_STEEL)
    ax.set_ylabel("Time (s)", fontsize=8.5, color=_STEEL)
    ax.legend(frameon=False, fontsize=9)
    plt.tight_layout()
    if path is not None:
        fig.savefig(path + "time.pdf", bbox_inches="tight")
    if show:
        plt.show()


# ---------------------------------------------------------------------------
# Minimum Circle (given data)
# ---------------------------------------------------------------------------

def plot_min_circle(mc_obj, data, path=None, show=False):
    """Plot minimum circles for a given dataset.

    Args:
        mc_obj: A ``MinimumCircle`` instance.
        data: Array-like of 2D points.
        path: Directory path to save PDF figures. None to skip saving.
        show: If True, display each figure interactively.
    """
    start_time = time.time()
    min_circle = mc_obj.minimum_circle(data)
    t_exact = time.time() - start_time

    start_time = time.time()
    min_circle_heuristic = mc_obj.minimum_circle_heuristic(data)
    t_heur = time.time() - start_time

    print(f"MinCircle  Center/Radius: {min_circle[0]} / {min_circle[1]}")
    print(f"Heuristic  Center/Radius: {min_circle_heuristic[0]} / {min_circle_heuristic[1]}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.patch.set_facecolor(_BG)

    for idx, (circ, label) in enumerate([
        (min_circle, "Exact"),
        (min_circle_heuristic, "Heuristic"),
    ]):
        ax = axes[idx]
        _style_ax(ax, label)
        ax.set_aspect("equal")
        ax.scatter(data[:, 0], data[:, 1], c=_INK, s=16,
                   zorder=3, edgecolors="white", linewidths=_PT_EDGE_W)
        circle_patch = plt.Circle(
            circ[0], circ[1],
            facecolor=(*_hex_to_rgb(_MIST), _CIRCLE_ALPHA),
            edgecolor=_CHARCOAL, linewidth=_LINE_W, zorder=1,
        )
        ax.add_patch(circle_patch)
        ax.plot(circ[0][0], circ[0][1], "+", color=_SILVER,
                markersize=7, markeredgewidth=1.0, zorder=4)
        margin = circ[1] * 0.15
        ax.set_xlim(circ[0][0] - circ[1] - margin, circ[0][0] + circ[1] + margin)
        ax.set_ylim(circ[0][1] - circ[1] - margin, circ[0][1] + circ[1] + margin)

    plt.tight_layout()
    if path is not None:
        fig.savefig(path + "plot_min_circle.pdf", bbox_inches="tight")
    if show:
        plt.show()

    # Runtime bar chart
    fig, ax = _new_fig(figsize=(4, 3.5))
    _style_ax(ax, "Runtime")
    ax.bar(["Exact", "Heuristic"], [t_exact, t_heur],
           color=[_CHARCOAL, _SILVER], width=0.40,
           edgecolor=_BG, linewidth=0.6)
    ax.set_ylabel("Time (s)", fontsize=8.5, color=_STEEL)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.4f"))
    plt.tight_layout()
    if path is not None:
        fig.savefig(path + "time.pdf", bbox_inches="tight")
    if show:
        plt.show()


# ---------------------------------------------------------------------------
# Triangulation
# ---------------------------------------------------------------------------

def plot_triangulation(tri_obj):
    """Plot the polygon and its triangulation diagonals.

    Args:
        tri_obj: A ``PolygonTriangulation`` instance.
    """
    fig, ax = _new_fig()
    _style_ax(ax, tri_obj.poly_name)
    ax.set_aspect("equal")

    xs = list(tri_obj.poly[:, 0]) + [tri_obj.poly[0, 0]]
    ys = list(tri_obj.poly[:, 1]) + [tri_obj.poly[0, 1]]
    ax.fill(xs, ys,
            facecolor=(*_hex_to_rgb(_MIST), _FILL_ALPHA),
            edgecolor=_INK, linewidth=_LINE_W, zorder=2)

    for diag in tri_obj.diagonals:
        ax.plot([diag[0][0], diag[1][0]], [diag[0][1], diag[1][1]],
                color=_SILVER, linewidth=_LINE_W_THIN,
                linestyle=(0, (3, 4)), zorder=3)

    ax.scatter(tri_obj.poly[:, 0], tri_obj.poly[:, 1],
               c=_INK, s=_PT_SIZE, zorder=4,
               edgecolors="white", linewidths=_PT_EDGE_W)

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Delaunay Triangulation
# ---------------------------------------------------------------------------

def plot_delaunay(dt_obj, title="Delaunay Triangulation", show_circumcircles=False):
    """Plot the Delaunay triangulation.

    Args:
        dt_obj: A ``DelaunayTriangulation`` instance.
        title: Title for the matplotlib figure.
        show_circumcircles: If True, draw dashed circumcircles for each triangle.
    """
    fig, ax = _new_fig()
    _style_ax(ax, title)
    ax.set_aspect("equal")

    triangles = dt_obj.get_triangles()
    edges = dt_obj.get_edges()

    # Light triangle fill
    for tri in triangles:
        xs = [p[0] for p in tri] + [tri[0][0]]
        ys = [p[1] for p in tri] + [tri[0][1]]
        ax.fill(xs, ys,
                facecolor=(*_hex_to_rgb(_MIST), 0.05),
                edgecolor="none", zorder=1)

    # Edges
    for edge in edges:
        ax.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]],
                color=_CHARCOAL, linewidth=_LINE_W, zorder=2)

    # Circumcircles
    if show_circumcircles:
        for circ in dt_obj.get_circumcircles():
            circle_patch = plt.Circle(
                circ[0], circ[1],
                facecolor="none",
                edgecolor=(*_hex_to_rgb(_SILVER), 0.4),
                linewidth=_LINE_W_THIN,
                linestyle=(0, (4, 5)), zorder=3,
            )
            ax.add_patch(circle_patch)

    # Points on top
    ax.scatter(
        dt_obj.points[:, 0], dt_obj.points[:, 1],
        c=_INK, s=_PT_SIZE, zorder=4,
        edgecolors="white", linewidths=_PT_EDGE_W,
    )

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Segment Intersections
# ---------------------------------------------------------------------------

def plot_intersections(si_obj, title="Segment Intersections"):
    """Plot segments and their intersection points.

    Args:
        si_obj: A ``SegmentIntersection`` instance.
        title: Title for the matplotlib figure.
    """
    fig, ax = _new_fig()
    _style_ax(ax, title)
    ax.set_aspect("equal")

    segs = si_obj.get_segments()
    intersections = si_obj.find_intersections()

    # Draw segments
    for seg in segs:
        ax.plot(
            [seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]],
            color=_CHARCOAL, linewidth=_LINE_W, zorder=2,
            solid_capstyle="round",
        )

    # Segment endpoints
    endpoints_x = [p[0] for seg in segs for p in seg]
    endpoints_y = [p[1] for seg in segs for p in seg]
    ax.scatter(endpoints_x, endpoints_y, c=_INK, s=_PT_SIZE, zorder=3,
               edgecolors="white", linewidths=_PT_EDGE_W)

    # Intersection points
    if intersections:
        ix = [p[0] for p in intersections]
        iy = [p[1] for p in intersections]
        ax.scatter(
            ix, iy, c=_INK, s=50, zorder=5,
            marker="o", edgecolors="white", linewidths=1.2,
        )

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Voronoi
# ---------------------------------------------------------------------------

def plot_voronoi(voronoi_obj, cells):
    """Plot the Voronoi diagram showing sites and cell edges.

    Args:
        voronoi_obj: A ``VoronoiDiagram`` instance.
        cells: The cell list returned by ``build_voronoi_diagram()``.
    """
    data = np.array(voronoi_obj.data)

    fig, ax = _new_fig(figsize=(5.5, 5.5))
    _style_ax(ax, "Voronoi Diagram")

    for c in cells:
        for line in c[1]:
            ax.plot([line[0][0], line[1][0]], [line[0][1], line[1][1]],
                    linewidth=_LINE_W_THIN,
                    color=(*_hex_to_rgb(_CHARCOAL), 0.45),
                    solid_capstyle="round", zorder=2)

    ax.scatter(data[:, 0], data[:, 1], c=_INK, s=_PT_SIZE, zorder=4,
               edgecolors="white", linewidths=_PT_EDGE_W)

    pad_x = (data[:, 0].max() - data[:, 0].min()) * 0.12
    pad_y = (data[:, 1].max() - data[:, 1].min()) * 0.12
    ax.set_xlim(data[:, 0].min() - pad_x, data[:, 0].max() + pad_x)
    ax.set_ylim(data[:, 1].min() - pad_y, data[:, 1].max() + pad_y)

    plt.tight_layout()
    plt.show()
