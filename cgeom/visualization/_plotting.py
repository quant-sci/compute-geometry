import time

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ---------------------------------------------------------------------------
# Grays palette — sample 4 evenly-spaced stops for discrete use
# ---------------------------------------------------------------------------
_CMAP = plt.cm.Greys
_C0 = _CMAP(0.85)  # near-black
_C1 = _CMAP(0.60)  # dark gray
_C2 = _CMAP(0.40)  # medium gray
_C3 = _CMAP(0.20)  # light gray
_GRID = "#e8e8e8"   # very faint grid

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _style_ax(ax, title=None):
    """Apply a minimal look to an axes object."""
    if title:
        ax.set_title(title, fontsize=12, fontweight="500", pad=10)
    ax.set_facecolor("white")
    ax.grid(False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("bottom", "left"):
        ax.spines[spine].set_color(_GRID)
        ax.spines[spine].set_linewidth(0.7)
    ax.tick_params(colors="#555555", labelsize=8, length=3, width=0.6)


def _new_fig(figsize=(6, 5)):
    """Create a figure with a white background."""
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
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

    ch = hull_obj.convex_hull()
    xs = [p[0] for p in ch] + [ch[0][0]]
    ys = [p[1] for p in ch] + [ch[0][1]]
    ax.fill(xs, ys, facecolor=(*_C2[:3], 0.12), edgecolor=_C1, linewidth=1.6, zorder=2)

    ax.scatter(
        hull_obj.points[:, 0], hull_obj.points[:, 1],
        c=_C0, s=30, zorder=3, linewidths=0,
    )
    ax.scatter(
        [p[0] for p in ch], [p[1] for p in ch],
        c=[_C1], s=44, zorder=4, edgecolors="white", linewidths=0.5,
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
        fig.patch.set_facecolor("white")

        for idx, (circ, label) in enumerate([
            (min_circle, f"Exact (n={num_points})"),
            (min_circle_heuristic, f"Heuristic (n={num_points})"),
        ]):
            ax = axes[idx]
            _style_ax(ax, label)
            ax.set_aspect(1)
            ax.scatter(x, y, c=[_C0], s=16, zorder=3, linewidths=0)
            circle_patch = plt.Circle(
                circ[0], circ[1],
                facecolor=(*_C2[:3], 0.10), edgecolor=_C1, linewidth=1.4,
                zorder=1,
            )
            ax.add_patch(circle_patch)
            ax.plot(circ[0][0], circ[0][1], '+', color=_C3,
                    markersize=8, markeredgewidth=1.5, zorder=4)
            ax.set_xlim(-10, 10)
            ax.set_ylim(-10, 10)

        plt.tight_layout()
        if path is not None:
            fig.savefig(path + f"plot_min_circle_{num_points}_points.pdf",
                        bbox_inches="tight")
        if show:
            plt.show()

    # Runtime comparison
    fig, ax = _new_fig(figsize=(6, 4))
    _style_ax(ax, "Runtime")
    ax.plot(sizes, time_min_circle, color=_C1, linewidth=1.6,
            marker="o", markersize=4, label="Exact")
    ax.plot(sizes, time_heuristic, color=_C2, linewidth=1.6,
            marker="s", markersize=4, label="Heuristic")
    ax.set_xlabel("Input size", fontsize=10)
    ax.set_ylabel("Time (s)", fontsize=10)
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
    fig.patch.set_facecolor("white")

    for idx, (circ, label) in enumerate([
        (min_circle, "Exact"),
        (min_circle_heuristic, "Heuristic"),
    ]):
        ax = axes[idx]
        _style_ax(ax, label)
        ax.set_aspect(1)
        ax.scatter(data[:, 0], data[:, 1], c=[_C0], s=16,
                   zorder=3, linewidths=0)
        circle_patch = plt.Circle(
            circ[0], circ[1],
            facecolor=(*_C2[:3], 0.10), edgecolor=_C1, linewidth=1.4,
            zorder=1,
        )
        ax.add_patch(circle_patch)
        ax.plot(circ[0][0], circ[0][1], '+', color=_C3,
                markersize=8, markeredgewidth=1.5, zorder=4)
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
           color=[_C1, _C2], width=0.45, edgecolor="white", linewidth=0.6)
    ax.set_ylabel("Time (s)", fontsize=10)
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
    ax.set_aspect(1)

    xs = list(tri_obj.poly[:, 0]) + [tri_obj.poly[0, 0]]
    ys = list(tri_obj.poly[:, 1]) + [tri_obj.poly[0, 1]]
    ax.fill(xs, ys, facecolor=(*_C0[:3], 0.08), edgecolor=_C0, linewidth=1.4, zorder=2)

    for diag in tri_obj.diagonals:
        ax.plot([diag[0][0], diag[1][0]], [diag[0][1], diag[1][1]],
                color=_C1, linewidth=1, linestyle="--", zorder=3)

    ax.scatter(tri_obj.poly[:, 0], tri_obj.poly[:, 1],
               c=[_C2], s=32, zorder=4, edgecolors="white", linewidths=0.4)

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Voronoi
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
    ax.set_aspect(1)

    triangles = dt_obj.get_triangles()
    edges = dt_obj.get_edges()

    # Light triangle fill
    for tri in triangles:
        xs = [p[0] for p in tri] + [tri[0][0]]
        ys = [p[1] for p in tri] + [tri[0][1]]
        ax.fill(xs, ys, facecolor=(*_C2[:3], 0.06), edgecolor="none", zorder=1)

    # Edges
    for edge in edges:
        ax.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]],
                color=_C1, linewidth=1.2, zorder=2)

    # Circumcircles
    if show_circumcircles:
        for circ in dt_obj.get_circumcircles():
            circle_patch = plt.Circle(
                circ[0], circ[1],
                facecolor="none", edgecolor=(*_C1[:3], 0.6),
                linewidth=1.0, linestyle="--", zorder=3,
            )
            ax.add_patch(circle_patch)

    # Points on top
    ax.scatter(
        dt_obj.points[:, 0], dt_obj.points[:, 1],
        c=[_C0], s=36, zorder=4, edgecolors="white", linewidths=0.5,
    )

    plt.tight_layout()
    plt.show()


def plot_voronoi(voronoi_obj, cells):
    """Plot the Voronoi diagram showing sites and cell edges.

    Args:
        voronoi_obj: A ``VoronoiDiagram`` instance.
        cells: The cell list returned by ``build_voronoi_diagram()``.
    """
    data = np.array(voronoi_obj.data)

    fig, ax = _new_fig(figsize=(6, 6))
    _style_ax(ax, "Voronoi Diagram")

    for c in cells:
        for line in c[1]:
            ax.plot([line[0][0], line[1][0]], [line[0][1], line[1][1]],
                    linewidth=0.8, color=(*_C0[:3], 0.5), zorder=2)

    ax.scatter(data[:, 0], data[:, 1], c=[_C1], s=36, zorder=4,
               edgecolors="white", linewidths=0.4)

    pad_x = (data[:, 0].max() - data[:, 0].min()) * 0.15
    pad_y = (data[:, 1].max() - data[:, 1].min()) * 0.15
    ax.set_xlim(data[:, 0].min() - pad_x, data[:, 0].max() + pad_x)
    ax.set_ylim(data[:, 1].min() - pad_y, data[:, 1].max() + pad_y)

    plt.tight_layout()
    plt.show()
