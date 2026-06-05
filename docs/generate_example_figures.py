"""Regenerate the figures embedded in ``docs/examples.md``.

Each visualization helper ends in ``plt.show()``; here we patch ``show`` to save
the current figure instead, so the rendered output matches the documented code.

Run from the repository root::

    .venv/bin/python docs/generate_example_figures.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from cgeom.algorithms import (  # noqa: E402
    ConvexHull,
    DelaunayTriangulation,
    MinimumCircle,
    PolygonTriangulation,
    SegmentIntersection,
    VoronoiDiagram,
)
from cgeom.visualization import (  # noqa: E402
    plot_convex_hull,
    plot_delaunay,
    plot_intersections,
    plot_min_circle,
    plot_triangulation,
    plot_voronoi,
)

OUT_DIR = "public/examples"
os.makedirs(OUT_DIR, exist_ok=True)

_saved = set()
_current = {"name": None}


def _save_show():
    """Stand-in for ``plt.show`` that saves the first figure per example."""
    name = _current["name"]
    if name and name not in _saved:
        fig = plt.gcf()
        fig.savefig(
            os.path.join(OUT_DIR, f"{name}.png"),
            dpi=150,
            facecolor=fig.get_facecolor(),
            bbox_inches="tight",
        )
        _saved.add(name)
    plt.close("all")


plt.show = _save_show


def render(name, fn):
    _current["name"] = name
    fn()
    print(f"saved {OUT_DIR}/{name}.png")


# --- Convex hull -----------------------------------------------------------
hull = ConvexHull(
    [
        (326, 237),
        (373, 209),
        (378, 265),
        (443, 241),
        (396, 231),
        (416, 270),
        (361, 335),
        (324, 297),
    ]
)
render("convex_hull", lambda: plot_convex_hull(hull))

# --- Minimum enclosing circle ---------------------------------------------
mc_points = np.array([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
render("min_circle", lambda: plot_min_circle(MinimumCircle(), mc_points, show=True))

# --- Delaunay triangulation ------------------------------------------------
dt = DelaunayTriangulation([(0, 0), (4, 0), (4, 4), (0, 4), (2, 2)])
dt.triangulate()
render("delaunay", lambda: plot_delaunay(dt, show_circumcircles=True))

# --- Voronoi diagram -------------------------------------------------------
voronoi = VoronoiDiagram(np.loadtxt("examples/points1.txt"))
cells = voronoi.build_voronoi_diagram()
render("voronoi", lambda: plot_voronoi(voronoi, cells))

# --- Polygon triangulation -------------------------------------------------
pt = PolygonTriangulation([[0, 0], [4, 0], [4, 4], [2, 2], [0, 4]])
render("triangulation", lambda: plot_triangulation(pt))

# --- Segment intersection --------------------------------------------------
si = SegmentIntersection(
    [[[0, 0], [4, 4]], [[0, 4], [4, 0]], [[0, 2], [4, 2]], [[2, 0], [2, 4]]]
)
render("intersection", lambda: plot_intersections(si))
