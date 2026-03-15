import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from cgeom.algorithms import VoronoiDiagram


class TestVoronoiTwoPoints:
    """Two points produce two cells split by a perpendicular bisector."""

    def setup_method(self):
        self.data = [[0, 0], [4, 2]]
        self.vd = VoronoiDiagram(self.data)
        self.cells = self.vd.build_voronoi_diagram()

    def test_cell_count(self):
        assert len(self.cells) == 2

    def test_each_cell_has_site(self):
        sites = {(c[0][0], c[0][1]) for c in self.cells}
        assert (0, 0) in sites
        assert (4, 2) in sites

    def test_each_cell_has_edges(self):
        for cell in self.cells:
            assert len(cell[1]) >= 1


class TestVoronoiThreePoints:
    """Three non-collinear points from known-good example data."""

    def setup_method(self):
        # Subset from examples/points1.txt — first two have different y-coords
        self.data = [[326, 237], [373, 209], [378, 265]]
        self.vd = VoronoiDiagram(self.data)
        self.cells = self.vd.build_voronoi_diagram()

    def test_cell_count(self):
        assert len(self.cells) == 3

    def test_sites_match_input(self):
        sites = sorted([(c[0][0], c[0][1]) for c in self.cells])
        expected = sorted([(326, 237), (373, 209), (378, 265)])
        for s, e in zip(sites, expected):
            assert s[0] == pytest.approx(e[0])
            assert s[1] == pytest.approx(e[1])


class TestVoronoiSixPoints:
    """Six points from example data — structural checks."""

    def setup_method(self):
        self.data = [
            [326, 237], [373, 209], [378, 265],
            [443, 241], [396, 231], [416, 270],
        ]
        self.vd = VoronoiDiagram(self.data)
        self.cells = self.vd.build_voronoi_diagram()

    def test_cell_count(self):
        assert len(self.cells) == 6

    def test_all_cells_have_edges(self):
        for cell in self.cells:
            assert len(cell[1]) >= 1

    def test_sites_are_unique(self):
        sites = [(round(c[0][0], 4), round(c[0][1], 4)) for c in self.cells]
        assert len(set(sites)) == 6


class TestVoronoiSearchForCell:
    """search_for_cell returns the index of the nearest site."""

    def setup_method(self):
        self.data = [[326, 237], [373, 209], [378, 265]]
        self.vd = VoronoiDiagram(self.data)
        self.cells = self.vd.build_voronoi_diagram()

    def test_point_near_first_site(self):
        idx = self.vd.search_for_cell([327, 238])
        site = self.cells[idx][0]
        assert site[0] == pytest.approx(326)
        assert site[1] == pytest.approx(237)

    def test_point_near_second_site(self):
        idx = self.vd.search_for_cell([374, 210])
        site = self.cells[idx][0]
        assert site[0] == pytest.approx(373)
        assert site[1] == pytest.approx(209)


class TestVoronoiFullExample:
    """Full example dataset from points1.txt (20 points)."""

    def setup_method(self):
        self.data = [
            [326, 237], [373, 209], [378, 265], [443, 241], [396, 231],
            [416, 270], [361, 335], [324, 297], [400, 306], [454, 315],
            [489, 285], [488, 234], [443, 185], [421, 137], [380, 169],
            [315, 160], [297, 204], [267, 248], [265, 344], [342, 263],
        ]
        self.vd = VoronoiDiagram(self.data)
        self.cells = self.vd.build_voronoi_diagram()

    def test_cell_count(self):
        assert len(self.cells) == 20

    def test_all_cells_have_edges(self):
        for cell in self.cells:
            assert len(cell[1]) >= 1
