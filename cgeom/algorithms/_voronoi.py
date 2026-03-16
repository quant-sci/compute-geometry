import numpy as np

class VoronoiDiagram:
    """Construct a Voronoi diagram for a set of 2D sites using incremental insertion.

    Attributes:
        data (np.ndarray): Array of site coordinates with shape (n, 2).
        cells (list): List of Voronoi cells, each containing a site and its edges.
    """

    def __init__(self, data):
        """Initialize the Voronoi diagram with a set of 2D sites.

        Args:
            data: Collection of 2D points as a list, tuple, or numpy array
                with shape (n, 2). The first two points must have different
                y-coordinates.

        Raises:
            pydantic.ValidationError: If fewer than 2 points, first two share
                a y-coordinate, non-numeric, or wrong shape.
        """
        from cgeom.elements.models import VoronoiDiagramInput
        validated = VoronoiDiagramInput(points=data)
        self.data = np.array(validated.points)
        self.cells = []

    def search_for_cell(self, point):
        """Find the index of the Voronoi cell whose site is nearest to *point*.

        Args:
            point: A 2D coordinate [x, y].

        Returns:
            int: Index of the nearest cell in ``self.cells``.
        """
        dist = float("inf")
        for i, cell in enumerate(self.cells):
            aux = np.sqrt((point[0] - cell[0][0])**2 + (point[1] - cell[0][1])**2)
            if aux < dist:
                dist = aux
                pos = i
        return pos

    def search_for_cell_by_edge(self, edge_ref, origin_cell):
        """Find a neighbouring cell that shares the given edge.

        Args:
            edge_ref: Reference edge as ``[[x1, y1, ...], [x2, y2, ...]]``.
            origin_cell: Index of the cell to exclude from the search.

        Returns:
            int: Index of the neighbouring cell that contains *edge_ref*.
        """
        precision = 0.8
        for i, cell in enumerate(self.cells):
            for j, edge in enumerate(cell[1]):
                if i != origin_cell:
                    if (
                        abs(edge[0][0] - edge_ref[0][0]) <= precision
                        and abs(edge[0][1] - edge_ref[0][1]) <= precision
                        and abs(edge[1][0] - edge_ref[1][0]) <= precision
                        and abs(edge[1][1] - edge_ref[1][1]) <= precision
                    ):
                        output = i

                    if (
                        abs(edge[1][0] - edge_ref[0][0]) <= precision
                        and abs(edge[1][1] - edge_ref[0][1]) <= precision
                        and abs(edge[0][0] - edge_ref[1][0]) <= precision
                        and abs(edge[0][1] - edge_ref[1][1]) <= precision
                    ):
                        output = i
        return output

    def adapt_cell(self, visited_cell, edge_added, point_added):
        """Update an existing cell after a new site is inserted.

        Clips or removes edges that are closer to *point_added* than to the
        cell's own site, and appends the new bisector edge.

        Args:
            visited_cell: Index of the cell to update.
            edge_added: The new bisector edge to insert.
            point_added: The newly added site coordinates.
        """
        drop = []
        for i, edge in enumerate(self.cells[visited_cell][1]):
            A = edge[0]
            B = edge[1]
            C = edge_added[0]
            D_0 = round(A[0] * B[1] + A[1] * C[0] + B[0] * C[1] - C[0] * B[1] - A[0] * C[1] - A[1] * B[0], 4)
            C = edge_added[1]
            D_1 = round(A[0] * B[1] + A[1] * C[0] + B[0] * C[1] - C[0] * B[1] - A[0] * C[1] - A[1] * B[0], 4)

            # Determine precision based on the magnitude of D
            if max(abs(D_0), abs(D_1)) >= 2000000:
                precision = 1.5
            elif max(abs(D_0), abs(D_1)) < 2000000 and max(abs(D_0), abs(D_1)) > 50000:
                precision = 0.8
            elif max(abs(D_0), abs(D_1)) <= 50000 and max(abs(D_0), abs(D_1)) > 1000:
                precision = 0.5
            elif max(abs(D_0), abs(D_1)) <= 1000:
                precision = 0.01

            if abs(D_0) <= precision or abs(D_1) <= precision:
                if abs(D_0) <= precision:
                    p = edge_added[0]
                else:
                    p = edge_added[1]

                dist_1 = np.sqrt((point_added[0] - A[0]) ** 2 + (point_added[1] - A[1]) ** 2)
                dist_2 = np.sqrt((self.cells[visited_cell][0][0] - A[0]) ** 2 + (self.cells[visited_cell][0][1] - A[1]) ** 2)

                if dist_1 < dist_2:
                    self.cells[visited_cell][1][i][0] = [p[0], p[1], 1][:]
                else:
                    self.cells[visited_cell][1][i][1] = [p[0], p[1], 1][:]

            if abs(D_0) > precision and abs(D_1) > precision:
                dist_1 = np.sqrt((point_added[0] - A[0]) ** 2 + (point_added[1] - A[1]) ** 2)
                dist_2 = np.sqrt((self.cells[visited_cell][0][0] - A[0]) ** 2 + (self.cells[visited_cell][0][1] - A[1]) ** 2)

                if dist_1 < dist_2:
                    drop.append(edge)

        self.cells[visited_cell][1].append(edge_added[:])

        for element in drop:
            self.cells[visited_cell][1].remove(element)

    def make_edge(self, cell, p):
        """Compute the bisector edge between a cell's site and a new point.

        The bisector is intersected with the cell's existing edges to
        determine its extent.

        Args:
            cell: Index of the existing cell.
            p: The new point [x, y].

        Returns:
            tuple: ``(edge_endpoints, neighbour_indices)`` where
                *edge_endpoints* is a list of intersection points and
                *neighbour_indices* lists the adjacent cells crossed.
        """
        C = self.cells[cell][0]
        M = [round((C[0] + p[0]) / 2, 4), round((C[1] + p[1]) / 2, 4)]
        V_1 = [C[0] - M[0], C[1] - M[1]]
        if V_1[1] == 0:
            V_2 = [0, 1]
        else:
            V_2 = [1, -V_1[0] / V_1[1]]
        output, NEIGHBOURS = [], []

        for edge in self.cells[cell][1]:
            A = edge[0]
            B = edge[1]

            if A[0] == B[0]:
                A[0] = A[0] + 0.1
                B[0] = B[0] - 0.1

            if A[1] == B[1]:
                A[1] = A[1] + 0.1
                B[1] = B[1] - 0.1

            V_3 = [B[0] - A[0], B[1] - A[1]]
            D = round((-V_3[0] * V_2[1] + V_3[1] * V_2[0]), 4)

            if abs(D) > 0.1:
                D_t = (M[1] - A[1]) * V_3[0] - (M[0] - A[0]) * V_3[1]
                t = D_t / D
                INTERSECT = [round(M[0] + t * V_2[0], 4), round(M[1] + t * V_2[1], 4), 1]
                IS_VALID = False

                if A[2] == 1 and B[2] == 1:
                    if INTERSECT[0] < max(A[0], B[0]) and INTERSECT[0] > min(A[0], B[0]):
                        IS_VALID = True

                if A[2] == 0 and B[2] == 1:
                    if A[0] > B[0] and INTERSECT[0] > B[0]:
                        IS_VALID = True
                    if A[0] < B[0] and INTERSECT[0] < B[0]:
                        IS_VALID = True

                if A[2] == 1 and B[2] == 0:
                    if A[0] > B[0] and INTERSECT[0] < A[0]:
                        IS_VALID = True
                    if A[0] < B[0] and INTERSECT[0] > A[0]:
                        IS_VALID = True

                if A[2] == 0 and B[2] == 0:
                    IS_VALID = True

                if IS_VALID:
                    T = t
                    output.append(INTERSECT)
                    NEIGHBOURS.append(self.search_for_cell_by_edge(edge, cell))

        if len(output) == 1:
            if self.search_for_cell(M) == cell:
                inf_point = [round(M[0] - 10 * T * V_2[0], 4), round(M[1] - 10 * T * V_2[1], 4), 0]
                i = 2
                while np.sqrt((self.cells[cell][0][0] - inf_point[0]) ** 2 + (self.cells[cell][0][1] - inf_point[0]) ** 2) < 1000:
                    inf_point = [round(M[0] - i * 10 * T * V_2[0], 4), round(M[1] - i * 10 * T * V_2[1], 4), 0]
                    i += 1

            else:
                inf_point = [round(M[0] + 10 * T * V_2[0], 4), round(M[1] + 10 * T * V_2[1], 4), 0]
                i = 2
                while (
                    np.sqrt((self.cells[cell][0][0] - inf_point[0]) ** 2 + (self.cells[cell][0][1] - inf_point[0]) ** 2)
                    < 1000
                    or self.search_for_cell(inf_point) != cell
                ):
                    inf_point = [round(M[0] + i * 10 * T * V_2[0], 4), round(M[1] + i * 10 * T * V_2[1], 4), 0]
                    i += 1

            output.append(inf_point)

        return output, NEIGHBOURS

    def construct_cell(self, cell, p):
        """Construct a new Voronoi cell for point *p* by traversing neighbours.

        Args:
            cell: Index of the starting cell (the one containing *p*).
            p: The new site coordinates.

        Returns:
            list: A cell entry ``[site, edges]`` for the new site.
        """
        edges, VISITED = [], []
        edge, NEIGHBOURS = self.make_edge(cell, p)
        VISITED.append(cell)
        edges.append(edge)

        for N in NEIGHBOURS:
            NEXT_CELL = N
            if NEXT_CELL not in VISITED:
                FLAG = True
            while FLAG:
                edge, NEIGHBOURS = self.make_edge(NEXT_CELL, p)
                edges.append(edge)
                VISITED.append(NEXT_CELL)
                FLAG = False
                if NEIGHBOURS != []:
                    for NEIGHBOUR in NEIGHBOURS:
                        if NEIGHBOUR not in VISITED:
                            FLAG = True
                            NEXT_CELL = NEIGHBOUR

        for i, VISIT in enumerate(VISITED):
            self.adapt_cell(VISIT, edges[i], p)

        return [p, edges]

    def build_voronoi_diagram(self):
        """Build the full Voronoi diagram by incrementally inserting each site.

        Returns:
            list: The list of Voronoi cells (also stored in ``self.cells``).
        """
        mid = [(self.data[1][0] + self.data[0][0]) / 2, (self.data[1][1] + self.data[0][1]) / 2]
        V = [self.data[0][0] - mid[0], self.data[0][1] - mid[1]]
        self.cells.append(
            [
                self.data[0],
                [
                    [
                        [round(mid[0] + 1000, 4), round(mid[1] - (V[0] / V[1]) * 1000, 4), 0],
                        [round(mid[0] - 1000, 4), round(mid[1] + (V[0] / V[1]) * 1000, 4), 0],
                    ]
                ],
            ]
        )
        self.cells.append(
            [
                self.data[1],
                [
                    [
                        [round(mid[0] + 1000, 4), round(mid[1] - (V[0] / V[1]) * 1000, 4), 0],
                        [round(mid[0] - 1000, 4), round(mid[1] + (V[0] / V[1]) * 1000, 4), 0],
                    ]
                ],
            ]
        )

        for i in range(2, len(self.data)):
            cell = self.search_for_cell(self.data[i])
            self.cells.append(self.construct_cell(cell, self.data[i]))

        return self.cells

    def plot_voronoi(self, cells):
        """Deprecated: use cgeom.visualization.plot_voronoi() instead."""
        import warnings
        warnings.warn(
            "VoronoiDiagram.plot_voronoi() is deprecated. "
            "Use cgeom.visualization.plot_voronoi(voronoi_obj, cells) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        from cgeom.visualization import plot_voronoi
        plot_voronoi(self, cells)
