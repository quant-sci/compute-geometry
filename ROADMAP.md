# Roadmap

Development roadmap for `compute-geometry`.

## v0.1.2 — Foundation (Testing & Code Quality)

- [x] Unit tests for all 4 algorithms (ConvexHull, MinimumCircle, PolygonTriangulation, VoronoiDiagram)
- [x] Input validation: type checks, degenerate cases (collinear points, duplicates, insufficient points)
- [x] Docstrings for all public classes and methods
- [x] Separate visualization from algorithm logic
- [x] Implement the `elements` module (Point, Line, Segment, Polygon, Circle primitives)

## v0.1.3 — Algorithm Expansion

- [ ] Delaunay triangulation
- [ ] Line segment intersection (Bentley-Ottmann)
- [ ] Closest pair of points
- [ ] Point-in-polygon testing
- [ ] Polygon area and centroid computation
- [ ] Boolean operations on polygons (union, intersection, difference)

## v0.1.4 — Performance & Usability

- [ ] NumPy-optimized computation paths
- [ ] Consistent API across all algorithms
- [ ] Type hints throughout
- [ ] Comprehensive documentation site (Sphinx)
- [ ] More example scripts
- [ ] Change the plots design for a better experience

## v0.1.5 — 3D Meshes

- [ ] 3D primitives: Point3D, Vector3D, Plane, Triangle3D
- [ ] Mesh data structure (half-edge or indexed face set)
- [ ] STL / OBJ / PLY file I/O
- [ ] Mesh normals computation (face and vertex normals)
- [ ] Mesh surface area and volume
- [ ] 3D convex hull
- [ ] Mesh boolean operations (union, intersection, difference)
- [ ] 3D visualization with matplotlib 3D projections

## v0.1.6 — Dev Ready

- [ ] Stable public API
- [ ] Robust numerical handling (floating-point edge cases)
- [ ] Benchmarks against scipy.spatial
- [ ] Full CI: coverage reporting, type checking

## v0.1.7 — ...
