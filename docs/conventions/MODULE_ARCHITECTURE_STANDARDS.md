# Module architecture

**Audience: internal.** Down-only import layers, extras, and house-rule splits.

↑ [conventions](README.md)

Top-level packages under `src/implicit_geometry_kernel/` mirror `docs/`, `tests/`,
and `examples/` where those exist.

## Import layers (down-only)

```
contract ← sdf / integration / frames ← structure ← catalogue ← section ← blade ← extract
```

`sdf` imports nothing from `section` or `blade`. Enforced by
`tests/test_import_boundary.py`.

## Extras (not layers)

Geometry modules must not import these. Core-install (`pip install .`) still
imports every non-plotting module because extras are lazy at the call site.

| Extra | Install | Lives in | Role |
|---|---|---|---|
| `[plots]` | matplotlib | `plotting/` | figures |
| `[mesh]` | gmsh | `mesh/gmsh/` | unadvertised convenience; consumers own real Gmsh jobs |
| `[cad]` | cadquery | `blade/cad.py` | STEP writer |
| `[splines]` | scipy | `blade/airfoil.py` CubicSpline only | dense LE resample; falls back to the polyline |
| `[reference]` | pyyaml | tests and example 08 | WindIO YAML; `blade/windio.py` takes mappings |
| profiling | none | `profiling.py` | timers; not on root `__all__` |

`mesh` may import `blade`, `extract`, `contract`, `geos`, and `profiling`.
Kernel polygon fill lives in `mesh/polygon.py` (Qhull Voronoi + numba clip;
faces have more than four nodes). Numpy Steiner Delaunay is `mesh/delaunay.py`.
The Gmsh backend is `mesh/gmsh/` (`shell.py` mapped strips, `plane.py` planar
cells). That extra is leftover convenience. The kernel station/CSG mesh is
`integration/quadtree` (`build_shared_quadtree`).

## House-rule splits

At most three levels: package / subpackage / module.py. A former oversized
module may become a subpackage (`blade/topology/`, `sdf/kernels/`, `sdf/nodes/`,
`frames/tensors/`, `integration/quadtree/`). Do not nest further.

Modules ≤ 400 lines; packages ≤ 12 modules. Split inside the owning folder.

## Dependencies

Runtime core: numpy, numba (numba pulls llvmlite). Offset and overlay need
**GEOS** (`libgeos_c`): system/conda (`libgeos-dev`, conda-forge `geos`) or
`GEOS_LIBRARY_PATH`. There is no Shapely dependency.

Constructive polygons go through `PolygonOps` (`offset_inward`, overlay,
`make_valid`, `rotate`, `envelope`). Algorithms are GEOS (JTS C++) via
`implicit_geometry_kernel.geos` (ctypes, WKB, kernel `Polygon`). Do not
reimplement buffer or overlay in numpy. `tests/section/test_offset_import_boundary.py`
forbids `import shapely`. Public products (`ShellMidlineStrip`, `CellWallLoop`,
properties) are numpy. A C++ port of the kernel links the same GEOS C functions.

Grey extra `[contours]`: `sdf.kernels.grid` may use **contourpy** for debug
level-sets if installed. It is not a core dependency. There is no matplotlib
fallback (geometry must not import matplotlib).

Span `TableLaw` is linear `numpy` interpolation, not SciPy.

## Frames

Twist, NACA pivot, and B/S maps live in `frames`. Do not reimplement 2×2 twist
outside that package.
