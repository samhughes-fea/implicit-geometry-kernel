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
| `[mesh]` | gmsh | `mesh/` | convenience shell mesh; consumers still own real Gmsh jobs |
| `[cad]` | cadquery | `blade/cad.py` | STEP writer |
| `[splines]` | scipy | `blade/airfoil.py` CubicSpline only | dense LE resample; falls back to the polyline |
| `[reference]` | pyyaml | tests and example 08 | WindIO YAML; `blade/windio.py` takes mappings |
| profiling | none | `profiling.py` | timers; not on root `__all__` |

`mesh` may import `blade`, `extract`, `contract`, and `profiling` only.

## House-rule splits

At most three levels: package / subpackage / module.py. A former oversized
module may become a subpackage (`blade/topology/`, `sdf/kernels/`, `sdf/nodes/`,
`frames/tensors/`, `integration/quadtree/`). Do not nest further.

Modules ≤ 400 lines; packages ≤ 12 modules. Split inside the owning folder.

## Dependencies

Runtime core: numpy, shapely, numba (numba pulls llvmlite).

Shapely is the C++-port risk: it wraps GEOS (`buffer`, overlay, `make_valid`).
The only `import shapely` in `src/` is `section/boundary/polygons.py`. Callers use
`PolygonOps` (`offset_inward`, overlay, `make_valid`, `rotate`, `envelope`).
`KernelPolygonOps` supplies numpy affine maps; overlay still uses GEOS until a
numpy offseter matches the oracle tests. Public products (`ShellMidlineStrip`,
`CellWallLoop`, properties) are numpy; Shapely `Polygon` is constructive scratch,
not a hand-off type. A C++ port copies SDF kernels and this protocol, not Shapely.
Golden numbers may move if GEOS does — changelog that. Do not add Shapely imports
outside the gateway (`tests/section/test_offset_import_boundary.py`).

Grey: `sdf.kernels.grid` may use **contourpy** for debug level-sets if installed.
It is not a core dependency. There is no matplotlib fallback (geometry must not
import matplotlib).

Span `TableLaw` PCHIP is in-house numpy, not SciPy.

## Frames

Twist, NACA pivot, and B/S maps live in `frames`. Do not reimplement 2×2 twist
outside that package.
