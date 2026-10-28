# sdf — signed-distance shape kernel

↑ [`section_geometry`](../README.md)

Domain-agnostic 2-D SDF library: closed-form primitives, CSG (`min` / `max`),
aerofoil construction, optional expression IR, and dense sampling via
`SDFGrid`. Blade region recipes live in [`../sections/`](../sections/README.md).

## Table of contents

1. [`AirfoilSDF`](#airfoilsdf)
2. [Organisation](#organisation)
3. [API](#api)
4. [Usage in the package](#usage-in-the-package)
5. [See also](#see-also)

---

## `AirfoilSDF`

The usual station outline. Construct from a NACA code, a vertex array, or a
Selig `.dat` file. The instance is callable: `phi = airfoil(x, y)` with
`phi < 0` inside.

**Contract (SG-008).** The discrete contour polygon *is* the section geometry.
Offsets, midlines, junctions, and solids are exact with respect to that
polyline — not with respect to an analytic aerofoil curve. Default
`n_points=200` is a sampling choice that sets how closely a NACA generator
matches its closed form (~62 µm LE deviation on NACA 2412 at 1 m chord);
raising `n_points` refines the polygon, it does not “fix” a package defect.
Micron-level junction figures elsewhere in the package are relative to this
polygon. See [`../docs/bug_log.md`](../docs/bug_log.md) SG-008.

```python
from blade_structure_utilities.geometry_utilities.section_geometry import AirfoilSDF

airfoil = AirfoilSDF.from_naca("2412", chord=1.0, n_points=200)
# airfoil.chord, airfoil.vertices — reference chord [m] and contour
# airfoil(x, y) → signed distance
```

| Constructor | Input |
|-------------|--------|
| `from_naca(code, *, n_points=200, chord=1.0, closed_te=True)` | NACA 4-digit code |
| `from_array(xy, chord=None)` | `(N, 2)` contour |
| `from_dat(path, chord=None)` | Selig-format file |
| `AirfoilSDF(vertices, chord=None)` | Direct vertices |

Helpers include camber line and thickness distribution (used by layout for
optional max-thickness cap anchors).

## Organisation

| File | Role |
|------|------|
| `primitives.py` | Closed-form 2-D SDFs |
| `csg.py` | Boolean / offset / shell operators |
| `ir.py` | Expression AST and evaluation |
| `airfoil.py` | `AirfoilSDF` |
| `naca.py` | NACA vertex generators |
| `grid.py` | Dense sampling |
| `transforms.py` | Rotate / translate / scale |
| `_polygon_numba.py` | Fast polygon edge-loop kernel |
| `cache.py`, `profiling.py` | Eval cache and timing |

## API

| Name | File | Role |
|------|------|------|
| `AirfoilSDF` | `airfoil.py` | Aerofoil SDF; prefer `from_naca(...)` |
| `SDFGrid` | `grid.py` | Sampled field |
| CSG helpers | `csg.py` | `union`, `intersect`, `difference`, `offset`, `slab`, … |
| Primitives | `primitives.py` | Capsule, circle, half-plane, … |

## Usage in the package

`sections/` composes these callables into named blade regions. Property and
midline paths consume the resulting host, not the raw primitives.

## See also

- [`../sections/`](../sections/README.md)
