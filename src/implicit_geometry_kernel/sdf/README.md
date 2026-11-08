# sdf

2-D SDF / CSG / feature field. Imports neither `section` nor `blade`.

**Import:** `implicit_geometry_kernel.sdf`

Tree: `kernels/` (arrays, `SDFGrid`), `nodes/` (`SDFNode` IR), `shapes/`
(`AirfoilSDF`), `solve.py` (`GeometrySolver`). Product integrals live in
`integration`; zero-set polylines in `extract`.

↑ [package](../README.md) · [docs](../../../../docs/sdf/README.md) ·
[examples](../../../../examples/sdf/README.md) · [tests](../../../../tests/sdf/README.md)

Prefer `Airfoil.to_sdf()` from `blade` for station outlines.
Barrel SDF names: [glossary](../../../../docs/glossary.md#sdf).
