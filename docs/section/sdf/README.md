# sdf (not a section module)

**Audience: internal.** Maintainer notes; not part of the library manual.

The shape kernel lives in the package `implicit_geometry_kernel.sdf`, not under
`section/`. Prefer `Airfoil.naca4(...).to_sdf()`; `AirfoilSDF` is the implicit
adapter (`sdf.shapes.airfoil`). CSG trees import from `sdf.nodes`.

See the kernel [README](../../README.md) and `src/implicit_geometry_kernel/sdf/`.
