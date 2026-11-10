# Signed-distance fields

↑ [docs](../README.md) · [src](../../src/implicit_geometry_kernel/sdf/README.md) ·
[examples](../../examples/sdf/README.md) · [tests](../../tests/sdf/README.md)

The `sdf` package supplies domain-agnostic two-dimensional signed-distance
fields, constructive solid geometry (CSG), feature tags, and an in-plane
Newton solver. It imports neither `section` nor `blade`.

**Import:** `implicit_geometry_kernel.sdf`

| Subpackage | Role |
|---|---|
| `sdf.kernels` | Array primitives, Numba polygon, functional CSG, `SDFGrid` |
| `sdf.nodes` | `SDFNode` expression tree (`Expr`, `Polygon`, `Union`, …) |
| `sdf.shapes` | One `AirfoilSDF`; NACA vertex generators (no blade types) |
| `sdf.solve` | `GeometrySolver` — in-plane Newton on fields |

Station callers obtain an aerofoil field from `Airfoil.to_sdf()`. That method
returns `sdf.shapes.airfoil.AirfoilSDF`. Do not add a second contour type.

## Seams

`SDFGrid` is an evaluation grid (sample, gradient, debug contour). Area and
moments of solids belong to [integration](../integration/README.md).
`SDFGrid.section_properties_fused` is debug, not a product (same policy as
`quadtree_mesh`).

`GeometrySolver` stays in `sdf.solve`. [extract](../extract/README.md) is a
higher package: polylines from zero sets and feature ids. It may import `sdf`;
`sdf` must not import `extract`.

[structure](../structure/README.md) assembles member fields from `sdf` nodes.

## Distance

Let \(\Omega\subset\mathbb{R}^2\) be a closed region with boundary \(\partial\Omega\).
The signed distance at a point \(\mathbf{x}\) is

$$
d(\mathbf{x})
= \mathrm{sign}(\mathbf{x};\Omega)\,
  \inf_{\mathbf{y}\in\partial\Omega}\lVert\mathbf{x}-\mathbf{y}\rVert,
\tag{1}
$$

where the sign is negative in the interior of \(\Omega\) (the kernel’s solid
convention) and positive in the exterior. Primitive fields (disk, box,
half-plane, segment, capsule, polygon) implement \((1)\) exactly or to the
polygon polyline.

## CSG

For fields \(d_A\) and \(d_B\),

$$
\begin{aligned}
d_{A\cup B} &= \min(d_A,d_B),\\
d_{A\cap B} &= \max(d_A,d_B),\\
d_{A\setminus B} &= \max(d_A,-d_B),\\
d_{\mathrm{offset}(\delta)} &= d - \delta.
\end{aligned}
\tag{2}
$$

Smooth blends replace \(\min\) and \(\max\) with \(C^1\) envelopes. An
expression tree (`Expr`, `Union`, …) is lowered to Numba kernels.

```mermaid
flowchart LR
  Prim[Primitives] --> CSG[Nodes]
  CSG --> Eval[NumbaKernels]
  Eval --> Sample[FieldSample]
  Sample --> Solver[GeometrySolver]
```

A two-and-a-half-dimensional feature field (`FieldSample`, `FeatureTable`,
`NO_FEATURE`, `TagNode`, `PriorityUnionNode`, `LoftedPolygonSDF`) tags walls
without becoming a Euclidean \(F(x,y,z)\). `LoftedPolygonSDF` caches a
uniform-grid AABB per station (`accelerate=True`, `accel_at` / `bvh`). Kernel barrel names:
[glossary](../glossary.md#sdf).

The page [section/sdf](../section/sdf/README.md) is a redirect: SDF is not a
section submodule.

## Example

[examples/sdf/01_signed_distance_csg](../../examples/sdf/01_signed_distance_csg/)
constructs a box beam from primitives. Gallery scripts live under
`examples/demos/section/`.

## See also

- Quadtree integrals: [integration](../integration/README.md)
- Member assembly on fields: [structure](../structure/README.md)
- Zero-level extraction: [extract](../extract/README.md)
