# bug_log — sdf and extract

**Audience: internal.** Maintainer notes; not part of the library manual.

↑ [`sdf`](README.md) · [R4](../r/R4_sdf_csg.md)

## Role

Package debug log for `implicit_geometry_kernel.sdf` and the zero-set extractors in
`implicit_geometry_kernel.extract`. Identifiers are `SD-nnn`. Cross-link
[`section/bug_log.md`](../section/bug_log.md) (`SG-`),
[`blade/bug_log.md`](../blade/bug_log.md) (`BG-`) and
[`frames/bug_log.md`](../frames/bug_log.md) (`CF-`) when a defect spans packages.

## Index

| ID | Title | Status |
|---|---|---|
| [SD-001](#sd-001) | `sdf_ellipse` is not the distance to the ellipse | **Open** |
| [SD-002](#sd-002) | `expr_bbox` of an offset of a disjoint intersection misses the region | **Open** |
| [SD-003](#sd-003) | `extract_zero_loops` drops two of four crossings in a cell | **Open** |
| [SD-004](#sd-004) | `LoftedPolygonSDF(accelerate=True)` builds an edge grid it never reads | **Open** |

## SD-001

**Symptom.** `sdf_ellipse(2, 0, rx=1, ry=0.5)` returns 2.06; the point is 1.0 from
the ellipse. Over random points in $[-1.6, 1.6]^2$ the value is wrong by up to 1.1.
The sign is right.

**Evidence.** Found while deriving [R4](../r/R4_sdf_csg.md), (R4.10).
`examples/sdf/03_distance_and_csg` check 1 compares it with the distance to the
boundary sampled at $4\times10^5$ points; Figure R4.1 maps the error.
`tests/section/test_primitives.py` checks only the sign and the circle case.

**Root cause.** The docstring calls it Quilez's method; the loop is a fixed-point
iteration on a unit direction $(t_x, t_y)$ whose update does not converge to the
nearest point, and four steps are run regardless. The nearest point is the root of
(R4.10).

**Fix.** Open. Solve (R4.10) for $s$ by bisection on its bracket (monotone $F$),
with the axis cases in closed form, as Eberly (2013). No kernel path calls
`sdf_ellipse`; `Airfoil.ellipse` builds a polygon.

## SD-002

**Symptom.** `Offset(Intersect(Box(0, 0, 1, 1), Box(2.5, 0, 1, 1)), 0.5).bbox()`
is $x\in[1, 2]$, $y\in[-1.5, -0.5]$; the region $\{\phi \le 0\}$ spans
$y\in[-1.43, 1.43]$.

**Evidence.** [R4](../r/R4_sdf_csg.md), Section 3.10; example 03, check 7.

**Root cause.** `_bbox_intersect` returns a degenerate box at $(\max x_0, \max y_0)$
when the boxes do not overlap. The intersection is empty, but its field is a bound
outside (R4.12), so $\{\max(\phi_A, \phi_B) \le \delta\} = A_\delta\cap B_\delta$ is
not empty when the shapes are closer than $2\delta$; the offset rule then grows a
point, not the true extent.

**Fix.** Open. Keep the per-axis interval $[\min(l, u), \max(l, u)]$ of an empty
overlap instead of collapsing it, so that growing it by $\delta$ covers
$(\square_A\oplus\delta)\cap(\square_B\oplus\delta)$.

## SD-003

**Symptom.** `extract_zero_loops` of two circles of radius 0.5 whose boundaries pass
0.02 apart, on a $48\times48$ grid, returns 13 loops instead of 2.

**Evidence.** [R4](../r/R4_sdf_csg.md), Section 3.11; example 03, check 7. One circle
extracts as one closed loop on $\phi = 0$ to $3\times10^{-14}$.

**Root cause.** A cell with four sign changes on its edges has two segments. The
loop keeps `hits[0]` and `hits[1]` only, pairing the first two edges in edge order,
which drops one segment and can join the wrong pair. The chain then breaks at each
such cell.

**Fix.** Open. Emit both segments, paired by the asymptotic decider (R4.31).

## SD-004

**Symptom.** `LoftedPolygonSDF(accelerate=True)` is documented as a uniform-grid
BVH for batched queries, but its evaluation time equals the unaccelerated path.

**Evidence.** [R4](../r/R4_sdf_csg.md), Section 6. `build_polygon_accel` fills
`offsets` and `flat` (the edge buckets); `sdf_polygon_bvh` reads only the edge
arrays `ax`, `ay`, `bx`, `by` and loops over every edge.

**Root cause.** The bucket grid was built but no nearest-edge search was written
against it.

**Fix.** Open. Either search the buckets outward from the query cell, stopping when
the ring distance exceeds the best edge distance, or drop the grid and describe the
flag as caching the station's edge arrays.
