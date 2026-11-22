# 0017. Surface derivatives and curvature are analytic, not finite differences

**Decision.** With smooth span laws ([0014](0014-smooth-authored-span-laws.md))
and one (t, z) patch per partition ([0015](0015-partition-mesher.md)), every
partition mid-surface is a smooth map x(t, z). The kernel evaluates its first and
second derivatives analytically, and from them the normal, the principal
curvatures and directions, and the gradient and Hessian of the signed distance
near the surface. Section curves get a smooth representation (a C2 spline through
given points, or the CST formula) instead of a polyline, with the same sampling
check as span laws. Junction nodes are put on both surfaces by a Newton solve on
the two stacked distance fields, and the check that they are on both is a mesh
invariant. Each patch is tagged FLAT, SINGLY, DOUBLY or UMBILIC by construction
where its construction decides it (a planar web is FLAT); such tags carry exact
zero curvature, and computed curvature must agree with the tag within a stated
tolerance. Issue #98.

**Why.** Earlier choices (polygon section distance, finite-difference gradients
in `sdf/solve.py`, "exact zeros are not meaningful") rested on a linear loft of
polylines, whose second derivatives are zero inside a bay and undefined at
stations. On smooth patches the derivatives exist everywhere, so finite
differences only add error and cost.

**Follows.** Polygon distance stays for section integrals and as a cross-check.
`sdf/` stays self-contained: the analytic evaluators live with the patch, and
`sdf/` receives callables. Curvature on the mesh is sampled from the patch, never
fitted from mesh nodes.
