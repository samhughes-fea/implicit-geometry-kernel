# 0014. Authored span laws may be smooth (PCHIP); tables stay linear; thickness is a field t = f(s, z)

**Decision.** Span laws can be interpolated with a shape-preserving C1 cubic
(PCHIP) as well as linearly. PCHIP is opt-in now. Tabulated inputs (windIO, CSV,
station lists read from a file) stay linear. At a later, flagged 0.9.x release
PCHIP becomes the default for authored parametric blades, with the golden numbers
re-baselined in that release. Thickness generalises from a gauge per run plus
ply-drop ramps to a field t = f(s, z) on each partition, sampled per mesh node.
Issue #95. This supersedes the step 8 row "Loft scalars: linear only; PCHIP
removed" in the plan doc.

**Why.** PCHIP was removed because sparse tables invented curvature between
knots. An authored law is the design itself, so its curvature is real. A linear
loft has a curvature kink at every station, which clashes with curvature on the
mesh and with curvature-based sizing. A field over (s, z) is the natural form of
`PlyDrop` and `StackThickness` and gives each node its own thickness.

**Follows.** PCHIP needs denser stations than a linear loft, so the kernel checks
the sampling. A law interpolated with PCHIP must pass a leave-one-out test: each
interior knot, dropped and predicted from its neighbours, must agree with its
given value within a stated tolerance. Inputs that fail are refused with the
station named, or flagged by `validate()`, and the fix is more stations, not a
looser tolerance. Linear laws are not checked this way.
