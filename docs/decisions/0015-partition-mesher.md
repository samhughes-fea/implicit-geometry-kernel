# 0015. The shell mesh is built per partition on a (t, z) patch from a shared skeleton

**Decision.** The shell mesh moves from ruled strips between span stations to one
structured patch per partition, with t in [0, 1] across the partition and z along
the span. Partition boundaries, web junctions, free edges and terminations form a
tagged skeleton. Each skeleton curve is discretised once and both neighbouring
patches use its node ids. The new mesher lands beside the ruled-strip mesher
(`shell_mesh(method="partition")`), is compared with it in CI on every catalogue
layout, and becomes the default at a flagged 0.9.x release; the old mesher is
kept for one release after that. Issue #96.

**Why.** Shared node ids by construction replace welding by exact coordinate,
which today needs after-the-fact checks (`near_miss`, `t_junction`). Most of the
event constructors (step, transition, fan) disappear, because a partition keeps
its node count along its span and changes only at authored terminations. One
(t, z) domain carries the thickness field (0014), sizing, curvature and, later,
higher-order nodes.

**Follows.** The product in [0005](0005-kernel-shell-mesh.md) is unchanged: an
all-quad conforming mesh on skin and webs that consumers use as is, and
triangles on ribs and the tip closure ([0007](0007-rib-cells-as-triangles.md)).
Only 0005's "Follows" line about mapped strips between span stations is
superseded. Terminations inside a partition (cap and web ends) need structured
transition bands; the mesh golden numbers change at the default switch.
