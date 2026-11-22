# 0016. Sub-partitions step properties, never geometry

**Decision.** Any run of the shell graph can be divided into authored
sub-partitions: spar cap, bond flange, web foot, wedge, promoted ply drop,
leading-edge protection, core extent. Each has its own thickness law and
eccentricity. Thickness and eccentricity step at sub-partition edges; nodes stay
on the one mid-surface. No per-partition reference surface (`partition_mid`) is
built. Issue #97.

**Why.** It generalises the cap partitions to every band a laminate changes over
and gives consumers one place to bind laminates. Stepped node geometry would need
connectors at every step for very little accuracy (see
[0004](0004-cap-on-skin-mid-surface.md)); the maintainer chose properties only on
2026-11-22.

**Follows.** [0004](0004-cap-on-skin-mid-surface.md) stands and is extended:
the shell graph keeps every sub-partition, not only caps. Spar caps are
re-expressed as sub-partitions first, with golden numbers unchanged. A wedge at
a non-normal web joint is one more sub-partition; its eccentricity is a vector
(normal, tangential) carried as data and applied downstream.
