# Decision records

Short records of choices the maintainer has made, so they are not reopened by
accident. Each says what was decided, why, and what follows. To change one, add a
new record that supersedes it; do not edit history.

| No. | Decision | Date |
|-----|----------|------|
| [0001](0001-geometry-only.md) | The kernel is geometry only | 2026-11-04 |
| [0002](0002-sdf-self-contained.md) | `sdf/` stays self-contained | 2026-11-04 |
| [0003](0003-layout-catalogue-front-door.md) | The layout catalogue is the front door; ribs are a spanwise axis | 2026-11-06 |
| [0004](0004-cap-on-skin-mid-surface.md) | Cap and skin share the skin mid-surface; eccentricity is applied downstream | 2026-11-10 |
| [0005](0005-kernel-shell-mesh.md) | The kernel ships the analysis shell mesh; consumers do not remesh | 2026-11-10 |
| [0006](0006-tip-closure-not-a-rib.md) | The tip closure is part of the skin shell, not a rib | 2026-11-10 |
| [0007](0007-rib-cells-as-triangles.md) | Rib and tip-closure cells are triangles | 2026-11-21 |
| [0008](0008-versioning.md) | Releases stay on 0.9.x; tags are the maintainer's | 2026-11-11 |
| [0009](0009-consumer-gaps-are-kernel-work.md) | Geometry a consumer computes on top of the kernel is kernel work | 2026-11-12 |
| [0010](0010-no-private-data.md) | No confidential blade data in this repository | 2026-11-05 |
| [0011](0011-tags.md) | Tags name releases only and are permanent | 2026-11-21 |
| [0012](0012-conventional-commits.md) | Commits follow Conventional Commits | 2026-11-21 |
| [0013](0013-release-workflow-tags.md) | The release workflow makes tags; supersedes "maintainer only" in 0008 and 0011 | 2026-11-22 |
| [0014](0014-smooth-authored-span-laws.md) | Authored span laws may be smooth (PCHIP) with a station-density check; tables stay linear; thickness is t = f(s, z) | 2026-11-22 |
| [0015](0015-partition-mesher.md) | The shell mesh is built per partition on a (t, z) patch from a shared skeleton | 2026-11-22 |
| [0016](0016-property-only-sub-partitions.md) | Sub-partitions step properties, never geometry | 2026-11-22 |
| [0017](0017-analytic-surface-derivatives.md) | Surface derivatives and curvature are analytic, not finite differences | 2026-11-22 |
