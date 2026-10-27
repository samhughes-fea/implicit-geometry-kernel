# Dual-rep acceptance — polygon vs SDF area

Two integration paths for region area / second moments:

| Path | Method |
|------|--------|
| **SDF / sampling** | Grid or quadtree occupancy of `φ ≤ 0` |
| **Constructive polygon** | Explicit boundary from the same recipe → shoelace |

SDF defines the solid. Mature labels (skin, webs, caps, cores, LE/TE inserts)
use the polygon path for integration by default. Residual φ remains available
when an insert recipe opts into `emission_maturity='residual_phi'`. This page
is the numeric parity check (polygon ≠ exact SDF zero set).

↑ [`section_geometry`](../../README.md) ·
[`../constructive/`](../../constructive/README.md) ·
[`../props_emit/`](../../props_emit/README.md) ·
[`../recipes/`](../../recipes/README.md)

## Table of contents

1. [Proof figures](#proof-figures)
2. [Label split (hybrid default)](#label-split-hybrid-default)
3. [Acceptance claims](#acceptance-claims)
4. [Caveats](#caveats)

---

## Proof figures

| Claim | Figure | Demo |
|-------|--------|------|
| Polygon boundary vs SDF contour gap is measurable | [`../visualisation/output/dual_rep_consistency/dual_rep_boundary_overlay.png`](../../visualisation/output/dual_rep_consistency/dual_rep_boundary_overlay.png) | [`dual_rep_consistency_demo.py`](../../visualisation/demos/dual_rep_consistency_demo.py) |
| LE/TE inserts on constructive path (opt residual for debug) | [`../visualisation/output/le_trailing_edge_insert_residual/inserts_hybrid_quadtree.png`](../../visualisation/output/le_trailing_edge_insert_residual/inserts_hybrid_quadtree.png) | [`le_trailing_edge_insert_residual_demo.py`](../../visualisation/demos/le_te_insert_residual_demo.py) |

Metrics JSON (regenerate with the demo):
[`../visualisation/output/dual_rep_consistency/dual_rep_consistency.json`](../../visualisation/output/dual_rep_consistency/dual_rep_consistency.json).

## Label split (hybrid default)

| Class | Examples | Path |
|-------|----------|------|
| Constructive | `skin_*_surface`, `shear_web_*`, `spar_cap_*_surface`, `sandwich_core_*`, `leading_edge_insert`, `trailing_edge_insert` | Polygon → shoelace |
| Residual | (empty by default; unknown labels) | Shared adaptive quadtree over the SDF |

Backend: `SECTION_GEOMETRY_PROPS_BACKEND=hybrid|quadtree|grid` (default **hybrid**).
`INSERTS_CONSTRUCTIVE_BY_DEFAULT` is **True**.

## Acceptance claims

| Check | Criterion | Evidence |
|-------|-----------|----------|
| Boundary gap | Max polygon↔SDF-contour gap | `dual_rep_consistency.json` → `max_boundary_gap_m` |
| Area parity (constructive labels) | Relative area error vs dense-grid SDF occupancy | same JSON `wall_solid_boundary_metrics.*.area_rel_err` |
| Insert constructive | Inserts on constructive labels under hybrid | `constructive_labels` includes LE/TE when present |
| Twist cases | Twisted solids in props / dual-rep demos | props dashboard + dual-rep demo |
| Overlap policy | Per-label moments independent (skin∩cap may double-count) | [`../constructive/README.md`](../../constructive/README.md) |

Exact rtol/atol for CI live in tests. Perf history:
[`PERFORMANCE_WORK_SUMMARY.md`](PERFORMANCE_WORK_SUMMARY.md).

<a id="caveats"></a>
## Caveats

- Constructive polygons are not Euclidean-exact SDF zero sets.
- Midlines are centrelines, not the solid boundary.
- No unique-material partition of the aerofoil disk.
