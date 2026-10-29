# section_geometry performance work — status and next steps

> **Status:** Historical performance handoff. Paths that cite deleted
> `engine/` / `interface/` / `core/` trees are **stale** — those homes are gone
> (see [G_ROW_OWNERSHIP_MATRIX](../reference/G_ROW_OWNERSHIP_MATRIX.md)). Prefer
> live domain homes (`sdf/`, `sections/`, `props_emit/`, …) when following links.

Phase timings for computational-efficiency work on `section_geometry`. Short
handoff:
[section_geometry_efficiency_review.md](../../../../../../docs/architecture/section_geometry_efficiency_review.md).

↑ [`docs index`](../README.md) · [`section_geometry`](../../README.md)

## Table of contents

1. [Goal](#goal)
2. [Architecture arc](#architecture-arc)
3. [What was done](#what-was-done)
4. [Measured results](#measured-results)
5. [How to reproduce profiling](#how-to-reproduce-profiling)
6. [Test coverage](#test-coverage)
7. [Completed since Phase 4](#completed-since-phase-4)
8. [Per-stage breakdown (many-station, post-Phase-9)](#per-stage-breakdown-many-station-post-phase-9)
9. [Phase 10](#phase-10)
10. [Remaining headroom](#remaining-headroom)
11. [Phase 11](#phase-11)
12. [Phase 12](#phase-12)
13. [Phase 13](#phase-13)
14. [Phase 14](#phase-14)
15. [Phase 15](#phase-15)
16. [Phase 16](#phase-16)
17. [Recommended next steps](#recommended-next-steps)

---

## Goal

Reduce wall-clock time for section property computation without changing numerical parity. The dominant cost was repeated `sdf_polygon` evaluation across ~8–16 labelled SDF components per station, originally on a dense uniform grid.

The default property backend is now **hybrid** (`SECTION_GEOMETRY_PROPS_BACKEND=hybrid`):
constructive wall polygons + shoelace for skin/webs/caps, with residual φ
quadtree for `core_*`, `leading_edge_insert`, and `trailing_edge_insert`. Full-label adaptive
quadtree remains available via `SECTION_GEOMETRY_PROPS_BACKEND=quadtree`.
The canonical benchmark is depth-7, 2-web NACA 2412 via the 2D-CN system layout
(same setup as [`baseline_ir_consolidation.json`](../../tests/baseline_ir_consolidation.json)).

See also [DUAL_REP_ACCEPTANCE.md](DUAL_REP_ACCEPTANCE.md).

## Architecture arc

```mermaid
flowchart LR
    uniformGrid["Uniform grid properties"]
    sharedQT["Shared adaptive quadtree"]
    irConsolidation["CSG IR consolidation"]
    perfLevers["Perf levers batch eval"]
    uniformGrid --> sharedQT --> irConsolidation --> perfLevers
```

## What was done

### Phase 1 — Adaptive quadtree property backend

- One shared quadtree refined near the **union** of all subcomponent boundaries (skin, webs, caps, cores).
- Per-label section properties integrated over the shared leaf partition (interior exact rectangles + boundary marching-squares clip).
- Wired as the default backend through `compute_section_properties` and orchestration.

| Area | Key files |
|------|-----------|
| Multi-field refine / shared mesh | `blade_structure_utilities/mesh_utilities/area_2d/quadtree.py` |
| Quadrature / center φ attach | `blade_structure_utilities/mesh_utilities/area_2d/quadrature.py` |
| Section adapter | `blade_structure_utilities/geometry_utilities/section_geometry/engine/area_mesh.py` |
| Backend wiring | `blade_structure_utilities/geometry_utilities/section_geometry/interface/props_compute.py`, `engine/props_backend.py` |

### Phase 2 — Profiling and orchestration

- Optional counters via `SECTION_GEOMETRY_PROFILE=1` (`blade_structure_utilities/geometry_utilities/section_geometry/sdf/profiling.py`).
- Per-station `SectionEvalCache` in `section_geometry_impl` (`blade_structure_precompute (archived stages path)`).
- Two-tier grids: 384×180 for properties reference, 768×330 for plot PNGs.
- Optional station parallelism via `ProcessPoolExecutor` when `section_solve_n_workers > 1`.

### Phase 3 — CSG IR consolidation

- Single CSG expression graph (`Expr`) per component; callable `__call__` derives from `eval_expr_batch`.
- `eval_expr_batch` / `compile_expr` subtree memoization keyed by `Expr` nodes.
- Component package split under `blade_structure_utilities/geometry_utilities/section_geometry/engine/sections/components/`.
- `SectionFieldEvaluator` bridges section labels to the IR on the quadtree hot path.
- Grid cache path uses IR by default (no feature flag).

| Area | Key files |
|------|-----------|
| IR primitives and eval | `blade_structure_utilities/geometry_utilities/section_geometry/sdf/ir.py` |
| Section evaluator | `blade_structure_utilities/geometry_utilities/section_geometry/engine/section_evaluator.py` |
| Components | `blade_structure_utilities/geometry_utilities/section_geometry/engine/sections/components/` |

### Phase 9 — Struct-of-arrays quadtree (latest)

After Phase 8, profiling showed the build/refine was dominated by **Python per-cell overhead**
(per-`QuadCell` objects, per-label `phi_fields` dicts, the `_phi_by_label_from_leaves` listcomp,
and ~900k throwaway `np.full(4, np.nan)` `.get` defaults) against a ~0.4 s genuine IR-eval floor.

1. **Quick allocation fix** (`_phi_by_label_from_leaves`, `quadtree.py`): hoisted the per-element
   `np.full(4, np.nan)` `.get` default to a shared `_NAN4` sentinel. Bit-identical;
   build/refine ~1002 → ~807 ms on its own.

2. **Struct-of-arrays refinement** (`_refine_soa`, `quadtree.py`). `build_shared_quadtree` now
   carries the refinement frontier as arrays — `bounds (M,4)`, `level (M,)`, per-label corner
   φ `(M,4)` — and expands refined cells into `4R` children by **array reshaping** (child-corner
   φ gathered from the parent corners + interior probes via `_CHILD_CORNER_RECIPE`, per-child
   support sentinel vectorized). No per-cell `QuadCell` objects or per-cell Python loops in the
   hot loop. The mesh stores `_leaf_bounds` / `_leaf_levels` / `phi_by_label` directly;
   `leaf_cells()` materialises `QuadCell`s lazily only for object-iterating consumers (contour,
   `summarize` fallback, examples). The quadrature hot path (`section_properties_from_quadtree`,
   `_attach_center_phi`, `_accumulate_moments_vectorized`) and `summarize_quadtree` read the arrays
   directly. The object-model path (`build_square_root_forest` + `refine_quadtree_multi` +
   `_subdivide_refined_reuse`) is retained as the equivalence reference.

   **Parity (reorder, not bit-identical):** the SoA leaf *set* and `quadtree_leaf_count` are
   **identical** to the object model (verified exact bounds equality across 2D-CN, twisted 2D-CN,
   3C-F), but leaves are emitted in refinement order, so per-label properties differ by float
   summation round-off (~1e-10, well within the `rtol=1e-6` regression bar). `phi_eval_count`
   (1,028,728) is unchanged; `sdf_polygon_count` 16 → 9 and `sdf_polygon_points` shift (batching)
   and were rebaselined.

   Canonical 2D-CN depth-7 (single station, untwisted): build/refine ~807 → **~478 ms**;
   `properties_s` ~1.06 → **~0.70 s** (and ~1.25 → ~0.70 s vs the Phase-8 baseline). Many-station
   bench (twisted mix): properties ~817 → **~461 ms/station**; export (recomputes properties
   internally) ~1119 → **~807 ms/station**. Locked by
   `tests/test_soa_matches_object_model.py`.

### Phase 8 — Subdivision φ-reuse + export memo sharing

Two memo-sharing wins, both **bit-identical** (verified: leaf count, per-label properties,
and per-field grid φ unchanged; only evaluation counts/wall drop).

1. **Quadtree subdivision reuses decision-phase φ** (`mesh/area_2d/quadtree.py`). The refine
   decision already evaluates each eligible cell's 4 corners + 5 interior probes
   (`_cell_interior_points`: center + 4 edge midpoints) — which are *exactly* the corners of
   the 4 children created by subdivision. Previously `_subdivide_cell_multi` re-evaluated all
   16 child corners from scratch, and the multi-sampler decision path also evaluated parent
   corners that were already in the per-cell cache. Now `refine_quadtree_multi`:
   - takes parent corners from the per-cell `phi_fields` cache (no corner re-eval), and
   - subdivides via `_subdivide_refined_reuse`, which gathers child-corner φ from the cached
     corners + probes (`_CHILD_CORNER_RECIPE`) with per-child support sentinel masking,
     instead of re-sampling.

   `_subdivide_cell_multi` is retained as the bit-for-bit reference (locked by
   `blade_structure_utilities/mesh_utilities/tests/test_quadtree_subdivide_reuse.py`).

   Canonical 2D-CN depth-7 (single station, untwisted):

   | Metric | Before | After |
   |--------|-------:|------:|
   | `properties_s` (min of repeats) | ~5.18 | **~1.25** (**~4.2×**) |
   | `phi_eval_count` | 3,216,536 | **1,028,728** |
   | `sdf_polygon_count` | 12,570 | **16** |
   | `quadtree_leaf_count` | 37,650 | 37,650 (unchanged) |

   `baseline_analytic_bounds.json` rebaselined for `phi_eval_count` / `sdf_polygon_count`
   (`sdf_polygon_points` too); leaf count and per-label properties are unchanged.

2. **Export evaluates all contour fields through one shared IR memo** (`interface/export.py`,
   `_shared_grid_phi`). The per-component `grid.eval` was ~100% of export cost (~3.2 s for the
   8 components + 2 skin boundaries on the 384×180 grid; contour extraction is ~10 ms). Every
   field embeds the same airfoil polygon subtree, so a single shared `eval_expr_batch` memo
   computes it once. `SectionComponent.__call__` *is* `eval_expr_batch(self.expr, …)`, so φ is
   bit-identical to per-field `grid.eval`. Twisted skin boundaries use the host's
   `_skin_*_boundary_expr` IR nodes (bit-identical to the `rotate_field` callables, and sharing
   the rotated-airfoil polygon already computed for the twisted components). Contour grid eval
   ~3216 → ~300 ms (~10×); locked by `tests/test_export_shared_phi.py`.

   Many-station bench (`--stations 4 --export`, 2D-CN depth-7, twisted mix): export ~1603 →
   ~1119 ms/station (remainder is the property recompute inside `SectionPropertiesReport`,
   itself ~4× cheaper after win 1).

### Phase 7 — Plot rendering from the shared quadtree mesh

The per-station PNG (`persist_pngs=True`, on by default in bulk runs) is now rendered directly
from the property backend's already-built quadtree mesh instead of a second dense 768×330
`SDFGrid` + `contourf`/`contour`.

1. **Mesh plumbing** — `compute_section_properties(..., with_mesh=True)` returns the
   `QuadTreeMesh`; `SectionPropertiesReport.mesh` exposes it (kept off `backend_meta`/`to_json`).
2. **Mesh→geometry helpers** (`mesh/area_2d/quadrature.py`) — `label_fill_polygons` (interior-leaf
   rectangles + boundary-leaf marching-squares clips; summed area matches the integrated section
   area to round-off) and `label_boundary_segments` (φ=0 segments). Reuse the `_BOUNDARY_RECIPES`
   edge-crossing machinery.
3. **`plot_section_from_mesh`** (`interface/plot.py`) — per-label `PolyCollection` fills
   (`rasterized=True`) + `LineCollection` boundary, same palette/legend/aspect as `plot_section`.
   The grid `plot_section` is kept as the grid-backend fallback.
4. **Orchestration** (`stages.py`, both plot sites) uses the renderer when `report.mesh` is
   available, skipping the 768×330 grid build.

Result: **plot ~4940 → ~285 ms/station (~17×)**; bulk per-station wall roughly halves.
Watertight by construction (leaf partition tiles the plane), so no T-junction stitching. Visual
parity confirmed for 2D-CN and 3C-F (twist + flapwise webs + discrete caps). Covered by
`tests/test_plot_from_mesh.py` (area tie-out + structural checks).

### Phase 6 — Stabilization + taxonomy generalization

1. **Quadtree leaf-cache fix** — `refine_quadtree` / `refine_quadtree_multi`
   (`mesh/area_2d/quadtree.py`) now invalidate the leaf cache after the refinement loop
   mutates the tree, so callers that don't explicitly refresh (`build_area_mesh_2d`) no longer
   see stale pre-refinement roots. Repairs `test_contour_from_circle_is_closed` and
   `test_quadtree_refines_near_zero_contour`; unblocks quadtree-derived contour work.
2. **Quadtree backend generalized to all 40 system types** — `SectionFieldEvaluator.from_section`
   now handles sections without an IR component graph (the `zero_web_catalogue` `0A`/`0B` layouts):
   each labelled callable is used via its own `.expr` when present (tight bounds, untwisted) or
   wrapped as an opaque `CallableField` (⇒ `expr_bbox` None ⇒ sampling fallback) instead of
   raising `KeyError`. All 38 multicell families (X=1..5, Y=A/B/C/D, Z=CN/F) already derive
   analytic bounds with zero sampling. Validated in `tests/test_support_bounds_taxonomy.py`
   (support coverage, zero-sampling, quadtree-vs-grid parity across representative keys).

### Phase 5 — Analytic IR bounds + boundary vectorization

1. **Analytic IR support bounds** — `expr_bbox` in `blade_structure_utilities/geometry_utilities/section_geometry/sdf/ir.py` derives a
   conservative AABB of `{expr <= 0}` directly from the IR graph (Polygon vertex extent,
   Circle/Box/Capsule extents, Offset/Shell/Slab expansion, Rotate corner rotation,
   Union/Intersect/Subtract recursion; `HalfPlane`/`CallableField` → unbounded). `_component_supports`
   (`engine/area_mesh.py`) now uses it instead of the 64×64 coarse grid, falling back to sampling
   only for labels with an unbounded/opaque leaf (none in the standard sections).
   Support computation: **148 ms → 0.18 ms per station (~848×)**.
2. **Vectorized boundary moments** — `_accumulate_moments_vectorized` (`mesh/area_2d/quadrature.py`)
   now buckets boundary leaves by their 4-bit marching-squares pattern and batches the
   shoelace moments (`_boundary_moments_batch` / `_shoelace_moments_batch`) instead of a
   per-cell Python loop. Bit-identical to the per-cell path (rtol 1e-12).
3. **Export contour dedup** — `export_section_json` (`interface/export.py`) extracts each
   component's φ=0 boundary once and shares it between the nested `geometry` block and the
   legacy flat `component_zero_contours` block.

Parity note: analytic bounds are conservative supersets of the true region, so no boundary or
interior cell is ever culled. Because the per-label support also gates *shared-mesh refinement*,
the tighter analytic boxes shift section properties by ≤~0.2% — **toward** the high-resolution
grid ground truth (e.g. core centroid error 1.2e-3 → 3.6e-5). This was an accepted, rebaselined
change; see `tests/baseline_analytic_bounds.json` and `tests/test_support_bounds_regression.py`.

Regression guard + workload bench: `tests/test_support_bounds_regression.py` locks the counters
and per-label properties (and asserts zero support-sampling);
`visualisation/demos/manystation_bench.py` reports per-stage wall time across N
stations (build / properties / plot / export).

### Phase 4 — Performance levers

Addressed the post-IR wall-time regression and unrealized cross-label sharing:

1. **Cached `Polygon.__hash__`** — precomputed in `__post_init__` so memo `get`/`set` is O(1) instead of O(N_vertices).
2. **`eval_labels` / `make_multi_sampler`** — one fresh memo per point batch, all labels evaluated together.
3. **`multi_sampler` on quadtree hot path** — optional batch evaluator threaded through `build_shared_quadtree` → `refine_quadtree_multi` → `_subdivide_cell_multi` and `_attach_center_phi`; inactive cells masked with sentinel after eval.
4. **`union_expr` built once** in `SectionFieldEvaluator.__init__`.
5. Removed fragile pointer-keyed batch memo (`_begin_batch` / `_pts_key`); single-label `make_sampler` uses per-call memo.

## Measured results

Benchmark: 2D-CN, NACA 2412 (`n_points=120`), `max_depth=7`, `max_cells=150000`, untwisted. Source: [`baseline_ir_consolidation.json`](../../tests/baseline_ir_consolidation.json).

| Stage | `properties_s` | `sdf_polygon_count` | `phi_eval_count` | `leaf_count` |
|-------|----------------|---------------------|------------------|--------------|
| Pre-IR (closure baseline) | 54.13 | 456,300 | 2,091,697 | 38,004 |
| Post-IR (regression) | 84.81 | 214,845 | 3,036,194 | 41,649 |
| After perf levers | 6.44 | 12,696 | 3,301,224 | 38,004 |
| Analytic bounds + boundary vec | ~5.2 | 12,570 | 3,216,536 | 37,650 |
| Phase 8 subdivision φ-reuse | ~1.25 | 16 | 1,028,728 | 37,650 |
| Phase 9 struct-of-arrays | ~0.70 | 9 | 1,028,728 | 37,650 |
| **Phase 10 centre-φ reuse** | **~0.70** | **9** | **983,656** | **37,650** |

The Phase-10 row is the rebaselined authority (`tests/baseline_analytic_bounds.json`):
leaf count is unchanged; `phi_eval_count` / `sdf_polygon_points` drop from skipping
re-evaluated early-stop centres; `properties_s` is the min of repeated runs (machine
dependent). Per-label properties stay within the existing `rtol=1e-6` bar.

**Interpretation**

- IR consolidation cut `sdf_polygon_count` ~53% by sharing `Polygon` / `Offset` subtrees across labels, but regressed wall time because every memo operation re-hashed the full vertex tuple.
- Perf levers fixed hashing and realized cross-label airfoil evaluation on the quadtree refine / center-attach paths.
- `phi_eval_count` rose slightly because `multi_sampler` evaluates all labels on full point batches (including culled cells); the net win comes from far fewer `sdf_polygon` calls.
- Leaf count returned to the pre-IR parity baseline (38,004).

The quadtree **properties** path at default depth is now ~0.70 s/station (~77× faster than the
original 54 s closure baseline). After the Phase-9 struct-of-arrays rewrite the ~0.70 s splits
~478 ms build/refine / ~188 ms centre-φ attach / ~20 ms quadrature, and the build/refine is now
close to its genuine IR-eval floor (~0.4 s); the residual Python overhead has been removed.

## How to reproduce profiling

From the repository root (PowerShell):

```powershell
$env:SECTION_GEOMETRY_PROFILE="1"
python -c "
from blade_structure_utilities.geometry_utilities.section_geometry.api.factory import build_section_view, resolve_system_type
from blade_structure_utilities.geometry_utilities.section_geometry.props_emit.props_compute import compute_section_properties
from blade_structure_utilities.geometry_utilities.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.airfoil import AirfoilSDF
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.profiling import reset_stats

airfoil = AirfoilSDF.from_naca('2412', n_points=120, chord=1.0)
section = build_section_view(airfoil, resolve_system_type('2D-CN'), twist_angle_rad=0.0)
spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=150000)
reset_stats()
props, meta = compute_section_properties(section, backend='quadtree', quadtree_spec=spec, airfoil_for_spec=airfoil)
print(meta['profiling'])
"
```

Depth / cost trade-offs: `blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/quadtree_convergence_study.py`.

Environment knobs:

| Variable | Default | Purpose |
|----------|---------|---------|
| `SECTION_GEOMETRY_PROPS_BACKEND` | `hybrid` | `hybrid` \| `quadtree` \| `grid` |
| `SECTION_GEOMETRY_QT_DEPTH` | `7` | Quadtree `max_depth` when not overridden |
| `SECTION_GEOMETRY_PROFILE` | off | Enable `sdf_polygon_count`, `properties_s`, etc. |

### Phase 10 — Dual-rep hybrid properties

- Recipe-faithful constructive polygons for laminates (`engine/constructive_polygons/`).
- Residual φ path retained `LEInsert` / `TEInsert` (first-class; less mature).
- `compute_section_properties(..., backend='hybrid')` merges shoelace + residual-label quadtree.
- Spike: `visualisation/demos/constructive_vs_quadtree_spike.py`; bench:
  `visualisation/demos/manystation_bench.py --compare`.
- Parity: `tests/test_hybrid_parity.py`, `tests/test_constructive_polygons.py`.

### Phase 11 — Constructive sandwich cores

- `core_*` moved onto constructive polygons: inner cavity − laminate union ∩ bay clip
  (chordwise or flapwise), matching SDF `SandwichCore`.
- Hybrid with no inserts skips residual φ entirely (`residual_quadtree_s ≈ 0`,
  `phi_eval_count = 0` for properties).
- LE/TE inserts remain residual φ when enabled.
- Measured (2D-CN depth-7, no inserts): hybrid properties ~6–7 ms/station vs
  ~94–136 ms full quadtree (**~15–19×**); constructive emit+shoelace ~8 ms.

## Test coverage

Validation suite at handoff (72 tests):

```text
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_csg_ir.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_multicell.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_props_backend_parity.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_eval_cache.py
blade_structure_utilities/mesh_utilities/tests/test_quadtree_quadrature.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_csg.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_structural_family.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_section_evaluator.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_expr_bbox.py
blade_structure_utilities/geometry_utilities/section_geometry/tests/test_support_bounds_regression.py
```

Phase 5 additions: `test_expr_bbox.py` (analytic bound primitives + conservative-superset
property check), `test_support_bounds_regression.py` (zero support-sampling + counter/property
guard vs `baseline_analytic_bounds.json`), and `test_boundary_moments_batch_matches_per_cell`
in `blade_structure_utilities/mesh_utilities/tests/test_quadtree_quadrature.py`.

Phase 8 additions: `blade_structure_utilities/mesh_utilities/tests/test_quadtree_subdivide_reuse.py`
(`_subdivide_refined_reuse` bit-identical to the `_subdivide_cell_multi` reference, with and
without supports) and `blade_structure_utilities/geometry_utilities/section_geometry/tests/test_export_shared_phi.py`
(`_shared_grid_phi` bit-identical to per-field `grid.eval`, untwisted and twisted).

Phase 9 additions: `blade_structure_utilities/geometry_utilities/section_geometry/tests/test_soa_matches_object_model.py`
(SoA leaf set + bounds exactly equal to the object-model reference; properties within rtol 1e-7,
untwisted and twisted).

## Completed since Phase 4

| Item | Status |
|------|--------|
| **B1 Cross-consumer eval cache** | Done — `plot_section` / `SectionPropertiesReport` accept `eval_cache` (`…/section_geometry/sdf/cache.py`). Note: keyed by `(id(grid), label)`, so the 768×330 plot grid and 384×180 props grid do not share entries. |
| **C1 Station parallelism** | Done — production uses `multiprocessing.Pool` in `orchestration/section_selection.py` via `BLADE_PARALLEL` / `BLADE_PROCESSES` (default `min(8, cpu_count()-1)`). The archived `section_solve_n_workers` / `stages.py` path is retired. |
| **A3 Fused grid properties** | Done — `SDFGrid.section_properties_fused` (single pass) in `blade_structure_utilities/geometry_utilities/section_geometry/sdf/grid.py`. |
| **A4 Contour backend** | Done — `SDFGrid.level_set` uses `contourpy.contour_generator`, matplotlib only as fallback. |
| **A5 Twist pre-transform** | Largely done — twist is a single `Rotate` IR node (`multicell.py`, `twist_component`); per-component `rotate_field` only wraps the skin-boundary helper SDFs. |
| **D Analytic IR support bounds** | Done — see Phase 5. |

## Per-stage breakdown (many-station, post-Phase-9)

`visualisation/demos/manystation_bench.py --stations 4 --plot --export` (2D-CN,
depth 7, twisted mix), ms/station:

| Stage | pre-Phase-8 | post-Phase-8 | post-Phase-9 | Notes |
|-------|------------:|-------------:|-------------:|-------|
| build | ~4 | ~4 | ~5 | section construction (negligible) |
| properties | ~4020 | ~820 | **~461** | quadtree backend; Phase-8 φ-reuse + Phase-9 SoA refinement |
| plot | ~285 | ~285 | ~279 | already render-from-mesh (Phase 7) |
| export | ~4480 | ~1120 | **~807** | Phase-8 shared-memo contour eval (~3216→~300 ms); remainder is the property recompute inside `SectionPropertiesReport`, now SoA |

The bulk per-station wall is now ~0.75 s (properties + plot). The build/refine is near its IR-eval
floor (~0.4 s); the dominant remaining levers are elsewhere (export property recompute, depth 8+).

<a id="phase-10"></a>
## Phase 10 — Centre-φ reuse + production mesh retain

1. **Centre-φ stash** — SoA `_refine_soa` (and object `refine_quadtree_multi`) stash
   interior probe[0] (cell centre) for early-stop leaves on `mesh._leaf_center_phi`.
   `_attach_center_phi` reuses finite stash entries and only evaluates NaN (typically
   max-depth interface leaves). Canonical 2D-CN depth-7: `phi_eval_count`
   1,028,728 → **983,656**; `sdf_polygon_points` 128,591 → **122,957**. Locked by
   `blade_structure_utilities/mesh_utilities/tests/test_center_phi_reuse.py` and rebaselined
   `baseline_analytic_bounds.json`.

2. **Production mesh retain** — `evaluate_station(..., retain_mesh=True)` returns the
   shared quadtree on `ShellStationEvaluateResult.quadtree_mesh`. Station sizing does a
   final property-only retain pass and stores `BladeSectionPropertyRow.quadtree_mesh`.
   Pool workers strip the mesh before IPC; cross-section viz prefers the retained mesh,
   re-attaches via `evaluate_station(retain_mesh=True)` on the parent when missing, and
   only then falls back to a local quadtree rebuild. Export examples pass
   `SectionPropertiesReport.from_precomputed` / `properties_report=`.

## Remaining headroom

### Medium leverage

| Item | Notes |
|------|-------|
| **A8 Medial filtering** | Removed — grid medial path deleted; constructive midlines are canonical. |

### Lower / conditional / examined-and-declined

- **Per-label quadrature batching** — `section_properties_all_labels` loops Python over labels,
  but measured at only ~20 ms (~1.6%) of the properties path; not worth the multi-label
  classifier complexity.
- **Airfoil / quadtree-spec memoization, web-anchor batching** — target the `build` stage
  (~4 ms/station); negligible, and `lru_cache` does not cross Pool workers.
- **float32 mesh storage** — memory only; the bottleneck is Python object overhead, not array
  size, so little wall-time benefit.
- **Eval-cache key `(id(grid), label)`** — leave as is. A geometry-hash key would cause *false*
  cross-station hits (same grid bbox/resolution, different section callables) ⇒ wrong φ;
  `id(grid)` is correct for the per-station design.
- **`multi_sampler` heuristic** — per-label gather fallback when supports are tiny.
- **Depth 8+** — cost grows sharply; tune `SECTION_GEOMETRY_QT_DEPTH` against parity budget.
- ~~**Export property recompute**~~ — API + production viz/export wiring done (Phase 9/10).
- ~~**Centre-φ attach**~~ — Done (Phase 10).
- ~~**Numba `sdf_polygon`**~~ — Done (Phase 11); now a **required** core dependency with no NumPy production fallback.

<a id="phase-11"></a>
## Phase 11 — Required Numba `sdf_polygon`

Profiling with `sdf_polygon_s` showed the polygon kernel was **~86%** of properties wall
on 2D-CN depth-7 without Numba (~0.58 s of ~0.68 s). An `@njit(cache=True)` edge-loop in
[`sdf/_polygon_numba.py`](../../sdf/_polygon_numba.py) is now the **sole** production path
(`numba>=0.59` is a core package dependency). `sdf_polygon_numpy` remains for parity tests
only. Parity: `blade_structure_utilities/tests/test_sdf_polygon_parity.py`.

<a id="phase-12"></a>
## Phase 12 — Topology substages + vectorized ring/ray kernels (E1–E3)

**E1 timers** on `build_shell_strip_topology` (when `SECTION_GEOMETRY_PROFILE=1`):
`topology_split` / `topology_cap_ray` / `topology_cap_web`, folded into worker `stage_timings`.

**E2 vectorization:** `nearest_on_closed_polyline_many`, NumPy edge loops in
`point_to_open_polyline_nearest` / `_ray_hits_open_polyline` / `skin_ring_vertex_tangents`.

**E3 DV reuse:** `StationLayoutCache` memos `resolve_station_layout` across spar/web
thickness trials; `GeometryResultCache` already skips midlines+topology when thicknesses
are unchanged (laminate/FI-only trials).

**Baseline (2D-CN, NACA2412 chord=1, midlines+properties, Numba on, post-E2):**

| Stage | s |
|-------|--:|
| assemble | 0.039 |
| properties | 0.007 |
| midlines | 0.001 |
| topology (total) | 0.012 |
| — topology_split | 0.001 |
| — topology_cap_ray | 0.003 |
| — topology_cap_web | 0.008 |

Post-E2, **assemble** dominates; topology is no longer the G7 ~50 ms hotspot.
**E4 (Shapely cavity/union reuse) skipped** — props/`constructive_s` are not the top stages.

<a id="phase-13"></a>
## Phase 13 — Assemble substages + batched cap feasibility (A1–A2)

**A1 timers** on `MultiCellSection` (via worker `stage_sink`):
`assemble_skin` / `assemble_anchors` / `assemble_webs` / `assemble_caps` /
`assemble_cores`.

**A1 blame (2D-CN, web_y pre-resolved):** `assemble_caps` ≈ 100% of assemble
(~27–30 ms); other assemble substages ≪ 1 ms.

**A2:** `cap_feasibility_clamp` now uses batched `inner_y_at_xs` with
`bisection_iters=0` (linear grid interp is enough for depth/camber clamp samples).
Family-D flapwise strip endpoints also batch via `inner_y_at_xs`.

**Baseline (2D-CN, NACA2412 chord=1, midlines+properties, Numba on, post-A2):**

| Stage | s |
|-------|--:|
| assemble | 0.005 |
| — assemble_caps | 0.004 |
| — assemble_skin/anchors/webs/cores | ~0 |
| properties | 0.008 |
| midlines | 0.001 |
| topology (total) | 0.008 |
| — topology_cap_web | 0.005 |
| — topology_cap_ray | 0.002 |
| — topology_split | 0.001 |

Assemble is no longer the wall (~8× down from Phase 12). Topology and properties
are now co-leading (~8 ms each). **Shapely cavity/`unary_union` reuse still not
justified** (props remain constructive/hybrid shoelace, not a multi-union spike).

<a id="phase-14"></a>
## Phase 14 — Cap–web topology vectorize + props stage keys (T1–T2)

**T1:** batched edge×edge intersections via `intersect_open_polylines_2d` in
[`ring_geometry.py`](../../midlines/ring_geometry.py), used by
[`insert_cap_web_interior_junctions`](../../topology/cap_web_junctions.py).
`topology_cap_web` ~0.004 → ~0.0005 s on 2D-CN.

**T2:** worker `stage_timings` now includes `properties_constructive` /
`properties_residual` folded from hybrid props meta.

**Baseline (2D-CN, NACA2412 chord=1, midlines+properties, Numba on, post-T1/T2):**

| Stage | s |
|-------|--:|
| properties | 0.006 |
| — properties_constructive | 0.006 |
| — properties_residual | ~0 |
| assemble | 0.004 |
| — assemble_caps | 0.004 |
| topology (total) | 0.003 |
| — topology_cap_ray | 0.002 |
| — topology_split | 0.001 |
| — topology_cap_web | 0.0005 |
| midlines | 0.001 |

**Properties constructive now leads** the midlines+properties path. Shapely
cavity/`unary_union` reuse (former E4) is **justified as the next efficiency
wave** if sizing remains wall-bound.

<a id="phase-15"></a>
## Phase 15 — Constructive emit substages + shared cavity / unions (S1–S2)

**S1 timers** on `emit_constructive_polygons` folded into worker `stage_timings`:
`emit_laminates` / `emit_unions` / `emit_cavity` / `emit_cores`.

**S2:**
- Single `unary_union(web_polys)` for both upper/lower cap differences
- One shared airfoil polygon + inner-cavity buffer reused by skin / webs / caps / cores
- One cavity−laminates difference shared across all core bays
- Same patterns on `emit_laminate_polygons` / `emit_core_polygons`

**Baseline (2D-CN, NACA2412 chord=1, midlines+properties, Numba on, post-S2):**

| Stage | s |
|-------|--:|
| properties | 0.007 |
| — properties_constructive | 0.006 |
| — emit_laminates | 0.002 |
| — emit_unions | 0.001 |
| — emit_cores | 0.001 |
| — emit_cavity | 0.0004 |
| — properties_residual | ~0 |
| assemble | 0.004 |
| topology | 0.003 |

Emit Booleans are no longer the bulk of `properties_constructive` (~4–5 ms emit
substages vs ~6 ms total including shoelace/filter). **Next candidate if sizing
still wall-bound:** laminate polygon construction cost inside `emit_laminates`
and/or shoelace in `section_properties_from_polygons` — not further union collapse.
Shapely remains.

<a id="phase-16"></a>
## Phase 16 — Shoelace vs laminate split + shared deeper offset (P1–P2)

**P1:** `constructive_s` split into `properties_emit` / `properties_filter` /
`properties_shoelace` (folded into worker `stage_timings`).

**P1 blame (2D-CN):** `emit_laminates` ~2.0 ms > `properties_shoelace` ~1.3 ms.

**P2:** share deeper airfoil offset (`buffer(-(skin+cap_h))`) across upper/lower
caps (and discrete bands on one surface) via `deeper_cavity` on
`emit_spar_cap_polygon`.

**Baseline (2D-CN, NACA2412 chord=1, midlines+properties, Numba on, post-P2):**

| Stage | s |
|-------|--:|
| properties_constructive | 0.006 |
| — properties_emit | 0.004 |
| — emit_laminates | 0.002 |
| — emit_unions / cores / cavity | ~0.002 |
| — properties_shoelace | 0.001 |
| assemble | 0.004 |
| topology | 0.003 |
| midlines | 0.001 |
| **main stages sum** | **~0.014** |

Further micro-opts on laminate stadiums/camber clips are sub-ms. **Efficiency
track parked** pending evidence that BSP sizing is still wall-bound on geometry
(cache misses). Prefer product/API/docs cleanup over more emit micro-opts.

## Recommended next steps

1. ~~Attack the export grid-backend cost~~ **Done (Phase 9).**
2. ~~Centre-φ reuse + production mesh retain~~ **Done (Phase 10).**
3. ~~Numba `sdf_polygon` (required)~~ **Done (Phase 11).**
4. ~~Topology substages + vectorize + layout cache~~ **Done (Phase 12).**
5. ~~Assemble substages + batched cap feasibility~~ **Done (Phase 13).**
6. ~~Cap–web topology vectorize + props stage keys~~ **Done (Phase 14).**
7. ~~Constructive emit substages + shared cavity/unions~~ **Done (Phase 15).**
8. ~~Shoelace vs laminate split + shared deeper offset~~ **Done (Phase 16).**
9. ~~Docs closeout (efficiency review index + README profile note)~~ **Done.**
10. **Efficiency track parked.** Ship Phases 12–16 as conventional-commit PRs when
    asked (split listed in
    [section_geometry_efficiency_review.md](../../../../../../docs/architecture/section_geometry_efficiency_review.md)).
11. Resume geometry perf only with BSP evidence that cache misses / rebuilds
    dominate sizing wall time — then profile hit rate, not emit Booleans.
