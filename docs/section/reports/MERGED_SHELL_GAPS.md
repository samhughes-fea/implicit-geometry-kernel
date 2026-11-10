# Merged shell graph — Pass 2 gap inventory

**Audience: internal.** Maintainer notes; not part of the library manual.

Evidence for G1–G11 and extract/TE/BVH. Waves 0–7 are shipped. G10/G11 are
closed for this kernel. Remaining park: TE classifier on OML only (WP6); BVH is
shipped. Geometry still required before a conformal shell-FE consumer can
delete its own kernel is a **separate** list (G-TE … G-ACC):
[`TIDALSHELL_GEOMETRY_REPLACEMENT_GAPS.md`](TIDALSHELL_GEOMETRY_REPLACEMENT_GAPS.md).
Characterization: [`tests/test_shell_graph_gaps.py`](../../../tests/test_shell_graph_gaps.py).

Fixture unless noted: 2D-F, `z = 4` m, NACA 0024→0018, chord 1.6 m, skin 3 mm, cap 8 mm, webs 2 mm at chord fractions `(0.18, 0.48)`.

↑ [`geometry_to_mesh`](../guides/geometry_to_mesh.md)

## Summary

| ID | Status | One-line hole |
|---|---|---|
| **G1** | Wave 3 shipped | `shell_strips_from_graph` = 8 runs; `cut.walls` still 10 (solid caps) |
| **G2** | Wave 1 shipped | `partition_edges` / `partition_edge_curve`; colouring unchanged |
| **G3** | Wave 2 shipped | `topology_diff` + `shared_loft_names`; no label lerp |
| **G4** | locked / closed | `T` step kept (8 mm); no geometric taper; Newton is not a spec at partition vertices |
| **G5** | Wave 0 shipped | Loft from graph; snap with SDF; joint `|φ| = 4.000 mm` is out of Pass 1 check |
| **G6** | Wave 4 shipped | `with_web_feet_on_loop`; loft re-resolves / drops webs on flapwise re-anchor |
| **G7** | Wave 5 shipped | `section_at` LRU 256 per `round(z, 9)`; `LoftedPolygonSDF` BVH on (`accelerate`, `accel_at`) |
| **G8** | Wave 6 shipped | `web_positions_from_windio_structure`; example 09 no longer hardcodes `(0.22, 0.50)` |
| **G9** | Wave 7 shipped | Cores / TE inserts stay constructive solids; not on `shell_graph` |
| **G10** | closed for this kernel | No `thin_wall` module; consumers use `mid_surface` / `shell_graph` |
| **G11** | closed for this kernel | No PCOMP / ply stagger; bind laminates via `offset_regions` / stacked `layers` |
| **WP6** | parked | named `feature_curve(geo, name_a, name_b, zs)` exists; TE classifier stays an OML helper; OML BVH shipped |

Implementation sequence is [Waves](#waves) below. **G4 is locked:** keep the `T` step; document Newton limits at partition vertices (no geometric taper).

Oracles (Wave 0): **loft strips** from `shell_graph` polylines; **in-plane snap** from `wall_mid_sdf` + `GeometrySolver`. Joint vertices are out of the Pass 1 `|φ| < 5e-4` check.

---

## G1 — One shell hand-off DTO

**Status:** Wave 3 shipped.

**Today:** [`build_shell_midline_strips`](../../../src/implicit_geometry_kernel/section/midlines/export.py) is skin+web (conformal mesh). Cap solids stay on [`cut.walls`](../../../src/implicit_geometry_kernel/blade/cut.py) with bonds; they are not a second mesh polyline. [`shell_strips_from_graph`](../../../src/implicit_geometry_kernel/blade/shell_graph.py) is the shell producer: one `ShellMidlineStrip` per graph run (`T = t_skin` on the loop, cap as partition `kind`, OFFSET on `offset_regions()`). [`ShellMidlineStrip`](../../../src/implicit_geometry_kernel/contract/shell_midline_strip.py) is **not** versioned.

**Evidence (2D-F, z=4):** 10 `SectionWall` rows vs 8 shell strips / `ShellRun`s.

**Out of repo:** adapters that counted cap rows on `result.midline_strips` should use `cut.walls` (solids) or `shell_strips_from_graph` (mesh).

**Non-goals:** do not version `ShellMidlineStrip`; do not delete caps from `cut.walls`.

---

## G2 — Feature lines at partition edges

**Status:** Wave 1 shipped (colouring in Pass 1; named edges on `ShellGraph`).

**Today:** `ShellRun.kind` is `cap` on stacked partitions; example 09 and `wall_meta` colour by that. Partition **edges** are only polyline endpoints. [`feature_curve`](../../../src/implicit_geometry_kernel/extract/curves.py) is `(field_a, field_b, zs, guess_xy)` with no names. [`station_network_from_meta`](../../../src/implicit_geometry_kernel/blade/topology/network.py) splits the mid loop by meta id and sets `thickness=0.0`.

**API:** `ShellGraph.partition_edges()` / `.edge(name_a, name_b)` and
`partition_edge_curve(geo, name_a, name_b, zs)` (omit if missing; no lerp).
`feature_curve` stays two SDFs + a guess. `station_network_from_meta` thickness
stays `0.0`.

**Non-goals:** do not change example 09 cap/bay colours.

---

## G3 — Spanwise topology events

**Status:** Wave 2 shipped.

**Today:** example 09 [`_wall_faces`](../../../examples/blade/08_reference_blades/reference_blades.py) lofts `a.keys() & b.keys()`. Missing names are omitted. No kernel `TopologyEvent`. Fixture `dropped_webs` at z=4 is empty. Cap runout / bay split would lerp silently if someone interpolated by index instead of name; the example does not lerp.

**API:** `topology_diff(a, b)` emits `cap_runout` / `dropped_web` / `bay_split`
(and onset / added_web / bay_merge). `shared_loft_names` is the loft rule.
Example 09 uses it. No silent label lerp.

**Non-goals:** no silent label lerp.

---

## G4 — Thickness jumps

**Status:** locked / closed (keep the step; Newton limits documented; no geometric taper).

**Today:** [`StackThickness`](../../../src/implicit_geometry_kernel/blade/topology/fields.py) assigns nearest loop-run `T` (step, not ramp). Not ply stagger.

**Evidence:** interior cap `T = 11` mm, interior skin `T = 3` mm, jump **8 mm** (`t_cap`). The two endpoints of one cap run evaluate to `{3 mm, 11 mm}` — nearest-run at the shared vertex is discrete. `GeometrySolver` default tolerance is `1e-6` m (`snap_wall_mid` uses `1e-8`); the jump is 8000× that.

**Decision (Wave 0, locked):** keep the nearest-run `T` step. Document Newton limits at partition vertices. No geometric taper.

**Non-goals:** ply-drop stagger; taper of `T`.

---

## G5 — SDF vs `mid_surface()` residual

**Status:** Wave 0 shipped (Pass 1 only checked midpoints).

**Today:** both oracles sit on the **skin mid**. `shell_graph` no longer moves cap-footprint vertices toward the stack NA. Cap-end `|φ|` matches run interiors. Constitutive eccentricity is `offset_regions()` (`t_cap/2` inward), not a geometric kink.

**Policy:** loft / conformal strips use `shell_graph` polylines; in-plane snap uses `wall_mid_sdf` + `GeometrySolver`. Characterization stays in [`tests/test_shell_graph_gaps.py`](../../../tests/test_shell_graph_gaps.py).

---

## G6 — Web feet on the stacked mid-surface

**Status:** Wave 4 shipped.

**Today:** [`ShellGraph.with_web_feet_on_loop`](../../../src/implicit_geometry_kernel/blade/cut.py) snaps web endpoints onto the stacked loop (applied by `shell_graph()`). 2D-F feet already shared vertices (gap **0**); unmerged skin walls remain ~`t_cap/2` away. Family D solids/recipes are unchanged.

**IEA 10 tip:** `SectionGeometryRunner.execute` (when layout is **not** a frozen sizing resolution) retries live `web_y` then drops webs that cannot re-anchor under twist. Passing `resolution=` still raises (SG-002).

---

## G7 — Live partition field without a full cut

**Status:** Wave 5 shipped.

**Today:** [`BladeGeometry.section_at`](../../../src/implicit_geometry_kernel/blade/loft/loft.py) caches cuts per `round(z, 9)` (LRU 256). [`CutWallMetaSDF`](../../../src/implicit_geometry_kernel/blade/topology/fields.py) and `StackThickness` keep derived LRU 256 on that same key. [`LoftedPolygonSDF`](../../../src/implicit_geometry_kernel/sdf/nodes/lofted.py) caches station polygons and, with ``accelerate=True`` (default), a uniform-grid AABB (`accel_at` / `bvh`). That is not the cut LRU.

---

## G8 — WindIO partition footprints

**Status:** Wave 6 shipped.

**Today:** [`web_positions_from_windio_structure`](../../../src/implicit_geometry_kernel/blade/windio.py) maps `structure.anchors` `nd_arc` (0/1 = TE, ~0.5 = LE) to span-median chord fractions. Example 09 `load_windio` passes that tuple into `BladeGeometry.web_positions` for any web count. [`cap_chord_windows_from_windio_structure`](../../../src/implicit_geometry_kernel/blade/windio.py) reads spar/cap layer windows; loft still places Family D caps from the webs.

Gauge tables still come from layer names.

---

## G9 — Cores / TE inserts / sandwich

**Status:** opt-in graph runs.

**Today:** sandwich cores and inserts remain **constructive solids** (`CoreWall.to_midline` / `InsertWall.to_midline` raise). `_mid_surface` still merges caps into skins + webs only. Default 2D-F `shell_graph` run names contain neither `core` nor `insert`; `surface` is `skin` / `web` (plus `rib` / `stringer` on faces and branches). `cut.walls` kinds stay `skin` / `cap` / `web`.

**Opt-in:** `shell_graph(cut, include_solids=True)` (also `BladeGeometry.shell_graph(z, include_solids=True)` and `SectionCut.shell_graph(include_solids=True)`) builds `ShellRun`s with `surface` `"core"` / `"insert"` from `to_boundary_polygon` rings, stored on `ShellGraph.solids` and included in `runs()`. Off by default so goldens and 2D-F run-name tests stay put.

**Decision:** keep them off the loop unless a loft consumer asks. Area/inertia keep the polygon path.

---

## G10 — Thin-wall vs shell consumers

**Status:** closed for this kernel.

**Today:** no `implicit_geometry_kernel.thin_wall` module (`test_g10_thin_walled_not_in_this_repo`). Handoff is `SectionCut.mid_surface()` / `shell_graph()` stacked `layers`. `cut.cells` are built from split **walls**, not from merged loop partitions.

**Not in this repo:** an external thin-wall consumer on `layers`. Verify that consumer separately.

---

## G11 — Combined vs tied laminates

**Status:** closed for this kernel.

**Today:** combined `T` is geometry. Laminate bind is `offset_regions()`, stacked `layers`, `shell_strips_from_graph`. No `PCOMP` string under `src/implicit_geometry_kernel` (`test_g11_no_pcomp_in_kernel_src`). ABD / ply stagger stay tidalshell / shell-fem.

---

## Former WP6 — feature-curve keys, TE on merged loop, BVH, IEA 10 tip

Parked, not a new wave. Named `feature_curve(geo, name_a, name_b, zs)` already
delegates to `partition_edge_curve`. The TE classifier remaining as an OML helper
is parked.

| Item | Evidence |
|---|---|
| Feature-curve keys | named form exists; Newton form stays two SDFs + a guess |
| TE on merged loop | `classify_trailing_edge(oml)` and on the stacked loop are both `sharp` on this NACA; the helper is still an **OML-ring** API |
| BVH on `LoftedPolygonSDF` | shipped: `accelerate=True`, `accel_at(z)`, `bvh` mapping (BG-015) |
| IEA 10 tip | G6: loft re-resolves / drops webs (Wave 4) |

---

## Waves

Opening a wave means implementing that row. **Waves 0–7 are shipped.** G10/G11 stay
closed for this kernel (consumer / PCOMP). No `ShellMidlineStrip` version, no `T`
taper. `LoftedPolygonSDF` BVH is shipped.

```text
Wave 0 (G5+G4) → Wave 1 (G2) → Wave 2 (G3) → Wave 3 (G1)
  → Wave 4 (G6) → Wave 5 (G7) → Wave 6 (G8) → Wave 7 (G9) → G10/G11 closed here
```

### Wave 0 — Oracle honesty (shipped)

**G5**

- **D0.1** Two oracles: loft / conformal strips use `shell_graph` polylines; in-plane snap uses `wall_mid_sdf` + `GeometrySolver`.
- **D0.2** Joint vertices are out of the Pass 1 `|φ| < 5e-4` check. Cap/skin ends: max `|φ|` = `t_cap/2` on this fixture (**4.000 mm**); cap midpoint ~0; snap residual ~4×10⁻⁷ m (98.6 % converged).
- **D0.3** Tests stay characterization (`test_shell_graph_gaps.py`). Do not fail midpoints at 4 mm.

**G4** (locked: document Newton limits, keep the step)

- **D0.4** `StackThickness` remains a step (11 mm vs 3 mm; jump 8 mm). Shared vertices may evaluate to either side (`{3, 11}` mm).
- **D0.5** Solver defaults (`tolerance=1e-6`, `snap_wall_mid` `1e-8`) are not a joint spec. Do not Newton-project along the thickness jump expecting the discrete graph vertex.
- **Non-goal:** geometric taper.

### Wave 1 — Partition edge keys (G2) (shipped)

Colouring is done. Missing mesh-facing **edges**.

- **D1.1** Table `(z, name_a, name_b) → xy` from adjacent `ShellRun` ends (cap edge stays a feature line).
- **D1.2** Optional: bind `feature_curve` to those names between stations (no spanwise topology yet).
- **D1.3** `station_network_from_meta` is not the source of thickness (`thickness=0.0` today).
- **Non-goal:** change example 09 colours.

### Wave 2 — Spanwise topology (G3) (shipped)

Depends on named runs (Wave 1).

- **D2.1** Event type: cap runout, dropped web, bay split. No silent label lerp.
- **D2.2** Loft rule: missing name → omit strip (already example 09 `keys() ∩ keys()`); promote that to a kernel helper.
- **D2.3** Tests with a dropped web / missing run. The 2D-F z=4 fixture has empty `dropped_webs`; use a dedicated station.

**Locked resampling (post-Wave 2, not G10):** `LockedLoftPlan` holds the same chord count on every run so a surviving name can loft through a topology event. `resample_locked` writes partition-edge vertices once so adjacent skin/cap strips share a seam. `shell_loft_zones` still records Wave 2 name-signature bands. Figure loft and `write_mid_surfaces` now join with `WallRegionGraph` / `map_vertex` / `loft_mapped_strips` (TE wrap and 1-station paint flicker are transition bands, not holes). A conformal mesher stays out of this repo.

### Wave 3 — One shell hand-off (G1) (shipped)

Depends on Waves 0–2.

- **D3.1** Shell consumers take `MidSurface` / `ShellGraph` (or a new DTO). Stop emitting extra cap strips for shells (10 `cut.walls` vs 8 `shell_graph` runs).
- **D3.2** `cut.walls` **keeps** caps (solids / bonds / cells). Do not collapse the two views.
- **D3.3** `ShellMidlineStrip` stays the frozen strip DTO until an explicit versioning WP. Adapter may live here or remain outside (`shell_mesh_inputs_from_geometry_result` is not in this repo).

### Wave 4 — Web feet and IEA 10 (G6) (shipped)

- **D4.1** 2D-F: feet already on the stacked loop (gap 0) via shared vertices; 4.00 mm from unmerged skin walls (`t_cap/2`). Decide if that is enough or if Family D needs an explicit re-host.
- **D4.2** IEA 10 tip (`z = 96.755` m): `flapwise web cannot re-anchor` — layout re-resolve or `web_y`, not a shell-graph merge. Separate from D4.1.

### Wave 5 — Cost (G7) (shipped)

- **G7:** `section_at` LRU 256 per `round(z, 9)`. Derived `CutWallMetaSDF` / `StackThickness` caches stay. `LoftedPolygonSDF` uses `accelerate` BVH (`accel_at` / `bvh`).

### Wave 6 — WindIO web footprints (G8) (shipped)

- Reader: `nd_arc` → x/c; median over the anchor grid (span-constant `web_positions`).
- Example 09 uses the reader; three-web blades get YAML fractions, not catalogue defaults.
- Cap windows are readable; caps on the loft still follow the layout between webs.

### Wave 7 — Cores / TE inserts (G9) (shipped)

- **Decision:** solids-only. Do not extend `_mid_surface` with core or insert runs.
- `CoreWall.to_midline` stays unimplemented. Properties keep the constructive polygon path.
- Re-open only if a shell consumer needs those names on the loop.

### Closed here / remaining park

- **G10** closed here (other repo owns the thin-wall adapter).
- **G11** closed here (PCOMP / ply ABD stay tidalshell / shell-fem).
- TE classifier on OML vs loop: both `sharp` on the NACA fixture — retarget only if a blunt merged loop appears.

---

## What this inventory does not do

G10/G11 stay out of `src/` (thin-wall consumer; PCOMP). No DTO version, no `T`
taper, no cores on `shell_graph` by default, no Gmsh driver, no tidalshell edits.
OML polygon BVH is in `LoftedPolygonSDF`.
