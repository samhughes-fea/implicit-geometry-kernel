# Merged shell graph — Pass 2 gap inventory

Evidence for G1–G11 and the parked extract/TE/BVH work. Waves 0–7 are shipped; G10/G11 stay parked.
Characterization: [`tests/test_shell_graph_gaps.py`](../../../tests/test_shell_graph_gaps.py).

Fixture unless noted: 2D-F, `z = 4` m, NACA 0024→0018, chord 1.6 m, skin 3 mm, cap 8 mm, webs 2 mm at chord fractions `(0.18, 0.48)`.

↑ [`geometry_to_mesh`](../guides/geometry_to_mesh.md)

## Summary

| ID | Status | One-line hole |
|---|---|---|
| **G1** | Wave 3 shipped | `shell_strips_from_graph` = 8 runs; `cut.walls` still 10 (solid caps) |
| **G2** | Wave 1 shipped | `partition_edges` / `partition_edge_curve`; colouring unchanged |
| **G3** | Wave 2 shipped | `topology_diff` + `shared_loft_names`; no label lerp |
| **G4** | Wave 0 shipped | `T` step kept (8 mm); Newton is not a spec at partition vertices |
| **G5** | Wave 0 shipped | Loft from graph; snap with SDF; joint `|φ| = 4.000 mm` is out of Pass 1 check |
| **G6** | Wave 4 shipped | `with_web_feet_on_loop`; loft re-resolves / drops webs on flapwise re-anchor |
| **G7** | Wave 5 shipped | `section_at` LRU 256 per `round(z, 9)`; no BVH |
| **G8** | Wave 6 shipped | `web_positions_from_windio_structure`; example 09 no longer hardcodes `(0.22, 0.50)` |
| **G9** | Wave 7 shipped | Cores / TE inserts stay constructive solids; not on `shell_graph` |
| **G10** | consumer-only | `ThinWalledSection.from_geometry` is not in this repo |
| **G11** | kernel-out-of-scope | No PCOMP / ply stagger in `src/implicit_geometry_kernel` |
| **WP6** | open | `classify_trailing_edge` is an OML helper; `feature_curve` unkeyed; no BVH |

Implementation sequence is [Waves](#waves) below. **G4 is locked:** keep the `T` step; document Newton limits at partition vertices (no geometric taper).

Oracles (Wave 0): **loft strips** from `shell_graph` polylines; **in-plane snap** from `wall_mid_sdf` + `GeometrySolver`. Joint vertices are out of the Pass 1 `|φ| < 5e-4` check.

---

## G1 — One shell hand-off DTO

**Status:** Wave 3 shipped.

**Today:** [`build_shell_midline_strips`](../../../src/implicit_geometry_kernel/section/midlines/shell_midline_export.py) still appends `kind="cap"` strips for **solids**. [`shell_strips_from_graph`](../../../src/implicit_geometry_kernel/blade/cut.py) is the shell producer: one `ShellMidlineStrip` per graph run (combined `T`, cap as partition `kind`). [`ShellMidlineStrip`](../../../src/implicit_geometry_kernel/contract/shell_midline_strip.py) is **not** versioned. `cut.walls` still lists solid caps.

**Evidence (2D-F, z=4):** 10 `SectionWall` rows vs 8 shell strips / `ShellRun`s.

**Out of repo:** `shell_mesh_inputs_from_geometry_result` still consumes the solids emit unless a downstream adapter is pointed at `shell_strips_from_graph`.

**Non-goals:** do not version `ShellMidlineStrip`; do not delete caps from `cut.walls`.

---

## G2 — Feature lines at partition edges

**Status:** Wave 1 shipped (colouring in Pass 1; named edges on `ShellGraph`).

**Today:** `ShellRun.kind` is `cap` on stacked partitions; example 09 and `wall_meta` colour by that. Partition **edges** are only polyline endpoints. [`feature_curve`](../../../src/implicit_geometry_kernel/extract/curves.py) is `(field_a, field_b, zs, guess_xy)` with no names. [`station_network_from_meta`](../../../src/implicit_geometry_kernel/blade/network.py) splits the mid loop by meta id and sets `thickness=0.0`.

**API:** `ShellGraph.partition_edges()` / `.edge(name_a, name_b)` and
`partition_edge_curve(geo, name_a, name_b, zs)` (omit if missing; no lerp).
`feature_curve` stays two SDFs + a guess. `station_network_from_meta` thickness
stays `0.0`.

**Non-goals:** do not change example 09 cap/bay colours.

---

## G3 — Spanwise topology events

**Status:** Wave 2 shipped.

**Today:** example 09 [`_wall_faces`](../../../examples/09_reference_blades/reference_blades.py) lofts `a.keys() & b.keys()`. Missing names are omitted. No kernel `TopologyEvent`. Fixture `dropped_webs` at z=4 is empty. Cap runout / bay split would lerp silently if someone interpolated by index instead of name; the example does not lerp.

**API:** `topology_diff(a, b)` emits `cap_runout` / `dropped_web` / `bay_split`
(and onset / added_web / bay_merge). `shared_loft_names` is the loft rule.
Example 09 uses it. No silent label lerp.

**Non-goals:** no silent label lerp.

---

## G4 — Thickness jumps

**Status:** Wave 0 shipped (keep the step; Newton limits documented).

**Today:** [`StackThickness`](../../../src/implicit_geometry_kernel/blade/fields.py) assigns nearest loop-run `T` (step, not ramp). Not ply stagger.

**Evidence:** interior cap `T = 11` mm, interior skin `T = 3` mm, jump **8 mm** (`t_cap`). The two endpoints of one cap run evaluate to `{3 mm, 11 mm}` — nearest-run at the shared vertex is discrete. `GeometrySolver` default tolerance is `1e-6` m (`snap_wall_mid` uses `1e-8`); the jump is 8000× that.

**Decision (Wave 0, locked):** keep the nearest-run `T` step. Document Newton limits at partition vertices. No geometric taper.

**Non-goals:** ply-drop stagger; taper of `T`.

---

## G5 — SDF vs `mid_surface()` residual

**Status:** Wave 0 shipped (Pass 1 only checked midpoints).

**Today:** two oracles — discrete `shell_graph` polylines vs `wall_mid_sdf` (OML polygon + stack `T`). Example 05 discrete stack pairing is **20.5 µm** (cap vs skin vertices along the bond, not this SDF).

**Evidence (2D-F, z=4):**

| Sample | max \|φ\| of `wall_mid_sdf` |
|---|---|
| Cap-run **midpoint** | ~0 (1e-17) |
| All **cap** vertices | **4.000 mm** (one partition end) |
| All **skin** vertices | **4.000 mm** |
| After `snap_wall_mid` on cap vertices | 4.35×10⁻⁷ m; 98.6 % converged |

The 4 mm spike is `t_cap/2` at the cap/skin joint — the polyline sits on the stacked mid, the SDF offset is still the skin mid (or the reverse) at that vertex. Pass 1 `|φ| < 5e-4` holds only away from the joint.

**Policy (Wave 0):** loft / conformal strips use `shell_graph` polylines; in-plane snap uses `wall_mid_sdf` + `GeometrySolver`. Joint vertices are **out of** the Pass 1 `|φ| < 5e-4` check. Do not treat 5e-4 as a joint-vertex spec. Characterization stays in [`tests/test_shell_graph_gaps.py`](../../../tests/test_shell_graph_gaps.py) — do not fail midpoints at 4 mm.

---

## G6 — Web feet on the stacked mid-surface

**Status:** Wave 4 shipped.

**Today:** [`ShellGraph.with_web_feet_on_loop`](../../../src/implicit_geometry_kernel/blade/cut.py) snaps web endpoints onto the stacked loop (applied by `shell_graph()`). 2D-F feet already shared vertices (gap **0**); unmerged skin walls remain ~`t_cap/2` away. Family D solids/recipes are unchanged.

**IEA 10 tip:** `SectionGeometryRunner.execute` (when layout is **not** a frozen sizing resolution) retries live `web_y` then drops webs that cannot re-anchor under twist. Passing `resolution=` still raises (SG-002).

---

## G7 — Live partition field without a full cut

**Status:** Wave 5 shipped.

**Today:** [`BladeGeometry.section_at`](../../../src/implicit_geometry_kernel/blade/loft.py) caches cuts per `round(z, 9)` (LRU 256). [`CutWallMetaSDF`](../../../src/implicit_geometry_kernel/blade/fields.py) and `StackThickness` keep derived LRU 256 on that same key. [`LoftedPolygonSDF`](../../../src/implicit_geometry_kernel/sdf/lofted.py) still caches polygons only — **no BVH**.

---

## G8 — WindIO partition footprints

**Status:** Wave 6 shipped.

**Today:** [`web_positions_from_windio_structure`](../../../src/implicit_geometry_kernel/blade/windio.py) maps `structure.anchors` `nd_arc` (0/1 = TE, ~0.5 = LE) to span-median chord fractions. Example 09 `load_windio` passes that tuple into `BladeGeometry.web_positions` for any web count. [`cap_chord_windows_from_windio_structure`](../../../src/implicit_geometry_kernel/blade/windio.py) reads spar/cap layer windows; loft still places Family D caps from the webs.

Gauge tables still come from layer names.

---

## G9 — Cores / TE inserts / sandwich

**Status:** Wave 7 shipped (solids-only).

**Today:** sandwich cores and `trailing_edge_insert` are **constructive solids** (`CoreWall.to_midline` is not implemented). `_mid_surface` merges caps into skins + webs only. 2D-F `shell_graph` run names contain neither `core` nor `insert`; `surface` is only `skin` / `web`. `cut.walls` kinds are `skin` / `cap` / `web`.

**Decision:** keep them off the loop. Area/inertia use the polygon path; shell loft uses the stacked skin+web graph. Put them on the graph only if a later consumer WP needs those runs.

---

## G10 — Thin-wall vs shell consumers

**Status:** consumer-only.

**Today:** `ThinWalledSection.from_geometry` is named in the README and `SectionCut.mid_surface` docstring; **no** `implicit_geometry_kernel.thin_wall` module. `cut.cells` are built from split **walls**, not from merged loop partitions.

**Missing work:** verify the external thin-wall consumer on `layers`; split vs merged cells. No patch in this repo.

---

## G11 — Combined vs tied laminates

**Status:** kernel-out-of-scope.

**Today:** combined `T` is geometry. No `PCOMP` string under `src/implicit_geometry_kernel`. Stays tidalshell / shell-fem.

---

## Former WP6 — feature-curve keys, TE on merged loop, BVH, IEA 10 tip

| Item | Evidence |
|---|---|
| Feature-curve keys | `feature_curve` has no partition names (G2) |
| TE on merged loop | `classify_trailing_edge(oml)` and on the stacked loop are both `sharp` on this NACA; the helper is still an **OML-ring** API |
| BVH on `LoftedPolygonSDF` | no `bvh` attribute (G7) |
| IEA 10 tip | G6: loft re-resolves / drops webs (Wave 4) |

---

## Waves

Opening a wave means implementing that row. **Waves 0–7 are shipped.** G10/G11 stay
parked. No `ShellMidlineStrip` version, no `T` taper, no BVH.

```text
Wave 0 (G5+G4) → Wave 1 (G2) → Wave 2 (G3) → Wave 3 (G1)
  → Wave 4 (G6) → Wave 5 (G7) → Wave 6 (G8) → Wave 7 (G9) → park G10/G11
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

### Wave 3 — One shell hand-off (G1) (shipped)

Depends on Waves 0–2.

- **D3.1** Shell consumers take `MidSurface` / `ShellGraph` (or a new DTO). Stop emitting extra cap strips for shells (10 `cut.walls` vs 8 `shell_graph` runs).
- **D3.2** `cut.walls` **keeps** caps (solids / bonds / cells). Do not collapse the two views.
- **D3.3** `ShellMidlineStrip` stays the frozen strip DTO until an explicit versioning WP. Adapter may live here or remain outside (`shell_mesh_inputs_from_geometry_result` is not in this repo).

### Wave 4 — Web feet and IEA 10 (G6) (shipped)

- **D4.1** 2D-F: feet already on the stacked loop (gap 0) via shared vertices; 4.00 mm from unmerged skin walls (`t_cap/2`). Decide if that is enough or if Family D needs an explicit re-host.
- **D4.2** IEA 10 tip (`z = 96.755` m): `flapwise web cannot re-anchor` — layout re-resolve or `web_y`, not a shell-graph merge. Separate from D4.1.

### Wave 5 — Cost (G7) (shipped)

- **G7:** `section_at` LRU 256 per `round(z, 9)`. Derived `CutWallMetaSDF` / `StackThickness` caches stay. **No** BVH on `LoftedPolygonSDF`.

### Wave 6 — WindIO web footprints (G8) (shipped)

- Reader: `nd_arc` → x/c; median over the anchor grid (span-constant `web_positions`).
- Example 09 uses the reader; three-web blades get YAML fractions, not catalogue defaults.
- Cap windows are readable; caps on the loft still follow the layout between webs.

### Wave 7 — Cores / TE inserts (G9) (shipped)

- **Decision:** solids-only. Do not extend `_mid_surface` with core or insert runs.
- `CoreWall.to_midline` stays unimplemented. Properties keep the constructive polygon path.
- Re-open only if a shell consumer needs those names on the loop.

### Parked

- **G10** thin-wall consumer (other repo).
- **G11** PCOMP (tidalshell / shell-fem).
- TE classifier on OML vs loop: both `sharp` on the NACA fixture — retarget only if a blunt merged loop appears.

---

## What this inventory does not do

This page does not close G10–G11. No DTO version, no taper, no cores on `shell_graph`, no BVH, no Gmsh, no PCOMP, no tidalshell edits — until the matching wave is opened.
