# bug_log — blade loft

**Audience: internal.** Maintainer notes; not part of the library manual.

↑ [`blade`](README.md) · [docs hub](../README.md)

Worker debug log for span loft (`BladeGeometry`, `shell_graph`, mapped strips).
Each entry records the **symptom**, **evidence**, **root cause**, and
**verification**. Entries stay after they close, per
[`DOCUMENTATION_STANDARDS`](../conventions/DOCUMENTATION_STANDARDS.md).

**Identifiers.** `BG-nnn` here. Independent of `SG-nnn`
([`section/bug_log.md`](../section/bug_log.md)) and `CF-nnn`
([`frames/bug_log.md`](../frames/bug_log.md)). Historical `BL-nnn` ids in closed
entries are from a retired repo-root log.

**Reference cases.** Example 09 catalogue keys on the synthetic 12 m planform
(`4D-F-PA-10R` for PA freeze). windIO: IEA 15 / 3.4 / 10 / 22, NREL 5 MW from
`reference_blades/` (example 08). Tests: [`tests/test_shell_topology.py`](../../tests/test_shell_topology.py),
[`tests/test_reference_transition.py`](../../tests/test_reference_transition.py),
[`tests/test_debug_mesh.py`](../../tests/test_debug_mesh.py).

## Index

| ID | Title | Status |
|---|---|---|
| [BG-001](#bg-001) | PA loft used catalogue LC fractions instead of the live intercept | **Fixed** |
| [BG-002](#bg-002) | Cylinder cluster radius 0.20 merged PA `web_2` | **Fixed** |
| [BG-003](#bg-003) | Pressure-skin holes at web onset/drop (loop name reuse) | **Fixed** |
| [BG-004](#bg-004) | VTK loft patches were grey in ParaView | **Fixed** |
| [BG-005](#bg-005) | IEA 15 TE wrap on/off wrinkled the loft | **Fixed** |
| [BG-006](#bg-006) | IEA 3.4 `fore_web` dropped at one station vs `cap_lo` | **Fixed** |
| [BG-007](#bg-007) | IEA 10 loft webs ≠ intercept; TE-oml ~0.41 | **Fixed** |
| [BG-008](#bg-008) | 5D CSG stations 05/06 lost suction TE (held wrap on a foil) | **Fixed** |
| [BG-009](#bg-009) | Family C conformal mesh had no per-web skin–sparcap patch | **Fixed** |
| [BG-010](#bg-010) | Cap loop branches partitioned rib mesh patches | **Fixed** |
| [BG-011](#bg-011) | Unconstrained rib fill sawtoothed skin/web edges | **Fixed** |
| [BG-012](#bg-012) | Example 09 at 12 m truncated the 20 m planform | **Fixed** |
| [BG-013](#bg-013) | Station CSG figures showed φ / OML disk, not wall CSG | **Fixed** |
| [BG-014](#bg-014) | Planform PNGs duplicated the legend in a grey caption | **Fixed** |
| [BG-015](#bg-015) | `LoftedPolygonSDF` BVH path disagreed with edge-loop φ | **Fixed** |
| [BG-016](#bg-016) | Cap loft identity flips `*_0` ↔ unindexed across web onset | **Fixed** |
| [BG-017](#bg-017) | Sparse PCHIP invented curvature; docs claimed it interpolated pitch | **Fixed** |
| [BG-018](#bg-018) | Bowtie loft quads on one span band | **Fixed** |
| [BG-019](#bg-019) | Live web anchors were not loop partition vertices | **Fixed** |
| [BG-020](#bg-020) | Jump-join of non-adjacent cap fragments opens the skin-mid loop | **Fixed** |
| [BG-021](#bg-021) | PA extra web sat two Δx past the last primary | **Fixed** |
| [BG-022](#bg-022) | `SectionCut.layout` became a sibling token after web drop | **Fixed** |
| [BG-023](#bg-023) | Example 09 mapped loft emitted collapsed or unmatched strip faces | **Fixed** |
| [BG-024](#bg-024) | Property CSG quadtree was rotated vs wall occupancy | **Fixed** |
| [BG-025](#bg-025) | Property CSG column 3 undersamples thin wall tubes | **Fixed** |
| [BG-026](#bg-026) | 5D occupancy QT hits `max_cells=150000` and truncates | **Fixed** |
| [BG-027](#bg-027) | Planform `axes.png` labelled cell-area locus as shear centre | **Fixed** |
| [BG-028](#bg-028) | Loft strips resampled the lower station: cracks on station lines | **Fixed** |
| [BG-029](#bg-029) | Rib/tip polygon fill left hanging nodes and 3–4-node cells | **Fixed** |
| [BG-030](#bg-030) | Web/cap start and end bands meshed with triangles and cracks | **Fixed** |
| [BG-031](#bg-031) | Mesh cap regions used one chord window for both sides | **Fixed** |
| [BG-032](#bg-032) | Cap edge passing a web foot between stations folded the strips | **Fixed** |
| [BG-033](#bg-033) | Side arcs walked through the other TE corner (IEA 15) | **Fixed** |
| [BG-034](#bg-034) | Sliver transition quads where node counts change | Resolved (K8) |
| [BG-035](#bg-035) | IEA 15 TE wrap folds at 15–22 % span | **Open** (SG-030) |
| [BG-036](#bg-036) | windIO webs kept root to tip without `start_nd_span` | **Fixed** (K9) |
| [BG-037](#bg-037) | Webs cross each other between stations (IEA 15, IEA 10) | **Fixed** (IEA 15 web anchors follow the placed arcs through the window end) |
| [BG-038](#bg-038) | windIO `nd_arc` read as chord position (IEA 15 webs 25–290 mm off) | **Fixed** (K9) |
| [BG-039](#bg-039) | Tip closure left out on a single-cell tip at a coarse `h` | **Fixed** (K9) |
| [BG-040](#bg-040) | `merge_cap_and_skin=True` drops the caps from the shell graph | **Open** |
| [BG-041](#bg-041) | windIO caps found by layer name only (IEA 10 gets none) | **Open** |
| [BG-042](#bg-042) | Station airfoil missed its own points at a blunt-TE corner (0.1 mm at 4.5 m chord) | **Fixed** (K9) |
| [BG-043](#bg-043) | Round trailing end bulged aft of its end point and was scrambled (0.1 mm at 5.4 m chord) | **Fixed** (K9) |
| [BG-044](#bg-044) | windIO web feet and cap edges interpolated in chord fraction between grid points (up to 10 mm) | **Fixed** (K9) |
| [BG-045](#bg-045) | A partial-span web started and ended at the nearest mesh station, not at its own ends | **Fixed** (K9) |

Keep-set / TE / VTK pass (2026-11-08): **BG-001**–**BG-007** (PA identity,
cylinder cluster, pressure-skin holes, ParaView RGB, IEA 15 wrap, IEA 3.4
web-vs-cap, IEA 10 keep-set and TE-oml). Freeze rules in the next section.

## Mapped loft freeze (2026-11-08)

Do not add per-blade branches (`if iea10`) in loft code.

| Rule | When |
|---|---|
| Loft webs = intercept `DroppedWeb.index` | `dropped_web_indices` is passed (`_template_from_cut`, examples 08/09) |
| No 0.08 cluster, no extra TE-plug omit | same keep-set path |
| `_cluster_web_members` (0.08, cylinder TE-plug) | only when `dropped_web_indices is None` |
| `_regularise_anchors`: `web_*` wins over nearby `cap_*` | always; do not merge two distinct webs |
| Placement graph law = section law | PA/LE freeze in metres (PA: uniform root-band Δx lattice, extras continue that lattice); LC is local chord fraction |
| Skin `t0` mould walk | PA onset/drop; not a cavity chord or collapsed TE |
| TE wrap includes the OML ridge | `_te_wrap_arc`; metrics use ridge / max chord-frac, not max blade-x |
| Wrap presence from skin-mid opening | Blunt vs sharp held on a loft `z` list; a single `section_at` / CSG frame uses the raw test (BG-005, BG-008) |
| Active web anchors are loop joints | Template rebuilds regions from `evaluate_anchors`; web mid ends share skin nodes (BG-019) |
| Bowtie loft quads | `_quad_bowtie` splits an opposing-normal quad into two triangles (BG-018); K7 removes their causes (BG-031–BG-033); IEA 15 TE wrap left (BG-035) |
| Rib rings = `n_live_webs + 1` | `rib_cell_polygons_from_samples` splits only at live web feet; `cap_lo` / `cap_hi` stay on the skin (BG-010) |
| Cut `layout` is the authored concept | Web drop lists `dropped_web_indices`; realised W/X are `station_webs` / `station_family` (BG-022) |

## BG-001

**Symptom.** On `4D-F-PA-10R`, `web_0` / `web_3` vanished mid-chord (ribs 5–6)
instead of walking to the LE/TE with the live intercept.

**Evidence.** Authored `WebMember.at(z)` followed catalogue LC fractions, not
the PA blade-x freeze. `test_09_pa_4d_loft_tracks_pitch_offset_not_catalogue`.

**Root cause.** Loft anchors used the layout token as a chord-fraction table
rather than the metre freeze (PA blade-x / LE offset).

**Fix.** PA is hard-coded blade-x about the pitch axis (uniform interiors of
the root band; extra webs continue that lattice). Labels keep authored
indices.

**Verification.** `test_09_pa_4d_loft_tracks_pitch_offset_not_catalogue`,
`test_09_pa_4d_loft_keeps_live_webs_on_cylinder`. **Fixed.**

## BG-002

**Symptom.** `web_2` missing through the cylinder root on the four-web PA lattice
(0.2 / 0.4 / 0.6 / 0.8).

**Evidence.** `_cluster_web_members` `min_sep` 0.20 merged neighbours;
`0.6 - 0.4` is slightly under 0.20 in float.

**Root cause.** Cylinder cluster radius was 0.20 chord, not `min_gap`.

**Fix.** Cluster radius stays `min_sep` (0.08 on the legacy path). TE-named plugs
still omit on a cylinder when `dropped_web_indices is None`.

**Verification.** `test_cluster_keeps_pa_lattice_on_cylinder`. **Fixed.**

## BG-003

**Symptom.** Two holes in the pressure skin at web-3 onset and just after drop
(VTK `bay_pr_*` missing).

**Evidence.** Loop names (`te_pr`, `bay_*`) reused across add/drop;
`loft_mapped_strips` parented the new panel to a collapsed TE or cavity chord.

**Root cause.** Skin loft at a topology event did not rule each t1 panel from
the t0 mould walk between mapped ends.

**Fix.** `_skin_parent` mould walk; dropped web anchors leave the loop.

**Verification.** `test_09_pa_4d_pressure_skin_covers_web_events`,
`test_09_pa_4d_web3_onset_and_drop_are_skin_ruled`. **Fixed.**

## BG-004

**Symptom.** Example PNG panel colours did not appear in ParaView (grey loft).

**Evidence.** VTK had `role` strings only; ParaView maps RGB from
`COLOR_SCALARS` (and “Map Scalars” off for PNG-true colour).

**Root cause.** `write_loft_vtk` omitted cell RGB.

**Fix.** Write `COLOR_SCALARS rgb 3` from the same `panel_colour` table; keep
`role` as the region name.

**Verification.** `test_write_loft_vtk_polydata_roles`. **Fixed.**

## BG-005

**Symptom.** IEA 15 TE gap at root/transition, then spanwise wrinkles
(wrap on/off, y-span ~0.04c). Example 08 dumps also lofted a zero-area `te`
wrap strip on every sharp station (~7520 degenerate faces = 160 bands × 47
chord steps; 09: 920 = 40 × 23).

**Evidence.** Sampling omitted wrap when suction/pressure keys snapped to the
ridge; a t/c-scaled inboard drop invented wrap on every sharp OML. Dump review
of skin-mid + web-mid OBJ/NPZ: almost all “deg” faces were role `te` with
coincident SS/PS TE keys.

**Root cause.** Wrap presence was a sampling threshold, not a loft-zone event
from the skin-mid opening. On a sharp TE the wrap chord SS→PS collapsed to a
line and was still lofted.

**Fix.** Classifier: blunt if TE opening `> 2 t_skin + 0.005 c` (skin-mid
offset of a point TE is ~2t; 0.005c is a finite OML edge). Else sharp. Hold
the boolean along a loft `z` list: leave blunt when opening `≤ ε`, enter when
`> ε`. A single `section_at` (including CSG frames) uses the raw test. Blunt
stations keep wrap on the real TE edge; sharp stations drop the wrap branch
and share one ridge key. Sampling also skips wrap when
`‖TE_ss − TE_ps‖ < 1e-6 c`. Wrap-off is only the collapsing patch at the zone
face (`loft_mapped_strips`). `_te_frac_drop` is superseded. `TE_ARC` remains
the web-exclusion band.

**Verification.** `test_ellipse_root_lofts_te_wrap`, `test_sharp_te_has_no_wrap_branch`,
`test_flatback_te_wrap_samples_skin_arc`, IEA 15
`test_windio_span_metrics` (`n_coverage_holes == 0`); regenerated 08/09
`samples.npz` show wrap-role `te` degenerates gone on sharp stations. **Fixed.**

Hold applied to independent CSG frames is [BG-008](#bg-008).

## BG-006

**Symptom.** IEA 3.4 `fore_web` missing from loft samples at z = 59.85 m
(neighbours present).

**Evidence.** `web_0` at 0.2488 vs `cap_lo` at 0.2469; `min_gap=0.002` dropped
the web so `place_station_keys` had no feet.

**Root cause.** Close-pair regularise kept `cap_*` (sorted first) and discarded
`web_*`.

**Fix.** Never discard a `web_*` for a cap; skip a cap that lands on an existing
web.

**Verification.** `test_regularise_keeps_web_over_nearby_cap`,
`test_iea34_fore_web_at_5985`. **Fixed.**

## BG-007

**Symptom.** IEA 10 midspan loft omitted `web_te` / `rear_web` while intercept
`n_int==3`. TE-oml ~0.41 (IEA 22 ~0.10).

**Evidence.** 0.08 cluster + TE-plug omit after `dropped_web_indices` filter.
TE sample used `argmax` blade-x; wrap could miss the OML ridge. Tip OML
`argmax x` ≠ chord-fraction TE.

**Root cause.** Extra keep-set gates; TE metric not locked to the ridge.

**Fix.** When `dropped_web_indices` is set, filter by index only. Wrap must
include `i_te`. `_te_xyz` nearest ridge / max chord-frac;
`span_metrics` uses `_oml_te_xyz`.

**Verification.** `test_iea10_loft_keepset_matches_intercept`,
`test_windio_span_metrics` IEA 10 `max_te_oml` / `max_te_chatter`,
`test_te_wrap_arc_includes_oml_ridge`. **Fixed.**

## BG-008

**Symptom.** `5D-F-PA-10R` CSG `05.png` / `06.png` (z = 3.2 m, 4.0 m): occupancy
and CSG solid were a foil; conformal nodes dropped the suction skin and grew a
TE spur past the OML. `07.png` (z = 4.8 m) was already a closed ridge.

**Evidence.** `te_06.png` wrap stub near x = 0.90 m; suction loop names missing
from `sample_region_polylines`. `figure_csg_mesh` passed `te_blunt_flags_for_cuts`
(loft hysteresis) into each independent frame. Leave-at-`0.5 ε` never fired
while opening sat near `2t`.

**Root cause.** Hold used on a single `section_at`. Farthest-pair TE corners in
the `0.005c` band could lie on one skin, so `_side_arcs` collapsed suction.

**Fix.** CSG uses the raw `te_is_blunt` / `te_blunt_of_cut`. Same-side corner
pair collapses to the ridge. Loft hold leaves blunt when opening `≤ ε`.

**Verification.** `test_09_5d_foil_stations_have_no_te_spur`; regenerated
`5D-F-PA-10R` `csg_mesh/05.png`–`07.png` and `te_05`–`te_07`. **Fixed.**

## BG-009

**Symptom.** Family C (`2C-F-LC-10R`) CSG solids showed a discrete pad on each
web, but the conformal-node panel (skin-shell midsurface mesh) had no
skin–sparcap patch at those pads. The strip between the webs was labelled as
the cap; the real pad halves were bays.

**Evidence.** Example 09 `2C-F-LC-10R` `csg_mesh/10.png` (z = 7.2 m): occupancy
and CSG solid had two short red pads; conformal nodes (`n_chord=8`) did not.
`WallRegionGraph.from_section_layout("2C-F-LC-0R")` authored `cap_su` as
`web_0 → web_1` and `bay_su_0` / `bay_su_1` as the outer pad halves.

**Root cause.** `_layout_parts` collapsed `per_web_pads` to one outer
`cap_lo`/`cap_hi` box. `_regions_from_ordered` then marked every web-to-web
span as a cap (the family D predicate).
`ShellTopologyTemplate.from_region_graph` rebuilt regions from that same
helper, so `sample_region_polylines` never saw per-web pad edges.

**Fix.** One chord window per web (`web_i ±` pad half-width, overlapping pads
merged). A loop branch is a cap only when its interval lies inside a window.
Family B (fixed band) and D (box between outer webs) unchanged. Template uses
`live_cap_windows`.

**Verification.** `test_family_c_regions_are_per_web_pads_not_a_box`,
`test_family_c_one_web_marks_both_pad_halves`,
`test_overlapping_family_c_pads_merge_into_one_window`,
`test_family_c_template_has_pad_patches_and_inter_pad_bay`,
`test_family_d_cap_still_spans_outer_webs`,
`test_family_b_cap_is_the_fixed_band`. **Fixed.**

## BG-010

**Symptom.** Rib faces in example 09 were extra-cut where caps painted the
skin. `2C-F-LC-10R` (two webs) showed five independently meshed rib patches
per station instead of three; `1B-F-LC-10R` (one web) showed three instead of
two. `0B-LC-0R` has no ribs — the brown root is `cap_pr` loft — but the same
seeder would have emitted three rings with zero webs.

**Evidence.** Cutaways
`examples/blade/09_blade_concepts/figures/2C-F-LC-10R/full_blade/cutaway.png`
and `0B-LC-0R/full_blade/cutaway.png`. Template suction branches for 2C were
`le, cap_lo, web_0, web_1, cap_hi, te` (five `LoopBranch` pairs). CSG
`build_cell_wall_loops` already ignores caps (`n_webs + 1`).

**Root cause.** `rib_cell_polygons_from_samples` treated every suction/pressure
`LoopBranch` pair as a rib patch. Those branches exist so the skin midsurface
can name constitutive OFFSET regions (`LaminateOffsetRegion`; xy stays the
skin). Cap anchors are not cell walls.

**Fix.** Group consecutive skin branches between `le` / live `web_*` / `te`
and concatenate cap-split arcs into one ring per void cell. Patch count is
`n_live_webs + 1`. Voronoi n-gon fill inside each patch is unchanged.

**Verification.** `test_rib_patches_split_only_at_live_webs` (2A/2C/2D → 3,
1B → 2), `test_rib_patches_follow_dropped_pa_webs`,
`test_rib_patches_zero_webs_one_ring`. Regenerated 09
`1B-F-LC-10R` and `2C-F-LC-10R` galleries. **Fixed.**

## BG-011

**Symptom.** Rib interiors meshed with unconstrained Delaunay or a vertex-0
fan did not keep the parent skin/web edges: loft-node rings became sawtooth
triangles, the LE C-cell could clip empty, and cutaways looked like solid fans
with no interior edges matching loft quads.

**Evidence.** 09 cutaways vs loft quads on the same walls. `fill_rib_cells`
lived in `mesh/delaunay.py` / ear-clip; matplotlib/STL fanned each n-gon from
vertex 0. C-cell tests needed a synthetic ring because LE ear-clip produced no
triangles. Cutaway captions: both LE/TE views use `skip_suction=True`.

**Root cause.** Steiner triangulation of a closed loft ring does not treat the
ring edges as hard constraints. A vertex-0 fan leaves a non-convex LE cell.

**Fix.** `mesh/polygon.py`: Qhull Voronoi of hex seeds, numba clip to the
parent ring, merge 3/4-gons so faces have more than four nodes. Examples
08/09 call `polygon.fill_rib_cells`. OBJ/VTK still dump the authored cell
n-gon. Rib stroke weights match loft quads (`MESH_EDGE`).

**Verification.** `tests/test_polygon_mesh.py` (len > 4, boundary coverage,
C-cell); `test_rib_triangles_use_only_ngon_vertices`,
`test_ear_clip_le_cell_stays_inside_ngon`; 5D and later 18-key 09 galleries.
**Fixed.**

## BG-012

**Symptom.** Capping example 09 at 12 m by dropping 20 m stations beyond that
span made the tip look truncated and killed PA web-drop events (those knots
sat past 12 m on the original planform). Root chord and wall gauges stayed at
the old 20 m scale.

**Evidence.** 5D cutaway chord envelope 1.08 → 1.44 → 0.45 m after the scale
fix, vs a chopped 20 m tip. Station tables printed `z` up to 12 m with
webs still dropping on PA keys.

**Root cause.** Span, chord, and gauges were not scaled; the blade was a
prefix of the 20 m knot list.

**Fix.** Scale the original 20 m stations by `12/20` (z, chord, wall
thicknesses `SKIN`/`CAP`/`WEB_T`). Catalogue keys and loft `N_SPAN` unchanged.
Shell-topology fixtures use the same scale.

**Verification.** `test_09_extra_keys_build` (`stations[-1].z == 12`);
`test_09_te_chatter_and_tip_span`; full 09 gallery on the scaled planform.
**Fixed.**

## BG-013

**Symptom.** Example 08/09 `csg_mesh` was a 1×2 of `wall_meta` **φ**
(`contourf` RdBu) and conformal nodes. That is a distance field, not the CSG
product: on IEA 15 at the root, |φ| spans ~1.3 m while the walls are
millimetres, so the panel is a rainbow with a thin φ=0 ring. A later
“CSG solid” panel filled the OML cavity (`oml_polyline`) with un-unioned
thickness bands, so the station still read as a foil silhouette rather than
wall solids.

**Evidence.** IEA 15 `figures/iea15/csg_mesh/01.png` (circular root, φ colourbar
to 1.28 m). `5D-F-PA-10R` `csg_mesh/10.png` (φ + nodes only); later
`csg_mesh/04.png` (uniform disk) vs `2D-F_labelled.png`
(`plot_section_hybrid` / `emit_constructive_polygons`). Occupancy used a coarse
regular grid; `wall_meta` leaf occupancy was unused.

**Root cause.** `figure_csg_mesh` (`examples/_shared/blade_viz.py`) sampled a
scalar field and plotted mesh nodes. `_paint_kind_solid` then painted the foil
interior, not `emit_constructive_polygons(cut.host)`.

**Fix.** 2×3: constructive polygons in blade \(xy\) (`place_host_xy`, including
cores) | `wall_meta` occupancy (quadtree depth 10 / 150k cells) | property
union; cell bays | mid-surface nodes | property occupancy by label. Per-bay
dashed loops use the bay hue. Property tree `max_depth` later set to 10 (same
as occupancy); filled-solid `thin_band=False` remains (BG-025).
Feature-ownership and constructive-band panels were tried on 5D and dropped.

**Verification.** `test_csg_constructive_polygons_are_wall_solids_not_oml_disk`,
`test_csg_property_quadtree_is_filled_solid_not_thin_band`; regenerated
`5D-F-PA-10R` `csg_mesh/01.png`–`16.png`; example 08/09 READMEs document the
2×3. **Fixed.**

## BG-014

**Symptom.** Planform `structure.png` / `axes.png` carried a grey sentence under
the axes repeating what the legend already named (web loci, ribs, LE/TE).

**Evidence.** `_finish_planform` (`blade_viz.py`) called `note()` / fig text
below the legend. Taxonomy Z pages showed the same caption under every
LC/LE/PA pair.

**Root cause.** Caption and legend both described the same artists.

**Fix.** Remove the bottom note. Keep the legend. Regenerated example 08 and 09
`planform/` PNGs and rebuilt `taxonomy.pdf`.

**Verification.** Planform PNGs have title + legend only; Z figures in
`docs/blade/concepts/taxonomy.pdf` match. **Fixed.**

## BG-015

**Symptom.** With `LoftedPolygonSDF(..., accelerate=True)` (WP6 optional AABB /
grid cache), batched OML queries disagreed with the unaccelerated edge-loop SDF
on the same polygon.

**Evidence.** `test_lofted_polygon_bvh_matches_unaccelerated_phi` on a 2D-F loft
at `z = 4` m: 34 / 40 query points mismatched; max absolute difference about
0.23 m (relative up to ~1.8). Unaccelerated path matched `sdf_polygon`.

**Root cause.** The first accelerator reimplemented nearest-edge distance with a
uniform-grid neighbourhood walk that could miss the true closest edge while the
winding-number sign still came from the full edge set, so φ magnitude was wrong.

**Fix.** Keep the acceleration structure as cached closed-ring edges (and an
AABB / grid inventory for later), but evaluate φ with the same Numba edge loop as
`sdf_polygon` (`_sdf_polygon_numba_core` via `sdf.kernels.polygon_bvh`).
`accelerate=False` still calls `sdf_polygon` directly. Files:
`sdf/kernels/polygon_bvh.py`, `sdf/nodes/lofted.py`.

**Verification.** `test_lofted_polygon_bvh_matches_unaccelerated_phi` (atol /
rtol `1e-12`). **Fixed.**

## BG-016

**Symptom.** Cap loft strips vanished or renamed across stations where a web
was added or dropped: `spar_cap_suction_surface_0` at the root became
`spar_cap_suction_surface` once a second web survived (~`z=2.4` m on IEA 10).
`shared_loft_names` then omitted the cap between those stations.

**Evidence.** Span scan of IEA 10 loop names; close-up windows that cross web
onset showed broken cap colour. Catalogue step-down rewrote
`continuous_box` → `per_web_pads`, and pad 0 was labelled `*_0`.

**Root cause.** Unstable loft keys: `per_web_pads` indexed every pad from
`_0`, and layout resolution / catalogue re-attach replaced a published
`chord_window` `fixed_band` with pads on web drop.

**Fix.** First pad uses the unindexed `spar_cap_{side}` label. Keep
`chord_window` through `_step_down_section_layout` and
`_reattach_catalogue_preset` (do not rewrite to `per_web_pads` when a
published window exists).

**Verification.** `test_iea10_cap_loft_key_stable_along_span`;
`test_one_web_chord_window_keeps_unindexed_cap_name`;
`test_step_down_keeps_published_chord_window`. **Fixed.**

## BG-017

**Symptom.** Sparse PCHIP on station chord / twist / offset invented curvature
between knots (offset “tent” became a round hump; mid-bay chord sagged). Docs
(`CHANGELOG`, `docs/REFACTOR_PLAN.md`, `BladeGeometry.span_kind`) said optional
monotone PCHIP interpolated **chord, twist, pitch axis and offset**. Setting
`span_kind="pchip"` still left pitch linear. Default loft and goldens were
already linear, so the mismatch looked like a smoothness defect rather than an
unwired interpolator.

**Evidence.** Example `05_pchip_span` / `loft_scalars.png`: three stations at
z = 0, 5, 10 m. Mid-bay chord z = 2.5 m: linear 1.600 m vs PCHIP 1.503 m;
pitch matched. `LoftCutMixin.interpolate`: chord / twist / station offset went
through `_span_scalar` (the only `TableLaw(..., kind="pchip")` site). Pitch was
always `self._pitch(z)` from `_as_pivot_fraction` / `_as_span_table`
(`np.interp`). Override `SpanTable`s ignored `span_kind`. Airfoil morph stayed
linear / smoothstep. `tests/golden/blade.json` did not use PCHIP.
`TableLaw` defaulted to `"pchip"` while `span_kind` defaulted to `"linear"`.

**Root cause.** PCHIP is C¹ Hermite: on a sparse knot set the slopes invent
geometry between stations. Ply-drop gauges are constant, then piecewise linear,
then constant — they do not need a cubic. Pitch was never passed through
`_span_scalar`. Docs copied the loft-scalar list onto `span_kind`.

**Fix.** Span interpolation is linear only (`TableLaw` / `np.interp`;
`span_kind` and `kind=` removed). Ply-drop gauges are plateaus with linear
ramps. Example `05_span_laws` shows table plateaus, `PlyDrop`, and station
chord / twist / offset. Example 09 regenerated on the linear loft. Historical
PCHIP wording stays in `CHANGELOG` and `docs/REFACTOR_PLAN.md` as shipped
history.

**Verification.** `loft_cut._span_scalar` is `np.interp`; package `__init__`
states scalars interpolate linearly; example 05 no longer claims a cubic pitch
law; `test_span_scalars_are_linear_between_stations`; 09 gallery on the linear
11-station planform. **Fixed.**

## BG-018

**Symptom.** One full station band of bay quads folded (opposing triangle
normals): IEA 22 ≈ 226 bowties (`bay_su/pr_0/1` 47 each + `te_pr`), IEA 3.4 ≈
199 (same four bays), IEA 10 ≈ 12 (`te_pr` only). 09 / NREL had none.

**Evidence.** 47 = `n_chord − 1` for one span interval of one region. Dump
`faces` / `roles` on example 08 `samples.npz`. Typical `map_vertex` /
signature flicker (one end maps to TE and the strip twists).

**Root cause.** `loft_mapped_strips` emitted every quad between matched samples
even when the two loft triangles had opposing normals.

**Fix.** `_quad_bowtie` skips (or collapses) strips whose triangle normals
oppose (`dot < −0.2 ‖n1‖ ‖n2‖`).

**Verification.** `test_bowtie_quad_is_not_lofted`; regenerated IEA 3.4 / 10 /
22 dumps with bowtie count 0. **Fixed.**

## BG-019

**Symptom.** IEA 10 dump alone had a web–skin gap (`max_web_skin ≈ 0.042 m` at
z ≈ 17.5 m; 294 misses on `rear_web` / `web_te` from z ≈ 1.2–89 m). Catalogue
09 and NREL / IEA 15 stayed at gap 0. All dumps had `max_loop_seam = 0`.

**Evidence.** Station names were only `le_*` / `te_*` while extra webs still
lofted. Web keys were projected onto the geometric loop, but `le_su` was
resampled LE→TE without holding those keys, so the nearest skin sample sat up
to ~one chord step away. Not a T-joint / web-foot issue: contract is SHARED
nodes (a shared-junction policy).

**Root cause.** Authored windIO regions stayed LE–TE; live web anchors were
not loop partition vertices. Distinct from [BG-007](#bg-007) (keep-set /
TE-oml).

**Fix.** `ShellTopologyTemplate.from_region_graph` rebuilds loop regions from
`evaluate_anchors` at `z`, so active `web_*` keys split the skin walk. Web mid
ends reuse those loop vertex ids; samples share nodes (gap 0).

**Verification.** `test_le_te_loop_still_holds_web_keys_on_skin`; regenerated
IEA 10 `samples.npz` (`max_web_skin = 0`, loop names include `bay_*` at web
keys). **Fixed.**

## BG-020

**Symptom.** Example 08 IEA 10 root (and 1C-F / one-web stations) lofted as an
open skin C: cap did not share a vertex with TE, web read as a filled plate,
`partition_edges` under-counted the ring.

**Evidence.** At IEA 10 `z=0`, skin walls included `bay_0`/`bay_1` bonded to
the cap, but `shell_graph().loop` was only `TE / cap / LE` per side. Cap SS
met LE and sat **0.18 m** from TE (`NEAR`). `closed_loop_chain()` returned
`None`. Regenerated `iea10/root.png` before the fix.

**Root cause.** `_runs_of` in `shell_graph.py` joined same-name cap fragments
even when endpoints did not meet (footprint across two web-bounded panels).
That concatenated a jump, dropped the unpainted host remnant from the loop,
and opened the cap–skin joint.

**Fix.** Join same-name arcs only when endpoints meet (correct orientation).
Non-adjacent arcs stay separate runs with unique names (`:bay_k` suffix).
Leftover skin stays on the loop. `ShellGraph.closed_loop_chain` exposes the
contract.

**Verification.** `assert_closed_conformal_loop` / `test_shell_graph.py` (2D-F,
1C-F, `cap_windows`); `test_reference_blades_closed_conformal_loop_at_landmarks`
(all five 08 blades); IEA 10 root `partition_edges` count 6 and chain ok;
regenerated `iea10/root.png` / `closeups.png`. **Fixed.**

## BG-021

**Symptom.** On `4D-F-PA-10R` planform structure, shear web 3 sat roughly
**two** root Δx past the last primary, not one — so the secondary did not match
the primary lattice spacing.

**Evidence.** Regenerated
`examples/blade/09_blade_concepts/figures/4D-F-PA-10R/planform/structure.png`
after the uniform-lattice rewrite. Root band interiors use
`x_i = x_LE + i·Δx` with `Δx = (x_TE − x_LE) / (n_primary + 1)`; the first
extra was coded as `x_TE + Δx`, which is already one Δx outside the last
primary, so the plotted gap was 2Δx.

**Root cause.** Secondary placement restarted from root TE instead of
continuing the same arithmetic lattice (`x_LE + (n_primary + k)·Δx`).

**Fix.** `pa_web_x_from_pitch` returns
`x_LE + i·Δx` for `i = 1 … n_primary + n_secondary` (first extra lands on
root TE). Example 09 adds `5D-F-PA-10R` for the second extra at two Δx.

**Verification.** `test_four_web_pa_last_slot_is_secondary`,
`test_five_web_pa_second_secondary_is_two_dx`,
`test_pa_five_web_freeze_matches_helper`; regenerated 4D/5D PA planforms.
**Fixed.**

## BG-022

**Symptom.** After LE/PA web drop, `SectionCut.layout` was rewritten to a
sibling catalogue token for the remaining webs (e.g. `5D-F-LE-0R` at the tip
became `1C-…`). Callers that keyed on the authored concept saw a different
layout than the blade they built.

**Evidence.** `5D-F-LE-0R` `section_at(10.0)`: remaining intercept is one web
(family C) and `dropped_web_indices` is non-empty, but the cut still named
the parent as `5D-F-LE-0R` after the fix. Before, cascade helpers minted a
new token. `tests/test_loft.py::test_cut_keeps_parent_concept_when_webs_drop`.

**Root cause.** Station assemble treated the live web count as a new concept.
`layout` was the realised W/X string, not the blade’s authored key.

**Fix.** `SectionCut.layout` stays the parent concept. Realised W and X are
`station_webs` / `station_family`. Cascade helpers remain assemble recipes
for that slice (`layout_resolution`, `system_layout`); loft copies the
authored token through (`cut.py`, `loft.py`, `loft_cut.py`).

**Verification.** `test_cut_keeps_parent_concept_when_webs_drop`; assemble
contract tests no longer expect a sibling token. PR #50 (`822e833`). **Fixed.**

## BG-023

**Symptom.** Example 09 3-D mapped-strip figures (`csg_mesh/`) produced
non-finite or collapsed faces: unmatched open polylines lofted vertex-to-vertex
without resampling, and birth fans from a point onto a strip collapsed to a
line. Matplotlib autoscale then failed or showed a degenerate axis.

**Evidence.** `5D-F-PA-10R` loft bands with web onset/drop; faces with NaN or
fewer than three distinct vertices. Tests:
`tests/test_shell_topology.py` (`_finite_loft_face`).

**Root cause.** `loft_mapped_strips` required equal vertex counts on both
station polylines and kept every emitted quad, including collapsed birth fans.

**Fix.** Arc-resample unmatched open strips to a common count; drop faces that
are non-finite or have fewer than three distinct vertices
(`map_loft._finite_loft_face` / `_append_face`). Example 09 runner timeout
raised so the gallery can finish after the extra resample.

**Verification.** Finite-face unit tests; regenerated 09 loft figures; PR #47
(`895dc64`). **Fixed.**

## BG-024

**Symptom.** Example 09 `csg_mesh` column 3 (property quadtree) showed a
diamond/rotated lattice; wall occupancy in the same blade \(xy\) window was
axis-aligned.

**Evidence.** `5D-F-PA-10R` `csg_mesh/04.png`: occupancy squares vs property
cells skewed by station twist.

**Root cause.** The property tree was grown in host NACA, then
`_map_naca_chord_xy` / `cut.to_blade` rotated every leaf into the plot.
Occupancy already refines on the twisted OML in blade \(xy\).

**Fix.** Grow the property tree on the occupancy OML box (`thin_band=False`).
Map sample points back to NACA for host \(\varphi\). Paint leaves in blade
\(xy\) with no twist of the lattice.

**Verification.** `test_csg_property_quadtree_is_filled_solid_not_thin_band`
(axis-aligned leaves on a twisted cut); regenerated `5D-F-PA-10R` `csg_mesh`.
**Fixed.**

## BG-025

**Symptom.** After the lattice was axis-aligned (BG-024) and `max_depth` was
raised to 10, column 3 still shows blocky OML and staircased webs; occupancy
on the same station hugs the wall tubes.

**Evidence.** `5D-F-PA-10R` `csg_mesh/04.png` and `08.png`: occupancy
`thin_band` tubes vs property union/by-label with coarse wall bands.

**Root cause.** Property spec keeps filled-solid refine (`thin_band=False`,
`refine_band ≈ 0.012c`). That ceiling stops a coarse cell whose probes sit in
the cavity from seeing a millimetre-scale web. Occupancy uses
`thin_band=True`, `refine_band = t/2`, `min_cell ≤ t/4`. Depth 10 does not
change the band law.

**Fix.** `_property_quadtree_spec` reuses `_csg_quadtree_spec` (`thin_band=True`,
`refine_band = t/2`, `min_cell ≤ t/4`) so wall-solid φ splits millimetre webs.
Cores share that tree. [`examples/_shared/blade_viz.py`](../../examples/_shared/blade_viz.py).

**Verification.** `test_csg_property_quadtree_uses_thin_band_wall_tubes`.
Regenerated example 09 `csg_mesh` (including `5D-F-PA-10R` `04.png`: occupancy
tubes and property union share thin-band refine).
**Fixed.**

## BG-026

**Symptom.** Example 09 `5D-F-PA-10R` occupancy (and possibly property) refine
stops with `_refine_soa: max_cells=150000 reached`; the mesh is truncated.

**Evidence.** Regen of `csg_mesh/04.png`–`06.png` and `09.png`–`10.png`:
212k / 295k / 330k / 156k / 176k leaves vs cap 150k. Foil stations 07–08 and
11–16 did not warn.

**Root cause.** Occupancy uses depth 10, `thin_band`, and 150k `max_cells` on
a five-web root/transition station. Truncation drops further splits.

**Fix.** CSG occupancy (and the property tree that reuses it) uses
`_CSG_QT_MAX_CELLS = 400000`, above the observed 330k five-web root. No
per-token branch. [`examples/_shared/blade_viz.py`](../../examples/_shared/blade_viz.py).

**Verification.** `test_csg_quadtree_max_cells_covers_five_web_root`.
Regenerated example 09 `5D-F-PA-10R` `csg_mesh/01.png`–`16.png` after
`max_cells=400000` (no truncation in the gallery pass).
**Fixed.**

## BG-027

**Symptom.** Example 08/09 `planform/axes.png` drew a purple dotted “shear centre”
next to the geometric centroid. The kernel has no materials or loads.

**Evidence.** [`5D-F-PA-10R/planform/axes.png`](../../examples/blade/09_blade_concepts/figures/5D-F-PA-10R/planform/axes.png)
legend: pitch axis, centroid, shear centre. Code: `_cell_centroid_x` in
[`examples/_shared/blade_viz.py`](../../examples/_shared/blade_viz.py) area-weighted
cell-polygon \(x\), plotted as `label="shear centre"`.

**Root cause.** A cell-bay geometric proxy was misnamed. Enclosed-cell centroids
are not a thin-wall shear centre.

**Fix.** Drop `_cell_centroid_x` and the shear series. `axes.png` keeps LE/TE,
pitch, and the thin-wall (midline × thickness) geometric centroid. The three
geometric-centroid routes (boundary integral, thin-wall strip, area integration)
are compared in [`examples/blade/14_geometric_centroid`](../../examples/blade/14_geometric_centroid/).

**Verification.** `test_sample_planform_includes_rib_lines` asserts `"shear" not
in data`. Regenerated example 09 `planform/axes.png`.
**Fixed.**

## BG-028

**Symptom.** Mapped loft strips cracked along interior station lines: the strip
below and the strip above a station met at different nodes.

**Evidence.** Edge-use count on the assembled faces: edges used once on interior
station lines (`tests/test_loft_conformity.py`).

**Root cause.** `_skin_parent` resampled the lower station's arc between the
mapped keys instead of reusing that station's own nodes.

**Fix.** `_chain_nodes_between` takes the t0 station's own nodes; where the two
stations' counts differ the band takes an all-quad transition (BG-030).

**Verification.** `test_strips_share_the_station_nodes` (commit 45af12a).
**Fixed.**

## BG-029

**Symptom.** Rib and tip-closure polygon fills left hanging nodes on the skin and
web edges of their ring, and some cells had 3 or 4 nodes.

**Evidence.** `tests/test_polygon_conforming.py`: interior nodes on ring edges;
area sums short of the ring area with the angle-sort clipper.

**Root cause.** Voronoi cells were clipped with a convex angle-sort clipper that
lost area and inserted new points on the ring.

**Fix.** Exact GEOS clip, cells snapped onto ring nodes (ring nodes returned bit
for bit), a tiling check, seed spacing grown until the cells tile, 3- and 4-node
cells merged into an edge neighbour. The numba clipper is deleted.

**Verification.** `tests/test_polygon_conforming.py` (commit dbb73d5).
**Fixed.**

## BG-030

**Symptom.** Where a web or cap starts or ends, loft strips fell back to
triangles and left cracks: at the TE ridge versus TE corner key, an upper-only TE
wrap that was never lofted, and a panel born from a point whose one-node side
was dropped.

**Evidence.** Edge-use and face-size counts on the five reference blades.

**Root cause.** Unequal node counts across an event band were zipped into
triangles; the TE ridge and corner were different keys; one-node sides were
skipped.

**Fix.** Odd node counts per region (an even number of edges), a mid-span row in
event bands, `_transition_quads` (coarse sides zipped, each triangle split into
three quads through real station nodes), the TE ridge and corner aliased, an
upper-only wrap lofted as a fan, a point-born panel meshed as an apex fan.

**Verification.** `test_09_pa_4d_web3_onset_and_drop_are_skin_ruled`,
`test_transition_quads_keep_both_stations_and_drop_nonfinite` (commit b0ab78c).
**Fixed.**

## BG-031

**Symptom.** Folded loft quads (split into triangles) near the root of IEA 22,
IEA 3.4 and NREL 5MW, mostly in pressure-side cap regions.

**Evidence.** IEA 22 at z = 13.8 m: the mesh put both cap edges at chord fraction
0.29 / 0.50 on both sides; the cut's suction cap was 0.23–0.44 and pressure cap
0.34–0.53. On the pressure side the `web_1` foot crossed the window edge and
`cap_pr_2` shrank to 2 mm and flipped.

**Root cause.** Cap keys came from the region graph's single chord window for
both sides. The section builds each side's cap where the layup puts it; on a
twisted section with upright webs the two sides differ.

**Fix.** Cap keys sit on the cut's own cap partition edges per side, matched by
lo/hi role (`_cap_keys_on_side`); templates order each side's regions by that
side's keys (`side_fracs`); per-web pads use the side's web foot.

**Verification.** `test_iea22_cap_keys_are_the_cut_cap_edges` (commit a4096d1).
Folded quads NREL 5MW / IEA 22 / IEA 3.4: 32 / 84 / 84 → 0.
**Fixed.**

## BG-032

**Symptom.** A crack and overlapping strips where a cap edge passes a web foot
between two stations (IEA 22 near the tip, IEA 3.4 at 7.5 % and 95 % span).

**Evidence.** IEA 22 z = 130.9 → 134.4 m: `cap_lo` moves from the LE side of the
`web_0` foot to its TE side; the t0 arc between them was covered by `le_su`,
`bay_su_0` (reversed) and `cap_su`.

**Root cause.** The region between two keys that swap order shrinks to nothing
inside the band and grows back reversed; a ruled strip cannot represent that.

**Fix.** `key_crossing_stations` bisects each band for the span where two keys
meet and the mesh adds a station there; met keys are merged into one point; a
t0 region whose keys meet on t1 dies to that point, and its neighbours take the
outer keys.

**Verification.** `test_iea34_cap_edge_passing_a_web_gets_a_station_and_stays_conforming`
(commit 82c4fc8).
**Fixed.**

## BG-033

**Symptom.** IEA 15 near 52 % span: pressure regions ten times their length and
folded quads; the pressure web feet sat on the suction skin.

**Evidence.** The pressure side arc's mean y was +0.26 (same as suction).

**Root cause.** `_side_arcs` walked each TE corner to the LE the shorter way
round, which passed the other TE corner when the LE was nearly opposite.

**Fix.** Each corner walks the way that reaches the LE before the other corner.

**Verification.** `tests/test_side_arcs.py` (fails on the old rule; commit
56fe1b6). IEA 15 folded quads 22 → 8.
**Fixed.**

## BG-034

**Symptom.** Sliver quads (a corner under 20° or over 160°) in bands where a
region's node count changes.

**Evidence.** NREL 5MW at `h = 0.25` m, 247 stations: 15 bad quads of 17 687,
all in the one band at 17.5–17.75 m where the TE wrap closes. There the webs
drop from 17 to 7 nodes over 0.25 m for the same 1.42 m length. The parametric
2D-F blade's 7 → 5 and 5 → 3 transitions had corners of 164–173°.

**Root cause.** Three parts.
1. One count per zone, from the zone's longest length: a zone that starts at
   the root keeps the root's fine count to its far end, then drops to the next
   zone's count in one band (17 → 7).
2. TE regions bunched their nodes at the ends even when sized from `h`, so the
   coarse nodes a transition zips were unevenly spaced.
3. The transition zip compared the next coarse node's arc fraction, so ties
   went to the first station: a 5 → 3 strip got a skewed middle triangle that
   the 3 → 5 strip did not (165° against 143° for the same shape).

**Fix.** `region_counts` sets each region's count from its own length at every
station, holding it while the element size stays in `[h/2, h]`. Regions
sized from `h` are evenly spaced. The zip advances the side whose next coarse
edge is centred earlier, so it is symmetric.

Quads with a corner under 20° or over 160°, `shell_mesh(geo, h=0.25)`:

| Blade | Before (zone counts) | After |
|---|---|---|
| Parametric 2D-F-LC-4R | 0 / 1 120 | 0 / 1 062 |
| NREL 5MW | 15 / 17 687 | 3 / 16 034 |
| IEA 15 MW | 100 / 41 413 | 3 / 35 964 |
| IEA 22 MW | 222 / 75 494 | 27 / 53 738 |
| IEA 3.4 MW | 92 / 15 817 | 33 / 15 340 |
| IEA 10 MW | 221 / 35 216 | 39 / 31 425 |

On IEA 22 every quad left is where a region ends at a point: the TE wrap
closing (100 m) and the caps and bays running out at the tip, 3 cm wide over
a 0.25 m band. The region itself is a wedge of a few degrees there, so its
fan is too. A transition needs a band about one to two fine
edges long; with `region_h` much finer than the span spacing, add `refine`
windows (test `test_region_counts.py`).
**Resolved** (K8).

## BG-035

**Symptom.** IEA 15 TE wrap quads fold between 15 % and 22 % span (8 split quads
at 41 stations, 2 bands at 161).

**Evidence.** Wrap length jumps 0.12 → 0.15 → 0.12 m between stations; at 16.3 %
span the raw wrap polyline runs down, back up and down again.

**Root cause.** The section's skin mid-line folds back on itself at the blunt TE
([SG-030](../section/bug_log.md#sg-030)); the loft is faithful to it.
**Open** (section geometry).

## BG-036

**Symptom.** windIO webs run root to tip unless `start_nd_span` is given.

**Evidence.** Found reading `web_span_windows_from_windio_structure` and the web
reader during K7: a web with no span keys gets the full span.

**Root cause.** The reader takes a web's span extent only from explicit span
keys, not from its `nd_arc` grid narrowed to the plane-intersection offset grid.

**Fix.** `windio_nd.web_span_extent`: the span window is the web's `nd_arc`
grid, narrowed to the offset grid of a plane-intersection anchor. Where no web
is live (the IEA 10 root) the layout steps down to a section without webs and
without web-bound caps (`loft_layout._layout_at`). **Fixed** (K9).

## BG-037

**Symptom.** Two webs swap chordwise order between stations: IEA 15 at 10 % span
(pressure side) and at the tip, IEA 10 at 2.5 %, 5 % and 50 % span.

**Evidence.** Per-side key order of `web_*` feet differs between consecutive
stations (station key scan, 41 stations).

**Root cause.** Not established. Webs cannot cross; this comes from the windIO
data or how the reader places webs. The mesh adds a station at the crossing
(BG-032), so the skin stays conforming, but the webs meet there.

**Update (K9).** Most of these came from the reader's arc mapping (BG-038) and
from webs kept outside their span window (BG-036). With arcs placed on the real
sections, the same 41-station scan finds no crossing on IEA 10, IEA 3.4, IEA 22
or NREL 5 MW, and one on IEA 15 near the web ends (z about 113.3 m of 117 m).

**Fix.** The section's webs keep placing the file's arc on the section at that
z through the end of the span window. The wall-region anchors were a chord
fraction table of the arc grid, which stops at the last sample (112.2 m here)
and holds that fraction. The cap edge, still placed on the live section, then
passes the frozen web anchor. The anchors now use the same arc placement as the
webs, so at 113.3 m the cap stays forward of the forward web and aft of the aft
web. **Fixed.**


## BG-038

**Symptom.** Webs and cap edges read from windIO sit at the wrong chord
position: IEA 15 webs 25 to 290 mm from the points the file gives.

**Evidence.** The OML point at a web's `nd_arc`, measured along the section's
own outline, against the kernel web's foot at the same span.

**Root cause.** `windio_nd` maps arc to chord with `x/c = 1 - 2 min(s, 1 - s)`:
arc length taken as chord position, leading edge at `s = 0.5`. The same value
was used for both feet of a web (one chord position) and for both sides' caps
(one window).

**Fix.** `blade/windio_arc.py` measures each arc on the section the blade has
at that span. Webs take two feet (`web_positions` pairs; the web mid-line is the
line through the two OML points) and each side keeps its own cap window
(`CapPlacementRule.from_side_windows`). `WallRegionGraph.from_windio_structure`
takes `geo=` to place its anchors the same way. **Fixed** (K9).

## BG-039

**Symptom.** IEA 10 at `h = 1` m meshes without a tip closure once its webs
stop before the tip (BG-036).

**Evidence.** The tip is one cell; `region_counts` gave its two loop regions
three nodes each, a 4-node ring, and `polygon_mesh_ring` fills only rings of
five nodes or more (it returned no faces, without an error).

**Root cause.** The count rule had a minimum per region but none per closure
cell. A cell bounded by a web always has six nodes or more; a cell with no web
is the bare loop.

**Fix.** `region_counts` closes a web-less loop on at least `MIN_CELL_NODES`
(5) nodes; the longest region takes the extra edges. **Fixed** (K9).

## BG-040

**Symptom.** With `merge_cap_and_skin=True` the shell graph has no cap
partitions, and one cap run meshes a degenerate triangle at `h = 0.25`.

**Evidence.** Seen building a parametric blade with stacked caps; with the
option off the caps are there and the mesh is clean.

**Root cause.** Not established. **Open.**

## BG-041

**Symptom.** The windIO reader finds no spar caps on IEA 10.

**Evidence.** Caps are taken from layers whose name contains `spar` or `cap`;
IEA 10 names its layers by keypoint pair and material.

**Root cause.** Identification by name. A structural rule (the uniaxial band
between the web feet) or an explicit list is needed. **Open.**

## BG-042

**Symptom.** A blade lofted through stations taken from another surface model
sits 0.1 mm off that surface at the trailing edge of a 4.5 m chord section,
however close the stations are.

**Evidence.** At a station, `Airfoil.surface` returned y/c = 0.003659 at the
last point before the trailing edge where the given point is 0.003701. The
section flares into a blunt trailing edge over the last 0.005 % of chord.

**Root cause.** `Airfoil._monotone` samples its arc-length spline at 4001
cosine-spaced points, fine at the leading edge and 0.04 % of chord apart at the
trailing edge, then interpolates linearly: a corner sharper than that spacing
is cut.

**Fix.** The given points are kept among the samples, so each surface passes
through them. The gap above falls to 0.001 mm. The golden numbers move by at
most 4e-8 (sections, unit chord) and 4e-7 (blade, sums of coordinates), with no
change of structure; re-baselined. **Fixed** (K9).

## BG-043

**Symptom.** A root transition section (80 % thick, round at both ends) lofted
from its own points sits 0.1 mm off them at the trailing end, at any station
spacing.

**Evidence.** The section's surfaces meet at (1, 0) with a vertical tangent.
The kernel's upper surface ended at x = 1.0000129 with y = 0.0028 instead of 0,
and its last samples alternated in y.

**Root cause.** Each surface is splined on its own with a free end. At a round
end the spline bulges 1e-5 of chord aft of the end point, y(x) is no longer
single-valued there, and `_monotone` sorts the samples by x, which interleaves
the two branches of the bulge.

**Fix.** A closed trailing end that a surface reaches steeper than 45 degrees is
round: its spline gets a vertical end tangent and the surface is cut at the end
point's x (under 1e-6 of chord is left aft of it). The gap falls to under
0.001 mm. A blunt base listed as points (the surface runs on down the
base to its mid-point: 11.5 mm off on a 5.4 m flatback section) is the same
defect at a corner: the corner is found by the jump in turning angle at the
last point (over 20 degrees and three times the turn before it), the spline
stops there and the base segment is straight. Golden numbers do not move (no golden section has a round trailing
end). **Fixed** (K9).

## BG-044

**Symptom.** Web feet placed from windIO arcs are exact at the file's grid
points and up to 10 mm off between them (61.5 m blade, 30 grid points).

**Evidence.** Distance from the outer-surface points at a web's arcs to the
kernel's web line: 0.007 mm at grid points, 3 to 10 mm half-way between.

**Root cause.** The first version of `windio_arc` placed the arcs at the grid
points and interpolated the resulting chord fractions in z. windIO interpolates
the arc; the chord fraction of a fixed arc is not linear in z where the section
changes.

**Fix.** `web_feet_fn_from_windio_structure` and
`cap_side_windows_fn_from_windio_structure` interpolate the arc to z and place
it on the section at z. The same check gives 0.09 mm at most. **Fixed** (K9).

## BG-045

**Symptom.** IEA 15 at `h = 1` m: the webs start at z = 12.0 m and end at
110.0 m in the mesh; the file starts them at 11.939 m and ends them at
109.837 m.

**Evidence.** `forced_stations()` was empty for the blade; the web nodes'
z-range was the first and last mesh station inside the web's window.

**Root cause.** Two things. The ends of a partial-span web were not forced
stations. And a span window compared z to 1e-12 while stations are keyed at
`round(z, 9)`, so a window end added as a station fell outside its own window
after rounding (by 2e-10 here).

**Fix.** `forced_stations()` includes the wall-region graph's web window ends,
and span windows include their ends to the station key resolution (1e-9).
**Fixed** (K9).
