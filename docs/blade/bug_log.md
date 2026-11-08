# bug_log — blade loft

↑ [`blade`](README.md) · [docs hub](../README.md)

Worker debug log for span loft (`BladeGeometry`, `shell_graph`, mapped strips).
Each entry records the **symptom**, **evidence**, **root cause**, and
**verification**. Entries stay after they close, per
[`DOCUMENTATION_STANDARDS`](../conventions/DOCUMENTATION_STANDARDS.md).

**Identifiers.** `BG-nnn` here. Independent of `SG-nnn`
([`section/bug_log.md`](../section/bug_log.md)) and repo-root `BL-nnn`.

**Reference cases.** Example 09 catalogue keys on the synthetic 20 m planform
(`4D-F-PA-10R` for PA lattice). windIO: IEA 15 / 3.4 / 10 / 22, NREL 5 MW from
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

## Mapped loft freeze (2026-11-08)

Do not add per-blade branches (`if iea10`) in loft code.

| Rule | When |
|---|---|
| Loft webs = intercept `DroppedWeb.index` | `dropped_web_indices` is passed (`_template_from_cut`, examples 08/09) |
| No 0.08 cluster, no extra TE-plug omit | same keep-set path |
| `_cluster_web_members` (0.08, cylinder TE-plug) | only when `dropped_web_indices is None` |
| `_regularise_anchors`: `web_*` wins over nearby `cap_*` | always; do not merge two distinct webs |
| Placement graph law = section law | PA/LE freeze in metres; LC is local chord fraction |
| Skin `t0` mould walk | PA onset/drop; not a cavity chord or collapsed TE |
| TE wrap includes the OML ridge | `_te_wrap_arc`; metrics use ridge / max chord-frac, not max blade-x |

## BG-001

**Symptom.** On `4D-F-PA-10R`, `web_0` / `web_3` vanished mid-chord (ribs 5–6)
instead of walking to the LE/TE with the live intercept.

**Evidence.** Authored `WebMember.at(z)` followed catalogue LC fractions, not
`pitch_axis + k·c_root/n`. `test_09_pa_4d_loft_tracks_pitch_offset_not_catalogue`.

**Root cause.** Loft anchors used the layout token as a chord-fraction table
rather than the metre freeze (PA lattice / LE offset).

**Fix.** PA is `k·c_root/n` (`k = 1…n`); last point on the root TE may omit until
a later station intercepts. Labels keep authored indices.

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
(wrap on/off, y-span ~0.04c).

**Evidence.** Wrap skipped when suction/pressure corners coincided or wrap
length exceeded a y-span switch; t/c-scaled drop was not always on.

**Root cause.** Threshold switch on y-span, not section t/c.

**Fix.** Wrap corners interpolate along the OML-TE ridge at `_te_frac_drop`
(t/c-scaled). Blunt roots still close; sharp TE keeps a short wrap.

**Verification.** `test_ellipse_root_lofts_te_wrap`, IEA 15
`test_windio_span_metrics` (`n_coverage_holes == 0`). **Fixed.**

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
