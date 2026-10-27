# bug_log — section_geometry

↑ [`docs index`](README.md) · [`section_geometry`](../README.md)

## Role

Worker debug log for this package. Each entry records the **symptom** as
observed, the **evidence** that pinned it, the **root cause**, and — once
fixed — the **verification**. Entries stay in the log after they close, per
[`DOCUMENTATION_STANDARDS`](../../../../../docs/conventions/documentation/DOCUMENTATION_STANDARDS.md);
superseded conclusions are struck through or narrated rather than deleted,
because several of these were mis-attributed on first pass and the correction is
the useful part.

**Identifiers.** Entries here are `SG-nnn`, local to this package. `BL-nnn`
refers to the repo-root [`docs/bug_log.md`](../../../../../docs/bug_log.md),
which tracks the section-selection sizing chain; the two numbering spaces are
independent.

**Reference case throughout:** NACA 2412, 1.0 m chord, 200 airfoil points,
6 mm skin / 14 mm spar cap / 8 mm shear web, layouts named by the
[system-type taxonomy](reference/system_type_xyz_taxonomy.md) (`2D-F`, `5D-CN`, …).
Deviations are stated per entry.

## Index

| ID | Title | Status |
|---|---|---|
| [SG-001](#sg-001) | Box-cap (family **D**) midline loses its curvature during topology preparation | **Fixed** (twice — the first fix was a better guess and was superseded) |
| [SG-002](#sg-002) | Shear-web midline kinks at the foot during junction preparation | **Fixed** at the root (emit extent) |
| [SG-003](#sg-003) | Band-cap (**B**) and discrete-cap (**C**) ends are relocated onto the skin | **Fixed** (no relocation at all) |
| [SG-004](#sg-004) | Cap–web junction has no shared node at a cap **end** | **Fixed** at emit |
| [SG-005](#sg-005) | Interior cap–web T-junctions have a node on the **cap** side only | **Fixed** — one solved crossing inserted into both strips |
| [SG-006](#sg-006) | `centerline_mesh` still draws the bounding-box filter SG-001 deleted | **Fixed** — plot mirrors the resampler; the constant is gone |
| [SG-007](#sg-007) | Absolute merge / coincide floors leave a 0.89 µm cap-node residual | **Fixed** — arc-fraction merge; coincide tol 1 nm |
| [SG-008](#sg-008) | Airfoil polygonisation is the package's real exactness floor, and is un-budgeted | **Closed** — polygon *is* the section; documented |
| [SG-009](#sg-009) | Junction detection reads sampled polylines, so exactness fights detectability | **Fixed** — post-resample web-axis × final-cap solve |
| [SG-010](#sg-010) | The surface dividing suction from pressure spar caps is a 33-point polyline | **Fixed** — `n_samples` derived from a 1 µm solid-clip budget |
| [SG-011](#sg-011) | Family-D continuous-box spar caps leave a bay gap against flapwise webs under twist | **Fixed** — CapWall flapwise strip + `c_station` SoT |
| [SG-012](#sg-012) | Constructive CapWall solids split suction/pressure on `y=0`, not cavity midline | **Fixed** — store `camber_xy` on CapWall recipes |
| [SG-013](#sg-013) | Family-C (and B) pads recenter onto `argmax(depth)` instead of the web / band center | **Fixed** — `extent_anchor="band_center"` symmetric half-width |
| [SG-014](#sg-014) | Cell wall loops pair the wrong web ends on B-frame sections and carry no wall membership; every consumer re-guesses it | **Fixed** — loops are faces of the labelled wall graph and carry `wall_labels` |
| [SG-015](#sg-015) | Closed-ring suction/pressure split and the side fallback read chord-first coordinates on B-frame strips | **Fixed** — transposed at the call; side from `X_B` |
| [SG-016](#sg-016) | 50 package tests fail at HEAD (`SparCap._camber_surface`, g-row hardening) | **Fixed** — wrap copies `_camber_surface`; barrel / comment gates cleared |
| [SG-017](#sg-017) | Caps are built from the skin but their bond to it is never emitted | **Fixed** — `bond_caps_to_skin` pairs every cap vertex with the skin; `ShellMidlineStrip.bonds` |
| [SG-018](#sg-018) | Residual nearest / argmin sites vs standing invariant (Exactness item 6 follow-on) | **Logged** (2026-09-26) — no relocation defect |

## Standing invariant (established 2026-08-30)

**A wall must emit its mid-surface at its true extent, so junctions
coincide by construction. Topology splits and indexes; it never relocates.**

SG-001 through SG-004 were one failure of this: a strip emitted short or in the
wrong place forces the topology pass to move it, moving it requires guessing a
direction or a tolerance, and each guess is a bug waiting to happen. The
heuristics that exist only to support relocation are the surface area for that
class of bug. Status after applying the invariant to both paths:

| Heuristic | Was | Now |
|---|---|---|
| web foot ← `point_to_open_polyline_nearest` | relocated the foot | **relocation gone** — web emits to the skin mid-surface; the call survives in `split_skin_at_junctions` but only to read an arc length |
| cap end ← `point_to_open_polyline_nearest` + `_CAP_ENDPOINT_ON_WEB_TOL_M` | relocated free ends onto the skin | **gone** — caps pass through `skin_split` untouched, and the constant and its last consumer have since been deleted (SG-004) |
| cap **end** ↔ web node ← mesher `ENDPOINT_TOL_FRAC` (6 mm at 1 m chord) | junction conformed on ~4 mm of margin | **gone** — web emits a vertex at each cap-end junction, coincident to 0.000 µm |
| cap **interior** ↔ web node ← mesher `ENDPOINT_TOL_FRAC` | junction conformed on ~4 mm of margin | **gone** — one solved axis×cap crossing inserted into both strips (SG-005, detection primitive SG-009) |
| cap resample ← `_SKIN_RAY_CAP_BBOX_MARGIN_M` | proximity box with a fixed margin | **gone** — skin vertices selected by declared suction/pressure side, in the resampler and in the plot that draws it; the constant is deleted (SG-001, SG-006) |

The merge tolerance that decides when two skin cuts are one cut is *indexing*,
not relocation, so it stays — but it is now derived from the arc being split
(`_MERGE_S_TOL_FRAC` × half-ring length) instead of being absolute. It was the
only tolerance in the junction path that did not scale: headroom against the
minimum web spacing (`min_web_gap_frac`, 2 % of chord) fell 200x → 48x → 10x
from a 1.0 m to a 0.05 m chord, and is now constant at ~200x. Behaviour at the
1 m reference chord is unchanged (99.93 µm vs the old fixed 100 µm).

### What a grep still finds (verified 2026-09-26, Exactness item 6 follow-on)

The tables above are about *relocation*, not about whether a symbol appears.
This section covers the **whole package plus its downstream callers**, including
`tests/` and `visualisation/` — unlike the tolerance audit below, which excludes
them. Nothing in this list is a defect except where noted, so do not "clean
them up" without reading the linked entry first. Full classification: [SG-018](#sg-018).

`point_to_open_polyline_nearest` — 4 live sites:

| Site | Still there because |
|---|---|
| `midlines/ring_geometry.py` | definition. |
| `midlines/skin_split.py`, in `split_skin_at_junctions` | **indexing only.** The web foot is already the junction; this reads the arc length at which to split the skin. It no longer writes to `arr[row]`. See SG-002. |
| `topology/shell_strip_topology.py` | re-export plus a `_point_to_open_polyline_nearest` compatibility alias. No call. |
| `visualisation/demos/shell_junction_topology_demo.py`, `visualisation/demos/topology_twist_stability_demo.py` | **measurement.** Both demos use it to *report* junction gaps, which is what produced the numbers in SG-002 and SG-005. |

It is **gone** from `topology/cap_web_junctions.py` and from `section_mesh`'s
`centerline_mesh.py`, neither of which imports it any more.

Other production `argmin` / nearest (2026-09-26): `cell_wall_loop_builders._node_of`
is a 1 nm coincidence **raise**; `cap_skin_bonds._foot` constructs SG-017 bonds
(does not relocate); LE `argmin` on `X_B` identifies the polygon LE vertex;
`nearest_on_closed_polyline*` is defined with **no** production callers in this
package. SDF primitive "nearest" is field math, not strip membership.

`_SKIN_RAY_CAP_BBOX_MARGIN_M` — **no live sites.** Deleted with SG-006; the only
surviving mention is a historical sentence in
`tests/test_cap_skin_ray_resample_fidelity.py`.

Deleted during the 2026-08-30 cleanup, 0 references each: `_CAP_ENDPOINT_ON_WEB_TOL_M`,
`_MERGE_S_TOL_M`, `_axis_polyline_crossing_nearest`, `required_cap_bbox_margin_m`,
`_SKIN_RAY_BBOX_MARGIN_SAFETY`, `_SKIN_RAY_CAP_BBOX_MARGIN_M`,
`ShellStripTopologySpec.skin_ray_cap_bbox_margin_m`, and the
`skin_ray_cap_bbox_margin_m` parameter of `shell_mesh_inputs_from_section`.

**Resolution / architecture (2026-08-30 follow-up):** SG-008 closed as
documentation — the discrete airfoil polygon *is* the section. SG-010 budgets
the solids camber clip under that contract. SG-007 drops the absolute merge
floor and tightens `vertex_coincide_tol_m` to 1 nm (the 0.89 µm residual was the
coincide guard skipping a true crossing next to a resampled vertex). SG-009 is
**Fixed**: interiors are solved post-resample as web axis × final cap
(`intersect_web_axis_open_polyline_2d`); cap ends stay at emit.

### One offset operator (hardened 2026-08-30)

Everything above rests on loci at different depths being **parallel by
construction**: the skin mid-surface at `t/2`, the inner cavity at `t`, the cap
slab at `t + h`, the cap mid-surface at `t + h/2`. That holds only while every
one of them comes from the same offset settings — a mitre join is the exact
equidistant point for a locally-convex ring, so interpolation along an offset
polyline stays on the locus, which a different join style would not give.

`buffer(-d, quad_segs=16, join_style=2, mitre_limit=5.0)` used to be hand-written
at **10 sites across 5 files**, so that contract was enforced by copy-paste: one
literal drifting would silently un-parallel two loci and reproduce the SG-001 /
SG-003 class with nothing to catch it. All ten now go through
`constructive/helpers.py::offset_polygon_inward`, which owns the three settings
as private constants. Verified separately that `mitre_limit = 5.0` never falls
back to a bevel (identical geometry against `mitre_limit = 1e6`) across
`n_points` 200/400, chords 1.0 / 0.241 m and depths 0.72–40 mm, so the limit is
not a hidden exactness cliff.

SG-001 through SG-004 were all found while building
[`visualisation/demos/geometry_to_mesh_flythrough.py`](../visualisation/demos/geometry_to_mesh_flythrough.py),
which draws raw emit and junction-prepared strips side by side at a junction
zoom — the two forms are hard to compare any other way.

## Tolerance audit (2026-08-30)

Sweep of every tolerance-like constant and inline numeric comparison in the
package (**excluding** `tests/` and `visualisation/`), classified by what the
number stands in for. A tolerance is only a defect when it substitutes for a
quantity that is knowable exactly.

**(a) Substituting for a known exact answer — the bug class.**

| Site | Verdict |
|---|---|
| `_SKIN_RAY_CAP_BBOX_MARGIN_M` proximity box | fixed — declared surface side (SG-001) |
| `_CAP_ENDPOINT_ON_WEB_TOL_M` deciding whether to relocate a cap end | fixed, then **deleted** — nothing relocates, so nothing decides (SG-003, SG-004) |
| web foot ← nearest-point projection | fixed — emit extent corrected (SG-002) |
| cap↔web node ← mesher `ENDPOINT_TOL_FRAC` | fixed — cap **ends** at emit (SG-004), cap **interiors** from one solved crossing (SG-005) |
| `_CAP_WEB_JUNCTION_MERGE_ARC_M` / 1 µm `vertex_coincide_tol_m` deciding a crossing is already noded | **fixed** — merge is arc-fraction; coincide default is 1 nm (SG-007) |

**(b) Correctly scaled design rules — not numerical tolerances at all.**

`min_web_gap_frac = 0.02` (webs closer than 2 % of chord collapse into one) and
`min_cavity_depth_m` are structural policy, chord-relative, and belong where
they are.

**(c) Legitimate numerics — correct as they stand.**

- `inner_y_at_x`: 240 samples + 10 bisections to locate inner-skin crossings.
  A polygon-offset crossing has no closed form, so root-finding *is* the exact
  answer available. Convergence, not classification.
- `abs(twist) > 1e-10`, `area > 1e-16`, `1e-9` ring-closure checks,
  `_MIN_SEGMENT_LENGTH_M = 1e-12`, `_SKIN_RAY_T_MIN_M = 1e-9` — floating-point
  is-zero and degenerate-geometry guards.
- `vertex_coincide_tol_m` in `insert_cap_web_interior_junctions` — asks whether a
  strip already has the node. Default was 1 µm and masked SG-007; now 1 nm so it
  is de-duplication of floating-point repeats, not classification of nearby
  geometry.

`_CAP_WEB_JUNCTION_MERGE_ARC_M` was listed here on the same grounds until SG-007
showed it is not: its absolute floor decides that two distinct crossings are one,
which *is* classification. Moved to (a), then fixed by deleting the absolute
floor the same way `_MERGE_S_TOL_M` was removed.

**(d) Scale-dependent but safe — fixed anyway.**

The skin-cut merge tolerance was the only one in the junction path expressed in
absolute metres while everything around it was chord-relative. It is now
`_MERGE_S_TOL_FRAC` × arc length. The absolute `_MERGE_S_TOL_M` constant that
used to sit beside it had no consumers once the derived form landed and has been
deleted; a caller that needs an absolute value still passes `merge_s_tol_m`
explicitly to `split_skin_at_junctions`.

### Cleared by measurement (2026-08-30)

Candidates that looked like sampling error and are not. Recorded so they are not
re-investigated: refining any of these buys nothing.

| Candidate | Measured |
|---|---|
| Cap midline sampling, `n_cap_samples = 80` | **0.00 µm** against a 5000-sample reference (`2D-F` suction cap). The loci are polygon offsets — piecewise linear — and interpolation along a segment stays on the locus, so `n` sets node density, not fidelity. Likewise `n_web_samples = 20` on a straight axis. |
| `inner_y_at_x`, 240 samples + 10 bisections | **<0.005 µm** against 20000 samples / 60 bisections, at x/c = 0.15–0.70. φ is piecewise linear along a vertical line, so the linear interpolation already *is* the exact crossing. |
| `inner_y_at_xs(..., bisection_iters=0)` in `cap_feasibility_clamp` | **0.00 µm** against 60 iterations at the same stations, for the same reason. The shortcut costs nothing, and its comment's "feasibility only needs clamp accuracy" caveat is unnecessary. |
| `mitre_limit = 5.0` silent bevel fallback | **Never triggers** — see "One offset operator" above. |
| Section-properties backend | Defaults to `hybrid`, i.e. constructive polygons integrated by shoelace, which is exact for a polygon, with the quadtree as cross-check. Area parity against dense-grid SDF occupancy is already instrumented in `docs/reports/DUAL_REP_ACCEPTANCE.md`. |

### The trap this class hides

Making one side of a junction exact can break the *other* side's detection.
Placing a web vertex precisely on a cap turns a polyline×polyline crossing into
a tangency, and `intersect_open_polylines_2d` then finds nothing. SG-005 records
the measured swap on `4D-F`: web gap 1570 µm → 0, cap gap 0 → 1056 µm. Junction
exactness has to be established for both participants in the same pass, from one
solved crossing, after both polylines are final. SG-009 removes the class of
trap for *detection* by solving interiors from the web axis rather than from
polyline×polyline crossings; the shared-node insert rule from SG-005 remains.

---

## SG-001

**Box-cap (family D) midline loses its curvature during topology preparation**

**Status:** Fixed — twice; the first fix was a better guess and was superseded
**Found:** 2026-08-30
**Area:** `midlines/cap_skin_ray_resample.py`

### Symptom

The family-**D** box-cap midline loses its curvature between the raw emit and
the junction-prepared strip. The mesh built on it is wrong, not just the picture.

### Evidence

`2D-F`, 6 mm skin / 14 mm cap: the cap goes 105 → **8** points, sagitta
10.139 → **6.121 mm**, and the middle ~80 % of the cap becomes one straight
chord. Mesh nodes sit **5.968 mm** off the true cap midline — **43 % of cap
thickness**. Refining does not help; the elements subdivide the wrong chord.
Only **6 of 24** skin vertices above the cap survive the bbox pre-filter, all
clustered at the cap ends.

### Root cause

`cap_skin_ray_resample.py` rejected skin vertices outside the cap bounding box
plus a **fixed** `_SKIN_RAY_CAP_BBOX_MARGIN_M = 5 mm`, but the skin midline sits
`(skin_t + cap_t)/2` from the cap midline — 10.0 mm at these gauges. The margin
was narrower than the thing it was filtering, so the rays had almost no origins.
Margin sweep: 5 mm → 8 pts / 5.971 mm error; 10 mm → 25 pts / 0.071 mm;
12 mm → 26 pts / 0.002 mm.

### Verification

First fixed by deriving the margin per cap from its own gauge
(`required_cap_bbox_margin_m` = `1.25 × (cap_t + skin_t)/2`): `2D-F` 26 pts,
sagitta 10.139 mm exact, mesh node error **0.0012 mm**.

That was a better *guess*, so it was superseded: the bounding box is gone
entirely, and cap resampling now selects skin vertices by the strips' **declared**
suction/pressure side. `required_cap_bbox_margin_m` and
`_SKIN_RAY_BBOX_MARGIN_SAFETY` have been deleted (0 references).

Pinned by `tests/test_cap_skin_ray_resample_fidelity.py`, which asserts the
*outcome* (prepared cap tracks raw cap) rather than the mechanism, and so
survived the rework unaltered.

### Cleanup, 2026-08-30

`ShellStripTopologySpec.skin_ray_cap_bbox_margin_m` and the matching
`shell_mesh_inputs_from_section(..., skin_ray_cap_bbox_margin_m=...)` parameter
outlived the mechanism they configured. `build_shell_strip_topology` never read
the field, so a caller could set the margin and get silence — the field's own
docstring admitted this, but the handoff signature did not. Both deleted; no
caller passed either. `_SKIN_RAY_CAP_BBOX_MARGIN_M` itself outlived them by one
consumer, a plot, and went with SG-006.

---

## SG-002

**Shear-web midline kinks at the foot during junction preparation**

**Status:** Fixed at the root (emit extent)
**Found:** 2026-08-30
**Area:** `sections/components/shear_web.py` — `ShearWeb.midline_polyline`

### Symptom

The projected segment that reaches the skin is not parallel to the rest of the
web.

### Evidence

The raw web is a capsule axis, straight to **0.000 µm**. After preparation the
last segment leans **1.03–3.43°** off the web axis at zero twist and
**3.92–8.45°** at 15° twist, putting the junction node **0.14–1.19 mm** off the
web axis on an 8 mm web (~15 % of web thickness).

Controls: a symmetric NACA 0012 gives identical error at both ends of a web
(2.14° / 0.297 mm), and `2D-F` ≡ `2D-CN` at zero twist, where flapwise and
chord-normal coincide.

### Root cause

**An emit-extent error, not a topology error.** `ShearWeb.midline_polyline`
interpolated between the two **inner-mould-line anchors** — the right locus for
building the *solid*, but half a skin gauge short of where the web's
**mid-surface** ends. Measured: the anchor sits at exactly 0.000 mm from the IML
and −3.000 mm from the skin mid-surface on a 6 mm gauge.

Because the strip arrived short, `split_skin_at_junctions` had to **relocate**
the foot, and relocating requires guessing a direction. It guessed perpendicular
to the skin (`point_to_open_polyline_nearest`), which equals the web axis only at
exactly 90° — never on a cambered or twisted section.

### Verification

`ShearWeb.midline_polyline` now emits the web mid-surface out to the skin
mid-surface, solving for the step along the web axis (`half_gauge / cos θ`;
measured 3.037–3.042 mm for a 3.000 mm half gauge, θ ≈ 9.6°) rather than
assuming it.

Raw web endpoints now land on the skin mid-surface to **0.000 µm** before any
topology pass, so `split_skin_at_junctions` no longer moves them: it only
*indexes*, finding the arc length at which to split the skin. Measured across
`2D-F`/`5D-F`/`3D-CN` at 0–15° twist: topology relocation **0.00 µm**,
foot-to-skin-node **0.00 µm**, kink **0.0000°**.

Pinned by `tests/test_web_foot_axis_snap.py` — 40 cases
(`test_prepared_web_stays_straight`, 5 layouts × 2 twists × 3 airfoils, plus
`test_web_foot_lands_exactly_on_a_skin_node`, 5 layouts × 2 twists).

~~An earlier fix made topology's *guess* correct, by intersecting the web axis
with the skin instead of projecting perpendicular to it.~~ Superseded by the
above: it removed the symptom, not the need to guess.

### Follow-on, 2026-09-27 (WP6 tip mono)

Confirmatory 2D-F e2e crashed on tip mono_after_redispatch: frozen thin-skin
`web_y` + raised skin (~7 mm) + tip twist made `flapwise_web_sframe_inputs` miss
B-frame IML intercepts and **rotate** the stale anchors; emit then returned a
blind `origin + t/2` foot (~3.5 mm off mid — SG-002 in `build_cell_wall_loops`).

Fix: refuse the rotate fallback (`InfeasibleSectionGeometry`); harden
`_extend_to_skin_mid_surface` (longer oracle / bidirectional SDF; no blind
offset); mono gauge projection keeps prior gauges when raise is infeasible.
Pinned by `test_tip_absolute_2a_f_loops_build_under_twist`,
`test_tip_fractional_2a_f_loops_build_under_twist`, and
`test_stale_frozen_web_y_under_twist_raises_infeasible_not_sg002_gap`.

### Cleanup, 2026-08-30

The superseded fix left its implementation behind —
`_axis_polyline_crossing_nearest` in `midlines/skin_split.py`, 45 lines,
defined and never called. Deleted. Its docstring was the clearest statement in
the package of *why* perpendicular projection is wrong, so that reasoning has
been folded into this entry rather than lost with it.

---

## SG-003

**Band-cap (B) and discrete-cap (C) ends are relocated onto the skin**

**Status:** Fixed — no relocation at all
**Found:** 2026-08-30
**Area:** `midlines/skin_split.py` — `split_skin_at_junctions`

### Symptom

Family-**B** and family-**C** cap midlines are displaced and their arc length
inflates during topology preparation. **Not** the same defect as SG-001 —
widening the bbox margin changes nothing here.

### Evidence

`2B-F` suction cap: `split_skin_at_junctions` moves the first point
`X_B 0.06373 → 0.07371` (**+9.98 mm**) and the last `0.06611 → 0.07611`
(**+10.0 mm**); arc **100.05 → 117.67 mm (+17.6 %)**. `2C-F` the same
(81.12 → 99.16 mm). `2D-F` is untouched (105 pts, arc identical).

Displacement ≈ `(skin_t + cap_t)/2` — the cap end is landing **on the skin
midline**. Ray-resample counters confirm it is not the margin: raising it
5 → 25 mm lifts `bbox_pass` 7 → 11 but `ray_hit` stays pinned at **7**, and the
endpoint move happens before the resampler runs. Residual error vs the raw cap
5.6–9.5 mm across gauges, unchanged by margin.

### Root cause

For each cap endpoint, if it was not within `_CAP_ENDPOINT_ON_WEB_TOL_M` of a
web polyline, the code snapped it onto the **nearest skin midline**. Family-D box
caps terminate on webs, so `on_web` was true and they were skipped — which is why
D never showed this. Families B and C terminate mid-bay with free ends, so both
ends got dragged ~`(skin_t + cap_t)/2` outward onto the skin, adding two diagonal
hooks. A free cap end should stay on the cap's own mid-surface; the snap is only
meaningful when the end really is at a junction.

### Verification

The answer was *no relocation at all*. Measurement settled it: a family-D cap end
is already at **0.000 mm** from the web mid-surface, and a B/C free end is
**38.6–148.6 mm** from the nearest web across 0–15° twist. The two populations
are four orders of magnitude apart, so there was never a tolerance to tune — only
a decision not to move anything.

`split_skin_at_junctions` now passes caps through untouched. Cap endpoint
displacement **10 mm → 0.0000 mm**, arc-length inflation **+17.6 % → 0.00 %**,
and B/C cap fidelity vs the raw curve **3.7–9.5 mm → 0.0000–0.0017 mm**.

B and C are now covered by `tests/test_cap_skin_ray_resample_fidelity.py`, which
also pins the invariant directly (`test_topology_never_relocates_a_cap_endpoint`).

Note that the SG-002 fix did **not** transfer: a cap end runs **0.00–0.95° from
parallel** to the skin, so "extend along its own axis to meet the skin" is
undefined for a cap — same call site, opposite correct answer.

---

## SG-004

**Cap–web junction has no shared node at a cap end**

**Status:** Fixed at emit
**Found:** 2026-08-30
**Area:** `midlines/shell_midline_export.py`

### Symptom

The cap end sat on the web mid-surface but between two web vertices. Latent, not
a live defect — output was correct.

### Evidence

Cap end to web **mid-surface** 0.000 mm, but to the nearest web **node**
0.624 mm (`2D-F`), 0.958 mm (`5D-CN`), **1.880 mm** (`3D-F`). The junction held
together only because `section_mesh` clusters endpoints within
`ENDPOINT_TOL_FRAC = 0.006` — **6 mm on a 1 m chord** — leaving `3D-F`
conforming on **4.1 mm of margin**, in a tolerance that scales with *chord*, not
with mesh density. Refining the web, shrinking the chord or raising the web count
erodes it silently, with no error raised.

### Root cause

The web is sampled independently of the caps, so nothing put a vertex where a cap
lands on it. `insert_cap_web_interior_junctions` patched this after the fact,
using `point_to_open_polyline_nearest` a fourth time.

### Verification

`midlines/shell_midline_export.py` (`_insert_web_cap_junction_nodes`) now inserts
a vertex into each web where a cap mid-surface meets it, at emit, where both are
in hand. Exact: a web midline is straight and the cap end already lies on it, so
the inserted vertex changes no geometry. Cap-end-to-web-node
**0.62–1.88 mm → 0.000 µm** across `2D-F`/`2D-CN`/`3D-F`/`4D-F`/`5D-CN`. The
on-axis membership test is a tolerance but not a fragile one — it separates
0.000 mm from 38.6–148.6 mm.

Only the **web** gets nodes at emit. Inserting into the cap there would be undone
by `resample_cap_interiors_from_skin_rays`, which rebuilds the cap interior
afterwards.

### Cleanup, 2026-08-30

Fixing the web side at emit made the cap side of the same junction unreachable.
`insert_cap_web_interior_junctions` carried a cap-endpoint-on-web branch that
inserted the nearest web-midline point for any cap end within
`_CAP_ENDPOINT_ON_WEB_TOL_M` of a web. Both of its guards now always reject:

- a family-D cap end is coincident with a web node to 0.000 µm, so the
  `vertex_coincide_tol_m = 1e-6` guard skips it;
- a family-B/C free end is 38.6–148.6 mm from the nearest web, so the
  `> _CAP_ENDPOINT_ON_WEB_TOL_M` guard skips it.

Confirmed by disabling the branch and diffing strip arrays across **108
configurations** (12 layouts × 3 twists {0°, 15°, 34.11°} × 3 chords
{1.0, 1.656, 0.241 m}): **byte-identical output, 0 strips differing**. The
branch and `_CAP_ENDPOINT_ON_WEB_TOL_M` are deleted, and
`cap_web_junctions.py` no longer imports `point_to_open_polyline_nearest`.

This is worth knowing before attacking SG-005: the earlier revision of this log
identified the deleted branch as SG-005's code, which it never was. SG-005 lives
in the **interior segment-intersection** path at the top of the same function.

**Follow-up worth checking:** three `section_shell_fea` cap-tie tests moved from
xfail to xpass (`test_untied_cap_s_and_q_diverge`,
`test_untied_cap_fi_diverges_so_the_test_discriminates`,
`test_connected_components_converge_regardless_of_cap_ties`), all tagged
[BL-014](../../../../../docs/bug_log.md#bl-014) in the repo-root log. They are
`strict=False` and documented as execution-context dependent. That is a signal
to run the **BL-014 triage track** (isolation-vs-suite polluter hunt — see
[BL-014 triage](../../../../../docs/bug_log.md#bl-014-triage)), not a claim the
tied-web Hashin FE defect is resolved.

---

## SG-005

**Interior cap–web T-junctions have a node on the cap side only**

**Status:** Fixed
**Found:** 2026-08-30
**Area:** `topology/cap_web_junctions.py` — `insert_cap_web_interior_junctions`

### Symptom

Where a web crosses a cap's **middle** rather than its end, only the cap carries
a node at the crossing.

### Evidence

Exact web-axis × cap-polyline crossings, measured after the full topology pass:

| Layout | Interior crossings | Cap-node gap | **Web-node gap** |
|---|---|---|---|
| `3D-F` | 2 | 0.0 µm | **481.8 µm** |
| `2C-F` | 2 | 0.0 µm | **624.0 µm** |
| `4D-F` | 4 | 0.0 µm | **1569.6 µm** |
| `5D-CN` | 6 | 0.0 µm | **1949.4 µm** |

Widened later to 9 layouts × 3 twists × 3 chords (81 configurations, 789 interior
crossings), the worst web-node gap was **4539.1 µm** — more than twice the worst
in the table above, which had sampled only the 1 m chord.

Conforming only through `section_mesh`'s `ENDPOINT_TOL_FRAC` (6 mm at 1 m
chord) — **4.05 mm of margin at worst**, in a chord-scaled tolerance. Affects
layouts where a web crosses a cap interior: family D with 3 or more webs, and
family C.

### Root cause

Webs and caps are sampled independently. `insert_cap_web_interior_junctions`
closes the cap side by intersecting the two polylines — verified: it inserts
exactly one cap vertex per interior crossing (`3D-F` +1 per cap, `4D-F` +2,
`5D-CN` +3, `2C-F` +1 per cap, `2B-F` +0, which has no interior crossings). But
nothing puts the matching node on the web.

### Why the obvious fix fails

An attempt to insert the web node at emit — alongside the cap-*end* junctions,
which `_insert_web_cap_junction_nodes` does fix — **made it worse and was
reverted**. Placing a web vertex exactly on the cap makes the two polylines
*touch* rather than cross, so `intersect_open_polylines_2d` stops finding the
junction and the cap side regresses. Measured on `4D-F`: web gap 1570 µm → 0,
cap gap 0 → 1056 µm — a clean swap of which side is broken.

The same mechanism is visible in correct output: re-running
`intersect_open_polylines_2d` on prepared strips finds **nothing**, because the
node the pass inserted has already turned every crossing into a tangency. That is
expected, and is not evidence the pass did nothing — but it does mean neither a
verification nor a test may use that detector.

### Fix

`insert_cap_web_interior_junctions` now returns modified **webs as well as caps**.
Each crossing is solved once, and the same point is inserted into both strips, so
the junction is conformal by construction rather than by the mesher clustering
two nearby nodes.

Three things make it exact:

- **Detection runs against the original polylines**, and every insertion is
  applied only afterwards. Because an inserted point already lies on the edge it
  is inserted into, neither strip's shape changes and no crossing is invalidated
  — so the tangency trap above never arises.
- `intersect_open_polylines_2d` now returns the **b-side edge index** as well as
  the a-side. A junction has two participants, and the caller placing a shared
  node needs to know which edge of each the crossing lands on. Recovering the
  second index afterwards would mean searching for the point with a tolerance,
  which is the failure mode this package keeps designing out.
- A crossing earns its web node only if the **cap actually carries it** — the
  pending list is committed against the cap's final polyline. A crossing dropped
  from the cap (slimmed away as too short) would otherwise get a web node alone,
  moving the defect rather than fixing it.

Measured across 9 layouts × 3 twists (0°, 15°, 34.11°) × 3 chords
(1.0, 1.656, 0.241 m) — 81 configurations, 789 interior crossings:
worst web-node gap **4539.1 µm → 0.0000 µm**, worst cap-node gap then still
at 0.8924 µm (SG-007 — absolute coincide / merge floors; fixed separately).

Pinned by `tests/test_cap_web_interior_junction_shared_node.py`, 40 cases. It
locates crossings from the **web axis**, deliberately not with
`intersect_open_polylines_2d`: a test built on that detector would pass vacuously
on a correct fix *and* on one that deleted every junction. Verified to
discriminate — against the pre-fix code all 20 web-side cases fail and all 20
cap-side cases pass, which is also independent confirmation that the cap side was
never the broken one.

**Later (SG-009):** production detection switched from
`intersect_open_polylines_2d` to `intersect_web_axis_open_polyline_2d`. The
detect-then-insert / pending-web half-junction rules above still apply; only the
crossing primitive changed. See [SG-009](#sg-009).

---

## SG-006

**`centerline_mesh` still draws the bounding-box filter SG-001 deleted**

**Status:** Fixed
**Found:** 2026-08-30
**Area:** `mesh_utilities/section_mesh/visualisation/centerline_mesh.py`

### Symptom

`_add_skin_inward_rays_class_c_near_caps` is documented as drawing "Class C ray
directions used by cap interior resampling", but it selects which skin vertices
get an arrow by reimplementing the cap bounding box plus
`_SKIN_RAY_CAP_BBOX_MARGIN_M` — the exact 5 mm proximity box that SG-001 proved
was the defect and removed from the resampler.

### Evidence

Skin vertices that would be drawn with an inward ray, old filter vs what the
resampler actually considers (NACA 2412, 1 m chord, zero twist):

| Layout | Skin vertices | Old bbox arrows | Declared-side arrows |
|---|---|---|---|
| `2D-F` | 187 | 16 | 187 |
| `3D-F` | 189 | 45 | 189 |
| `2C-F` | 187 | **6** | 187 |
| `5D-CN` | 193 | 65 | 193 |

`2C-F` drew **6 of 187** — 3 % of the ray origins the resampler uses. That is the
SG-001 failure mode itself, rendered as though it were the algorithm.

Geometry, strips and mesh were unaffected throughout — this code path only emits
arrows.

### Fix

`_add_skin_inward_rays_class_c_near_caps` now mirrors
`_resample_one_cap_interior_from_skin_rays`: a skin vertex contributes a ray when
its declared surface side matches a cap's, aimed inward at that cap's centroid,
with no proximity re-filter. It imports from `cap_skin_ray_resample` directly
rather than through `shell_strip_topology`, so the plot and the algorithm now
share a module as well as a rule. The quiver label changed from "Class C ray dir
(perp tangent, in cap bbox)" to "(perp tangent, declared side)".

With the last consumer gone, `_SKIN_RAY_CAP_BBOX_MARGIN_M` and its re-export from
`topology/shell_strip_topology.py` are deleted — 0 references in code.
Covered by `section_mesh/tests/test_centerline_mesh_plots.py`, which passes
unaltered: it asserts the figure renders, not what it filters.

---

## SG-007

**Absolute merge / coincide floors leave a 0.89 µm cap-node residual**

**Status:** Fixed — arc-fraction merge; coincide tol 1 nm
**Found:** 2026-08-30, while verifying SG-005
**Fixed:** 2026-08-30
**Area:** `topology/cap_web_junctions.py`

### Symptom

On a short / densely resampled cap, an interior web crossing can leave the cap
node **0.89 µm** off the web axis even after SG-005 inserts shared nodes.

### Evidence

Worst case over the SG-005 sweep (81 configurations, 789 interior crossings):
**0.8924 µm**, on `5C-CN` at zero twist, chord 0.241 m,
`spar_cap_suction_surface_1`. Measured identically before and after the SG-005
fix.

On that case the then-production polyline×polyline detector finds one crossing
whose solved point lies **0.89 µm** from an existing skin-ray-resampled cap
vertex. With `vertex_coincide_tol_m = 1e-6` the insert was skipped, so the cap
kept the nearby sample instead of the true junction. Separately, the absolute
merge floor `_CAP_WEB_JUNCTION_MERGE_ARC_M = 1e-6` could also collapse two
crossings closer than 1 µm on a short cap (same absolute-tolerance pattern as
the old `_MERGE_S_TOL_M`).

### Root cause

Absolute 1 µm floors in a path where everything else scales with the section:
`max(1e-7 * arc_total, 1e-6)` for merge, and `vertex_coincide_tol_m = 1e-6` for
“already noded”. On a 19 mm tip-scale cap the intended arc-fraction merge is
~2 nm, so the absolute floor dominated; the coincide guard was the direct cause
of the measured `5C-CN` residual.

### Fix

Two absolute floors were involved:

1. **`_CAP_WEB_JUNCTION_MERGE_ARC_M = 1e-6`** — deleted. Default merge is now
   `_CAP_WEB_JUNCTION_MERGE_ARC_FRAC * arc_total` (`1e-7`), matching the skin-cut
   `_MERGE_S_TOL_FRAC` pattern. Callers that need an absolute window still pass
   `merge_arc_m` explicitly.

2. **`vertex_coincide_tol_m = 1e-6`** — the residual on `5C-CN` was this, not
   merge de-dup. After skin-ray resample, a cap vertex already sat **0.89 µm**
   from the true web×cap crossing; the coincide guard treated that as “node
   already present” and skipped inserting the solved point, leaving the cap
   0.89 µm off the axis. Default lowered to **1 nm** (same order as the pass’s
   `min_seg` and the shared-node test).

Checked against the shell handoff: `section_mesh` rejects only segments shorter
than `_MIN_SEGMENT_LENGTH_M = 1e-12` m, and this pass already refuses insertions
that would create edges shorter than its own `min_seg` default (1 nm). A
sub-µm cap segment between a resampled vertex and the exact crossing is
therefore acceptable. Cap-side shared-node test tightened from 1 µm to 1 nm.

---

## SG-008

**Airfoil polygonisation is the package's real exactness floor, and is un-budgeted**

**Status:** Closed — decision: the polygon *is* the section (reading 1)
**Found:** 2026-08-30
**Closed:** 2026-08-30
**Area:** `sdf/airfoil.py` — `AirfoilSDF.from_naca(..., n_points=...)`, and every
call site that picks `n_points`

### Symptom

Junctions are now exact to 0.000 µm against a shape that is itself tens of
microns away from the aerofoil it claims to be. Nothing in the package measures,
pins or documents that gap, and `n_points` is chosen per call site by habit.

### Evidence

Maximum deviation of the airfoil polygon from the analytic NACA 2412 curve,
1 m chord, in the leading-edge region (x/c < 0.02) where curvature is highest:

| `n_points` | LE deviation |
|---|---|
| 100 | 246.8 µm |
| **200** (the usual default) | **62.1 µm** |
| 400 | 15.5 µm |
| 800 | 3.9 µm |

Clean second-order convergence — halving the step quarters the error, exactly as
a chord approximation of a curve should. So this is discretisation, not a defect
in the sampling.

For scale, 62 µm is **~1 % of a 6 mm skin gauge**, at the most curved part of the
skin, and roughly 70x the worst junction residual that remained while SG-007 was
open (0.89 µm).

Measurements away from the LE are dominated by a trailing-edge closure
difference between the SDF's airfoil and an open-TE analytic reference; that
part is a reference artefact, not a package error, and was excluded.

### Decision

**The polygon *is* the section.** The package's input contract is a discretised
aerofoil; everything downstream is exact with respect to it. The ~62 µm LE gap
at default `n_points=200` is how closely the NACA *generator* matches its closed
form — not a junction defect, and not something the topology path should chase.
Raising `n_points` refines the polygon when a caller wants a finer outline.

Documented in the package README Notes, `sdf/README.md`, and the
`AirfoilSDF.from_naca` docstring. Reading (2) — budget `n_points` against an
analytic LE error — is deliberately not taken; reverse that only with an
explicit product decision.

---

## SG-009

**Junction detection reads sampled polylines, so exactness fights detectability**

**Status:** Fixed — post-resample web-axis × final-cap solve
**Found:** 2026-08-30, generalising SG-005
**Fixed:** 2026-08-30
**Area:** `topology/cap_web_junctions.py`, `midlines/ring_geometry.py`

### Symptom

Making a junction exact destroys the evidence used to find it. Placing a shared
node at a crossing turns two crossing polylines into two *touching* polylines,
and `intersect_open_polylines_2d` — a strict-interior detector — then returns
nothing.

### Evidence

This has now bitten twice, in opposite directions:

- SG-005's first attempt inserted the web node at emit, which made the cap side
  undetectable and swapped the defect (`4D-F`: web 1570 µm → 0, cap 0 → 1056 µm).
- After SG-005 was fixed correctly, re-running the polyline×polyline detector over
  the prepared strips finds **zero** interior crossings across 81 configurations —
  correct behaviour that looks identical to "the pass did nothing".

The second point is the sharper problem: it means neither a verification script
nor a regression test may use polyline×polyline as the existence proof, because
such a test passes on a correct fix *and* on one that deleted every junction.

### Root cause

Interior junctions were recovered from **sampled polylines** via
`intersect_open_polylines_2d` rather than solved from the web axis against the
final (post-resample) cap. Cap interiors cannot close at emit because ray
resample rebuilds them.

### Fix

`insert_cap_web_interior_junctions` now solves interiors with
`intersect_web_axis_open_polyline_2d` (web axis × final cap), then inserts the
same point into both strips — same insertion / pending-web machinery as SG-005.
Cap-end junctions stay at emit (`_insert_web_cap_junction_nodes`).
`intersect_open_polylines_2d` remains for tools; it is not topology truth.

Design record:
[`reports/SG009_JUNCTION_PRIMITIVES_SPIKE.md`](reports/SG009_JUNCTION_PRIMITIVES_SPIKE.md).

### Verification

- `test_cap_web_interior_junction_shared_node.py` — LAYOUTS × twist × chord, 1 nm
  gaps on both strips; oracle is the same axis×cap primitive.
- Related web-foot / shell-strip / cap-resample fidelity suites green.
- BL-014 triage / tied-web Hashin FI (out of scope for this rewrite).
- Pre-sample analytic offset-ring × axis left out of scope.

---

## SG-010

**The surface dividing suction from pressure spar caps is a 33-point polyline**

**Status:** Fixed — `n_samples` derived from a 1 µm solid-clip budget
**Found:** 2026-08-30
**Fixed:** 2026-08-30
**Area:** `sections/cap_feasibility.py` → `constructive/emitters.py::_camber_clip`
and `sections/components/spar_cap.py::_camber_side_split` (consumed via
`outer_skin.py`)

### Symptom

Which material belongs to the suction cap and which to the pressure cap is
decided by a coarse sampled polyline standing in for the cavity midline.

### Evidence

`cap_feasibility_clamp` samples the cavity depth at `n_samples = 33` stations
across the cap band and returns those stations' midpoints as
`info["camber_xy"]`. That 33-point polyline is then used directly as a clip
boundary — by `_camber_clip` in the constructive path and by
`_camber_side_split` in the SDF path.

Deviation of the 33-segment polyline from a densely sampled cavity midline,
NACA 2412, 6 mm skin, band x/c 0.15–0.65:

| Chord | Deviation |
|---|---|
| 1.0 m | **18.56 µm** |
| 0.241 m | 4.47 µm |

The error is a chord approximation of a curve, so it is systematic rather than
oscillating: it biases area consistently from one cap to the other rather than
cancelling. On a 14 mm cap that is ~0.13 % of cap height.

Both representations consume the *same* 33 points, so the constructive and SDF
paths stay consistent with each other — this is not a dual-rep divergence. Both
are simply that far from the true midline.

Note what is **not** wrong here: the `bisection_iters=0` shortcut used to build
these points costs nothing (0.00 µm — see "Cleared by measurement"). The sampling
density is the whole of the error.

### Why it matters less than it looks, and still matters

The cap **mid-surfaces** exported for shell meshing come from
`_airfoil_offset_arc`, i.e. the offset ring, not from this polyline. So SG-010
does not touch the strip geometry that SG-001 to SG-005 were about. It affects
the cap **solids** — areas, and therefore section properties and mass
distribution.

### Fix

Under the SG-008 polygon contract, the solids clip should not be coarser than
the outline it inherits. Default `n_samples` is now derived from band width so
chord-approximation deviation stays within
`CAMBER_CLIP_DEV_BUDGET_M = 1e-6` (1 µm) — well below the ~62 µm LE polygon
floor. Callers may still pass an explicit `n_samples`. Feasibility clamping
behaviour is unchanged; only the clip polyline densifies.

## SG-011

**Family-D continuous-box spar caps leave a bay gap against flapwise webs under twist**

**Status:** Fixed
**Found:** 2026-09-06
**Fixed:** 2026-09-06
**Area:** `recipes/sdf_ir/spar_cap.py`, `constructive/emitters.py::emit_spar_cap_polygon`,
`sections/build_spar_caps.py` (`continuous_box`)

### Symptom

On `2D-F` with non-zero twist (e.g. section-selection e2e st018, twist ~11°),
constructive wall solids and shell midlines show spar caps floating inside the
box bay — several millimetres clear of one or both shear webs. At twist = 0°
the same gauges look flush.

### Evidence

NACA 2412, chord 1 m, skin 3 mm / spar 11.625 mm / web 4 mm, layout `2D-F`:

| Twist | Solid cap–web gap | Midline end → web mid-surface |
|---|---|---|
| 0° | 0.000 mm | 0.000 mm |
| ~11° (pre-fix) | ~5 mm (one side per cap) | ~8–9 mm |

CapWall already stored `flapwise_nx/ny/c_lo/c_hi`, but recipe emission ignored
them unless `_strip_clip` was attached. CapWall never passes `_strip_clip`, so
`SparCap` IR and `emit_spar_cap_polygon` fell back to `ChordBand(x_start, x_end)`.
At θ = 0, `n = (-1, 0)` so ChordBand coincides with web midplanes; under twist
the webs sit on slanted `c_station` planes and the chordwise band drifts.

Secondary: `continuous_box` recomputed strip bounds from camber at `x_cap`
instead of using assemble SoT `c_station[0]/[-1]` (what sandwich cores already use).

### Root cause

Flapwise strip parameters on CapWall were documentation for the Shapely path but
were not applied unless an explicit `_strip_clip` callable was present. Recipe-
driven solids, properties, and midlines therefore stayed chordwise under twist.

### Fix

1. `SparCap.__init__` derives `_strip_clip` from `flapwise_*` when missing and
   builds the IR clip from that strip (midlines inherit the same `_strip_clip`).
2. `emit_spar_cap_polygon` uses `FlapwiseStrip` whenever `flapwise_*` is set —
   same contract as sandwich cores; no longer requires `_strip_clip`.
3. `build_spar_caps` `continuous_box` sets `c_lo/c_hi` from `c_station[0]/[-1]`.

`export.py` `cap.difference(webs)` unchanged — it only trims overlap after the
strip is authored.

### Verification

- `tests/test_family_d_cap_web_flush_under_twist.py` — solids and midlines flush
  at 0° and 15° (`gap ≤ 1e-6` m).
- Section-selection e2e `2D-F` / `new_3.2m_rotor` / 8 stations with plots: cross-
  section PNGs show caps flush to both webs.

## SG-012

**Constructive CapWall solids split suction/pressure spar caps on `y=0`, not the cavity midline**

**Status:** Fixed
**Found:** 2026-09-06
**Fixed:** 2026-09-06
**Area:** `recipes/walls.py::CapWall`, `sections/build_spar_caps.py`,
section-selection viz `_rebuild_sized_section_for_viz`

### Symptom

On outboard `2D-F` stations (e.g. section-selection e2e st045 / st054 / tip
st063), suction and pressure spar caps appear to merge across the bay. The
shared boundary sits on a near-horizontal chord line; suction area is ~2.7×
pressure for equal recipe heights. Tip stations with `spar_cap_thickness`
run out (`None`) still plot thick caps.

### Evidence

NACA 2412, chord 0.123 m, skin 3 mm / spar 6.125 mm / web 4 mm, layout `2D-F`,
twist 0°:

| Path | Area ratio su/pr | Shared boundary |
|---|---|---|
| Pre-fix constructive solids | ≈ 2.70 | `y = 0` |
| Post-fix constructive solids | ≈ 1.00 | cavity midline (`camber_xy`) |

`cap_feasibility_clamp` already returned `info["camber_xy"]` (SG-010) and
`_project_cap` attached it to SDF `SparCap` components. CapWall recipes —
the SoT for `emit_wall_solid_polygons` / cross-section plots — never stored
those samples, so `emit_spar_cap_polygon` fell through to
`half_plane_clip(..., y=0)`.

Separately, viz `_rebuild_sized_section_for_viz` replaced a missing/zero spar
DV with `skin_thickness` while keeping the root `2D-F` layout, redrawing caps
after BL-006 run-out.

### Root cause

Camber split lived only on the SDF projection, not on the constructive CapWall
recipe consumed by plot emission. Tip viz invented a spar thickness instead of
stepping to the cap-less sibling layout.

### Fix

1. Store `camber_xy` on `CapWall`; `_cap_wall` freezes feasibility samples into
   the recipe for all placement kinds (`continuous_box`, `fixed_band`,
   `per_web_pads`).
2. CapWall `_component` / `to_boundary_polygon` default to the stored camber so
   constructive emit uses `_camber_clip`.
3. Viz: when the layout has caps but `spar_cap_thickness <= 0`, switch to
   `layout_without_spar_caps` instead of inventing a thickness.

### Verification

- `tests/test_cap_camber_split_constructive.py` — constructive area ratio in
  `[0.75, 1.35]` at 0° and 2.7°; untwisted interface tracks cavity midline.
- `tests/unit/test_retain_mesh.py::test_viz_rebuild_drops_caps_when_spar_dv_missing`
  — no `spar_cap_*` labels when the DV is omitted.

## SG-013

**Family-C (and B) pads recenter onto `argmax(depth)` instead of the web / band center**

**Status:** Fixed
**Found:** 2026-09-07
**Fixed:** 2026-09-07
**Area:** `sections/cap_feasibility.py`, `sections/build_spar_caps.py`

### Symptom

On cambered sections (e.g. section-selection e2e `1C-F` at st038 / z/L≈0.60),
family-C spar caps are visibly skewed about the shear web — more chordwise
extent LE-ward than TE-ward — despite being authored as `web_x ± dhw`.

### Evidence

`cap_feasibility_clamp(..., aft_half_depth_rule=False)` rebuilt the feasible
interval by growing from `nanargmax(depth)`. On NACA 2412 the deepest cavity
sample sits LE-ward of a mid/aft web, so the contiguous feasible island slides
off the web. Taxonomy already requires **Y = C** caps “centred on each web
chord station” (`reference/system_type_xyz_taxonomy.md`).

### Root cause

The `False` branch of `_contiguous_extent` was written as a deepest-point
island grower (legacy unused path) and then wired to B/C under BL-008 without
preserving the authored band midpoint. Family D’s `aft_half_depth_rule=True`
path never touches extent and was not at fault.

### Fix

1. `cap_feasibility_clamp` gains `extent_anchor` (`"band_center"` default,
   `"deepest"` legacy) and optional `center_x`.
2. `"band_center"` finds the largest half-width such that
   `[center − half_w, center + half_w]` is feasible everywhere — equal trim,
   midpoint fixed (SG-013). Collapsed bands scale heights to the worst sample
   in the authored band; a web that cannot carry any centred pad runs out
   (`pad_run_out`, height 0) instead of emitting a zero-width strip.
3. `build_spar_caps` passes `extent_anchor="band_center"` for B/C with
   `center_x` = pitch/max-thickness or each `xs[i]`. Multi-web C re-clamps
   extents after the shared min height, and skips run-out pads when taking that
   min / emitting CapWalls. Family D unchanged.

### Verification

- `tests/test_cap_feasibility_clamp.py` — 2412 band_center midpoint stays on
  `center_x`; deepest drifts LE-ward; D aft-rule bit-identical; B/C tip meshability.
- `tests/test_structural_family.py::test_structural_C_cap_recipes_centered_on_webs`
  — CapWall midpoints match web stations on 2412.

## SG-014

**Cell wall loops pair the wrong web ends on B-frame sections and carry no wall membership; every consumer re-guesses it**

**Status:** Fixed
**Found:** 2026-09-25 (4D-F `new_3.2m_rotor` re-qualification; repo [BL-040](../../../../../docs/bug_log.md#bl-040))
**Fixed:** 2026-09-25
**Area:** `cell_wall_loops/`, `section_contract/cell_wall_loop.py`; consumers in
`section_properties.reference_axes.wall_network` and `section_loads.thinwall`

### Symptom

4D-F station 21 (layout 4A-F, four webs, 0.2 m chord) never reached a
closed-section `GJ`: every trial logged `closed-section GJ unavailable (segment 0
is claimed by 3 cells; a wall can bound at most two); polar kept`, so the
torsion gate could not pass however the station was sized.

### Evidence

The failing station's midlines and loops were captured from the campaign
workers. The emitted cell areas summed to **1.71×** the skin-ring area
(0.01133 vs 0.00662 m²): bay 1 spanned `Y_B` −0.069…+0.023 (it had swallowed the
LE cell) and bays 2 and 3 both ran to the TE, so the TE skin sat in cells 2, 3
and 4. The same strips transposed to chord-first gave five correct LE→TE cells
tiling the ring exactly.

A worker sweep (1C, 2D, 3D, 4D in -F and -CN, 2B-F, 3C-F; twist 0–1.2 rad;
chord 1.0 and 0.2 m) found **36 of 80** configurations emitting cells that do not
tile the section — 14 raising "claimed by 3 cells" and **22 solving silently on
overlapping cells**, including 2D-F and 3D-F at moderate twist (area +18 % to
+100 %). Every multi-web layout was affected somewhere; one web never was.

Production audit (a hook on every closed-section solve of the 16-station
`new_3.2m_rotor` campaign at `ad560972`):

| Layout | Solves | On wrong cells | Raised | Σ cell area / ring |
|---|---|---|---|---|
| 2D-F | 113 | 32 | 0 | ≤ 1.28 |
| 3D-F | 121 | 48 | 0 | ≤ 1.67 |

Bredt `GJ` scales with the enclosed area squared, so on those solves the
torsional stiffness was overstated by up to ~2.8× without warning.

### Root cause

Two defects, one of them structural.

1. **The contract dropped what the producer knew.** `build_cell_wall_loops`
   assembled each cell from named pieces — the skin arcs `split_skin_at_junctions`
   cut at the web feet, and specific webs — then emitted only a polyline and an
   area. Every consumer re-derived which walls bound which cell from
   coordinates, each with its own tolerance:
   `section_properties.network_from_midlines_props` took the nearest segment
   midpoint within an absolute 5 cm (over *all* strips, spar caps included), and
   `section_loads` took segment midpoints within 0.25 × median segment length of
   the loop polyline. Both are guesses that fail as walls get close — the class
   of heuristic SG-001…SG-004 removed from the junction pass.
2. **The builder reasoned chord-first on canonical B.** It sorted webs LE→TE on
   `x`, chose each web's top on `y`, and chose the LE/TE skin arcs on `mean_x`.
   In B, `X_B` is flapwise and `Y_B` chordwise. A flapwise web has both ends at
   the same `Y_B`, so "top" was arbitrary per web; a bay joined web *a*'s top to
   web *b*'s top, and when those lay on opposite surfaces the shorter skin arc
   ran round the nose or the tail. When the midlines moved to B,
   `section_loads` gained an explicit `_b_to_chordwise_first` adapter; this
   builder did not, and its unit tests were written chord-first, so nothing
   caught it.

Nothing checked the result: no test that each web bounds two cells, each skin
arc one, and that the cells tile the ring. `section_properties`' log records
that "a segment appearing in two cells **is** a shared wall, so nothing is
re-derived from coordinates" — true of the solver, but the network it was
handed had been re-derived from coordinates.

### Fix

- **Contract.** `CellWallLoop.wall_labels`: the strip labels bounding the cell,
  in loop order. It **is** the membership, and consumers read it.
- **Builder** (`cell_wall_loop_builders.py`, rewritten). The section is a planar
  graph — nodes are the skin-arc junctions, edges the labelled strips — and the
  cells are its bounded faces, found by walking the ring and turning onto a web
  wherever one starts. No coordinate axis is read: web order comes from the
  walk, LE→TE numbering from the leading-edge arcs. A web end must **coincide**
  with a junction to 1 nm (measured ≤ 5×10⁻¹³ m; closest junction pair 4 mm), so
  a web stopping short is refused rather than extrapolated (standing
  invariant). It also refuses two webs on one junction, crossing webs, a skin
  that does not chain, and any result whose cell areas do not tile the ring to
  1×10⁻⁹. `cell_wall_web_attach.py` and the ray/arc helpers are deleted.
- **`section_properties`.** `network_from_midlines_props` builds cells from
  `wall_labels` and refuses a loop without them; the nearest-midpoint match is
  deleted. It now adapts `CellWallLoop` objects itself — before, a DTO passed
  straight in failed on `.get`, and `shear_centre_from_midlines` silently
  returned the elastic centroid, so the BIS `shear_centre_xy` has never been a
  multi-cell shear centre. That value now changes.
- **`section_loads`** — Channel C builds its own cells; logged and fixed in its
  own log.

### Verification

- `tests/test_cell_wall_loops.py` (rewritten in B): 1–5 flapwise webs give
  contiguous LE→TE cells with exact areas; web order and direction are
  irrelevant; membership is invariant under rotation; caps bound nothing;
  short, crossing and coincident webs are refused; the worker sweep (64
  layout × twist × chord cases) tiles every section.
- `section_properties/tests/test_closed_section_torsion.py`: unlabelled loops
  are refused; four flapwise webs on a 0.2 m chord solve with five cells.
- Station 21's captured inputs: five cells, areas 1.53 / 2.39 / 1.90 / 0.77 /
  0.027 ×10⁻³ m², summing to the ring.

## SG-015

**Closed-ring suction/pressure split and the side fallback read chord-first coordinates on B-frame strips**

**Status:** Fixed
**Found:** 2026-09-25 (while fixing SG-014)
**Fixed:** 2026-09-25
**Area:** `midlines/skin_split.py`

### Symptom

`split_skin_at_junctions`, given a single closed skin ring in B, split it at the
wrong points: rotating a box section by 0.3 rad moved the "leading-edge" arcs to
the trailing-edge web.

### Root cause

`_ensure_ss_ps_skins` passed B-frame points to `split_closed_ring_ss_ps`, which
takes the LE as `argmin x` and suction as the higher mean `y` — chord-first
conventions. `_normalise_side`'s geometric fallback likewise classified by
`mean(y) >= 0`, i.e. by LE/TE half rather than by surface. Production is not
affected: multi-web sections emit labelled suction/pressure halves, and
`OuterSkin` calls the splitter in its chord-first authoring frame. The path is
reached by closed or unlabelled skins (legacy inputs, tests).

### Fix

Transpose to chord-first around the `split_closed_ring_ss_ps` call and back; the
side fallback reads `X_B` (suction on `+X_B`).

### Verification

`test_cell_wall_loops.py::test_a_closed_b_frame_skin_is_split_into_suction_and_pressure_by_x_b`.

## SG-016

**50 package tests fail at HEAD (`SparCap._camber_surface`, g-row hardening)**

**Status:** Fixed
**Found:** 2026-09-25
**Fixed:** 2026-09-25
**Area:** `sections/build_spar_caps.py::_wrap_expr_cap`, `sections/cap_ply_snap.py`, package `__all__`

### Symptom
At `ad560972`, `test_spar_cap_midline_offset` failed 48 cases with
`AttributeError: 'SparCap' object has no attribute '_camber_surface'`, and
`test_g_row_hardening` failed on a `laminate` comment and barrel extras.

### Root cause
`_wrap_expr_cap` rebuilds Family-D caps via `SparCap.__new__` and copied
`_camber_xy` but not `_camber_surface`, which `midline_polyline` requires.
Separately, a production comment used the word `laminate`, and
`apply_spar_cap_runout` / `max_spar_from_cavity` were on `__all__` without being
on the G-row allowlist (matrix silent → deep-only).

### Fix
Copy `_camber_surface` in `_wrap_expr_cap`; reword the ply-snap comment; drop the
run-out helpers from the public barrel (lazy deep import remains).

### Verification
`test_spar_cap_midline_offset` + `test_g_row_hardening` — 83 passed. Package
graphify refreshed under `section_geometry/graphify-out/`.

## SG-017

**Caps are built from the skin but their bond to it is never emitted**

**Status:** Fixed
**Found:** 2026-09-25 (thin_wall_section TW-001; section_loads SL-012)
**Fixed:** 2026-09-25
**Area:** `topology/cap_skin_bonds.py` (new), `section_contract.ShellMidlineStrip.bonds` / `StripBond`

### Symptom

Every consumer that needs to know which skin a spar cap is attached to guessed
it: Channel C decided which skin segments a cap "stacks" on by a chordwise
bounding box and a flap tolerance of 15 % of the section's flap span
(section_loads SL-012). The exact thin-wall kernel cannot solve shear with a cap
at all without the attachment (TW-001): the cap's axial rate has to enter the
skin somewhere.

### Root cause

`cap_skin_ray_resample` builds each cap's interior vertices from inward rays cast
from skin vertices, so the correspondence exists while the cap is built — but the
hits are merged and filtered, and only the cap polyline is emitted. As with
SG-014, the producer knew the topology and the contract dropped it.

### Fix

`bond_caps_to_skin`, the last topology pass (nothing after it may move cap or
skin vertices):

- the skin arcs of the cap's declared surface are chained head to tail;
- each cap vertex is paired with its perpendicular foot on that chain, and the
  feet must advance monotonically;
- the cap-to-skin map is piecewise linear between those pairs; every skin vertex
  inside the footprint takes its cap partner **from the same map**. An earlier
  draft projected skin vertices back onto the cap independently; near a web
  junction (a 0.28 mm skin segment beside a 1.2 mm cap segment) the two
  projections are not mutual inverses and interleaved (3D-F and 3C-F at 0.2–1.0 m),
  which the one-map construction rules out;
- a foot within a quarter of the adjacent segment of an existing vertex is paired
  with that vertex — this chooses which points correspond, never where anything
  is; every other foot is inserted **on** its polyline, so no geometry moves;
- one `StripBond` per host skin arc, as coordinate pairs (they survive strip
  reversal in the mesh handoff).

### Verification

60 worker sections (1C-F/CN, 2D-F/CN, 3D-F, 3C-F, 4D-F/CN, 2B-F, 0B × twist 0 /
0.6 / 1.2 rad × chord 1.0 / 0.2 m): every cap bonded one-to-one, and the thin-wall
kernel integrates every unit load back to ≤ 6.4×10⁻¹¹ with the caps at their own
offset.

## SG-018

**Residual nearest / argmin sites vs standing invariant**

**Status:** Logged — 2026-09-26 (no relocation defect)
**Area:** package-wide production nearest / `argmin` (Exactness item 6 follow-on)

### Symptom
Programme item 6 asked for a review of remaining "nearest" sites against
"topology splits and indexes; it never relocates."

### Classification (production, excluding tests/demos unless they drive emit)

| Site | Verdict |
|---|---|
| `split_skin_at_junctions` + `point_to_open_polyline_nearest` | OK — indexes arc-s; does not write the web foot |
| `cell_wall_loop_builders._node_of` | OK — coincidence check; raises above `WEB_FOOT_COINCIDENCE_TOL_M` |
| `cap_skin_bonds._foot` | OK — SG-017 bond construction; inserts on the polyline, does not move existing vertices |
| `recipes/sdf_ir/outer_skin.py`, `sdf/airfoil.py`, `sections/clips.py` LE `argmin` | OK — identify LE vertex on the discrete airfoil polygon |
| `sections/cap_feasibility.py` `argmin` on `xs_q` | OK — seed a design sample, not topology membership |
| `nearest_on_closed_polyline*` | definition only; no production caller |
| `sdf/primitives.py` nearest box | SDF distance, not strip identity |

No site rewrites geometry or invents wall/strip membership when the producer
already knows it. Survivors stay as documented indexing / coincidence / emit
construction.

### Verification
Re-grep of `point_to_open_polyline_nearest`, `nearest_on_closed_polyline`, and
production `np.argmin` under `section_geometry/` (2026-09-26).

