# R10 · Feasibility, drops and failure modes

## 1. Problem

A station is described by an aerofoil, a layout ([R5](R5_section_layout.md)) and wall
thicknesses. Not every description can be built: a wall may be thinner than can be
made, a web may stand where the cavity has closed, and a pair of caps may be deeper
than the cavity between them. This document derives the rules the kernel applies
before and during assembly to decide what the station actually builds: the
thickness floors, the probe that finds a web's intercepts with the inner skin, the
three reasons a web is dropped, the step-down that follows a drop, the clamp that
fits the caps into the cavity, the cap run-out, and the conditions under which the
kernel raises `InfeasibleSectionGeometry` instead of building.

The working guide for callers is [feasibility and drops](../section/guides/feasibility_and_drops.md).

## 2. Definitions and assumptions

1. **Section frame.** Aerofoil coordinates $(x, y)$, leading edge at $x = 0$, chord
   $c$ along $x$ ([R1](R1_frames.md)). All rules below act before the twist turn.
2. **Inner skin.** The zero level of $\phi_i = \phi + t_s$, the aerofoil's signed
   distance offset inwards by the skin thickness $t_s$ ((R4.15) with $\delta = -t_s$). Both skins use the
   larger of the two thicknesses.
3. **Cavity depth.** At a chord station $x$ where the vertical line meets the inner
   skin, $y_t(x)$ and $y_b(x)$ are its highest and lowest crossings, and
   $d(x) = y_t(x) - y_b(x)$. Where the line misses the cavity, $d = 0$ (the kernel
   stores $-\infty$).
4. **Bounding box.** $[x_{\min}, x_{\max}] \times [y_{\min}, y_{\max}]$ of the
   aerofoil polygon; $T = y_{\max} - y_{\min}$ is its thickness and
   $W = x_{\max} - x_{\min}$ its chord extent.
5. **Caps.** Suction and pressure cap heights $h_u, h_l$, authored over a chordwise
   band $[x_\ell, x_r]$ with centre $x_c$.
6. **Assumption.** The section is a single simply connected aerofoil whose inner skin
   meets every vertical line at most twice inside the band of interest (true for the
   NACA and blended families of [R2](R2_aerofoils.md) away from a cusped trailing
   edge). Where it is not, $y_t$ and $y_b$ are still the outermost crossings.

## 3. Derivation

### 3.1 Thickness floors

Walls are not built thinner than a manufacturing minimum. Every thickness enters the
section builder through

\begin{equation*}
\begin{gathered}
t_s \leftarrow \max(t_s,\ t_s^{\min}), \quad
t_w \leftarrow \max(t_w,\ t_w^{\min}), \quad
h \leftarrow \max(h,\ h^{\min}),
\\
t_s^{\min} = 2.5\ \mathrm{mm}, \quad t_w^{\min} = 4\ \mathrm{mm}, \quad h^{\min} = 5\ \mathrm{mm}.
\end{gathered}
\tag{R10.1}
\end{equation*}

The floors are applied in `plan_section`, so they bind every section; a cap height
of zero is raised to $h^{\min}$, not removed (a capless station uses a family A
layout). The loft avoids them for small chords by building the cut at the scale $k$
of (R8.3). The web anchors of Section 3.2 are found by `resolve_station_layout`
with the skin thickness as given, before (R10.1) is applied: a skin under
$t_s^{\min}$ is probed at one depth and built at another (SG-040, Section 6).

### 3.2 The intercept probe

The crossings $y_t, y_b$ at a web station $x$ are found by the probe of (R5.10): it
samples $\phi_i(x, y)$ at $n = 240$ equally spaced $y$ in a window padded by $p$
above and below the bounding box,

\begin{equation*}
y_k = y_{\min} - p + k\,\Delta y, \qquad
\Delta y = \frac{T + 2p}{n - 1}, \qquad k = 0, \dots, n - 1,
\tag{R10.2}
\end{equation*}

and keeps the indices where the sign of $\phi_i$ changes. The window ends lie outside
the aerofoil, where $\phi_i > 0$, so a column has at least two sign changes exactly
when some sample lies inside the cavity, $\phi_i(x, y_k) < 0$. Inside the cavity
$y_b < y_k < y_t$, and the samples are $\Delta y$ apart, so

\begin{equation*}
\begin{gathered}
d(x) > \Delta y \ \Rightarrow\ \text{found},
\\
d(x) \le \Delta y\ \Rightarrow\ \text{found only if some } y_k \in (y_b, y_t).
\end{gathered}
\tag{R10.3}
\end{equation*}

Below $\Delta y$ the outcome depends on where the grid falls; a cavity that exists is
then reported as `no_inner_skin_intercept`. The pad decides $\Delta y$. The probe's
default window pads by $p = \max(0.05\,T, 0.01\ \mathrm{m})$; `resolve_web_anchors`
builds its own window with $p = \max(0.05\,W, 0.01\ \mathrm{m})$, the chord extent in
place of the thickness. On a 1 m NACA 0006 that is $\Delta y = 0.67$ mm against
0.33 mm, and cavities of 0.4 and 0.6 mm at $0.8c$ are dropped although they exist
(SG-038).

### 3.3 Crossing accuracy

Between the two bracketing samples $(y_k, \phi_k)$ and $(y_{k+1}, \phi_{k+1})$ the
probe bisects $m$ times and then takes the secant root

\begin{equation*}
y^\ast = y_a - \phi_a\,\frac{y_b' - y_a}{\phi_b' - \phi_a},
\qquad
|y^\ast - y_{\mathrm{true}}| \le 2^{-m}\,\Delta y,
\tag{R10.4}
\end{equation*}

where $[y_a, y_b']$ is the bracket left by the bisection. The bound holds because the
root stays in the bracket. The web anchors use $m = 10$; the cap clamp
(Section 3.6) uses $m = 0$. On a polygon aerofoil $\phi_i$ is the distance to one
edge throughout a short bracket, which is linear in $y$, so the secant root is exact
to round-off unless the bracket straddles a vertex's medial axis.

### 3.4 Drop reasons

`resolve_web_anchors` visits the requested web stations in increasing $x$. A web at
$x_j$ is dropped, with the first reason that applies, when

\begin{equation*}
\begin{gathered}
\text{(i) no intercept:}\quad \text{fewer than two sign changes by (R10.3)},
\\
\text{(ii) cavity too shallow:}\quad d(x_j) < d_{\min},
\\
\text{(iii) merged:}\quad x_j - x_{\mathrm{prev}} < 0.02\,c,
\end{gathered}
\tag{R10.5}
\end{equation*}

where $x_{\mathrm{prev}}$ is the last web kept (not the last requested) and $d_{\min}$
is `min_cavity_depth_m`, zero unless the caller sets it. Reason (ii) is a structural
test that (i) cannot make: as the chord shrinks along the span, a web anchored in
metres moves aft in $x/c$ into a cavity that still has two crossings but is too
shallow to carry a cap. A dropped web is not moved; the legacy walk towards the
leading edge (`allow_le_ward_migration`) is off by default.

### 3.5 Step-down

With $N$ webs requested and $N'$ kept, the effective layout is the nominal one with
the dropped webs removed:

\begin{equation*}
\text{layout}' =
\begin{cases}
\text{raise} & N' = 0 < N,
\\
\text{family C, per-web pads} & N' = 1,\ \text{nominal family D},
\\
\text{nominal family, }N'\text{ webs} & \text{otherwise},
\end{cases}
\tag{R10.6}
\end{equation*}

unless the family D layout has a published chord window, and the catalogue key is
re-attached for $N'$ webs. This is the ladder of (R5.3)
applied to the survivors: D falls to C at one web because a continuous box needs two,
and the ladder floors at one web because a section with no web has no shear path.
`StationLayoutResolution.stepped_down` is true when the web count, family or cap
placement changed. Assembly never repeats the drop: given web positions without
frozen anchors, it raises if any web would be dropped, so the step-down is always
the caller's visible decision.

### 3.6 The cap allowance

Upper and lower caps must leave a void between them. With a void of at least
$v = \max(1\ \mathrm{mm}, 0.15\,d)$, the room for the two caps together is

\begin{equation*}
\begin{gathered}
a(x) = d(x) - \max\big(v_0,\ \alpha\, d(x)\big),
\qquad v_0 = 1\ \mathrm{mm},\quad \alpha = 0.15,
\\
\text{feasible at } x \iff h_u + h_l \le a(x).
\end{gathered}
\tag{R10.7}
\end{equation*}

The two terms are equal at $d = v_0/\alpha = 6.67$ mm; deeper cavities keep 15 % of
their depth as void. The depth is sampled at $n$ equally spaced stations of the band,
with $n$ chosen so that the polyline through the cavity mid-points, which the caps
are later split along, deviates from the true mid-line by at most $\varepsilon =
1\ \mu$m. The deviation was measured as $18.56\ \mu$m at $n = 33$ on a 0.5 m band
(NACA 2412, 1 m chord, 6 mm skin) and taken to scale as $C B / n^2$ for a band of
width $B$, so

\begin{equation*}
n = \max\Big(33,\ \Big\lceil \sqrt{C B/\varepsilon}\, \Big\rceil\Big),
\qquad
C = \frac{18.56\times10^{-6} \cdot 33^2}{0.5} = 0.0404.
\tag{R10.8}
\end{equation*}

This is a fitted law, not a derived one: a chord-to-arc error scales as
$\kappa B^2/n^2$, and the fit absorbs $\kappa B$ for sections like the one measured.

### 3.7 The pad extent

For families B and C (`extent_anchor="band_center"`) the cap keeps its centre $x_c$
and is trimmed equally on both sides. The kept half-width is the largest $w$ for which
every sample in the interval is feasible:

\begin{equation*}
\begin{gathered}
w^\ast = \max\big\{\, w \le w_{\max} :\ h_u + h_l \le a(x_q)\ \ \forall\, x_q \in [x_c - w, x_c + w] \,\big\},
\\
w_{\max} = \min\big(x_c - x_\ell,\ x_r - x_c,\ \tfrac12 (x_r - x_\ell)\big),
\end{gathered}
\tag{R10.9}
\end{equation*}

found by 48 bisection steps on $w$. The extent never grows past the authored band and
never slides its centre, so a pad stays centred on its web.

### 3.8 Height scaling and run-out

The heights are scaled down in two cases. If the caps fit nowhere,
$h_u + h_l > \max_q a(x_q)$, both are scaled by one factor so that their sum fits the
central half of the band, or half the best allowance if that is larger:

\begin{equation*}
s_1 = \min\Big(1,\ \frac{\max\big(\min_{q \in Q_{1/2}} a_q,\ \tfrac12 \max_q a_q\big)}{h_u + h_l}\Big),
\tag{R10.10}
\end{equation*}

where $Q_{1/2}$ is the middle half of the sample indices. If the extent (R10.9) is
then narrower than half the authored band, the caps are thinned rather than left as a
needle: for families B and C both heights are scaled so their sum fits the whole band,

\begin{equation*}
s_2 = \frac{\min_q a_q}{h_u + h_l} \quad \text{when } \min_q a_q > 0,
\tag{R10.11}
\end{equation*}

and the extent is recomputed. When some sample of the band has $a_q \le 0$, for
example where the band reaches past the end of the cavity, (R10.11) does not apply.
If the extent is then empty, or the centre sample is infeasible, the pad runs out:

\begin{equation*}
\begin{gathered}
w^\ast < 10^{-12}\ \text{m} \ \ \text{or}\ \ h_u + h_l > a(x_c)
\\
\Rightarrow\ h_u = h_l = 0,\quad x_\ell' = x_r' = x_c .
\end{gathered}
\tag{R10.12}
\end{equation*}

For per-web pads the heights are then made common, the smallest over the pads that
did not run out, and the pads are clamped again with that height. Every height is
finally snapped down to a whole number of plies when a ply pitch is given; a height
under one ply becomes zero.

### 3.9 Family D

A continuous box spans the webs, so its extent is not trimmed. Each cap is instead
limited to half the depth at the aft end of the band, the aft web, where an aerofoil
aft of maximum thickness is shallowest:

\begin{equation*}
h_u \leftarrow \min\big(h_u,\ \tfrac12\, d(x_r)\big),
\qquad
h_l \leftarrow \min\big(h_l,\ \tfrac12\, d(x_r)\big).
\tag{R10.13}
\end{equation*}

The caps can then at most meet at the cavity mid-line at the aft web; (R10.13) keeps
no void, unlike (R10.7). If the box begins forward of maximum thickness the forward
web can be the shallower one; the rule does not look there.

### 3.10 The cap run-out gate

`apply_spar_cap_runout` replaces a capped layout by the same layout without caps
(R5.3) when a minimum-gauge cap cannot be built. Its first test bounds each cap by

\begin{equation*}
h \le \tfrac12\big(h_{\mathrm{loc}} - t_s - g\big),
\qquad g = 2\ \mathrm{mm},
\tag{R10.14}
\end{equation*}

which is the cavity room per cap when $h_{\mathrm{loc}}$ is the depth between the
two skin mid-surfaces, $h_{\mathrm{loc}} = d + t_s$; it then becomes $(d - g)/2$. Its
second test runs the clamp of Sections 3.6 to 3.9 at minimum gauge. The two void rules
differ: the gate keeps $g = 2$ mm and the clamp keeps $\max(1\ \mathrm{mm}, 0.15\,d)$,
equal at $d = g/\alpha = 13.3$ mm, so the gate is stricter in shallower cavities and
more lenient in deeper ones; and $h_{\mathrm{loc}}$ is not defined in the code
(SG-039). The gate has no caller in the library; consumers call it.

### 3.11 Failure modes that raise

Everything above degrades the station. The cases below cannot be degraded and raise
`InfeasibleSectionGeometry`, a `ValueError`:

| Condition | Where | Note |
|---|---|---|
| No requested web anchors, $N' = 0 < N$ | `resolve_station_layout`, `catalogue.system_layout` | (R10.6) |
| Assembly asked to build webs that would drop | `section.assembly.plan._anchor_webs` | Resolve first |
| Frozen anchor outboard of the skin mid-surface | `ShearWebComponent._extend_to_skin_mid_surface` | Anchors resolved for a thinner skin (SG-031, SG-040) |
| Web axis does not reach the skin mid-surface | `ShearWebComponent` mid-line | |
| Inward offset empty, not a polygon, or under three vertices | `section.boundary.polygons` | Skin deeper than the section |
| Flapwise web cannot re-anchor under twist | `web_layout.flapwise_web_sframe_inputs` | |
| Flapwise web within 1° of the chord | `plan._check_flapwise_webs_cross_chord` | SG-033 |
| Cap placement missing a side | `section.assembly.members` | |
| Continuous box with fewer than two webs | `plan_section`, `build_spar_caps` | Unreachable after (R10.6) |
| Offset branch not ascending in $x$ | `section.assembly.clips` | |
| Two-foot web feet closer than $4\,t_s$ | `blade.loft.loft_cut` | (R8.5) |

## 4. Algorithm

At one station (`SectionGeometry.resolve`, then `execute`):

1. Probe each requested web station by (R10.2)–(R10.4) and drop webs by (R10.5).
2. Step the layout down by (R10.6), or raise if no web is left.
3. Apply the floors (R10.1) and plan the walls; refuse a flapwise web along the chord.
4. Build the webs on the frozen anchors of step 1.
5. Clamp the caps: (R10.7)–(R10.12) for families B and C, (R10.13) for family D;
   snap to the ply pitch.
6. Optionally, before step 3, a consumer applies the run-out gate (R10.14).

## 5. Implementation

| Result | Function (`implicit_geometry_kernel.`) |
|---|---|
| (R10.1) | `section.layout.thickness_limits.clamp_skin_thickness_m`, `clamp_web_thickness_m`, `clamp_spar_thickness_m`; applied in `section.assembly.plan.plan_section` |
| (R10.2)–(R10.4) | `section.assembly.web_layout.inner_y_at_x`, `inner_y_at_xs`, `_default_y_search`, `_interp_crossing_1d` |
| (R10.5) | `section.assembly.web_layout.resolve_web_anchors`; `section.spec.results.DroppedWeb` |
| (R10.6) | `section.spec.layout_resolution.resolve_station_layout`, `_step_down_section_layout`, `_reattach_catalogue_preset`, `StationLayoutResolution` |
| (R10.7)–(R10.12) | `section.assembly.cap_feasibility.cap_feasibility_clamp`, `default_camber_clip_n_samples`, `_symmetric_extent_about_center`, `_scale_heights_to_central_half`; `section.assembly.build_spar_caps_pads.fill_per_web_pads`, `snap_cap_heights` |
| (R10.13) | `cap_feasibility_clamp(..., aft_half_depth_rule=True)` from `section.assembly.build_spar_caps._caps_continuous_box` |
| (R10.14) | `section.layout.spar_cap_runout.max_spar_from_cavity`, `apply_spar_cap_runout`, `_placement_buildable_at_gauge` |

## 6. Limits

1. **Shallow cavities (SG-038).** `resolve_web_anchors` pads its probe window by 5 %
   of the chord, not of the thickness, so a cavity shallower than (R10.2)'s
   $\Delta y$ may be dropped although it exists; on thin tip sections $\Delta y$ is
   about 0.7 mm per metre of chord.
2. **Floors after the probe (SG-040).** Anchors are resolved at the given skin and
   built at the floored skin of (R10.1). Below $t_s^{\min}/2 = 1.25$ mm the frozen
   anchor lies outboard of the built skin's mid-surface and the build raises; between
   1.25 and 2.5 mm the drop decisions are made for the wrong skin.
3. **Zero caps.** A cap height of zero is floored to 5 mm, not removed (R10.1).
4. **Two void rules (SG-039).** The run-out gate (R10.14) and the clamp (R10.7) keep
   different voids, and the gate's $h_{\mathrm{loc}}$ is undocumented.
5. **Run-out gate geometry (SG-041).** The gate tests a fixed band centred at $0.5c$
   and per-web pads of half-width $0.04c$; the build centres a fixed band at the
   pitch axis or maximum thickness (or its published window) and takes the pad
   half-width from the section options.
6. **Fitted sample law.** (R10.8) is fitted on one section; for a much more curved
   cavity mid-line the 1 µm budget is not guaranteed.
7. **Family D void.** (R10.13) keeps no void at the aft web and reads only the aft
   end of the band.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R10.1) | Asked 2, 1, 0 mm; built 2.5, 4, 5 mm. SG-040: a 1.0 and 1.2 mm skin raise with the anchor $t_s^{\min}/2 - t_s$ outboard; 1.3 mm builds | [section example 06](../../examples/section/06_feasibility_checks/README.md), check 1 |
| (R10.2), (R10.3) | Six cavities from 0.2 to 2 mm at $0.8c$ of NACA 0006: every outcome equals "a sample lies in the cavity"; SG-038 printed | example 06, check 2 |
| (R10.4) | 18 columns: error $4\times10^{-18}$ m against the bound $9\times10^{-7}$ m | example 06, check 3 |
| (R10.5) | One web of each drop reason, in order | example 06, check 4; `tests/section/test_assemble_web_layout_contract.py`, `tests/section/test_spanwise_web_anchoring.py` |
| (R10.6) | 2D-F to 1C-F, 3D-F to 2D-F, no web raises | example 06, check 5; `tests/section/test_layout_step_down_cascade.py` |
| (R10.7), (R10.9) | Full, trimmed and scaled pads; the trimmed pad is symmetric and the next sample infeasible | example 06, check 6 |
| (R10.8) | Sample counts equal the formula for four band widths | example 06, check 7 |
| (R10.10)–(R10.12) | Scale factor equal to (R10.10); a pad with a non-positive band minimum runs out | example 06, check 8 |
| (R10.13) | Suction cap 50 mm clamped to $d(x_r)/2$, pressure cap kept, extent unchanged | example 06, check 9 |
| (R10.14) | Gate equals $(d - g)/2$; equal to the clamp at 13.3 mm; SG-041 printed | example 06, check 10; `tests/section/test_spar_cap_runout.py` |
| Section 3.11 | Flapwise web at 90° twist; assembly refusing a drop | example 06, check 11 |

```{figure} ../../examples/section/06_feasibility_checks/figures/feasibility_checks.png
:alt: Left, the cavity depth and cap allowance along a NACA 2412 with four cap pads; middle, cavities at 0.8c of a NACA 0006 kept or dropped against the probe spacing; right, the per-cap ceiling of the clamp and of the run-out gate against cavity depth
:width: 100%

**Figure R10.1.** (a) Cavity depth $d(x)$ and allowance $a(x)$ (R10.7) of a NACA 2412,
1 m chord, 6 mm skin, with four pads of half-width $0.04c$: authored band dotted,
kept extent solid at the height $h_u + h_l$ it keeps. (b) Cavities at $0.8c$ of a
NACA 0006 kept or dropped by the probe (R10.3), with the spacing $\Delta y$ of the
chord pad and of the thickness pad (SG-038). (c) Per-cap ceiling of the clamp
(R10.7) and of the run-out gate (R10.14) (SG-039). Drawn by
`examples/section/06_feasibility_checks/feasibility_checks.py`.
```

These checks move to V1 when it is written.

## 8. References

The intercept probe is the one of [R5](R5_section_layout.md), (R5.10); the offset
field is (R4.15) of [R4](R4_sdf_csg.md); the similarity scaling that keeps loft cuts
above the floors is (R8.3) of [R8](R8_loft_members.md). Full bibliography entries are in
[R0](R0_notation.md#bibliography).
