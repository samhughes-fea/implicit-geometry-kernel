# R5 · Section layout: skin, caps, webs and the WX-Y-Z-NR catalogue

## 1. Problem

A layout key such as `2D-F-LC-0R` names a structural concept. At each station the
kernel turns that key, a scaled aerofoil and three gauges (skin $t_s$, cap $h$,
web $t_w$) into the walls of the section: a skin, spar caps and shear webs, each a
field in the sense of [R4](R4_sdf_csg.md). This document derives the grammar and
counts of the catalogue, the step-down rule that maps a station that cannot carry
every web onto a smaller layout, the field of each wall, the web axes for the two
orientations, the cap bands of the four cap families, and the web stations along
the span for the three anchorings. The mid-surfaces and cells cut from these walls
are in [R6](R6_midsurfaces_cells.md); cap and web feasibility checks are in R10.

## 2. Definitions and assumptions

1. **Chord frame.** Section coordinates $(x, y)$ in metres with the leading edge at
   the origin, $x$ along the chord towards the trailing edge, before twist. The
   aerofoil $\Omega_A$ is the polygon of [R2](R2_aerofoils.md) scaled by the chord
   $c$; its field is $\phi_A$, (R4.9).
2. **Twist.** The section is turned by the angle $\theta$ passed as `twist_rad`: a
   chord-frame point $\mathbf{p}$ lands at $\mathbf{Q}(\theta)\mathbf{p}$ with
   $\mathbf{Q}$ of (R4.6). The blade frame of this document is that turned frame;
   its relation to the rotor frames is [R1](R1_frames.md).
3. **Gauges.** Skin thickness $t_s$ (the larger of the two sides when they differ),
   cap height $h$ (per side), web thickness $t_w$, all constant within a section.
4. **Mid-line of the aerofoil** $\bar y(x)$: the vertical mid-line of the polygon
   (`AirfoilSDF.camber_line`), the discrete form of (R2.19).
5. **Web stations** $x_1 < \dots < x_N$, in metres from the leading edge.
6. **Assumption A1.** $t_s$, $h$ and $t_w$ are small against the local aerofoil
   thickness and nose radius, so each wall's inner face is reached along the outer
   face's normal (Section 3.4 of R4). Where they are not, R10 thins or drops the
   wall.

## 3. Derivation

### 3.1 The catalogue

A key has five parts, `[W][X]-[Y]-[Z]-[R]`:

\begin{equation*}
\begin{gathered}
\texttt{key} = \underbrace{W X}_{\text{core}}\,\texttt{-}\,\underbrace{Y}_{\texttt{CN}\,|\,\texttt{F}}\,\texttt{-}\,\underbrace{Z}_{\texttt{LC}\,|\,\texttt{LE}\,|\,\texttt{PA}}\,\texttt{-}\,\underbrace{n\texttt{R}}_{n \ge 0},
\\
W\in\{0,\dots,5\},\qquad X\in\{\texttt{A},\texttt{B},\texttt{C},\texttt{D}\} .
\end{gathered}
\tag{R5.1}
\end{equation*}

$W$ is the number of shear webs, $X$ the cap family (A no caps, B one band, C a pad
on each web, D a continuous box between the outer webs), $Y$ the web orientation
(chord-normal or flapwise), $Z$ the spanwise anchoring of the webs (Section 3.6)
and $n$ the number of interior ribs. Zero-web cores have no orientation (`0A`,
`0B`). The anchoring and the rib count are required: `2D-F` and `2D-F-LC` are
refused.

Family D needs a fore and an aft web, so `1D` does not exist, and families C and D
need webs, so `0C` and `0D` do not. Counting the cores and then the keys with no
ribs,

\begin{equation*}
\begin{gathered}
N_{\mathrm{core}} = \underbrace{2}_{\texttt{0A},\,\texttt{0B}} + \sum_{W=1}^{5}\big(4 - [W = 1]\big)\times 2 = 2 + 38 = 40,
\\
N_{\mathrm{key}}^{0\mathrm{R}} = 3\,N_{\mathrm{core}} = 120,
\end{gathered}
\tag{R5.2}
\end{equation*}

the 40 entries of `SYSTEM_TYPE_LAYOUT_KEYS` (each with `-LC-0R`) and the 120 of
`SYSTEM_TYPE_KEYS`. Other rib counts resolve but are not enumerated. A core with $W$
webs has $W + 1$ cells.

### 3.2 Step-down and run-out

Webs are anchored along the span (Section 3.6), so outboard a station may be able
to carry only $n < W$ of them. The kernel then builds the station with a smaller
core, never a larger one. With the family order $\mathrm{D}\succ\mathrm{C}\succ
\mathrm{B}\succ\mathrm{A}$, the core for $n$ webs is

\begin{equation*}
\begin{gathered}
X'(n) = \max_{\succ}\big\{X' \preceq X :\ nX'\texttt{-}Y \in \mathcal{C},\ \operatorname{caps}(X') = \operatorname{caps}(X)\big\},
\\
n = W - 1, \dots, 1,
\end{gathered}
\tag{R5.3}
\end{equation*}

Here $\mathcal{C}$ is the set of cores of (R5.2) and the maximum is taken in the
family order, so the most elaborate admissible family wins. The orientation,
anchoring and rib count are copied from the parent. For A, B and C
this keeps the family down to one web; for D it gives D down to two webs and C at
one, since `1D` does not exist and C keeps the caps. The ladder floors at one web: a
station that cannot anchor a single web raises rather than becoming `0A` or `0B`,
because a section with no web has no shear path. Losing the caps is a separate
event (`layout_without_spar_caps`): a cavity too thin for a minimum-gauge cap runs
the caps out to family A with the same webs, anchoring and ribs.

### 3.3 Default web stations

A catalogue core places its webs at fractions $f_k$ of the chord from the leading
edge:

\begin{equation*}
(f_k) =
\begin{cases}
(0.50), & W = 1,\\
(0.15,\ 0.50), & W = 2,\\
\big(k/(W + 1)\big)_{k=1}^{W}, & W \ge 3 .
\end{cases}
\tag{R5.4}
\end{equation*}

The two-web pair is a fixed legacy choice; three or more webs divide the chord
evenly. A layout from a reference blade or from windIO gives its own stations
instead.

### 3.4 Wall fields

**Skin.** The skin is the slab of depth $t_s$ inside the outer mould line, (R4.18):

\begin{equation*}
\phi_{\mathrm{skin}} = \max\big(\phi_A,\ -\phi_A - t_s\big) .
\tag{R5.5}
\end{equation*}

The outer shape is held fixed: a thicker skin grows inwards. (R5.5) is exact inside
the wall where A1 holds, and its mid-surface is the inward offset of $\partial
\Omega_A$ by $t_s/2$. Other reference faces (`inner_mould_line`, `midsurface`) move
the faces but keep the solid authored from the outer mould line.

**Suction and pressure sides.** Each skin, and each cap, is split at the mid-line
$\bar y$ by the side fields

\begin{equation*}
\begin{gathered}
\sigma_{\mathrm{SS}}(x, y) = \bar y(x) - y,
\qquad
\sigma_{\mathrm{PS}}(x, y) = y - \bar y(x),
\\
\phi_{\mathrm{skin,SS}} = \max\big(\phi_{\mathrm{skin}},\ \sigma_{\mathrm{SS}}\big) .
\end{gathered}
\tag{R5.6}
\end{equation*}

Splitting at $\bar y$ rather than at the chord $y = 0$ keeps material of one skin
from crossing to the other on a cambered section. $\sigma$ is a vertical offset,
not a distance: its gradient has length $\sqrt{1 + \bar y'^2}$, so the side cut
of (R5.6) is $\sqrt{1 + \bar y'^2}$-Lipschitz rather than 1-Lipschitz (1.005 for a
mid-line slope of 0.1).

**Caps.** The inner mould line is the erosion of the aerofoil by $t_s$, (R4.16),
and a cap is the slab of depth $h$ inside it, on its own side, within its band, with
the webs cut out:

\begin{equation*}
\begin{gathered}
\phi_I = \phi_A + t_s,
\\
\phi_{\mathrm{cap}} = \max\Big(\phi_I,\ -\phi_I - h,\ \sigma_{\mathrm{side}},\ \beta(x, y),\ -\min_k \phi_{\mathrm{web},k}\Big) .
\end{gathered}
\tag{R5.7}
\end{equation*}

By (R4.12) and (R4.13) the cap field is exact inside the cap away from its band
edges and the webs. The band field $\beta$ is a pair of half-planes, either chordwise
or, when every web is flapwise, normal to the webs:

\begin{equation*}
\begin{gathered}
\beta_{\mathrm{chord}} = \max\big(x_0 - x,\ x - x_1\big),
\\
\beta_{\mathrm{strip}} = \max\big(c_{\mathrm{lo}} - \mathbf{n}\cdot\mathbf{p},\ \mathbf{n}\cdot\mathbf{p} - c_{\mathrm{hi}}\big),
\end{gathered}
\tag{R5.8}
\end{equation*}

with $\mathbf{n}$ the unit normal of the webs (R5.11) and $c_{\mathrm{lo}},
c_{\mathrm{hi}}$ the values of $\mathbf{n}\cdot\mathbf{p}$ on the two bounding web
axes. Both are exact distances inside the band (R4.12). The cap's mid-surface is the
inward offset of $\partial\Omega_A$ by $t_s + h/2$ within the band.

**Webs.** A web is a strip of width $t_w$ about its axis from $\mathbf{a}$ (suction
foot) to $\mathbf{b}$ (pressure foot), clipped to the cavity. The kernel builds the
strip as a capsule (R4.7) whose axis is extended by $e$ beyond both feet, so that its
round ends fall outside the cavity:

\begin{equation*}
\begin{gathered}
\mathbf{u} = \frac{\mathbf{a} - \mathbf{b}}{\lVert\mathbf{a} - \mathbf{b}\rVert},
\qquad
e = 2t_w + t_s,
\\
\phi_{\mathrm{web}} = \max\Big(u\big(\mathbf{p};\ \mathbf{a} + e\mathbf{u},\ \mathbf{b} - e\mathbf{u}\big) - \tfrac{t_w}{2},\ \phi_I\Big) .
\end{gathered}
\tag{R5.9}
\end{equation*}

Inside the cavity and along the web, $u - t_w/2$ is the distance to the strip's
sides, so (R5.9) is exact there. The round end stays outside the cavity when the
inner skin crosses the strip within $e$ of the foot; for an inner surface at angle
$\gamma$ to the web's normal that needs $\tfrac{t_w}{2}\tan\gamma \le e$, that is
$\tan\gamma \le 4 + 2t_s/t_w$ ($\gamma \lesssim 76°$), which holds everywhere except
very close to the nose.

### 3.5 Web axes

**Feet.** A web at station $x_k$ stands on the inner mould line. Its feet are the
outermost roots of $\phi_I$ on the line through $x_k$:

\begin{equation*}
\begin{gathered}
\phi_I(x_k, y) = 0,
\qquad
y_{\mathrm{top}} = \max\{y\},\quad y_{\mathrm{bot}} = \min\{y\},
\\
\mathbf{a} = (x_k, y_{\mathrm{top}}),\ \mathbf{b} = (x_k, y_{\mathrm{bot}}) .
\end{gathered}
\tag{R5.10}
\end{equation*}

The roots are bracketed on 240 samples of $y$ and refined by ten bisections and a
final secant step. A chord-normal (CN) web keeps this vertical axis: it is normal to
the chord in the untwisted section and turns with it.

**Flapwise webs.** A flapwise (F) web stays vertical in the twisted section, so in
the chord frame its axis leans by the twist. Its direction and unit normal in the
chord frame are $\mathbf{Q}(-\theta)$ applied to the twisted frame's vertical and
horizontal:

\begin{equation*}
\mathbf{u}_{\mathrm{F}} = (\sin\theta,\ \cos\theta),
\qquad
\mathbf{n} = (-\cos\theta,\ \sin\theta) .
\tag{R5.11}
\end{equation*}

The kernel finds its feet in the twisted frame: it places the web at
$x' = x_k\cos\theta - \bar y_k\sin\theta$, with $\bar y_k$ the mid-height of the CN
feet, solves (R5.10) on the vertical line $x = x'$ of the turned inner mould line,
and turns the two feet back by $\mathbf{Q}(-\theta)$. As $|\theta|$ approaches
$90°$ a flapwise web lies along the chord and cannot cross from skin to skin; the
kernel refuses it when

\begin{equation*}
|\cos\theta| < \sin 1° .
\tag{R5.12}
\end{equation*}

### 3.6 Cap bands

The band $[x_0, x_1]$ of (R5.8) depends on the cap family.

**Family B**, one band. Its centre is either the pitching anchor, a fraction
$f_p$ of the chord from the leading edge, or the station of maximum vertical
thickness; its half-width is $w_B$:

\begin{equation*}
\begin{gathered}
x_c = \begin{cases} f_p\,c, & \text{pitching},\\ \arg\max_x\big(y_u(x) - y_l(x)\big), & \text{max thickness},\end{cases}
\\
[x_0, x_1] = [x_c - w_B,\ x_c + w_B],
\qquad
f_p = \tfrac13,\quad w_B = 0.05\,c .
\end{gathered}
\tag{R5.13}
\end{equation*}

The catalogue's B layouts use the maximum-thickness anchor; the maximum is taken on
200 evenly spaced samples of the polygon's vertical thickness. A layout read from a
blade definition gives the band directly as chord fractions
`((lo, hi), ...)`, with $[x_0, x_1] = [\mathrm{lo}\,c, \mathrm{hi}\,c]$, and may give
each side its own band.

**Family C**, a pad on each web:

\begin{equation*}
[x_0, x_1]_k = [x_k - w_C,\ x_k + w_C],
\qquad
w_C = 0.04\,c .
\tag{R5.14}
\end{equation*}

**Family D**, a continuous box between the outer webs:

\begin{equation*}
[x_0, x_1] = [x_1,\ x_N],
\qquad
N \ge 2 ,
\tag{R5.15}
\end{equation*}

with the strip of (R5.8) between the first and last web axes when the webs are
flapwise. In every family the band is then clamped against the cavity
(`cap_feasibility_clamp`, R10) and the cap heights snapped to the ply pitch when one
is given.

### 3.7 Web stations along the span

A web is held along the span in one of three ways. With $c(z)$ the chord, $c_0$
the root chord, $p(z)$ the pitch-axis fraction and $f_k$ the authored fraction:

**Local chord (LC)**: the web keeps its chord fraction, sliding with the chord,

\begin{equation*}
x_k(z) = f_k\,c(z) .
\tag{R5.16}
\end{equation*}

**Leading edge (LE)**: the web keeps its distance from the leading edge, set at
the root,

\begin{equation*}
x_k(z) = f_k\,c_0 .
\tag{R5.17}
\end{equation*}

**Pitch axis (PA)**: the web keeps its distance $\xi_k$ from the pitch axis, so it
runs parallel to it on the planform. The first $\min(N, 3)$ webs divide the root
band evenly; further webs continue the same lattice, the first of them at the root
trailing edge. With $\xi_{\mathrm{LE},0} = -p_0 c_0$ and $\xi_{\mathrm{TE},0} = (1 -
p_0) c_0$,

\begin{equation*}
\begin{gathered}
\Delta\xi = \frac{\xi_{\mathrm{TE},0} - \xi_{\mathrm{LE},0}}{\min(N, 3) + 1},
\qquad
\xi_k = \xi_{\mathrm{LE},0} + k\,\Delta\xi,
\\
x_k(z) = \big(p(z) + \xi_k / c(z)\big)\,c(z) = p(z)\,c(z) + \xi_k .
\end{gathered}
\tag{R5.18}
\end{equation*}

An LE or PA web whose station leaves the inner mould line outboard (where the chord
has shrunk past it) cannot anchor there; the step-down of (R5.3) then removes it.

## 4. Algorithm

From a key to the walls at one station (`multicell_from_layout`, `MultiCellSection`):

1. Resolve the key by (R5.1) to a `SystemLayoutSpec`; read $W$, $X$, $Y$, $Z$, $n$.
2. Place the web stations: (R5.4), or the layout's own, mapped to this station by
   (R5.16)–(R5.18).
3. Find the feet (R5.10), or (R5.11) for flapwise webs; refuse (R5.12). If not every
   web can anchor, step the core down by (R5.3) (`resolve_station_layout`).
4. Build the skin (R5.5), split it by (R5.6).
5. Build the webs (R5.9).
6. Build the cap bands (R5.13)–(R5.15), clamp them to the cavity, and build the caps
   (R5.7) with the webs cut out.
7. Turn every wall by $\theta$ (R4.22).

## 5. Implementation

| Result | Function (`implicit_geometry_kernel` unless stated) |
|---|---|
| (R5.1), (R5.2) | `catalogue.system_layout.resolve_system_type`, `_LAYOUT_REGISTRY`, `SYSTEM_TYPE_LAYOUT_KEYS`, `SYSTEM_TYPE_KEYS`; `blade.concepts.BladeConcept.from_key` |
| (R5.3) | `catalogue.system_layout.effective_layout_for_web_count`, `step_down_cascade`, `family_min_webs`, `layout_without_spar_caps` |
| (R5.4) | `catalogue.system_layout._default_web_fracs` |
| (R5.5), (R5.6) | `section.walls.outer_skin.OuterSkin`; `section.assembly.clip_algebra.camber_side_field` |
| (R5.7), (R5.8) | `section.walls.spar_cap.SparCap`; `section.assembly.clips.chordwise_clip_expr`, `_parallel_web_strip_clip`; `section.assembly.build_spar_caps` (web cut) |
| (R5.9) | `section.walls.shear_web.ShearWebComponent` |
| (R5.10) | `section.assembly.web_layout.inner_y_at_x`, `inner_y_at_xs`, `resolve_web_anchors` |
| (R5.11), (R5.12) | `section.assembly.clips.flapwise_web_normal`; `section.assembly.web_layout.flapwise_web_sframe_inputs`; `section.assembly.plan._check_flapwise_webs_cross_chord` |
| (R5.13)–(R5.15) | `section.layout.cap_placement.CapPlacementRule`; `section.assembly.build_spar_caps._caps_fixed_band`, `_caps_continuous_box`, `build_spar_caps_pads.fill_per_web_pads`; `section.layout.section_axes.pitch_axis_x_from_le`, `max_thickness_chord_x` |
| (R5.16)–(R5.18) | `blade.concepts.web_anchoring.lc_web_fractions`, `le_web_offsets_m`, `pa_web_x_from_pitch`, `pa_regions`, `chord_fraction_from_pitch` |
| Algorithm | `section.assembly.from_layout.multicell_from_layout`, `section.assembly.plan.plan_section`, `section.assembly.multicell.MultiCellSection` |

## 6. Limits

1. **Thick walls.** (R5.5), (R5.7) and (R5.9) are exact only under A1. Near a sharp
   trailing edge or a small nose the inner face does not exist at the full gauge;
   R10 covers how the kernel thins, clamps or drops the wall.
2. **The side split** (R5.6) is vertical, not normal to the mid-line, and its field
   is not 1-Lipschitz; a field built with it is a bound only up to the factor
   $\sqrt{1 + \bar y'^2}$.
3. **The maximum-thickness anchor** (R5.13) is sampled on 200 points: on a NACA 2418
   of 2 m chord it is 1.4 mm (0.07 % of chord) aft of the fine-grid maximum.
4. **Web ends near the nose.** The capsule extension of (R5.9) leaves a rounded web
   end inside the cavity if the inner skin meets the web at more than about $76°$
   from the web's normal.
5. **Default stations** (R5.4) are catalogue defaults for figures and goldens, not
   design values; a real blade supplies its own.
6. **Feet** (R5.10) take the outermost two roots; a station where the line through
   $x_k$ crosses the inner mould line more than twice (a re-entrant section) is
   outside the derivation.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R5.1), (R5.2) | 40 cores and 120 keys, equal to the set built from the rules; six malformed keys refused | [section example 04](../../examples/section/04_layout_construction/README.md), check 1 |
| (R5.3) | Ladders from 5A, 5B, 5C, 5D; cap run-out keeps webs, anchoring and ribs | example 04, check 2 |
| (R5.4) | Default stations for one to five webs | example 04, check 3 |
| (R5.5), (R5.7), (R5.9) | Skin, cap and web fields along exact normals against $\max(-d, d - t_s)$, $\max(t_s - d, d - t_s - h)$ and $|s| - t_w/2$: $6.2\times10^{-17}$; web feet on $\phi_I = 0$ to $3.5\times10^{-18}$ | example 04, check 4 |
| (R5.11) | Flapwise axes at the twist from the chord normal for $-10°$ to $30°$: $3.1\times10^{-16}$ rad | example 04, check 5 |
| (R5.13)–(R5.15) | Bands of D, C and B equal the rule | example 04, check 6 |
| (R5.16)–(R5.18) | LC, LE and PA stations on a tapered planform; PA lattice | example 04, check 7 |
| (R5.2), all walls | Counts of webs, caps and cells of all 40 cores; cells tile the skin mid-surface | [section example 01](../../examples/section/01_layout_catalogue/README.md) |
| (R5.1)–(R5.3) | Registry, step-down and run-out | `tests/catalogue/test_registry.py`, `tests/section/test_layout_step_down_cascade.py` |

```{figure} ../../examples/section/04_layout_construction/figures/layout_construction.png
:alt: Left, the walls of a 2D-F section on a NACA 2418 in the chord frame; right, web stations for LC, LE and PA anchoring on a tapered planform
:width: 100%

**Figure R5.1.** Left: the walls of `2D-F` on a NACA 2418 of 2 m chord at 15° twist,
in the chord frame: the skin (R5.5) split at the mid-line (R5.6), the box caps
(R5.7) between the outer webs (R5.15), and the flapwise webs (R5.9), leaning by the
twist (R5.11). Right: two webs on a planform tapering from 2.0 to 0.8 m, held by
local chord (R5.16), from the leading edge (R5.17) and from the pitch axis (R5.18);
the second LE and PA webs leave the section near the tip. Drawn by
`examples/section/04_layout_construction/layout_construction.py`.
```

These checks move to [V1](../v/index.md) when it is written.

## 8. References

The walls are built from the fields and operations of [R4](R4_sdf_csg.md); the
aerofoil and its mid-line from [R2](R2_aerofoils.md). The catalogue is the kernel's
own; full bibliography entries are in [R0](R0_notation.md#bibliography).
