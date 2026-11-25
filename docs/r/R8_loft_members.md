# R8 · Spanwise members and the loft

## 1. Problem

A `BladeGeometry` holds stations, a layout and span laws. From them it builds the
section at any span position $z$ in the blade frame, the members that run along the
span or sit at one station, and smooth surfaces through the walls of successive
sections. This document derives the placement of a section in the blade frame and
its inverse, the scaling that builds thin walls, the tilted and anchored webs, the
cap edges given on the outer surface, the stations the members force, the locked
resampling that lofts named walls through topology changes, the partition patch
with its curvature and projection, and the vertex-normal offset of the skin.

The scalar laws (chord, twist, pitch axis, offset) and the aerofoil morph are
derived in [R3](R3_span_laws.md); the layout and the web stations along the span
in [R5](R5_section_layout.md), Section 3.7; the mid-surfaces and cells of one
section in [R6](R6_midsurfaces_cells.md).

## 2. Definitions and assumptions

1. **Cut frame.** A `SectionCut` stores $(x, y, z)$: $x$ along the chord towards
   the trailing edge (TE), $y$ towards the suction side, origin on the pitch axis,
   $z$ along the span. This is the swapped B frame of (R1.11).
2. **Unit section.** $(\xi, \eta)$ are chord-normalised aerofoil coordinates: the
   leading edge (LE) at $(0, 0)$, the TE at $\xi = 1$. $\eta_u(\xi)$ and
   $\eta_l(\xi)$ are the suction and pressure surfaces, sampled at the same $\xi_j$.
3. **Station laws.** At $z$: chord $c$, twist $\tau$ (degrees, nose up positive
   under `twist_sense="nose_up"`; `"nose_down"` negates the input), pitch-axis
   fraction $p$ and in-plane offset $\mathbf{o} = (o_x, o_y)$, all from R3.
4. **Thicknesses.** Skin $t_s$, cap height $h$ and web thickness $t_w$ at $z$, from
   the thickness laws; manufacturing minima $t_s^{\min} = 2.5$ mm,
   $t_w^{\min} = 4$ mm and $h^{\min} = 5$ mm (`section.layout.thickness_limits`).
5. **Run.** A named open polyline of the shell graph of a cut ([R6](R6_midsurfaces_cells.md),
   Section 3.8): a skin partition of the mid ring between two partition points, or
   a web. Its normalised arc length is $t \in [0, 1]$.
6. **Signature.** The ordered names of the loftable runs (skin, web, core and insert
   surfaces) of a cut. Ribs and stringers are not in it.

## 3. Derivation

### 3.1 Placement in the blade frame

A unit-section point is scaled by the chord, moved so that the pitch axis is the
origin, turned by the twist (R1.12) and shifted by the offset:

\begin{equation*}
\begin{pmatrix} x \\ y \end{pmatrix}
= \mathbf{Q}(\tau)\, c \begin{pmatrix} \xi - p \\ \eta \end{pmatrix} + \mathbf{o},
\qquad
\mathbf{Q}(\tau) = \begin{pmatrix} \cos\tau & \sin\tau \\ -\sin\tau & \cos\tau \end{pmatrix},
\tag{R8.1}
\end{equation*}

with $z$ unchanged. The outer mould ring is (R8.1) applied to the suction surface
from TE to LE followed by the pressure surface from LE to TE; a repeated TE point is
dropped once and the ring is closed by repeating its first point. Because
$\mathbf{Q}$ is a rotation (R1.5), the ring's area is $c^2$ times the unit
section's. Solving (R8.1) for $(\xi, \eta)$ with $\mathbf{Q}^{-1} = \mathbf{Q}^{\mathsf T}$:

\begin{equation*}
\begin{pmatrix} \xi \\ \eta \end{pmatrix}
= \frac{1}{c}\,\mathbf{Q}(\tau)^{\mathsf T}\big(\mathbf{x} - \mathbf{o}\big)
+ \begin{pmatrix} p \\ 0 \end{pmatrix}.
\tag{R8.2}
\end{equation*}

### 3.2 Thin walls by similarity

The section builder refuses walls thinner than the minima of Definition 4. A cut
whose walls are thinner is built at a larger scale. Every length of a section built
from a chord $c$ and wall thicknesses $t_i$ is a homogeneous function of degree one
in $(c, t_i)$ and the metre-valued layout inputs, so building with all of them
multiplied by $k$ and dividing the result by $k$ reproduces the section exactly. The
kernel takes

\begin{equation*}
\begin{gathered}
k = \max\Big(1,\ \max_i \frac{t_i^{\min}}{t_i}\Big)\,(1 + 10^{-9}) \quad (k > 1),
\\
\mathbf{x} = \mathbf{x}^{(k)}/k, \qquad A = A^{(k)}/k^2, \qquad t = t^{(k)}/k,
\end{gathered}
\tag{R8.3}
\end{equation*}

where the factor $1 + 10^{-9}$ keeps the scaled walls strictly above the minima. The
chord, thicknesses, web offsets and two-foot web axes are scaled; metre-valued
fields of a user `SectionLayout` are not, and land at $1/k$ of their value
([BG-058](../blade/bug_log.md#bg-058)).

### 3.3 Tilted webs

A web given by one chord fraction $f$ and a tilt $\beta$ (`web_tilt_deg`) is the line
through a pivot $(f, \eta_0)$ at angle $\beta$ from the thickness direction, in the
unit section. With $\theta = -\beta$, its points satisfy

\begin{equation*}
r(\xi, \eta) = (\xi - f)\cos\theta - (\eta - \eta_0)\sin\theta = 0,
\qquad
\mathbf{d} = (-\sin\beta,\ \cos\beta),
\tag{R8.4}
\end{equation*}

where $\mathbf{d}$ is the direction from the pressure end towards the suction end.
A positive tilt therefore puts the pressure foot aft of the suction foot
([BG-061](../blade/bug_log.md#bg-061)). The pivot is on the chord line,
$\eta_0 = 0$ (`web_tilt_anchor="chord"`), or on one surface, $\eta_0 = \eta_u(f)$ or
$\eta_l(f)$, which keeps that foot at $f$. Since (R8.1) scales both axes by $c$ and
rotates, the angle is the same in the cut.

A foot is where (R8.4) meets a surface polyline. On a segment $j$ where $r$ changes
sign, linear interpolation of $r$ gives the crossing

\begin{equation*}
u = \frac{r_j}{r_j - r_{j+1}},
\qquad
\xi_{\mathrm{foot}} = \xi_j + u\,(\xi_{j+1} - \xi_j),
\tag{R8.5}
\end{equation*}

and of several crossings the highest is taken on the suction side and the lowest on
the pressure side. If the line meets no segment, the foot stays at $f$.

A two-foot web (`web_positions` as pairs) is the line through the outer-mould points
at its two feet. Its axis points are moved $2t_s$ inward from each foot along the
line, inside the skin mid-surface for webs within about $75^\circ$ of the skin
normal, and the section builder extends the web from them to the mid ring
([R6](R6_midsurfaces_cells.md), (R6.13)). Feet closer than $4t_s$ are refused.

### 3.4 Anchors in metres

Webs frozen in metres, stringers given by `position_m`, and cap edges given by
`cap_x_from_pitch` are converted to chord fractions at each $z$. A distance $\xi_k$
from the pitch axis, or $x_k$ from the LE, gives

\begin{equation*}
f_k(z) = p(z) + \frac{\xi_k}{c(z)} \quad\text{(pitch axis)},
\qquad
f_k(z) = \frac{x_k}{c(z)} \quad\text{(leading edge)},
\tag{R8.6}
\end{equation*}

which is (R5.18) and (R5.17) solved for the fraction. The same law gives the web
identity used by the loft (`_web_chord_frac_law`), so a run keeps its name while
its fraction changes along the span.

### 3.5 Cap edges given on the outer surface

With `cap_edge_reference="outer"`, a cap edge at fraction $f$ is measured on the
outer mould line. The cap mid-surface lies $t_s + h/2$ inside it
([R6](R6_midsurfaces_cells.md), (R6.11)). The kernel moves the outer point
$\mathbf{x}_A = c\,(f, \eta(f))$ along the unit normal $\boldsymbol\nu$ of the
surface and keeps the chordwise coordinate:

\begin{equation*}
f' = \frac{1}{c}\Big(\mathbf{x}_A + \big(t_s + \tfrac{h}{2}\big)\,\boldsymbol\nu\Big)\cdot\mathbf{e}_\xi ,
\tag{R8.7}
\end{equation*}

where the tangent is the central difference over $f \pm 10^{-3}$, and
$\boldsymbol\nu$ is the normal to it that points towards $(c/2, 0)$. The shift is
$f' - f = (t_s + h/2)\,\nu_\xi / c$, of the order of $(t_s + h/2)\,\eta'(f)/c$. The
point $(c/2, 0)$ is inside the section only where the surface does not cross the
chord line; on a pressure surface above the chord line $\boldsymbol\nu$ points out
([BG-059](../blade/bug_log.md#bg-059)).

### 3.6 Members and forced stations

A rib or bulkhead (`TransverseMember`) is a plane member at one station $z_r$; the
tip closure is the face at the last station, carries the skin thickness unless given
its own, and is not a rib. A stringer exists on the cut where $z$ lies in its
`extent`; its position and depth follow their laws or (R8.6). Webs may start or end
inside the span (R5, Section 3.7). The stations where a member starts, ends or sits,

\begin{equation*}
\begin{gathered}
F = \{z_0, z_1 \text{ of each partial web}\} \cap (z_{\mathrm{root}}, z_{\mathrm{tip}})
\ \cup\ \{z_r \text{ of each rib}\}
\\
\cup\ \{\text{ends of each stringer extent}\}
\ \cup\ \{z_0, z_1 \text{ of each ply drop}\},
\end{gathered}
\tag{R8.8}
\end{equation*}

are returned sorted and unique by `forced_stations()`. A spanwise mesh that
contains them changes topology only at a station. The ply-drop ramps of (R3.19) end
at $z_0 - \ell$ and $z_1 + \ell$; those kinks of the thickness are not in $F$.

### 3.7 Loft zones and locked resampling

Lofting a run between two cuts joins its samples by index, so the run must exist on
both cuts and carry the same number of samples on each. Let $\sigma(z)$ be the
signature of the cut at $z$. A zone is a maximal run of consecutive stations with

\begin{equation*}
\sigma(z_i) = \sigma(z_{i+1}) = \dots = \sigma(z_j),
\tag{R8.9}
\end{equation*}

and a name present on only part of the stations starts or ends where it appears or
disappears: it is never interpolated across a missing station. Each run is
resampled to $n$ points equally spaced in arc length. With cumulative lengths $s_m$
at its vertices and total $L$, sample $i$ is at $\bar s_i = iL/(n-1)$, between the
vertices $m$ and $m+1$ with $s_m \le \bar s_i < s_{m+1}$:

\begin{equation*}
\mathbf{x}(\bar s_i) = \mathbf{x}_m + \frac{\bar s_i - s_m}{s_{m+1} - s_m}\,(\mathbf{x}_{m+1} - \mathbf{x}_m),
\qquad i = 0, \dots, n-1 .
\tag{R8.10}
\end{equation*}

Both ends are kept exactly, and adjacent skin partitions get their shared partition
point as an end. A web foot that lies inside a skin run must stay a vertex of the
run (SG-035). For pins at vertex indices $j_1 < \dots < j_M$, the output index of
pin $m$ is the rounded arc-length position, held inside the run and in order:

\begin{equation*}
\begin{gathered}
k_m = \min\big(\max(\rho_m,\ k_{m-1} + 1),\ n - 2 - (M - m)\big),
\\
\rho_m = \operatorname{round}\Big(\frac{s_{j_m}}{L}(n-1)\Big),
\qquad k_0 = 0,
\end{gathered}
\tag{R8.11}
\end{equation*}

and each piece between consecutive pins is resampled by (R8.10) with $k_{m+1} - k_m + 1$
points. This needs $n - 2 \ge M$. The web ends are then snapped to the nearest skin
sample, which is the pin itself; a gap larger than $10^{-9}$ of the section's size
is an error.

### 3.8 The partition patch

`PartitionPatch` rebuilds one run as a surface $\mathbf{x}(t, z)$. Across the run,
the points of each station are splined by their normalised arc length $t_j$ with a
not-a-knot cubic spline. Its knot second derivatives $M_j$ solve, with $h_j = t_{j+1} - t_j$
and $\delta_j = (\mathbf{x}_{j+1} - \mathbf{x}_j)/h_j$,

\begin{equation*}
\begin{gathered}
h_{j-1}M_{j-1} + 2(h_{j-1} + h_j)M_j + h_j M_{j+1} = 6(\delta_j - \delta_{j-1}),
\\
h_1 M_0 - (h_0 + h_1)M_1 + h_0 M_2 = 0,
\end{gathered}
\tag{R8.12}
\end{equation*}

for the interior knots $0 < j < n-1$, and the mirror image of the last row at the
other end. The end rows say that the
third derivative is continuous at the second and second-last knots, so the end
curvature is not forced to zero. Three knots give the parabola,
$M = 2(\delta_1 - \delta_0)/(h_0 + h_1)$, and two a line. On $[t_j, t_{j+1}]$, with
$a = (t_{j+1} - t)/h_j$ and $b = 1 - a$,

\begin{equation*}
\mathbf{x}(t) = a\mathbf{x}_j + b\mathbf{x}_{j+1}
+ \big((a^3 - a)M_j + (b^3 - b)M_{j+1}\big)\frac{h_j^2}{6},
\tag{R8.13}
\end{equation*}

whose derivatives in $t$ follow by differentiating $a$ and $b$. A cubic is
reproduced exactly. Each station's spline is sampled at $n_t$ equally spaced
$t_i$; along the span, each row $\mathbf{x}(t_i, \cdot)$ is a PCHIP law in $z$
(R3.3)–(R3.9), density-checked by (R3.12), or linear (a ruled bay). At any $z$ the
rows give a column $\mathbf{c}(z) = (\mathbf{x}(t_i, z))_i$, which is splined in $t$
again by (R8.12)–(R8.13). The spline is linear in its data, so it commutes with
differentiation in $z$:

\begin{equation*}
\partial_z^m \mathcal{S}[\mathbf{c}]\,(t) = \mathcal{S}[\partial_z^m \mathbf{c}]\,(t),
\qquad m = 0, 1, 2,
\tag{R8.14}
\end{equation*}

where $\mathcal{S}$ is the spline operator. $\mathbf{x}_z$, $\mathbf{x}_{tz}$ and
$\mathbf{x}_{zz}$ are therefore exact derivatives of the patch, with no finite
differences.

The unit normal is $\mathbf{n} = \mathbf{x}_t \times \mathbf{x}_z / |\mathbf{x}_t \times \mathbf{x}_z|$.
With the first and second fundamental forms {cite:p}`docarmo2016`,

\begin{equation*}
\begin{gathered}
\mathrm{I} = \begin{pmatrix} \mathbf{x}_t\cdot\mathbf{x}_t & \mathbf{x}_t\cdot\mathbf{x}_z \\
\mathbf{x}_t\cdot\mathbf{x}_z & \mathbf{x}_z\cdot\mathbf{x}_z \end{pmatrix},
\qquad
\mathrm{II} = \begin{pmatrix} \mathbf{x}_{tt}\cdot\mathbf{n} & \mathbf{x}_{tz}\cdot\mathbf{n} \\
\mathbf{x}_{tz}\cdot\mathbf{n} & \mathbf{x}_{zz}\cdot\mathbf{n} \end{pmatrix},
\\
\mathrm{II}\,\mathbf{v} = \kappa\,\mathrm{I}\,\mathbf{v},
\end{gathered}
\tag{R8.15}
\end{equation*}

the principal curvatures $\kappa_1 \le \kappa_2$ are the eigenvalues of
$\mathrm{I}^{-1}\mathrm{II}$, positive where the surface bends towards
$\mathbf{n}$, and the principal directions are $v_t\mathbf{x}_t + v_z\mathbf{x}_z$.
Curvatures below `curvature_tol` are set to zero and name the point flat, singly or
doubly curved, or umbilic. On a cylinder of mid radius $r_m$, $|\kappa| = 1/r_m$
across the run and 0 along the span.

The closest point to $\mathbf{p}$ minimises $\tfrac12|\mathbf{x}(t, z) - \mathbf{p}|^2$.
With $\mathbf{r} = \mathbf{x} - \mathbf{p}$, Newton's method on its gradient takes
the step $\Delta$ from

\begin{equation*}
\begin{pmatrix}
\mathbf{x}_t\cdot\mathbf{x}_t + \mathbf{r}\cdot\mathbf{x}_{tt} & \mathbf{x}_t\cdot\mathbf{x}_z + \mathbf{r}\cdot\mathbf{x}_{tz} \\
\mathbf{x}_t\cdot\mathbf{x}_z + \mathbf{r}\cdot\mathbf{x}_{tz} & \mathbf{x}_z\cdot\mathbf{x}_z + \mathbf{r}\cdot\mathbf{x}_{zz}
\end{pmatrix}\Delta
= -\begin{pmatrix} \mathbf{r}\cdot\mathbf{x}_t \\ \mathbf{r}\cdot\mathbf{x}_z \end{pmatrix},
\tag{R8.16}
\end{equation*}

clipped to $[0, 1] \times [z_0, z_{\mathrm{end}}]$. At convergence $\mathbf{r}$ is
normal to the patch, the signed distance is $(\mathbf{p} - \mathbf{x})\cdot\mathbf{n}$,
and its gradient at $\mathbf{p}$ is $\mathbf{n}$.

### 3.9 The vertex-normal skin offset

`offset_rings` moves every vertex of a placed outer mould ring a distance $d$ along
its vertex normal, the normalised sum of the unit normals $\boldsymbol\nu_{i-1,i}$
and $\boldsymbol\nu_{i,i+1}$ of its two edges, oriented by the vector to the vertex
mean $\bar{\mathbf{x}}$:

\begin{equation*}
\begin{gathered}
\mathbf{n}_i = \frac{\boldsymbol\nu_{i-1,i} + \boldsymbol\nu_{i,i+1}}{|\boldsymbol\nu_{i-1,i} + \boldsymbol\nu_{i,i+1}|},
\qquad
s_i = \operatorname{sign}\big((\bar{\mathbf{x}} - \mathbf{x}_i)\cdot\mathbf{n}_i\big),
\\
\mathbf{x}_i' = \mathbf{x}_i + d\,s_i\,\mathbf{n}_i .
\end{gathered}
\tag{R8.17}
\end{equation*}

$\mathbf{n}_i$ bisects the turning angle $\alpha_i$ at the vertex, so it makes the
angle $\alpha_i/2$ with each edge normal. The edge between two moved vertices is
therefore displaced by

\begin{equation*}
d_{\mathrm{edge}} = d\cos\frac{\alpha_i}{2} = d\Big(1 - \frac{\alpha_i^2}{8} + O(\alpha_i^4)\Big),
\tag{R8.18}
\end{equation*}

$d\cos(\pi/N)$ on a regular $N$-gon, against $d$ for an exact offset (R6.2): the
error is second order in the turning angle, so it vanishes as the sampling is
refined except at true corners such as a sharp TE. The orientation $s_i$ is right
only while $\bar{\mathbf{x}}$ lies on the inner side of every vertex's tangent;
on thin cambered sections it does not ([BG-060](../blade/bug_log.md#bg-060)).

## 4. Algorithm

`BladeGeometry.section_at(z)`:

1. Round $z$ to $10^{-9}$, clamp it to the blade, and return the cached cut if there
   is one (least recently used, 256 cuts).
2. Interpolate the station laws and the aerofoil at $z$ (R3), and the layout at $z$:
   web fractions or two-foot webs, metre anchors by (R8.6), cap windows (converted by
   (R8.7) when they are given on the outer surface).
3. Find $k$ by (R8.3). For tilted webs find the feet by (R8.4)–(R8.5); for two-foot
   webs set the axis points of Section 3.3 at scale $k$.
4. Build the section at chord $kc$ with thicknesses $kt_i$ and web offsets $kx_k$,
   with midlines and cell loops (R5, R6).
5. Bond caps to the skin and split skins at cap footprints; place every midline and
   cell loop by $1/k$ and (R8.1); store wall thicknesses $t/k$ and cell areas $A/k^2$.
6. Attach the ribs, stringers and tip closure that exist at $z$ (Section 3.6), the
   dropped webs, and the realised web count and family. With
   `merge_cap_and_skin=True`, replace the stacked walls by one wall per shell-graph run.

`locked_mid_surface_patches(graphs, n)`: resample each cut's runs by (R8.10)–(R8.11)
with its own plan, share partition points, snap web ends to skin samples, and stack
each name's samples over each maximal set of consecutive stations that carry it
(R8.9).

## 5. Implementation

| Result | Function |
|---|---|
| (R8.1) | `blade.loft.span_place.place_oml`, `place_chord_normal`; `BladeGeometry.oml_polyline`, `oml_polylines`; `blade.cut._twist_offset` |
| (R8.2) | `blade.loft.span_place.to_chord_normal` |
| (R8.3) | `blade.loft.loft_cut._min_thickness_scale`, `_scale_gauge`; `place` in `LoftCutMixin.section_at` |
| (R8.4), (R8.5) | `blade.loft.loft_cut._tilted_feet`, `_surface_crossing`; two-foot axes `_with_web_axes` |
| (R8.6) | `blade.concepts.web_anchoring.chord_fraction_from_pitch`; `LoftLayoutMixin._web_chord_frac_law`, `_stringer_on_chord`, `_cap_metres_at` |
| (R8.7) | `LoftLayoutMixin._outer_fraction_to_midplane`, `_cap_windows_from_outer_surface` |
| (R8.8) | `BladeGeometry.forced_stations`; `_transverse_at`, `_stringers_at`, `_tip_closure_at`; `blade.loft.loft._split_tip_closure` |
| (R8.9) | `blade.loft.locked_loft.loft_signature`, `shell_loft_zones`, `locked_mid_surface_patches` |
| (R8.10) | `blade.loft.span_place.resample_open_many`; `blade.loft.locked_loft.arc_resample_open` |
| (R8.11) | `blade.loft.locked_loft._resample_with_pins`, `resample_locked` |
| (R8.12), (R8.13) | `blade.patch._spline_second`, `_spline_eval` |
| (R8.14)–(R8.16) | `blade.patch.PartitionPatch.derivatives`, `at`, `closest` |
| (R8.17) | `blade.loft.span_place.offset_rings` (used by `mesh.shell._stacked_shell`) |

## 6. Limits

1. **BG-058 (open).** Metre-valued fields of a user `SectionLayout`
   (`web_stations_m`, a fixed cap half-width) are not scaled by $k$, so they shrink
   by $1/k$ on a cut whose walls are thinner than the minima.
2. **BG-059 (open).** Cap edges given on the outer surface move outward on a
   pressure surface above the chord line. The error is
   $2(t_s + h/2)|\nu_\xi|/c$ in chord fraction, $2\times10^{-3}$ on the example.
3. **BG-060 (open).** The vertex-normal offset (R8.17) moves vertices outward on
   thin cambered sections. It is not yet on a production path.
4. **BG-061 (open).** The tilt sense is stated wrongly in a comment; (R8.4) gives
   the sense the code implements.
5. **Tilted feet** (R8.5) are on the linear interpolant of the surface samples;
   a line that misses a surface leaves the foot at $f$ without a warning.
6. **Outer-surface tangent** (R8.7) is a central difference of step $10^{-3}$ on
   the linear interpolant; its error is of the order of the sample spacing.
7. **Forced stations** (R8.8) omit the outer ends of ply-drop ramps.
8. **Locked resampling** (R8.11) moves a pin to the rounded index, so the samples
   either side of it are spaced unevenly by up to half a sample.
9. **Partition patch.** The run's points set its accuracy: on a circle the
   curvature is within about 0.1 % of $1/r_m$ away from the ends; a vertex the
   section inserted near an end (a cap edge or web foot) can move the end curvature
   by several per cent. The run must have the same orientation at every station.
   A PCHIP row is C1, so $\mathbf{x}_{zz}$ jumps at stations; a linear row has
   $\mathbf{x}_{zz} = 0$ inside a bay.
10. **Newton projection** (R8.16) converges from a start near the foot; the
    Hessian is singular where $\mathbf{p}$ is at a centre of curvature.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R8.1), (R8.2) | Placed OML equals the hand map; inverse and re-placement: $4\times10^{-16}$ m; area ratio $c^2$ | [blade example 19](../../examples/blade/19_loft_checks/README.md), check 1 |
| (R8.3) | A 0.4 m chord with a 1 mm skin ($k = 2.5$) equals 0.4 times the 1 m cut: $3\times10^{-14}$; BG-058 printed | example 19, check 2 |
| (R8.4), (R8.5) | Both feet on the tilted line: $5\times10^{-17}$; one-surface anchors keep their foot | example 19, check 3; `tests/test_web_tilt.py`, `tests/test_skeleton.py` |
| (R8.6) | LE-frozen webs at 0.3 and 0.7 m on a chord tapering 2 m to 1 m: $10^{-16}$ m | example 19, check 4 |
| (R8.7) | Suction side within $3\times10^{-5}$ of the normal offset; BG-059 printed | example 19, check 5 |
| (R8.8) | Rib, stringer and ply-drop ends collected | example 19, check 6 |
| (R8.9)–(R8.11) | 13 runs of a 3D-F cut at 41 points; web ends on skin samples; pins at the rounded indices | example 19, check 7; `tests/test_locked_loft.py`; [blade example 12](../../examples/blade/12_locked_loft_zones/README.md) |
| (R8.12)–(R8.16) | Spline exact on a cubic; cylinder curvature within $9\times10^{-4}$; projection recovers a 0.03 m offset | example 19, check 8; `tests/test_partition_patch.py` |
| (R8.17), (R8.18) | Edge displacement $d\cos(\pi/N)$ to $10^{-16}$; BG-060 counts printed | example 19, check 9 |

```{figure} ../../examples/blade/19_loft_checks/figures/loft_checks.png
:alt: Top left, a unit aerofoil, the same scaled about the pitch axis, and twisted and offset; top right, three tilted webs through their pivots on a symmetric section; bottom left, the patch curvature on a cylinder against the exact value; bottom right, a NACA 4406 offset by vertex normals, with the vertices moved outward marked
:width: 100%

**Figure R8.1.** Top left: placement (R8.1) of a NACA 2418 with $c = 1.7$ m,
$p = 0.3$, $\tau = 12^\circ$ and an offset. Top right: webs tilted about $(f, 0)$
(R8.4). Bottom left: $|\kappa|\,r_m$ on the skin patch of a cylinder (R8.15).
Bottom right: the vertex-normal offset (R8.17) of a NACA 4406; the vertices moved
outward are BG-060. Drawn by `examples/blade/19_loft_checks/loft_checks.py`.
```

These checks move to V1 when it is written.

## 8. References

The placement uses the frames of [R1](R1_frames.md), the span laws of
[R3](R3_span_laws.md) and the layout of [R5](R5_section_layout.md). The fundamental
forms and the shape operator are as in {cite:t}`docarmo2016`; the not-a-knot spline
is that of {cite:t}`deboor2001`. Full bibliography entries are in
[R0](R0_notation.md#bibliography).
