# R2 · Aerofoil parametrisation and families

## 1. Problem

Each station of a blade carries an aerofoil: a closed curve in chord-normalised
coordinates that the loft scales by the chord, places at the pitch axis and twists
([R1](R1_frames.md)). This document derives how `Airfoil` stores that curve, how
each family it builds is constructed (NACA 4-digit, 5-digit and 6-series, CST,
the ellipse), and what each shape operator (thickness scaling, flatback, blend)
does to the surfaces. How stations are blended along the span is in
[R3](R3_span_laws.md).

## 2. Definitions and assumptions

1. **Chord-normalised coordinates** $(x, y) = (x_c/c,\ y_t/c)$, with $x_c, y_t$ the
   aerofoil coordinates of [R1](R1_frames.md): $x$ runs from the leading edge (LE,
   $x=0$) to the trailing edge (TE, $x=1$) and $y$ is positive towards the suction
   side.
2. **Upper and lower surfaces** $y_u(x)$ and $y_l(x)$: the parts of the contour
   from the LE to the TE on the suction and pressure sides. The LE is the contour
   point farthest from the mid-point of the two TE points.
3. **Camber and thickness, NACA sense.** A mean (camber) line $y_c(x)$ and a
   half-thickness $h(x)$ laid off normal to it {cite:p}`abbott1959`. The thickness
   ratio $t$ is the largest value of $2h$.
4. **Camber and thickness, kernel sense.** For the shape operators of Section 3.6,
   the vertical mid-line $\bar y(x) = \tfrac12(y_u + y_l)$ and the vertical thickness
   $\delta(x) = y_u - y_l$. These differ from item 3 by terms of order
   $h\,(\mathrm{d}y_c/\mathrm{d}x)^2$ (Section 3.6).
5. **Assumption A1.** Each surface is a single-valued function $y(x)$ on
   $0 \le x \le 1$. It holds everywhere except within about $10^{-4}$ of the LE,
   where the tangent is vertical (Section 6).
6. **Selig order** lists the contour TE → upper surface → LE → lower surface → TE;
   **Lednicer order** lists the upper surface LE → TE, then the lower surface LE →
   TE {cite:p}`selig2026`.

## 3. Derivation

### 3.1 Storage of a surface

Given points of one surface $\mathbf{p}_0, \dots, \mathbf{p}_N$ from the LE to the
TE, let $u$ be the cumulative chord length,

\begin{equation*}
u_0 = 0,\qquad u_{i} = u_{i-1} + \lVert \mathbf{p}_i - \mathbf{p}_{i-1} \rVert ,
\tag{R2.1}
\end{equation*}

and fit two cubic splines $x(u)$ and $y(u)$ through the points (not-a-knot ends).
Parametrising by $u$ rather than by $x$ keeps the curve smooth through the LE, where
$\mathrm{d}y/\mathrm{d}x$ is unbounded. The splines are sampled at 4001 values of
$u$ clustered towards the LE, $u = u_N(1 - \cos\tfrac{\pi}{2}\varsigma)$ for
$\varsigma\in[0,1]$, together with the given $u_i$ so that the surface passes
through every given point; the samples are then sorted by $x$. A query at a
chord fraction interpolates these samples linearly after mapping the fraction onto
the stored extent $[x_0, x_1]$:

\begin{equation*}
y_u(x) = \mathcal{L}\big[\,x_0 + x\,(x_1 - x_0)\,\big],
\tag{R2.2}
\end{equation*}

where $\mathcal{L}$ is piecewise-linear interpolation of the sorted samples, and
likewise for $y_l$. A contour whose LE sits slightly ahead of $x = 0$ (a cambered
NACA section, Section 3.2) is therefore renormalised to unit chord. Two TE forms get
special treatment: a *corner* (a last segment that turns more than 20° and more than
three times the turn before it) ends the spline one point early and keeps a straight
last segment, and a *round* closed end (steeper than 45°) gets a vertical end
tangent. Both stop the spline bulging aft of the TE, where (A1) would fail.

### 3.2 NACA 4-digit sections

The code $MPTT$ gives the maximum camber, its chordwise position and the thickness
ratio {cite:p}`abbott1959`:

\begin{equation*}
m = \frac{M}{100},\qquad p = \frac{P}{10},\qquad t = \frac{TT}{100}.
\tag{R2.3}
\end{equation*}

Points are spaced by the cosine rule, which clusters them at both ends where the
curvature is largest:

\begin{equation*}
x_j = \tfrac12\left(1 - \cos\beta_j\right),\qquad \beta_j = \frac{j\pi}{n-1},\quad j = 0,\dots,n-1 .
\tag{R2.4}
\end{equation*}

**Thickness.** The half-thickness is

\begin{equation*}
\begin{gathered}
h(x) = \frac{t}{0.2}\left(a_0\sqrt{x} + a_1 x + a_2 x^2 + a_3 x^3 + a_4 x^4\right),\\
(a_0,\dots,a_3) = (0.2969,\,-0.1260,\,-0.3516,\,0.2843),
\end{gathered}
\tag{R2.5}
\end{equation*}

with the published $a_4 = -0.1015$, which leaves a TE of thickness
$2h(1) = 10t\,(a_0+a_1+a_2+a_3+a_4) = 0.021\,t$. Setting $h(1)=0$
and solving (R2.5) for $a_4$ gives the closed-TE coefficient the kernel uses by
default:

\begin{equation*}
a_4 = -(a_0 + a_1 + a_2 + a_3) = -0.1036 .
\tag{R2.6}
\end{equation*}

The $\sqrt{x}$ term sets the LE radius: near $x = 0$, $h \approx 5t\,a_0\sqrt{x}$,
which is the parabola $x = h^2/(5ta_0)^2$ with vertex radius
$\rho_{LE} = \tfrac12 (5ta_0)^2 = 1.1019\,t^2$, the value given by Abbott and von Doenhoff.

**Camber.** The mean line is two parabolas that meet with zero slope at $x = p$:

\begin{equation*}
y_c(x) =
\begin{cases}
\dfrac{m}{p^2}\left(2px - x^2\right), & 0 \le x \le p,\\[2ex]
\dfrac{m}{(1-p)^2}\left(1 - 2p + 2px - x^2\right), & p < x \le 1 .
\end{cases}
\tag{R2.7}
\end{equation*}

Both branches give $y_c(p) = m$, $y_c(0) = y_c(1) = 0$. Differentiating,

\begin{equation*}
\frac{\mathrm{d}y_c}{\mathrm{d}x} =
\begin{cases}
\dfrac{2m}{p^2}(p - x), & x \le p,\\[2ex]
\dfrac{2m}{(1-p)^2}(p - x), & x > p,
\end{cases}
\qquad
\phi = \arctan\frac{\mathrm{d}y_c}{\mathrm{d}x},
\tag{R2.8}
\end{equation*}

which vanishes at $x = p$ from both sides, so the mean line is $C^1$ and its
maximum is $m$ at $p$. For $m = 0$ or $p = 0$ the section is symmetric,
$y_c \equiv 0$.

**Surfaces.** The half-thickness is laid off along the unit normal of the mean
line, $(-\sin\phi, \cos\phi)$:

\begin{equation*}
\begin{aligned}
x_u &= x - h\sin\phi, & y_u &= y_c + h\cos\phi,\\
x_l &= x + h\sin\phi, & y_l &= y_c - h\cos\phi .
\end{aligned}
\tag{R2.9}
\end{equation*}

Because $\phi > 0$ near the LE of a cambered section, $x_u < 0$ there: the
surface reaches slightly ahead of $x = 0$, which (R2.2) absorbs by renormalising.
With a closed TE the two last points are set to their mid-point so that the contour
closes exactly.

### 3.3 NACA 5-digit sections

The code $LPQTT$ gives the design lift coefficient, the position of maximum camber,
a reflex flag and the thickness ratio {cite:p}`abbott1959`:

\begin{equation*}
C_{l,i} = \frac{3L}{20},\qquad p = \frac{P}{20},\qquad
Q \in \{0, 1\},\qquad t = \frac{TT}{100},
\tag{R2.10}
\end{equation*}

where $Q = 1$ marks a reflexed mean line. The thickness is (R2.5). The normal ($Q = 0$) mean line is a cubic ahead of a
join $x = r$ and a straight line behind it:

\begin{equation*}
y_c(x) =
\begin{cases}
\dfrac{k_1}{6}\left(x^3 - 3rx^2 + r^2(3 - r)\,x\right), & 0 \le x \le r,\\[2ex]
\dfrac{k_1 r^3}{6}\,(1 - x), & r < x \le 1 .
\end{cases}
\tag{R2.11}
\end{equation*}

At $x = r$ the first branch equals $\tfrac{k_1}{6}(r^3 - r^4)$ and has slope
$\tfrac{k_1}{6}(3r^2 - 6r^2 + 3r^2 - r^3) = -\tfrac{k_1 r^3}{6}$, the value and slope
of the second branch, so the mean line is $C^1$. Maximum camber is where the
slope of the first branch vanishes. Setting $\mathrm{d}y_c/\mathrm{d}x = 0$ at
$x = p$ gives the equation for the join:

\begin{equation*}
3p^2 - 6pr + r^2(3 - r) = 0 ,
\tag{R2.12}
\end{equation*}

whose root in $(p, 1)$ reproduces the tabulated joins to $3\times10^{-4}$
(for example $p = 0.15 \Rightarrow r = 0.2027$, tabulated 0.2025). The factor $k_1$
fixes the design lift. Thin-aerofoil theory is linear in the camber, so $k_1$ is
the tabulated value for $C_{l,i} = 0.3$ scaled by $C_{l,i}/0.3$; for NACA 23012,
$r = 0.2025$ and $k_1 = 15.957$, which gives a maximum camber of 0.0184 at
$x = 0.15$. The surfaces follow from (R2.9).

**Kernel status.** `naca_five_digit_vertices` does not evaluate (R2.10)–(R2.12) as
stated: it takes the two digits $PQ$ together as $r = PQ/200$ and reads $k_1$ and
a scale from a table that is not Abbott and von Doenhoff's. For NACA 23012 its mean
line peaks at 0.0017 at $x = 0.11$ instead of 0.0184 at 0.15. This is logged as
[BG-055](../blade/bug_log.md#bg-055) (open). The thickness, and hence $t$, is
correct.

### 3.4 NACA 6-series sections

A 6-series section is defined by a conformal-mapping thickness form and a
uniform-load mean line {cite:p}`abbott1959`, which the kernel does not construct.
Instead `naca6(family, t_percent)` holds one tabulated reference section per family
(NACA 63-415 and NACA 64₁-212, unit chord, from the UIUC database
{cite:p}`selig2026`) and scales its ordinates:

\begin{equation*}
(x, y) = \left(x_{\mathrm{ref}},\ \frac{t}{t_{\mathrm{ref}}}\,y_{\mathrm{ref}}\right),
\qquad t_{\mathrm{ref}} = \max y_{\mathrm{ref}} - \min y_{\mathrm{ref}} .
\tag{R2.13}
\end{equation*}

Two consequences follow from (R2.13). Because the whole ordinate is scaled, the
camber of the reference section scales with $t$ too: `naca6(64, 18)` has a
maximum camber of 0.0166, 1.5 times that of the 64₁-212 it came from. And because
$t_{\mathrm{ref}}$ takes the upper maximum and lower minimum at different $x$, it
overestimates the thickness ratio slightly, so the result is about $3.5\times10^{-4}$
thinner than asked. The design-lift digit is ignored. The first is logged as
[BG-056](../blade/bug_log.md#bg-056) (open).

### 3.5 Class-shape transformation (CST)

Kulfan's CST writes each surface as a class function, which fixes the LE and TE
form, times a shape function {cite:p}`kulfan2008`. With a round nose and a sharp
TE the class function is

\begin{equation*}
C(x) = x^{1/2}(1 - x),
\tag{R2.14}
\end{equation*}

and the shape function is a Bernstein polynomial of degree $n$ with weights
$A_0,\dots,A_n$:

\begin{equation*}
S(x) = \sum_{i=0}^{n} A_i \binom{n}{i} x^i (1 - x)^{n - i} .
\tag{R2.15}
\end{equation*}

The kernel's surfaces, with a TE thickness $\Delta$ shared equally between them,
are

\begin{equation*}
y_u = C(x)\,S_u(x) + \tfrac12 \Delta x,\qquad
y_l = -\left(C(x)\,S_l(x) + \tfrac12 \Delta x\right).
\tag{R2.16}
\end{equation*}

The minus sign in $y_l$ means *positive lower weights put the lower surface below the
chord*. Kulfan writes $y_l = C S_l + \Delta_l x$ with negative lower weights, so
weights taken from the literature need their lower-surface signs reversed. Near
the LE, $S \to A_0$ and $y \approx A_0\sqrt{x}$, the parabola of Section 3.2 with
vertex radius

\begin{equation*}
\rho_{LE} = \tfrac12 A_0^2 ,
\tag{R2.17}
\end{equation*}

so $A_0$ alone sets the nose radius. Equal upper and lower weights give a symmetric
section, $y_l = -y_u$. At $x = 1$, $C = 0$ and $y_u - y_l = \Delta$.

### 3.6 Ellipse and shape operators

**Ellipse.** A root section of thickness ratio $t$ (a circle for $t = 1$) is

\begin{equation*}
x = \tfrac12\left(1 + \cos\varphi\right),\qquad y = \tfrac12 t \sin\varphi,
\qquad 0 \le \varphi \le 2\pi,
\tag{R2.18}
\end{equation*}

which has unit chord and maximum thickness $t$ at $x = \tfrac12$.

**Mid-line and thickness.** The operators below act on the vertical split of
item 4, Section 2:

\begin{equation*}
\bar y(x) = \tfrac12\left(y_u + y_l\right),\qquad \delta(x) = y_u - y_l .
\tag{R2.19}
\end{equation*}

For a NACA section, (R2.9) puts the upper point that lies above $x$ at the mean-line
station $x + h\sin\phi$. Expanding $y_c$ there to first order in $h$,
$y_u(x) = y_c + h(\tan\phi\sin\phi + \cos\phi) = y_c + h\sec\phi$, and likewise
$y_l(x) = y_c - h\sec\phi$, so $\bar y = y_c$ and
$\delta = 2h\sec\phi = 2h\,(1 + \tfrac12\phi^2 + \dots)$: the vertical and normal
thicknesses agree to second order in the camber slope, the vertical one being
larger. `Airfoil.thickness` returns the analytic $2h$
for NACA 4-digit sections and $\delta$ for all others.

**Thickness scaling** (`with_thickness`). The section is rescaled about its mid-line so
that its largest vertical thickness becomes $t^\ast$:

\begin{equation*}
y_{u,l} = \bar y \pm \tfrac12\,\delta\,\frac{t^\ast}{\max\delta} .
\tag{R2.20}
\end{equation*}

**Flatback** (`with_te_thickness`). A TE thickness $\Delta^\ast$ is added by
ramping both surfaces linearly from the LE:

\begin{equation*}
y_u \leftarrow y_u + \tfrac12 x\,\left(\Delta^\ast - \delta(1)\right),\qquad
y_l \leftarrow y_l - \tfrac12 x\,\left(\Delta^\ast - \delta(1)\right),
\tag{R2.21}
\end{equation*}

so the LE is unchanged and $\delta(1) = \Delta^\ast$ exactly.

**Blend** (`blend`, `cylinder_to_airfoil`). Two sections are combined surface by
surface at the same $x$:

\begin{equation*}
y_{u,l}(x; s) = (1 - s)\,y^{a}_{u,l}(x) + s\,y^{b}_{u,l}(x),\qquad 0 \le s \le 1 .
\tag{R2.22}
\end{equation*}

`cylinder_to_airfoil` blends from the ellipse (R2.18) whose $t$ equals the
aerofoil's, so at $s = 0$ it is an ellipse, and a circle only when $t = 1$.
`interpolate_by_thickness(a, b, t^*)` sets $s = (t^\ast - t_a)/(t_b - t_a)$ in (R2.22)
and then applies (R2.20) with $t^\ast$. The same pair of steps is the airfoil morph
along the span, derived in [R3](R3_span_laws.md).

## 4. Algorithm

Building an `Airfoil` from points:

1. Take the TE as the mid-point of the first and last points, and the LE as the
   point farthest from it (Section 2, item 2).
2. Split the contour there; call the side with the larger mean $y$ the upper
   surface. A contour given in reverse order is therefore read correctly.
3. For each side, drop repeated points, classify the TE end (corner, round or
   free, Section 3.1), fit the arc-length splines (R2.1), and sample, sort and
   store them.
4. Record the extent $[x_0, x_1]$ over both sides for (R2.2).

`from_selig` and `from_lednicer` parse `.dat` text into points in this form;
Lednicer files are reordered into Selig order first.

## 5. Implementation

| Result | Function (`implicit_geometry_kernel` unless stated) |
|---|---|
| (R2.1), (R2.2), Section 4 | `Airfoil.__init__`, `Airfoil._monotone`, `Airfoil.surface` |
| (R2.3)–(R2.9) | `Airfoil.naca4`; `sdf.shapes.naca.naca_four_digit_vertices`, `naca_four_digit_thickness` |
| (R2.10)–(R2.12) | `Airfoil.naca5`; `sdf.shapes.naca.naca_five_digit_vertices` (BG-055) |
| (R2.13) | `Airfoil.naca6`; `sdf.shapes.naca.naca_six_series_vertices` (BG-056) |
| (R2.14)–(R2.16) | `Airfoil.cst` |
| (R2.18) | `Airfoil.ellipse` |
| (R2.19) | `Airfoil.thickness`, `Airfoil.max_thickness` |
| (R2.20) | `Airfoil.with_thickness` |
| (R2.21) | `Airfoil.with_te_thickness` |
| (R2.22) | `Airfoil.blend`, `Airfoil.cylinder_to_airfoil`, `Airfoil.interpolate_by_thickness` |
| Readers | `Airfoil.from_selig`, `Airfoil.from_lednicer` |

## 6. Limits

1. **The leading edge in $y(x)$ form.** Assumption A1 fails within about $10^{-4}$
   of the LE, where the surface is vertical and, for a cambered NACA section, runs
   ahead of $x = 0$. There the sorted samples are not single-valued and
   `surface` can be off by up to $2.5\times10^{-3}$ of chord (NACA 2412, 81 points a
   side, against 2001 points). Beyond $x = 0.01$ the same comparison agrees to
   $3\times10^{-6}$. Section geometry built from the contour itself, rather than
   from `surface` queries at the nose, is not affected.
2. **5-digit mean lines** are wrong (BG-055), and reflexed sections ($Q = 1$) cannot
   be requested from a code: the reflex branch of `naca_five_digit_vertices` is
   reached only by its `reflex` argument.
3. **6-series sections** are scaled tabulated sections, not the 6-series
   construction, and their camber follows the thickness (BG-056). Only families 63
   and 64 are available.
4. **Vertical, not normal, thickness** in the shape operators (R2.19)–(R2.22). The
   difference is second order in the camber slope (Section 3.6).
5. **CST** has the fixed class exponents $(\tfrac12, 1)$ of (R2.14); sharp-nosed or
   other class functions are not available.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R2.6) | $-(a_0+a_1+a_2+a_3) = -0.1036$; open TE $2h(1) = 0.021\,t$ | [blade example 17](../../examples/blade/17_aerofoil_construction/README.md), check 1 |
| (R2.7), (R2.8) | NACA 2412: maximum camber 0.02000 at $x = 0.40$, `max_thickness` 0.12001 | example 17, check 2 |
| (R2.9) | Surfaces are the mean line offset by $h$ along its normal, to $1.1\times10^{-16}$ | example 17, check 3 |
| (R2.16) | `Airfoil.cst` against (R2.16) to $1.5\times10^{-5}$ of chord; TE gap exact; equal weights symmetric | example 17, check 4 |
| (R2.12) | Join $r$ against the tabulated values to $3\times10^{-4}$ | example 17, check 5 |
| (R2.20), (R2.21) | Scaled $t/c$ to $10^{-3}$ relative; flatback gap; symmetric CST | `tests/test_airfoil_families.py` |
| Readers, Section 4 | Selig and Lednicer files read back to the section they were written from | `tests/test_airfoil_families.py` |

```{figure} ../../examples/blade/17_aerofoil_construction/figures/aerofoil_construction.png
:alt: Top, a NACA 4418 built from its camber line and thickness; bottom, the families the kernel builds and the NACA 23012 mean line against the kernel's
:width: 100%

**Figure R2.1.** Top: NACA 4418 from (R2.5), (R2.7) and (R2.9); the grey segments
are the half-thickness laid off along the mean-line normal. Bottom: one section of
each family, chord-normalised, and the NACA 23012 mean line from (R2.11) (dashed)
against the kernel's (dotted, BG-055). Drawn by
`examples/blade/17_aerofoil_construction/aerofoil_construction.py`.
```

These checks move to [V1](../v/index.md) when it is written.

## 8. References

{cite:t}`abbott1959` for the NACA families, {cite:t}`kulfan2008` for CST, and
{cite:t}`selig2026` for the 6-series coordinates and the `.dat` formats. Full
entries are in [R0](R0_notation.md#bibliography).
