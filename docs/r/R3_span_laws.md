# R3 · Planform and span laws

## 1. Problem

A blade is given at a few spanwise stations $z_0 < z_1 < \dots < z_N$: at each, a
chord, a twist, a pitch-axis fraction, an offset and an aerofoil. Between stations,
and for every member dimension that varies along the span (thicknesses, web and cap
positions), the kernel needs a value at any $z$. This document derives the three
interpolants it uses (constant-ended linear tables, the shape-preserving cubic
PCHIP, and the airfoil morph), the test that decides whether a set of stations is
dense enough for PCHIP, and the composite laws built from them. The frames in which
the interpolated section is then placed are in [R1](R1_frames.md), and the aerofoil
operators in [R2](R2_aerofoils.md).

## 2. Definitions and assumptions

1. **Knots** $(z_k, v_k)$, $k = 0,\dots,N$, with $z_k$ strictly increasing: the
   stations, or the rows of a `(z, values)` table.
2. **Bay** $k$: the interval $[z_k, z_{k+1}]$, of length $h_k = z_{k+1} - z_k$,
   with secant slope $\delta_k = (v_{k+1} - v_k)/h_k$.
3. **Span law** $v(z)$: a scalar function of the span (`SpanLaw`). Outside
   $[z_0, z_N]$ every law holds its end value.
4. **Planform scalars**: chord $c(z)$, twist $\tau(z)$, pitch-axis fraction $f(z)$
   and offset $(o_x, o_y)(z)$ of `BladeGeometry`.
5. **Assumption A1.** The stations bound the blade: `BladeGeometry` refuses
   $z < z_0$ or $z > z_N$, so the hold-end rule of item 3 applies only to member
   laws, never to the planform.

## 3. Derivation

### 3.1 Bay and weight

The bay that contains $z$ and the position within it are

\begin{equation*}
k = \min\big(\max\{\,j : z_j \le z\,\},\ N - 1\big),
\qquad
w = \frac{z - z_k}{h_k} \in [0, 1] .
\tag{R3.1}
\end{equation*}

Every interpolant below is a function of $k$ and $w$ alone.

### 3.2 Linear tables

The linear law (`TableLaw`, and a `(z, values)` table with `span_interp="linear"`)
is

\begin{equation*}
\begin{gathered}
v(z) = (1 - w)\,v_k + w\,v_{k+1},\\
v(z) = v_0 \ (z < z_0),\qquad v(z) = v_N \ (z > z_N).
\end{gathered}
\tag{R3.2}
\end{equation*}

It reproduces the knots, never leaves $[\min(v_k, v_{k+1}), \max(v_k, v_{k+1})]$
within a bay, and keeps a run of equal knots flat. Its slope is $\delta_k$ in bay
$k$, so it jumps by $\delta_k - \delta_{k-1}$ at each interior knot. Neighbouring
equal knots give the ply-drop shape: constant, a linear ramp, constant.

### 3.3 Cubic Hermite form

PCHIP is a cubic on each bay fixed by the two end values and two end slopes
$d_k, d_{k+1}$. With the Hermite basis

\begin{equation*}
\begin{aligned}
H_{00} &= 2w^3 - 3w^2 + 1, & H_{10} &= w^3 - 2w^2 + w,\\
H_{01} &= -2w^3 + 3w^2, & H_{11} &= w^3 - w^2,
\end{aligned}
\tag{R3.3}
\end{equation*}

the law on bay $k$ is

\begin{equation*}
v(z) = H_{00}(w)\,v_k + H_{10}(w)\,h_k d_k + H_{01}(w)\,v_{k+1} + H_{11}(w)\,h_k d_{k+1} .
\tag{R3.4}
\end{equation*}

Since $\mathrm{d}w/\mathrm{d}z = 1/h_k$, differentiating (R3.3) term by term gives
the slope and curvature that `PchipLaw.derivative` returns:

\begin{equation*}
\begin{aligned}
v'(z) = \frac{1}{h_k}\Big[&(6w^2 - 6w)(v_k - v_{k+1})\\
&+ (3w^2 - 4w + 1)\,h_k d_k + (3w^2 - 2w)\,h_k d_{k+1}\Big],
\end{aligned}
\tag{R3.5}
\end{equation*}

\begin{equation*}
v''(z) = \frac{(12w - 6)(v_k - v_{k+1}) + (6w - 4)\,h_k d_k + (6w - 2)\,h_k d_{k+1}}{h_k^2} .
\tag{R3.6}
\end{equation*}

At $w = 0$, (R3.4) and (R3.5) give $v_k$ and $d_k$; at $w = 1$, $v_{k+1}$ and
$d_{k+1}$. The value and slope at a knot are the same from both bays, so the law is
$C^1$. Evaluating (R3.6) at $w = 1$ of bay $k-1$ and $w = 0$ of bay $k$ gives two
different curvatures in general, so the law is not $C^2$: its curvature jumps at
the knots.

### 3.4 Interior slopes

The slopes make the cubic *shape-preserving*. Where the data turn (an extremum) or
stop (a plateau), the slope is zero; elsewhere it is a weighted harmonic mean of the
two neighbouring secants {cite:p}`fritsch1984`:

\begin{equation*}
d_k =
\begin{cases}
\dfrac{w_1 + w_2}{w_1/\delta_{k-1} + w_2/\delta_k}, & \delta_{k-1}\,\delta_k > 0,\\[2ex]
0, & \text{otherwise},
\end{cases}
\tag{R3.7}
\end{equation*}

with weights $w_1 = 2h_k + h_{k-1}$ and $w_2 = h_k + 2h_{k-1}$.

**Bound.** Take $\delta_{k-1}, \delta_k > 0$ (the negative case is the mirror image).
Dropping one positive term from the denominator of (R3.7) can only increase the
quotient, so $d_k < \delta_{k-1}(1 + w_2/w_1)$ and $d_k < \delta_k(1 + w_1/w_2)$.
Both weight ratios lie strictly between $\tfrac12$ and 2 for any positive $h$, hence

\begin{equation*}
0 < d_k < 3\min(\delta_{k-1}, \delta_k) .
\tag{R3.8}
\end{equation*}

### 3.5 End slopes

At the ends there is one secant on each side to use. The slope is the one-sided
three-point estimate, limited so the end bay keeps the sign of its secant and does
not exceed three times it {cite:p}`moler2004`:

\begin{equation*}
\begin{gathered}
\hat d_0 = \frac{(2h_0 + h_1)\,\delta_0 - h_0\,\delta_1}{h_0 + h_1},
\\
d_0 =
\begin{cases}
0, & \operatorname{sign}\hat d_0 \ne \operatorname{sign}\delta_0,\\
3\delta_0, & \operatorname{sign}\delta_0 \ne \operatorname{sign}\delta_1 \text{ and } |\hat d_0| > 3|\delta_0|,\\
\hat d_0, & \text{otherwise},
\end{cases}
\end{gathered}
\tag{R3.9}
\end{equation*}

and symmetrically at $z_N$. With two knots $d_0 = d_1 = \delta_0$. When $\delta_0$
and $\delta_1$ have the same sign, $\hat d_0 < (2h_0 + h_1)\delta_0/(h_0+h_1) < 2\delta_0$
for $\delta_0 > 0$, so the cap of (R3.9) is needed only when they differ. Either way
$0 \le d_0/\delta_0 \le 3$.

**Limiting case: straight data.** If every $\delta_k = \delta$, (R3.7) gives
$d_k = \delta(w_1+w_2)/(w_1+w_2) = \delta$ and (R3.9) gives
$\hat d_0 = \delta$. Substituting $d_k = d_{k+1} = \delta$ in (R3.4), and
$v_{k+1} = v_k + h_k\delta$, the cubic and quadratic terms cancel:

\begin{equation*}
v(z) = v_k + w\,h_k\,\delta = v_k + (z - z_k)\,\delta .
\tag{R3.10}
\end{equation*}

PCHIP through collinear knots is the line itself, at any spacing.

### 3.6 Monotonicity and no overshoot

On bay $k$ let $\alpha = d_k/\delta_k$ and $\beta = d_{k+1}/\delta_k$. Fritsch and
Carlson show that the cubic (R3.4) is monotone on the bay if $\alpha, \beta \ge 0$
and either $\alpha + \beta \le 2$ or $(\alpha, \beta)$ lies inside the ellipse
$(\alpha-1)^2 + (\alpha-1)(\beta-1) + (\beta-1)^2 = 3(\alpha + \beta - 2)$, and
that this region contains the square $[0, 3]^2$, whose corners $(0,3)$, $(3,0)$ and
$(3,3)$ lie on the ellipse {cite:p}`fritsch1980`. (R3.8) at knots $k$ and $k+1$, and
(R3.9) at the ends, put both ratios in $[0, 3]$, so

\begin{equation*}
\begin{gathered}
\delta_k \ne 0:\ \ v \text{ monotone on bay } k,\ \
\min(v_k, v_{k+1}) \le v \le \max(v_k, v_{k+1});
\\
\delta_k = 0:\ \ v \equiv v_k .
\end{gathered}
\tag{R3.11}
\end{equation*}

The second case follows because $\delta_k = 0$ sets $d_k = d_{k+1} = 0$ by (R3.7),
and (R3.4) with equal end values and zero end slopes is the constant. No bay
overshoots its knots, and plateaus stay flat. The law is
`PchipLaw`; its slopes agree with SciPy's `PchipInterpolator`, which implements the
same rules.

### 3.7 Whether the knots are dense enough

A linear loft between sparse stations is wrong only between them; a cubic one can
also invent curvature that no station supports. Decision
[0014](../decisions/0014-smooth-authored-span-laws.md) therefore admits a PCHIP law
only if each interior knot can be predicted from the others. Let $v^{(-i)}$ be the
PCHIP law through every knot except $i$. The leave-one-out error and its tolerance
are

\begin{equation*}
\begin{gathered}
e_i = \left|\,v^{(-i)}(z_i) - v_i\,\right| \le \varepsilon,\qquad 1 \le i \le N - 1,\\
\varepsilon = \varepsilon_a + \varepsilon_r\left(\max_k v_k - \min_k v_k\right),
\end{gathered}
\tag{R3.12}
\end{equation*}

with defaults $\varepsilon_a = 0$, $\varepsilon_r = 5\times10^{-3}$. End knots
are not tested. The report also lists the largest $|v - v_{\text{linear}}|$ in each
bay: the geometry PCHIP adds over a linear loft.

**At a smooth extremum.** Near a maximum, let the knots sample
$v = v_i - \tfrac12\kappa(z - z_i)^2$ at a uniform spacing $h$, with $\kappa > 0$
the curvature. Dropping knot $i$ leaves $z_i \pm h$ with equal values. The secants
next to them have opposite signs, so (R3.7) sets both slopes to zero, and the bay
between them is the constant $v_i - \tfrac12\kappa h^2$ by (R3.11). Hence

\begin{equation*}
e_i = \tfrac12\,\kappa\,h^2,
\qquad\text{and the check passes only for}\qquad
h \le \sqrt{2\varepsilon/\kappa} .
\tag{R3.13}
\end{equation*}

The flattening that keeps PCHIP from overshooting makes it a second-order predictor
at an extremum, so laws with an interior maximum (chord is the usual one) need much
denser stations near it than elsewhere. For a chord law of range 1.5 m peaking with
$\kappa = 0.083\ \mathrm{m^{-1}}$, (R3.13) asks for $h \le 0.42$ m around the peak;
away from it, 1 m spacing passes.

### 3.8 Which law each quantity uses

`BladeGeometry(span_interp=...)` builds the chord, twist, offset and pitch-axis
fraction from the stations (or from a `(z, values)` table) as follows:

* `"linear"`: (R3.2) for all of them, the loft before 0.9.6. Tabulated inputs
  (windIO, CSV) and `from_section` use it.
* `"pchip"`: (R3.4)–(R3.9) for all of them. A law that fails (R3.12) raises
  `UnderSampledSpanLaw`, or with `span_sampling="flag"` is listed by `validate()`.
* `"auto"` (default): PCHIP for each law that passes (R3.12) and linear for each that
  fails; `BladeGeometry.span_law_interp` records the choice per quantity.

A constant offset or pitch-axis fraction stays constant. Callables pass through
unchanged. The aerofoil is not a scalar law: it is morphed between the two stations
of its bay (Section 3.9), with the weight $w$ of (R3.1) whatever the scalars use.
On the `from_section` path the stations already carry the consumer's blend, and
surfaces and scalars are interpolated linearly in $w$.

### 3.9 Airfoil morph

The section at $z$ in bay $k$ blends the upper and lower surfaces of the bounding
station aerofoils $a$ and $b$ at the same chord fractions, by one of three rules.

**Linear.** (R2.22) with $s = w$.

**Smoothstep.** (R2.22) with

\begin{equation*}
s = S(w) = w^2(3 - 2w),
\tag{R3.14}
\end{equation*}

for which $S(0) = 0$, $S(1) = 1$, $S'(0) = S'(1) = 0$ and $S(\tfrac12) = \tfrac12$.
The zero end slopes make the section change slowly near both stations; it is used
between a root circle and an aerofoil, whose shapes differ most.

**Thickness.** The target thickness ratio is linear in $w$, and the surfaces are
blended and then rescaled to it. With $t_a, t_b$ the largest vertical thicknesses
(R2.19) of the parents,

\begin{equation*}
t^\ast = (1 - w)\,t_a + w\,t_b,\qquad
s = \frac{t^\ast - t_a}{t_b - t_a} = w \quad (t_a \ne t_b),
\tag{R3.15}
\end{equation*}

so the blend weight of `interpolate_by_thickness` is exactly the span weight, and
the rule is (R2.22) at $s = w$ followed by (R2.20) with

\begin{equation*}
\lambda = \frac{t^\ast}{\max_x \delta_w(x)},\qquad
\delta_w = (1 - w)\,\delta_a + w\,\delta_b .
\tag{R3.16}
\end{equation*}

Because a maximum of a sum is at most the sum of the maxima,
$\max \delta_w \le (1-w)\max\delta_a + w\max\delta_b = t^\ast$, so

\begin{equation*}
\lambda \ge 1,
\tag{R3.17}
\end{equation*}

with equality when both parents peak at the same $x$. The rescale restores the
thickness that a plain blend loses when the two thickest points are at different
chord fractions.

**Choice of rule.** With `airfoil_morph="auto"`, a parent counts as near-circular
if its thickness ratio is at least 0.6 or its name contains "ellipse". Smoothstep is
used when exactly one parent is near-circular, linear when both are, and thickness
when neither is. The other settings force one rule on every bay.

**Trailing edge.** After any rule, the TE gap is reset to the linear blend of the
parents' gaps about the blended mid-point:

\begin{equation*}
\delta(1) = (1 - w)\,\delta_a(1) + w\,\delta_b(1) .
\tag{R3.18}
\end{equation*}

### 3.10 Composite laws and thickness fields

`ScaledLaw`, `SumLaw` and `ProductLaw` evaluate $\lambda v(z) + \mu$,
$v_1(z) + v_2(z)$ and $v_1(z)\,v_2(z)$ pointwise. A ply drop adds a thickness
$\Delta t$ inside a window $[z_0, z_1]\times[s_0, s_1]$ of span and chord fraction,
ramped linearly over a taper length $\ell > 0$ on every side ($\ell = 0$ gives a
step):

\begin{equation*}
\begin{gathered}
t_{\mathrm{drop}}(z, s) = \Delta t\,\min\big(\rho(z; z_0, z_1),\ \rho(s; s_0, s_1)\big),
\\
\rho(x; a, b) = \operatorname{clip}\!\left(\min\!\left(\frac{x - a + \ell}{\ell}, \frac{b + \ell - x}{\ell}\right), 0, 1\right).
\end{gathered}
\tag{R3.19}
\end{equation*}

`GridThickness` holds a thickness table over the run coordinate $t\in[0,1]$ across
a partition and the span; each row is a span law (linear or PCHIP, checked by
(R3.12)), and rows are blended linearly in $t$.

## 4. Algorithm

`BladeGeometry.interpolate(z)`:

1. Refuse $z$ outside $[z_0, z_N]$ (A1); find $k$ and $w$ by (R3.1).
2. If the blade came from `from_section`, blend the two stored samples linearly in
   $w$ and stop.
3. Otherwise morph the two station aerofoils by Section 3.9 on the cosine grid
   (R2.4) of `n_airfoil_points` points.
4. Evaluate chord, twist and offset from their span laws (Section 3.8) and the
   pitch-axis fraction from its law.
5. Hand the result to the section build, which scales, places and twists it
   ([R1](R1_frames.md), (R1.12)).

`interpolate_many` does the same for many $z$, sampling each pair of parents once per
bay.

## 5. Implementation

| Result | Function (`implicit_geometry_kernel` unless stated) |
|---|---|
| (R3.1) | `BladeGeometry.interpolate`, `interpolate_many` (`blade.loft.loft_cut`) |
| (R3.2) | `structure.TableLaw`; `blade.loft.loft_tables._as_span_table` (linear branch) |
| (R3.3)–(R3.6) | `structure.laws._hermite`; `PchipLaw.__call__`, `PchipLaw.derivative` |
| (R3.7)–(R3.9) | `structure.laws.pchip_slopes`, `_pchip_end` |
| (R3.12) | `structure.check_pchip_sampling`, `SamplingReport`, `UnderSampledSpanLaw` |
| Section 3.8 | `BladeGeometry.__init__` (`span_interp`, `span_sampling`, `span_rtol`), `BladeGeometry._pchip_factory`, `span_law_interp` |
| (R3.14)–(R3.18) | `blade.loft.loft_tables._morph_surfaces`, `_morph_surfaces_many`, `_smoothstep`, `_near_circle` |
| (R3.19) | `PlyDrop.extra`, `PlyDrop.extra_many` |
| Section 3.10 | `structure.ScaledLaw`, `SumLaw`, `ProductLaw`, `GridThickness` |

## 6. Limits

1. **Stations near an extremum.** By (R3.13), the default tolerance refuses a
   smooth maximum unless the stations near it are closer than
   $\sqrt{2\varepsilon/\kappa}$. Under `span_interp="auto"` most chord laws given at
   ordinary station spacing therefore stay linear. This is the intended effect of
   decision 0014: add stations near the maximum, or accept the linear loft.
2. **The aerofoil is not smooth across stations.** The morph is linear (or
   smoothstep) in $w$ within each bay, so surfaces are only $C^0$ in $z$ at a
   station, even when the planform scalars are PCHIP.
3. **Trailing-edge point under smoothstep.** (R3.18) resets the last point with $w$,
   while the rest of the section uses $S(w)$. Between a closed circle and a flatback
   of gap $\Delta$, the last point therefore sits $(w - S(w))\Delta$ away from where
   the blend puts its neighbours: $0.09\Delta$ at $w = 0.25$. Logged as
   [BG-052](../blade/bug_log.md#bg-052) (open).
4. **The near-circle test** reads the aerofoil's name. A thick section named
   otherwise is classified by its thickness ratio alone.
5. **Member laws hold their end values** outside their knots (item 3, Section 2);
   a law given over part of the span is constant beyond it, not zero.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R3.4)–(R3.9) | `PchipLaw` equals SciPy's `PchipInterpolator` in value and slope to $7.5\times10^{-14}$ over 300 random knot sets | [blade example 18](../../examples/blade/18_pchip_span_laws/README.md), check 1 |
| (R3.10) | Collinear, unequally spaced knots: distance from the line $1.7\times10^{-16}$ | example 18, check 2 |
| (R3.8), (R3.11) | 300 monotone knot sets with plateaus: largest slope ratio 2.86, no decrease beyond $1.8\times10^{-15}$, no bay overshoot beyond $8.9\times10^{-16}$ | example 18, check 3 |
| (R3.6) | Slope continuous at every knot to $10^{-9}$; curvature jumps up to $0.86\ \mathrm{m^{-1}}$ | example 18, check 4 |
| (R3.12), (R3.13) | Leave-one-out error at a smooth maximum equals $\tfrac12\kappa h^2$ to 1.5 % for 21–81 stations; 41 stations refused, 61 accepted | example 18, check 5 |
| (R3.15)–(R3.18) | $s = w$ exactly; $\lambda \ge 1$; peak $t/c$ equals $t^\ast$ to $3\times10^{-17}$; TE gap linear in $w$ for all three rules | example 18, check 6 |
| Section 3.8 | `auto`, `pchip` and `linear` choices; refusals; golden numbers | `tests/test_pchip_span_laws.py`, `tests/golden/` |
| (R3.18) | A blunt TE gap survives the loft | `tests/test_loft.py::test_interpolate_preserves_blunt_te_gap` |

```{figure} ../../examples/blade/18_pchip_span_laws/figures/pchip_span_laws.png
:alt: A chord law through ten stations, linear, PCHIP and natural spline; its slope; and a circle morphed into a flatback aerofoil
:width: 100%

**Figure R3.1.** Top: a chord law through ten stations by (R3.2), by PCHIP
(R3.4)–(R3.9) and by a natural cubic spline, which overshoots the plateau and the
root. Bottom left: the slope, continuous for PCHIP and stepped for the linear table.
Bottom right: the smoothstep morph (R3.14) from a circle to a NACA 2418 flatback.
Drawn by `examples/blade/18_pchip_span_laws/pchip_span_laws.py`.
```

These checks move to [V1](../v/index.md) when it is written.

## 8. References

{cite:t}`fritsch1980` for the monotonicity region, {cite:t}`fritsch1984` for the
harmonic-mean slopes, and {cite:t}`moler2004` for the end slopes. Full entries are
in [R0](R0_notation.md#bibliography).
