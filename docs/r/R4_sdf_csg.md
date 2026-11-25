# R4 · Signed distance fields and constructive solid geometry

## 1. Problem

Every region the kernel builds (an aerofoil, a skin wall, a spar cap, a web, a
cavity) is held as a scalar field $\phi(\mathbf{x})$ whose sign says whether a point
$\mathbf{x}$ in the section plane is inside, and whose magnitude says, exactly or as
a bound, how far the point is from the region's boundary. This document derives the
fields of the primitive shapes, the operations that combine them (union,
intersection, difference, offset, shell, slab, smooth blends, transforms), what each
operation preserves of the distance, the integer feature ids carried alongside, the
two-and-a-half-dimensional field used along the span, the Newton solves that put
points on a zero set, the bounding boxes used for culling, and the extraction of
zero loops. The section integrals over these regions are in R7, and the layouts
built from them in [R5](index.md).

## 2. Definitions and assumptions

1. **Region** $\Omega \subset \mathbb{R}^2$: a closed set with boundary
   $\partial\Omega$ and interior $\Omega^\circ$. Its **complement** is
   $\Omega^c = \overline{\mathbb{R}^2 \setminus \Omega}$.
2. **Distance** to a set $S$: $\operatorname{dist}(\mathbf{x}, S) =
   \inf_{\mathbf{y}\in S}\lVert\mathbf{x} - \mathbf{y}\rVert$.
3. **Signed distance** $d_\Omega$ of (R4.1): negative inside, positive outside
   (the kernel's solid convention).
4. **Field** $\phi$: any function the kernel evaluates for a region. It is
   **exact** at $\mathbf{x}$ when $\phi(\mathbf{x}) = d_\Omega(\mathbf{x})$, and a
   **bound** when $\operatorname{sign}\phi = \operatorname{sign} d_\Omega$ and
   $|\phi| \le |d_\Omega|$.
5. **Medial axis** $\mathcal{M}_\Omega$: the points with more than one nearest
   boundary point. $d_\Omega$ is differentiable off $\mathcal{M}_\Omega \cup
   \partial\Omega$.
6. **Assumption A1.** Fields are evaluated in the section plane, in metres, at
   fixed $z$. Section 3.8 states what holds across $z$.
7. **Assumption A2.** Regions are regular: $\Omega = \overline{\Omega^\circ}$. A
   boolean operation that leaves a set of zero area (two boxes that touch along an
   edge) is outside the derivation.

## 3. Derivation

### 3.1 The signed distance

The signed distance to $\Omega$ is

\begin{equation*}
d_\Omega(\mathbf{x}) =
\begin{cases}
-\operatorname{dist}(\mathbf{x}, \partial\Omega), & \mathbf{x}\in\Omega,\\
\phantom{-}\operatorname{dist}(\mathbf{x}, \partial\Omega), & \mathbf{x}\notin\Omega .
\end{cases}
\tag{R4.1}
\end{equation*}

Outside $\Omega$, the nearest point of $\Omega$ lies on its boundary, so
$d_\Omega = \operatorname{dist}(\mathbf{x}, \Omega)$ there; inside, $d_\Omega =
-\operatorname{dist}(\mathbf{x}, \Omega^c)$. The triangle inequality applied to the
nearest boundary point of either of two points gives

\begin{equation*}
|d_\Omega(\mathbf{x}) - d_\Omega(\mathbf{y})| \le \lVert\mathbf{x} - \mathbf{y}\rVert,
\qquad
\lVert\nabla d_\Omega\rVert = 1 \ \text{off}\ \mathcal{M}_\Omega\cup\partial\Omega .
\tag{R4.2}
\end{equation*}

The second statement (the eikonal equation) holds because $d_\Omega$ grows at unit
rate along the segment from the nearest boundary point and the gradient is that
segment's direction {cite:p}`osher2003`.

### 3.2 Primitive fields

Each primitive is written in its own frame and evaluated exactly.

**Disk** of centre $\mathbf{c}$ and radius $r$ (`sdf_circle`, `Circle`): the
nearest boundary point lies on the ray from $\mathbf{c}$ through $\mathbf{x}$, so

\begin{equation*}
\phi = \lVert\mathbf{x} - \mathbf{c}\rVert - r .
\tag{R4.3}
\end{equation*}

**Box** of centre $\mathbf{c}$ and half-sizes $\mathbf{b} = (b_x, b_y)$
(`sdf_box`, `Box`). By the box's two mirror symmetries it is enough to take the
first quadrant, $\mathbf{q} = |\mathbf{x} - \mathbf{c}| - \mathbf{b}$
component by component. With $q_x > 0$ and $q_y > 0$ the nearest point is the
corner, at distance $\lVert\mathbf{q}\rVert$; with only $q_x > 0$ it is on the
vertical side, at $q_x$; inside, it is on the nearer side, at
$\min(-q_x, -q_y)$. The three cases combine into

\begin{equation*}
\phi = \lVert\max(\mathbf{q}, \mathbf{0})\rVert + \min\big(\max(q_x, q_y), 0\big) .
\tag{R4.4}
\end{equation*}

**Rounded box** of corner radius $r \le \min(b_x, b_y)$ (`sdf_rounded_box`,
`RoundedBox`). It is the box of half-sizes $\mathbf{b} - r$ dilated by $r$. That box
is convex, and the dilation of a convex set by $r$ subtracts $r$ from its signed
distance everywhere (Section 3.4), so

\begin{equation*}
\phi = \phi_{\mathrm{box}}(\mathbf{x};\ \mathbf{c},\ \mathbf{b} - r) - r .
\tag{R4.5}
\end{equation*}

The kernel clips $r$ to $\min(b_x, b_y)$ before using (R4.5).

**Oriented box** turned by $\alpha$ anticlockwise (`sdf_oriented_box`,
`OrientedBox`). A rotation keeps distances, so the field is (R4.4) evaluated at the
point turned back into the box frame (Section 3.6):

\begin{equation*}
\phi(\mathbf{x}) = \phi_{\mathrm{box}}\big(\mathbf{Q}(-\alpha)(\mathbf{x} - \mathbf{c}) + \mathbf{c}\big),
\qquad
\mathbf{Q}(\alpha) = \begin{pmatrix}\cos\alpha & -\sin\alpha\\ \sin\alpha & \cos\alpha\end{pmatrix}.
\tag{R4.6}
\end{equation*}

**Segment and capsule** from $\mathbf{a}$ to $\mathbf{b}$ (`sdf_segment`,
`sdf_capsule`, `Capsule`). The nearest point of the segment to $\mathbf{x}$ is the
orthogonal projection, clamped to the segment's ends:

\begin{equation*}
\begin{gathered}
t^\ast = \operatorname{clip}\!\left(\frac{(\mathbf{x} - \mathbf{a})\cdot(\mathbf{b} - \mathbf{a})}{\lVert\mathbf{b} - \mathbf{a}\rVert^2},\ 0,\ 1\right),
\qquad
u(\mathbf{x}) = \big\lVert\mathbf{x} - \mathbf{a} - t^\ast(\mathbf{b} - \mathbf{a})\big\rVert,
\\
\phi_{\mathrm{capsule}} = u - r .
\end{gathered}
\tag{R4.7}
\end{equation*}

$u$ is unsigned. The capsule is the segment dilated by $r$, a convex set, so
$u - r$ is exact on both sides.

**Half-plane** $\{\mathbf{n}\cdot\mathbf{x} + d_0 \le 0\}$ with outward normal
$\mathbf{n}$ (`sdf_half_plane`, `HalfPlane`): the distance to a line is the
normalised residual of its equation,

\begin{equation*}
\phi = \frac{\mathbf{n}\cdot\mathbf{x} + d_0}{\lVert\mathbf{n}\rVert} .
\tag{R4.8}
\end{equation*}

**Polygon** with vertices $\mathbf{v}_0,\dots,\mathbf{v}_{n-1}$ and edges
$e_i = [\mathbf{v}_i, \mathbf{v}_{i+1}]$, the last wrapping to $\mathbf{v}_0$
(`sdf_polygon`, `Polygon`). The boundary is the union of the edges, so the
unsigned distance is the least of the edge distances $u_i$ of (R4.7). The sign comes
from the winding number of the boundary about $\mathbf{x}$, counted on the ray from
$\mathbf{x}$ in the $+x$ direction {cite:p}`hormann2001`: an edge crossing the ray
upward with $\mathbf{x}$ to its left adds one, an edge crossing downward with
$\mathbf{x}$ to its right subtracts one. With
$c_i = (\mathbf{v}_{i+1} - \mathbf{v}_i)\times(\mathbf{x} - \mathbf{v}_i)$, the
cross product's $z$ component,

\begin{equation*}
\begin{gathered}
w(\mathbf{x}) = \sum_{i} \Big( [\,v_{i,y} \le y < v_{i+1,y},\ c_i > 0\,] - [\,v_{i+1,y} \le y < v_{i,y},\ c_i < 0\,] \Big),
\\
\phi = \begin{cases} -\min_i u_i, & w \ne 0,\\ \phantom{-}\min_i u_i, & w = 0, \end{cases}
\end{gathered}
\tag{R4.9}
\end{equation*}

where $[\cdot]$ is one when its conditions hold and zero otherwise. The half-open
test on $y$ counts a ray through a vertex once. For a simple polygon $w$ is $\pm 1$
inside and zero outside, so (R4.9) is exact; for a self-intersecting one it is the
non-zero rule. The polygon is the section's geometry: an aerofoil given as a
polyline is exact with respect to those vertices, and its distance to the analytic
aerofoil is the polygonisation error of [R2](R2_aerofoils.md).

**Ellipse** with semi-axes $a \ge b$ along $x$ and $y$ (`sdf_ellipse`). By symmetry
take $\mathbf{x} = (x, y)$ with $x, y \ge 0$. The nearest point $\mathbf{p}$ of the
ellipse makes $\mathbf{x} - \mathbf{p}$ parallel to the normal there; a Lagrange
multiplier $s$ on the constraint $(p_x/a)^2 + (p_y/b)^2 = 1$ gives
$\mathbf{p} = \big(a^2 x/(a^2 + s),\ b^2 y/(b^2 + s)\big)$, and the constraint becomes
one equation for $s$ {cite:p}`eberly2013`:

\begin{equation*}
\begin{gathered}
F(s) = \left(\frac{a x}{a^2 + s}\right)^{2} + \left(\frac{b y}{b^2 + s}\right)^{2} - 1 = 0,
\qquad s > -b^2,
\\
\phi = \pm\lVert\mathbf{x} - \mathbf{p}(s)\rVert .
\end{gathered}
\tag{R4.10}
\end{equation*}

$F$ decreases monotonically on $s > -b^2$, so the root is unique and bracketed by
$[\,b y - b^2,\ \lVert(a x, b y)\rVert - b^2\,]$ for $x, y > 0$; the axis cases
$x = 0$ or $y = 0$ are closed-form: on the minor axis the nearest point is
$(0, b)$; on the major axis it is $\big(a^2x/(a^2 - b^2),\ b\sqrt{1 - p_x^2/a^2}\big)$
for $x < (a^2 - b^2)/a$ and $(a, 0)$ beyond. The sign is that of
$(x/a)^2 + (y/b)^2 - 1$. The kernel bisects $F$ on that bracket to machine
precision (fixed in SD-001; it used to run a fixed-point iteration that did not
converge).

### 3.3 Union, intersection and difference

**Union** (`union`, `Union`). Outside $A \cup B$, the nearest point of the union is
the nearer of the nearest points of $A$ and of $B$:
$\operatorname{dist}(\mathbf{x}, A\cup B) = \min\big(\operatorname{dist}(\mathbf{x}, A),
\operatorname{dist}(\mathbf{x}, B)\big)$. Inside, take $\mathbf{x} \in A$: the
boundary of the union lies in $A^c$ (where $A$'s interior ends), so its distance is
at least $|d_A|$, and likewise $|d_B|$ when $\mathbf{x}\in B$. Hence

\begin{equation*}
\phi_{A\cup B} = \min(\phi_A, \phi_B):
\quad \text{exact outside } A\cup B,
\quad \text{a bound inside}.
\tag{R4.11}
\end{equation*}

The bound is strict where the two regions overlap: a point inside both sees the
nearer of two boundaries that are partly buried in the other region.

**Intersection** (`intersect`, `Intersect`). Since $A\cap B = (A^c\cup B^c)^c$ and
$d_{\Omega^c} = -d_\Omega$, negating (R4.11) applied to the complements gives

\begin{equation*}
\phi_{A\cap B} = \max(\phi_A, \phi_B):
\quad \text{exact inside } A\cap B,
\quad \text{a bound outside}.
\tag{R4.12}
\end{equation*}

**Difference** (`subtract`, `Subtract`). $A\setminus B = A\cap B^c$, so by
(R4.12)

\begin{equation*}
\phi_{A\setminus B} = \max(\phi_A, -\phi_B):
\quad \text{exact inside } A\setminus B,
\quad \text{a bound outside}.
\tag{R4.13}
\end{equation*}

**Trees.** The maximum and minimum of 1-Lipschitz functions are 1-Lipschitz, and so
is a negation; each of (R4.11)–(R4.13) gives the correct sign and never
overestimates $|d|$ when its arguments do not. By induction over a tree of these
operations and the offsets of Section 3.4,

\begin{equation*}
|\phi(\mathbf{x}) - \phi(\mathbf{y})| \le \lVert\mathbf{x} - \mathbf{y}\rVert,
\qquad
\operatorname{sign}\phi = \operatorname{sign} d_\Omega,
\qquad
|\phi| \le |d_\Omega| .
\tag{R4.14}
\end{equation*}

This is what culling relies on: a cell of half-diagonal $\rho$ whose centre has
$|\phi| > \rho$ contains no boundary, whether $\phi$ is exact there or a bound.

### 3.4 Offsets, shells and slabs

**Dilation** by $\delta > 0$ (`offset`, `Offset`). The dilated set is the Minkowski
sum $\Omega_\delta = \Omega\oplus B_\delta$, all points within $\delta$ of $\Omega$.
A point outside it is at distance $\operatorname{dist}(\mathbf{x},\Omega) - \delta$
from it: the nearest point of $\Omega_\delta$ lies $\delta$ along the segment to the
nearest point of $\Omega$. Hence

\begin{equation*}
\phi_{\Omega_\delta} = \phi_\Omega - \delta,
\qquad \text{exact outside } \Omega_\delta \text{ when } \phi_\Omega \text{ is exact outside } \Omega .
\tag{R4.15}
\end{equation*}

For a convex $\Omega$ the same subtraction is exact inside $\Omega_\delta$ as well:
the dilation of a convex set is convex and its nearest boundary point lies on the
same normal. This is the step used for the rounded box (R4.5) and the capsule
(R4.7).

**Erosion** by $\delta > 0$ (`Offset` with a negative amount) is the dilation of
the complement, $\Omega_{-\delta} = (\Omega^c\oplus B_\delta)^c$. Applying (R4.15)
to $\Omega^c$, with $d_{\Omega^c} = -d_\Omega$,

\begin{equation*}
\phi_{\Omega_{-\delta}} = \phi_\Omega + \delta,
\qquad \text{exact inside } \Omega_{-\delta} \text{ when } \phi_\Omega \text{ is exact inside } \Omega .
\tag{R4.16}
\end{equation*}

**Shell** of thickness $t$ about $\partial\Omega$ (`shell`, `Shell`): the band
$\{|d_\Omega| \le t/2\}$, the dilation by $t/2$ of the boundary curve, whose
unsigned distance is $|d_\Omega|$. By (R4.15),

\begin{equation*}
\phi_{\mathrm{shell}} = |\phi_\Omega| - \tfrac{t}{2},
\qquad \text{exact outside the band}.
\tag{R4.17}
\end{equation*}

**Slab** of depth $D$ inside the boundary (`Slab`) and **mid-surface**
(`FieldOffsetNode`, `MidSurfaceNode`). The wall between the outer mould line and a
depth $D$ is $\Omega\cap(\Omega_{-D})^c$; its field by (R4.12), and the field of
the cavity $\Omega_{-D}$ by (R4.16), are

\begin{equation*}
\begin{gathered}
\phi_{\mathrm{slab}} = \max\big(\phi_\Omega,\ -\phi_\Omega - D\big),
\qquad
\phi_{\mathrm{cav}} = \phi_\Omega + D,
\\
\phi_{\mathrm{mid}} = \tfrac12\big(\phi_\Omega + \phi_{\mathrm{cav}}\big) = \phi_\Omega + \tfrac{D}{2},
\qquad
D = \phi_{\mathrm{cav}} - \phi_\Omega .
\end{gathered}
\tag{R4.18}
\end{equation*}

The zero set of $\phi_{\mathrm{mid}}$ is the wall's mid-surface, at depth $D/2$, and
`MidSurfaceNode.thickness` returns the last identity. Inside the wall, (R4.12) is exact only
where both of its arguments are; $\phi_\Omega$ is, but $\phi_\Omega + D$ is exact only
inside the cavity (R4.16). The slab field is therefore exact in the wall where the
inner face is reached along the outer face's normal (depth $D$ less than the distance
to $\Omega$'s medial axis), and a bound elsewhere. Where the thickness is a
field $D(\mathbf{x})$ rather than a constant, $\phi_\Omega + D(\mathbf{x})$ is
$(1 + \lVert\nabla D\rVert)$-Lipschitz: still the right zero set, but no longer a
distance.

### 3.5 Smooth blends

The polynomial smooth minimum {cite:p}`quilez2013` of $a = \phi_A$ and
$b = \phi_B$ with blend width $k$ (`smooth_union`) is

\begin{equation*}
\begin{gathered}
h = \operatorname{clip}\!\left(\tfrac12 + \frac{b - a}{2k},\ 0,\ 1\right),
\\
\operatorname{smin}_k(a, b) = h\,a + (1 - h)\,b - k\,h(1 - h) .
\end{gathered}
\tag{R4.19}
\end{equation*}

For $|a - b| \ge k$, $h$ is 0 or 1 and (R4.19) is $\min(a, b)$. Inside the band,
write $u = b - a$ so that $\operatorname{smin}_k = b - h u - k h (1 - h)$ with
$\mathrm{d}h/\mathrm{d}u = 1/(2k)$. Differentiating, the $u$ terms cancel and

\begin{equation*}
\begin{gathered}
\frac{\partial\,\operatorname{smin}_k}{\partial a} = h,
\qquad
\frac{\partial\,\operatorname{smin}_k}{\partial b} = 1 - h,
\\
\nabla\operatorname{smin}_k = h\,\nabla\phi_A + (1 - h)\,\nabla\phi_B .
\end{gathered}
\tag{R4.20}
\end{equation*}

Three consequences follow. The partial derivatives are continuous at $|u| = k$
(they reach 1 and 0, the derivatives of $\min$), so the blend is $C^1$; the second
derivative jumps from $-1/(2k)$ to zero there, so it is not $C^2$. The gradient is
a convex combination of two vectors of length at most one, so the blend is
1-Lipschitz. And $\min - \operatorname{smin}_k = k h(1 - h) - |u|\min(h, 1-h)$ peaks
at $u = 0$, where it is $k/4$: the blend never exceeds $\min$, so the blended region
contains $A\cup B$ and adds a fillet whose depth along the bisector is $k/4$ in
field units.

Smooth intersection and difference (`smooth_intersect`, `smooth_subtract`) follow by
the complement rule, as (R4.12) and (R4.13) follow from (R4.11):

\begin{equation*}
\operatorname{smax}_k(a, b) = -\operatorname{smin}_k(-a, -b),
\qquad
\phi_{A\setminus_k B} = \operatorname{smax}_k(\phi_A, -\phi_B) .
\tag{R4.21}
\end{equation*}

The linear `blend` $(1 - t)\phi_A + t\phi_B$ is a field between two shapes; it is
1-Lipschitz but not a distance, and its zero set is not the morph of
[R3](R3_span_laws.md), which acts on aerofoil coordinates.

### 3.6 Transforms

A transform moves the geometry by evaluating the field at the point carried back
into the shape's frame (pull-back). For a rigid motion $T(\mathbf{x}) =
\mathbf{Q}\mathbf{x} + \mathbf{t}$ and a uniform scaling by $\sigma > 0$ about
$\mathbf{c}$, distances are kept, or scaled by $\sigma$:

\begin{equation*}
\phi_{T\Omega}(\mathbf{x}) = \phi_\Omega\big(T^{-1}\mathbf{x}\big),
\qquad
\phi_{\sigma\Omega}(\mathbf{x}) = \sigma\,\phi_\Omega\!\left(\frac{\mathbf{x} - \mathbf{c}}{\sigma} + \mathbf{c}\right).
\tag{R4.22}
\end{equation*}

`rotate_field`, `translate_field`, `mirror_field_x`, `mirror_field_y`, `Rotate` and
`scale_field` with equal factors are (R4.22). `SDFFrame` applies its operations to
the geometry in the order they are added, so it pulls the query point back through
them in reverse. An anisotropic scaling $\mathbf{S} = \operatorname{diag}(\sigma_x,
\sigma_y)$ does not keep distances; the kernel multiplies by the smaller factor,

\begin{equation*}
\phi(\mathbf{x}) = \min(\sigma_x, \sigma_y)\;\phi_\Omega\big(\mathbf{S}^{-1}(\mathbf{x} - \mathbf{c}) + \mathbf{c}\big),
\tag{R4.23}
\end{equation*}

which is 1-Lipschitz, because $\lVert\mathbf{S}^{-1}\rVert = 1/\min(\sigma_x,
\sigma_y)$, and has the right sign: a bound in the sense of item 4.

### 3.7 Feature ids

A node returns a `FieldSample`: the field $\phi$ and an integer id $m$ naming the
member that owns the point (`NO_FEATURE` = −1 for none). The combinators pass the id
of the child that sets the value:

\begin{equation*}
\begin{gathered}
m_{A\cup B} = \begin{cases} m_A, & \phi_A \le \phi_B,\\ m_B, & \text{otherwise},\end{cases}
\qquad
m_{A\cap B} = \begin{cases} m_A, & \phi_A \ge \phi_B,\\ m_B, & \text{otherwise},\end{cases}
\\
m_{\mathrm{priority}} = \begin{cases} m_{j^\ast},\ j^\ast = \arg\max\{\pi_j : \phi_j < 0\}, & \text{inside some child},\\ m_{\arg\min_j \phi_j}, & \text{otherwise}.\end{cases}
\end{gathered}
\tag{R4.24}
\end{equation*}

Offsets, shells, slabs, rotations, trims and differences keep the id of their base;
`TagNode` stamps a constant; `MaxNode` takes it from a named child when given. The
priority $\pi$ of `PriorityUnionNode` comes from `KIND_PRIORITY` (core 0, skin 1,
web and stringer 2, cap 3, transverse 4, bond 5), so a cap overlapping a skin owns
the overlap. The field of the priority union is still $\min_j\phi_j$, (R4.11): the
priority decides only the id.

### 3.8 The field along the span

A blade member is a region $\Omega(z)$ in each section plane. Its field is the
distance in that plane, at each point's own $z$:

\begin{equation*}
\phi(x, y, z) = d_{\Omega(z)}(x, y),
\qquad
|\phi(\mathbf{x}, z) - \phi(\mathbf{y}, z)| \le \lVert\mathbf{x} - \mathbf{y}\rVert .
\tag{R4.25}
\end{equation*}

This is not a distance in three dimensions: $\partial\phi/\partial z$ follows the
planform and is not bounded by one. `LoftedPolygonSDF` evaluates (R4.9) on the
station polygon at the point's $z$; it does not interpolate distances between
stations, because a weighted mean of two distance fields is not the distance to the
weighted shape (and, under twist, misplaces a thin trailing edge). The span
combinators act plane by plane:

\begin{equation*}
\begin{gathered}
\phi_{\mathrm{law}}(\mathbf{x}, z) = \phi(\mathbf{x}, z) - \sigma\,\ell(z),
\\
\phi_{[z_0, z_1]}(\mathbf{x}, z) = \begin{cases}\phi(\mathbf{x}, z), & z_0 \le z \le z_1,\\ +\infty, & \text{otherwise}.\end{cases}
\end{gathered}
\tag{R4.26}
\end{equation*}

`LawOffsetNode` is (R4.15) or (R4.16) with an amount that is constant in each plane,
so its exactness statements hold plane by plane. `SpanIntervalNode` returns
$+\infty$ outside its interval, the field of the empty set, which leaves unions and
priorities unchanged. `ComplementNode` is $-\phi$.

### 3.9 Newton solves on fields

**Projection** onto $\phi = 0$ (`GeometrySolver.project`). Newton's method on the
scalar equation along the gradient is

\begin{equation*}
\mathbf{x}_{n+1} = \mathbf{x}_n - \phi(\mathbf{x}_n)\,\frac{\nabla\phi(\mathbf{x}_n)}{\lVert\nabla\phi(\mathbf{x}_n)\rVert^2} .
\tag{R4.27}
\end{equation*}

For an exact field off the medial axis, $\nabla\phi$ is the unit vector from the
nearest boundary point by (R4.2), so one step lands on that point. The kernel
estimates the gradient by central differences of step $\varepsilon$ (error
$O(\varepsilon^2)$), limits the step to `max_step`, and accepts a step only if
$|\phi|$ falls.

**Crossing** of two zero sets (`GeometrySolver.intersect`, `feature_curve`). Newton's
method on $(\phi_a, \phi_b) = \mathbf{0}$ in the plane solves, at each step,

\begin{equation*}
\begin{pmatrix}\partial_x\phi_a & \partial_y\phi_a\\ \partial_x\phi_b & \partial_y\phi_b\end{pmatrix}
\Delta\mathbf{x} = -\begin{pmatrix}\phi_a\\ \phi_b\end{pmatrix},
\qquad
\det = \nabla\phi_a\times\nabla\phi_b .
\tag{R4.28}
\end{equation*}

For exact fields the determinant is the sine of the angle between the two curves, so
the iteration converges quadratically where they cross transversally and is refused
(`|det| <= 1e-12`) where they are tangent. `zero_along` solves $\phi(\mathbf{p} +
t\mathbf{u}) = 0$ along a fixed direction with a fixed slope (the chord method),
which converges linearly.

### 3.10 Bounding boxes

`expr_bbox` returns a box that contains the region $\{\phi \le 0\}$, built from the
tree without evaluating it. For leaves it is the shape's box (a polygon's vertex
box; a rounded box's box grown by $r$, larger than needed). For the combinators,
with $\square$ the box of a child,

\begin{equation*}
\begin{gathered}
\square_{A\cup B} = \operatorname{hull}(\square_A, \square_B),
\qquad
\square_{A\cap B} = \square_A\cap\square_B,
\qquad
\square_{A\setminus B} = \square_A,
\\
\square_{\Omega_\delta} = \square_\Omega \oplus \max(\delta, 0),
\qquad
\square_{\mathrm{shell}} = \square_\Omega\oplus \tfrac{t}{2},
\\
\square_{\mathbf{Q}\Omega} = \operatorname{box}(\mathbf{Q}\,\text{corners of }\square_\Omega) .
\end{gathered}
\tag{R4.29}
\end{equation*}

The growth rule for an offset needs $\{\phi \le \delta\}$ to lie within $\delta$ of
$\{\phi \le 0\}$. That holds when $\phi$ is exact outside the region, but an
intersection is only a bound there (R4.12): $\{\max(\phi_A, \phi_B) \le \delta\} =
A_\delta\cap B_\delta$, which is not empty when $A$ and $B$ are disjoint but closer
than $2\delta$. The right box for that case is
$(\square_A\oplus\delta)\cap(\square_B\oplus\delta)$. The kernel keeps an empty
overlap as the gap interval $[\min(l, u), \max(l, u)]$ on each axis, which grown by
$\delta$ covers that box (fixed in SD-002; it used to collapse the overlap to a
point).

### 3.11 Zero loops

`extract_zero_loops` samples $\phi$ on an $n\times n$ grid and marches the squares
{cite:p}`lorensen1987`. On a cell edge from $\mathbf{p}_0$ to $\mathbf{p}_1$ whose
end values change sign, the crossing is placed by linear interpolation,

\begin{equation*}
\mathbf{x} = \mathbf{p}_0 + \frac{\phi_0}{\phi_0 - \phi_1}(\mathbf{p}_1 - \mathbf{p}_0),
\tag{R4.30}
\end{equation*}

exact for a field that is linear along the edge and in error by $O(h^2\kappa)$ for a
boundary of curvature $\kappa$ and cell size $h$. Each chained loop is then
projected onto $\phi = 0$ by (R4.27). A cell whose four corners alternate in sign
has four crossings and two segments, and the pairing is ambiguous; the bilinear
interpolant resolves it {cite:p}`nielson1991` by the sign of its saddle value

\begin{equation*}
\phi_{\mathrm{s}} = \frac{\phi_{00}\phi_{11} - \phi_{10}\phi_{01}}{\phi_{00} + \phi_{11} - \phi_{10} - \phi_{01}},
\tag{R4.31}
\end{equation*}

where $\phi_{ij}$ are the corner values. When $\phi_{\mathrm{s}}$ has the sign of
$\phi_{00}$ the corners of that sign are joined through the cell centre, otherwise
they are separated. A cell with four crossings also occurs when two boundaries pass
through it. The kernel emits both segments of such a cell, paired by (R4.31)
(fixed in SD-003; it used to keep the first two crossings only, and the loops broke
there).

## 4. Algorithm

Evaluating a tree on $N$ points (`SDFNode.evaluate`):

1. Check that the points are an $(N, 3)$ array; open an `EvalCache` for them.
2. Visit the tree depth first. A node already in the cache returns its stored
   `FieldSample`, so a shape shared by several members is evaluated once per batch.
3. A leaf evaluates its primitive (R4.3)–(R4.10) on $(x, y)$, or, for a lofted
   polygon, groups the points by $z$ and evaluates (R4.9) on each station polygon.
4. A combinator combines its children's fields by (R4.11)–(R4.26) and their ids by
   (R4.24).

`eval_expr_batch` with a memo does the same for `Expr` trees on coordinate arrays,
keyed by the (hashable) subtree, so a memo can be shared across trees.

## 5. Implementation

| Result | Function (`implicit_geometry_kernel.sdf` unless stated) |
|---|---|
| (R4.3)–(R4.8) | `kernels.primitives.sdf_circle`, `sdf_box`, `sdf_rounded_box`, `sdf_oriented_box`, `sdf_segment`, `sdf_capsule`, `sdf_half_plane`; leaves `Circle`, `Box`, `RoundedBox`, `OrientedBox`, `Capsule`, `HalfPlane` |
| (R4.9) | `kernels.polygon_numba._sdf_polygon_numba_core` (via `sdf_polygon`, `Polygon`, `AirfoilSDF`); unsigned form `_sdf_polyline_unsigned_numba_core` (`PolylineBandSDF`) |
| (R4.10) | `kernels.primitives.sdf_ellipse` (SD-001) |
| (R4.11)–(R4.13) | `kernels.csg.union`, `intersect`, `subtract`, `union_all`, `intersect_all`; nodes `Union`, `Intersect`, `Subtract`, `union_expr`, `intersect_expr`, `subtract_expr`; `nodes.span.MaxNode`, `nodes.tagged.TrimNode` |
| (R4.15)–(R4.18) | `kernels.csg.offset`, `shell`; nodes `Offset`, `Shell`, `Slab`; `nodes.tagged.FieldOffsetNode`, `MidSurfaceNode`, `PolylineBandSDF` |
| (R4.19)–(R4.21) | `kernels.csg.smooth_union`, `smooth_intersect`, `smooth_subtract`, `blend` |
| (R4.22), (R4.23) | `kernels.transforms.rotate_field`, `translate_field`, `scale_field`, `mirror_field_x`, `mirror_field_y`, `SDFFrame`; node `Rotate` |
| (R4.24) | `nodes.base.FieldSample`, `FeatureTable`, `KIND_PRIORITY`; `nodes.tagged.TagNode`, `PriorityUnionNode` |
| (R4.25), (R4.26) | `nodes.lofted.LoftedPolygonSDF`; `nodes.span.LawOffsetNode`, `SpanIntervalNode`, `ComplementNode` |
| (R4.27), (R4.28) | `solve.GeometrySolver.project`, `intersect`, `zero_along`, `gradient`; `extract.feature_curve` |
| (R4.29) | `nodes.bbox.expr_bbox` (via `Expr.bbox`) |
| (R4.30), (R4.31) | `extract.marching.extract_zero_loops` (SD-003) |
| Algorithm, Section 4 | `nodes.base.SDFNode.evaluate`, `EvalCache`; `nodes.csg.eval_expr_batch`, `compile_expr` |

## 6. Limits

1. **Bounds, not distances.** By (R4.11)–(R4.13) a CSG field is a distance only on
   one side of each operation. Anything that reads $|\phi|$ as a distance (a wall
   thickness, a clearance, a step size) is exact only there; culling by (R4.14) is
   safe everywhere.
2. **Offsets of non-convex shapes.** (R4.15) and (R4.16) are exact on one side only.
   Inside a dilation of a concave shape, or outside an erosion, $\phi - \delta$ is a
   bound.
3. **Varying thickness.** A wall offset by a thickness field (Section 3.4) has the
   right zero set but is not a distance.
4. **Smooth blends** move the zero set by up to $k/4$ along the bisector of the two
   boundaries, (R4.19); they are $C^1$, not $C^2$.
5. **Empty intersections.** The box of an empty intersection (SD-002, fixed) is
   the gap between the two boxes, not empty, so culling by it is conservative.
6. **Zero loops.** Two boundaries closer than one grid cell in a cell with only two
   crossings still merge into one loop; (R4.31) resolves the four-crossing case
   only (SD-003, fixed).
7. **Lofted polygons.** `LoftedPolygonSDF(accelerate=True)` caches each station
   polygon's edge arrays; every query still visits every edge, so the cost is
   linear in the number of edges (SD-004, fixed by describing it so).
8. **Along the span** (R4.25) is a distance in each plane only. A step in $z$ is not
   bounded by $|\phi|$.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R4.3)–(R4.9) | Circle, box, rounded box, oriented box, capsule and a star polygon against the distance to a boundary sampled at $4\times10^5$ points: largest difference $4.6\times10^{-6}$, the sampling spacing; half-plane exact | [sdf example 03](../../examples/sdf/03_distance_and_csg/README.md), check 1 |
| (R4.10) | Ellipse against the sampled distance: $1.2\times10^{-7}$, the sampling spacing | example 03, check 1; Figure R4.1; `tests/sdf/test_sd_fixes.py` |
| (R4.11)–(R4.13) | Union, intersection and difference of a disk and a box against the true boundary: exact side to $2.8\times10^{-6}$; on the bound side $|\phi|/|d| \le 1$ and signs agree | example 03, check 2 |
| (R4.14) | A min/max/offset tree is 1-Lipschitz over $2\times10^4$ random pairs | example 03, check 2 |
| (R4.15)–(R4.17) | Dilation, erosion and shell of an L-shaped polygon against independently built level sets | example 03, check 3 |
| (R4.19)–(R4.21) | Smooth union equals min beyond the band, sits $k/4$ below it at $a = b$, is 1-Lipschitz, and its slope is continuous at $|a-b| = k$ | example 03, check 4 |
| (R4.22), (R4.23) | Rotate-then-translate equals the oriented box; uniform scaling exact; anisotropic scaling 1-Lipschitz with correct sign | example 03, check 5 |
| (R4.27), (R4.28) | One projection step lands on a circle; the circle-line crossing converges quadratically | example 03, check 6; Figure R4.1 |
| (R4.29)–(R4.31) | The box of an offset of two disjoint boxes contains its region; two circles 0.05 apart extract as two loops; a single circle as one loop on $\phi = 0$ | example 03, check 7; `tests/sdf/test_sd_fixes.py` |
| (R4.24) | Ids of union, intersection, priority union and tags | `tests/sdf/test_field_meta.py` |
| (R4.25) | Lofted polygon field against the station polygon; Newton on the lofted field | `tests/sdf/test_lofted_solve.py` |
| (R4.9) | Primitive values, signs and polygon winding | `tests/section/test_primitives.py` |

```{figure} ../../examples/sdf/03_distance_and_csg/figures/distance_and_csg.png
:alt: Four panels: the shortfall of a union field inside the union; the smooth union against min; the ellipse field's error; Newton residuals for a circle-line crossing
:width: 100%

**Figure R4.1.** Top left: union of a disk and a box by (R4.11), with field contours;
the shading is $|d| - |\phi|$, zero outside the union and positive inside, where the
field is a bound. Top right: the smooth union (R4.19) against $\min$, $k/4$ below it
at $a = b$. Bottom left: the error of `sdf_ellipse` against the distance to a
densely sampled ellipse, at the sampling spacing. Bottom right: residuals of the two-field Newton solve (R4.28). Drawn by
`examples/sdf/03_distance_and_csg/distance_and_csg.py`.
```

These checks move to [V1](../v/index.md) when it is written.

## 8. References

{cite:t}`osher2003` for distance fields and the eikonal equation,
{cite:t}`hormann2001` for the winding number, {cite:t}`eberly2013` for the distance
to an ellipse, {cite:t}`quilez2013` for the polynomial smooth minimum,
{cite:t}`lorensen1987` for marching cells and {cite:t}`nielson1991` for the
asymptotic decider. Full entries are in [R0](R0_notation.md#bibliography).
