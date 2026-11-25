# R6 · Mid-surfaces, cells and the thin-wall idealisation

## 1. Problem

A shell model of a section does not use the wall solids of [R5](R5_section_layout.md)
directly. It uses a mid-surface for each wall, with a thickness, and the closed cells
the mid-surfaces enclose. This document derives:

- the skin mid-surface as an offset locus of the aerofoil, and the length and area
  laws of that offset;
- the identity that makes "thickness times mid-surface length" the exact wall area,
  and the error the same idealisation makes in the second moment;
- where the skin mid-surface ends at a sharp trailing edge;
- the cap mid-surface and the eccentricity carried when cap and skin share the skin
  mid-surface;
- the web feet, and how far a web mid-surface reaches past the inner skin;
- the wall graph whose bounded faces are the cells, with their number and the area
  identity the kernel checks.

Section properties from the solids are in R7; the shell mesh built on these
mid-surfaces is in R9.

## 2. Definitions and assumptions

1. **Aerofoil polygon** $\Omega_A$ with $N$ vertices in counter-clockwise order,
   perimeter $L$ and area $A$, in the chord frame of R5 (leading edge at the origin).
   At vertex $i$ the boundary turns by the exterior angle $\alpha_i \in (-\pi, \pi)$,
   positive where the polygon is convex. A closed simple polygon has
   $\sum_i \alpha_i = 2\pi$.
2. **Reflex vertices.** A vertex with $\alpha_i < 0$. The kernel's NACA polygons have
   none (example 05, check 2). Assumption A2 says the polygon has none; Section 6
   covers what changes where it has some.
3. **Sharp trailing edge.** The polygon closes at the trailing edge with included angle
   $\beta = \pi - \alpha_{\mathrm{TE}}$.
4. **Gauges.** As in R5: skin $t_s$, cap height $h$, web thickness $t_w$, each
   constant within a section.
5. **Offset locus.** For a depth $d > 0$, $\Omega_A(d)$ is the set of points of
   $\Omega_A$ at distance at least $d$ from $\partial\Omega_A$. Its boundary is the
   locus at depth $d$; its perimeter and area are $L(d)$ and $A(d)$.
6. **Smooth wall.** For the second moment, a mid-curve $\mathbf{r}(s)$ of arc length
   $s$ and length $L_m$, with unit tangent $\mathbf{T}$, unit normal $\boldsymbol\nu$
   to the left of $\mathbf{T}$, and signed curvature $\kappa$, with
   $\mathrm{d}\boldsymbol\nu/\mathrm{d}s = -\kappa\mathbf{T}$. On a counter-clockwise
   convex curve, $\boldsymbol\nu$ points inward and $\kappa > 0$.
7. **Assumption A3.** The wall is thin against its own radius of curvature:
   $t\,|\kappa| < 2$ everywhere, so the wall band does not fold over itself.

## 3. Derivation

### 3.1 The skin mid-surface

The skin is authored by one of four reference faces. Each one fixes the depth of the
skin mid-surface below the outer mould line. With $t$ the skin thickness and $D_I$ the
depth of the inner mould line where one is given:

\begin{equation*}
d_m =
\begin{cases}
t/2 & \text{outer mould line,}\\
D_I - t/2 & \text{inner mould line } (D_I = t \text{ if not given}),\\
D_I & \text{mid-surface } (t/2 \text{ if not given}),\\
D_I + t/2 & \text{skin bond line.}
\end{cases}
\tag{R6.1}
\end{equation*}

The skin mid-surface is the locus at that depth:

\begin{equation*}
M = \partial\,\Omega_A(d_m),
\qquad
\Omega_A(d) = \{\mathbf{p}\in\Omega_A : \operatorname{dist}(\mathbf{p}, \partial\Omega_A) \ge d\} .
\tag{R6.2}
\end{equation*}

$\Omega_A(d)$ is the erosion of $\Omega_A$ by a disc of radius $d$. At a convex vertex,
the two offset edges meet at a corner. That corner is at distance $d$ from both
edges, so it lies on $M$. At a reflex vertex, the exact locus is an arc of radius $d$
about the vertex. Under A2 there are no reflex vertices, so the polygon whose edges are
the original edges moved inward by $d$, and that meet at their intersections (a mitre
join), is exactly $\Omega_A(d)$. The kernel computes it as an inward buffer with mitre
joins.

### 3.2 Length and area of the offset

Move edge $i$, from vertex $i$ to vertex $i+1$, inward by $d$. At vertex $i$ the two
moved edges meet on the bisector of the vertex. Each edge loses
$d\tan(\alpha_i/2)$ there. So edge $i$ has length

\begin{equation*}
\ell_i(d) = \ell_i - d\left[\tan\frac{\alpha_i}{2} + \tan\frac{\alpha_{i+1}}{2}\right].
\end{equation*}

Summing over the edges, each vertex is counted twice:

\begin{equation*}
L(d) = L - 2dK,
\qquad
K = \sum_{i=1}^{N} \tan\frac{\alpha_i}{2} .
\tag{R6.3}
\end{equation*}

An eroded point lies inside the locus at every smaller depth, so
$\mathrm{d}A/\mathrm{d}d = -L(d)$. Integrating (R6.3) from zero:

\begin{equation*}
A(d) = A - dL + d^2 K .
\tag{R6.4}
\end{equation*}

Both laws hold while every edge keeps a positive length. Edge $i$ vanishes when
$\ell_i(d) = 0$, so (R6.3) and (R6.4) hold for $d \le d^*$, where

\begin{equation*}
d^* = \min_i \frac{\ell_i}{\tan(\alpha_i/2) + \tan(\alpha_{i+1}/2)} .
\tag{R6.5}
\end{equation*}

For a smooth convex curve $\tan(\alpha_i/2) \to \alpha_i/2$, which gives $K \to \pi$
and the Steiner formula $A(d) = A - dL + \pi d^2$. A sharp trailing edge adds
$\tan(\alpha_{\mathrm{TE}}/2) - \alpha_{\mathrm{TE}}/2$ to $K$. On the NACA 2418 of
example 05 that is 3.27, and $K = 6.41$.

### 3.3 The band identity

The skin solid of thickness $t$ is $\Omega_A \setminus \Omega_A(t)$. Its area, by
(R6.4) with $d = t$, compared with $t$ times the mid-surface length from (R6.3) with
$d = t/2$:

\begin{equation*}
\begin{gathered}
A - A(t) = tL - t^2K,
\qquad
t\,L(t/2) = tL - t^2K,
\\
A_{\text{skin}} = t\,L_m \quad\text{exactly, for } t \le d^* .
\end{gathered}
\tag{R6.6}
\end{equation*}

The thin-wall area, thickness times mid-surface length, is therefore not an
approximation for a polygonal skin. It stops being exact only when an edge vanishes
before the inner face is reached ($t > d^*$). On example 05's polygon $d^*$ is 2.74 mm,
because the trailing edge edges are 12.5 mm long. The relative error is then
$3\times10^{-6}$ at $t = 4$ mm and $1.5\times10^{-5}$ at 8 mm.

### 3.4 The end of the ring at a sharp trailing edge

The two edges at the trailing edge meet at angle $\beta$. Their offsets meet on the
bisector, at distance $r$ from the trailing edge with $r\sin(\beta/2) = d$:

\begin{equation*}
r_{\mathrm{TE}}(d) = \frac{d}{\sin(\beta/2)}, \qquad d \le d^* .
\tag{R6.7}
\end{equation*}

The skin mid-surface therefore ends $r_{\mathrm{TE}}(d_m)$ short of the trailing edge.
No mid-surface of constant thickness exists closer than that. The kernel does not trim
the skin here: the offset itself ends there. With $\beta = 24.5°$ and $d = 1$ mm,
$r = 4.72$ mm.

### 3.5 The second moment of a thin wall

Parametrise the wall band by arc length $s$ along the mid-curve and the offset $\eta$
along $\boldsymbol\nu$, with $|\eta| \le t/2$: $\mathbf{x} = \mathbf{r}(s) +
\eta\boldsymbol\nu(s)$. By Definition 6, $\partial\mathbf{x}/\partial s =
(1 - \eta\kappa)\mathbf{T}$ and $\partial\mathbf{x}/\partial\eta = \boldsymbol\nu$.
The area element is

\begin{equation*}
\mathrm{d}A = (1 - \eta\kappa)\,\mathrm{d}\eta\,\mathrm{d}s ,
\tag{R6.8}
\end{equation*}

and A3 keeps $1 - \eta\kappa > 0$. The area is $\int_0^{L_m}\!\int (1 -
\eta\kappa)\,\mathrm{d}\eta\,\mathrm{d}s = tL_m$, because the $\eta$ term is odd. This
is (R6.6) for a smooth wall.

The second moment about the axis $y = 0$ is
$\iint (y + \eta\nu_y)^2 (1 - \eta\kappa)\,\mathrm{d}\eta\,\mathrm{d}s$. Expand the
product and keep the even powers of $\eta$. With
$\int\eta^2\,\mathrm{d}\eta = t^3/12$:

\begin{equation*}
I_{xx} = t\int_0^{L_m} y^2\,\mathrm{d}s
+ \frac{t^3}{12}\int_0^{L_m}\left(\nu_y^2 - 2y\,\nu_y\,\kappa\right)\mathrm{d}s .
\tag{R6.9}
\end{equation*}

This is exact for any band that satisfies A3, because the integrand is a polynomial
in $\eta$. The thin-wall second moment keeps only the first term. Its relative error
is

\begin{equation*}
\frac{I_{\text{thin}}}{I_{xx}} - 1
= -\frac{t^2}{12}\,\frac{\int(\nu_y^2 - 2y\nu_y\kappa)\,\mathrm{d}s}{\int y^2\,\mathrm{d}s}
+ O(t^4)
= O\!\left(\frac{t^2}{\rho^2}\right),
\tag{R6.10}
\end{equation*}

where $\rho$ is the smaller of the radius of gyration of the mid-curve about the axis
and its radius of curvature.

On a circle of radius $R$ with inward normal, $\nu_y = -\sin\vartheta$ and
$\kappa = 1/R$, so the integrand is $3\sin^2\vartheta$. (R6.9) gives
$\pi R^3 t + \pi R t^3/4$. That equals the annulus value
$\tfrac{\pi}{4}[(R + t/2)^4 - (R - t/2)^4]$ exactly. The thin-wall error is then
$-t^2/(4R^2)/(1 + t^2/(4R^2))$.

### 3.6 Caps on the skin mid-surface

A cap is the slab of the inner mould line from depth $t_s$ to $t_s + h$ (R5.7). Its own
mid-surface is the offset locus (R6.2) at

\begin{equation*}
d_c = t_s + \frac{h}{2} .
\tag{R6.11}
\end{equation*}

The kernel samples the cap mid-surface at that depth when a cap is kept as its own
solid layer. A shell model instead keeps cap and skin as one laminate on the skin
mid-surface (decision [0004](../decisions/0004-cap-on-skin-mid-surface.md)). That
laminate has thickness $t_s + h$, from depth 0 to depth $t_s + h$, so its own
mid-plane is at depth $(t_s + h)/2$. Its offset from the skin mid-surface, positive
towards the cavity, is

\begin{equation*}
e = \frac{t_s + h}{2} - \frac{t_s}{2} = \frac{h}{2} .
\tag{R6.12}
\end{equation*}

The kernel leaves the nodes on the skin mid-surface, so they do not step at a cap
edge, and carries $e$ on every element over the cap. A consumer applies the offset in
its constitutive law. For the laminate stiffness $(\mathbf{A}, \mathbf{B},
\mathbf{D})$ about its own mid-plane, moving the reference surface by $e$ gives
$\mathbf{B}' = \mathbf{B} + e\mathbf{A}$ and
$\mathbf{D}' = \mathbf{D} + 2e\mathbf{B} + e^2\mathbf{A}$, with $e$ measured along the
laminate's own thickness coordinate.

### 3.7 Web feet

The web solid ends on the inner mould line (R5.9). Its mid-surface is the segment of
the web axis $\mathbf{p}_k + \lambda\mathbf{u}$ between the two points where the axis
meets the skin mid-surface:

\begin{equation*}
\mathbf{f}_k^{\pm} = \mathbf{p}_k + \lambda^{\pm}\mathbf{u},
\qquad \mathbf{f}_k^{\pm} \in M ,
\tag{R6.13}
\end{equation*}

where $\lambda^{\pm}$ are the roots nearest the web on each side. So the web and skin
mid-surfaces meet at a shared point, and the junction is a node of the shell.

Let $\gamma$ be the angle between the web axis and the skin normal at the foot. The
skin is locally straight across its own thickness under A3. From the inner mould line
to the skin mid-surface, the axis then covers a normal depth of $t_s/2$ over a length

\begin{equation*}
\delta = \frac{t_s}{2\cos\gamma} + O(\kappa t_s^2) .
\tag{R6.14}
\end{equation*}

That extension is the usual shell-junction idealisation. Over the length $\delta$ the
web band and the skin band overlap, so a thin-wall integral over all walls counts an
area of about $t_w\delta$ twice at each foot. Relative to the web area, that is
$t_s/(2H\cos\gamma)$ per foot for a web of height $H$.

### 3.8 The wall graph and the cells

Split the skin mid-surface at the $2n$ web feet, and at the leading and trailing edge
points that separate the suction and pressure sides. This gives a planar graph:

- the nodes are the split points;
- the edges are the skin arcs around the ring and the $n$ webs across it.

With $n$ webs whose feet are distinct:

\begin{equation*}
V = 2n + 2, \qquad E = (2n + 2) + n = 3n + 2 .
\tag{R6.15}
\end{equation*}

The graph is connected, and the webs do not cross. Euler's formula for a connected
planar graph, $V - E + F = 2$, then counts the faces. One face is the exterior.

\begin{equation*}
F = E - V + 2 = n + 2,
\qquad
\text{cells} = F - 1 = n + 1 .
\tag{R6.16}
\end{equation*}

The skin arcs bound the exterior on one side. So every skin arc bounds exactly one
cell, and every web bounds exactly two {cite:p}`deberg2008`.

The shoelace area of a closed polyline is a sum over its edges,
$A = \tfrac12\sum (x_j y_{j+1} - x_{j+1} y_j)$. A term changes sign when its edge is
reversed. Each cell boundary walks its skin arcs once, in the direction of the ring,
and each web is walked once in each direction by the two cells it separates. Summing
over the cells, the web terms cancel, and the skin terms add up to the ring:

\begin{equation*}
\sum_{c=1}^{n+1} A_c = A_M ,
\tag{R6.17}
\end{equation*}

with each $A_c$ of the same sign as $A_M$. This holds to rounding for any planar
subdivision, so the kernel uses it as a check: relative tolerance $10^{-9}$, not a
geometric tolerance.

## 4. Algorithm

The cells of one station (`build_cell_wall_loops`):

1. Take the skin and web mid-surfaces. Split the skin at the web feet (R6.13) if it is
   not already split, labelling each arc by side (suction, pressure) and position
   (leading edge, bay $i$, trailing edge).
2. Chain the arcs head to tail around the ring. Refuse if two consecutive arcs are more
   than $10^{-4}$ m apart, or the ring does not close.
3. Attach each web to the two nodes its ends coincide with. The match tolerance is
   1 nm, because a foot is a split point by construction. Refuse a web end that is off
   the ring, two webs on one node, or a web that starts and ends at the same node.
4. Refuse webs that cross or touch: any contact is a crossing, since no two webs share
   a node.
5. Walk the faces. Start at an unused arc and follow the ring; at a node where a web
   starts, take the web to its other node and continue along the ring from there. Stop
   when the start node is reached. Each arc is used once, and each web twice.
6. Order the cells from the leading edge cell (the one holding the leading edge arcs)
   across one shared web at a time.
7. Check that there are $n + 1$ cells (R6.16), that every cell has the ring's
   orientation, and the area identity (R6.17).

## 5. Implementation

| Result | Function (`implicit_geometry_kernel` unless stated) |
|---|---|
| (R6.1) | `section.layout.anchoring.airfoil_inward_depth_for_midsurface` |
| (R6.2)–(R6.4), (R6.7) | `section.boundary.airfoil_offset_ring` (mitre joins, `OFFSET_JOIN_STYLE`, `OFFSET_MITRE_LIMIT` 5); `section.walls.outer_skin._closed_mid_ring`, `OuterSkin.midline_polyline` |
| (R6.6), (R6.9) | Consumers integrate the strips with `ShellMidlineStrip.thickness_m`; the kernel's own solid properties are R7 |
| (R6.11) | `section.assembly.clips._airfoil_offset_arc`, `section.walls.spar_cap.SparCap.midline_polyline`, `section.midlines.export.build_solid_cap_strips` |
| (R6.12) | `blade.shell_graph_build` (`ShellRun.eccentricity_m`, `requires_eccentricity`); `structure.subpartition` for sub-partitions |
| (R6.13) | `section.midlines.junction_primitives.solve_web_foot_on_skin_mid`, `intersect_axis_offset_ring`; `section.walls.shear_web.ShearWebComponent._extend_to_skin_mid_surface` |
| Algorithm step 1 | `section.midlines.skin_split.split_skin_at_junctions` |
| (R6.15)–(R6.17), algorithm | `section.cells.builders.build_cell_wall_loops`, `_chain_ring`, `_attach_webs`, `_check_webs_do_not_cross`, `_faces`, `_order_le_to_te`; `section.cells.geom.shoelace_area` |

## 6. Limits

1. **Vanishing edges.** (R6.3)–(R6.7) hold for $d \le d^*$ (R6.5). Above $d^*$, the
   offset is still the exact locus, but the polynomial laws pick up the lost edges.
   The band identity (R6.6) is then off by $3\times10^{-6}$ at 4 mm and
   $1.5\times10^{-5}$ at 8 mm on a 2 m NACA 2418.
2. **Reflex vertices.** A mitre join at a reflex vertex lies $d/\cos(\alpha_i/2)$ from
   the vertex, not $d$: an error of about $d\alpha_i^2/8$, small where the polygon is
   densely sampled. The NACA polygons of R2 have no reflex vertex; an imported
   aerofoil with a concave aft pressure side and coarse sampling would need round
   joins at those vertices.
3. **Thin-wall second moment.** Dropping the $t^3$ term of (R6.9) has relative error
   $O(t^2/\rho^2)$ (R6.10): $3\times10^{-4}$ for an 8 mm skin on a 2 m chord. Thick
   root sections need the full (R6.9), or the solid properties of R7.
4. **Junction overlap.** (R6.14): the web mid-surface runs $t_s/(2\cos\gamma)$ into the
   skin band at each foot. A consumer that integrates every strip counts that overlap
   twice.
5. **Eccentricity** (R6.12) is a constant offset per element. It assumes the cap and
   the skin are bonded with no gap and share a normal, which holds under A1 of R5.
6. **Cells.** Webs must run from one skin side to the other. A web with both feet on
   one side, or webs that cross, are refused rather than cut into more cells.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R6.1) | Depths for the four reference faces | [section example 05](../../examples/section/05_midsurfaces_cells/README.md), check 1 |
| (R6.2) | Skin mid ring at distance $t/2$ from the outer mould line: $8\times10^{-17}$ m; no reflex vertex | example 05, check 2 |
| (R6.3)–(R6.6) | $L(d)$ exact and $A(d)$ to $7\times10^{-15}$ for $d \le d^*$; band identity to $10^{-12}$ for $t \le d^*$ | example 05, check 3 |
| (R6.7) | Trailing edge retreat $d/\sin(\beta/2)$ to $10^{-9}$ relative | example 05, check 4 |
| (R6.9), (R6.10) | Circle exact; aerofoil thin-wall error with log-log slope 2.01 | example 05, check 5 |
| (R6.11), (R6.12) | Cap mid at $t_s + h/2$ to $4\times10^{-17}$ m; cap runs carry $t_s$ and $e = h/2$ | example 05, check 6; `tests/test_shell_graph.py` |
| (R6.13), (R6.14) | Feet on the mid ring to $6\times10^{-17}$ m; reach $t_s/(2\cos\gamma)$ to $10^{-3}$ relative, $\gamma$ up to 22° | example 05, check 7; `tests/section/test_junction_primitives_web_foot.py`, `tests/section/test_web_foot_axis_snap.py` |
| (R6.15)–(R6.17) | $V$, $E$ and $n + 1$ cells for $n = 0$ to 5; $\sum A_c / A_M - 1$ below $2\times10^{-14}$; webs bound two cells, skin arcs one | example 05, check 8; [section example 01](../../examples/section/01_layout_catalogue/README.md) (all 40 cores) |

```{figure} ../../examples/section/05_midsurfaces_cells/figures/midsurfaces_cells.png
:alt: Top, the four cells of a 3D-F section bounded by the skin mid ring and three web mid-surfaces; bottom, the trailing edge offsets, the thin-wall second moment error against skin thickness, and a web foot on the skin mid ring
:width: 100%

**Figure R6.1.** Top: `3D-F` on a NACA 2418 of 2 m chord at 15° twist, in the blade
frame. The skin mid ring is split at the six web feet, and the three webs divide it into
$n + 1 = 4$ cells (R6.16), whose areas sum to the ring's (R6.17). Bottom left: at the
trailing edge, the mid ring and the inner mould line end short of the edge by
$d/\sin(\beta/2)$ (R6.7). Bottom centre: the thin-wall second moment error grows as
$t^2$ (R6.10). Bottom right: a web foot on the skin mid ring (R6.13); the web
mid-surface runs past the inner mould line by $t_s/(2\cos\gamma)$ (R6.14). Drawn by
`examples/section/05_midsurfaces_cells/midsurfaces_cells.py`.
```

These checks move to [V1](../v/index.md) when it is written.

## 8. References

The offset and the fields are those of [R4](R4_sdf_csg.md), and the walls those of
[R5](R5_section_layout.md). The thin-wall idealisation follows
{cite:t}`megson2017`; the face walk and Euler's formula follow {cite:t}`deberg2008`.
Full bibliography entries are in [R0](R0_notation.md#bibliography).
