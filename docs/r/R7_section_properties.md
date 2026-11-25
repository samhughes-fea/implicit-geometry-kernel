# R7 · Section properties by shoelace and quadtree integration

## 1. Problem

Each wall solid of a section (R5) has an area $A$, a centroid $(\bar X, \bar Y)$ and
second moments $I_{xx}$, $I_{yy}$, $I_{xy}$, as defined in (R1.17). The kernel
computes them in two ways:

- **Shoelace.** Green's theorem on the solid's boundary polygon. This is exact for
  the polygon.
- **Quadtree.** Quadrature of the solid's field $\phi$ over an adaptive quadtree.
  The error converges as the square of the smallest leaf.

The default "hybrid" backend uses the shoelace for the wall solids, whose polygons
the kernel builds, and the quadtree for regions that are known only as fields. This
document derives both rules and their errors, the leaf classes and the refinement
rule of the quadtree, and the reduction to the centroid.

## 2. Definitions and assumptions

1. **Frame.** Coordinates $(X, Y)$ in the station frame S of [R1](R1_frames.md),
   with its origin on the pitch axis. In the code these are written $(x, y)$.
2. **Polygon.** A simple polygon with vertices $\mathbf{p}_i = (x_i, y_i)$,
   $i = 1, \dots, n$, closed by $\mathbf{p}_{n+1} = \mathbf{p}_1$, counter-clockwise.
   Its edge cross products are $c_i = x_i y_{i+1} - x_{i+1} y_i$.
3. **Raw moments** about the origin:
   $S_x = \int X\,\mathrm{d}A$, $S_y = \int Y\,\mathrm{d}A$,
   $S_{xx} = \int X^2\,\mathrm{d}A$, $S_{yy} = \int Y^2\,\mathrm{d}A$,
   $S_{xy} = \int XY\,\mathrm{d}A$. In the notation of (R1.17), $I_{xx} = S_{yy}$,
   $I_{yy} = S_{xx}$ and $I_{xy} = S_{xy}$ about the origin.
4. **Field.** $\phi$ is negative inside the solid and 1-Lipschitz, as in
   [R4](R4_sdf_csg.md).
5. **Quadtree.** Square leaves of side $w$ and diagonal $D = \sqrt2\,w$. The
   smallest side is $h$ (`min_cell_size`). Each leaf stores $\phi$ at its four
   corners, numbered counter-clockwise from the lower left, and at its centre.

## 3. Derivation

### 3.1 The shoelace rule

Join the origin to every edge of the polygon. Triangle $T_i = (\mathbf{0},
\mathbf{p}_i, \mathbf{p}_{i+1})$ has signed area $c_i/2$. It is positive where the
edge runs counter-clockwise about the origin and negative otherwise. Every point
of the polygon is covered with net weight one, so any integral over the polygon is
the signed sum over the triangles. For the area:

\begin{equation*}
A = \frac12 \sum_{i=1}^{n} c_i .
\tag{R7.1}
\end{equation*}

A polynomial integrated over a triangle with vertices $\mathbf{v}_0, \mathbf{v}_1,
\mathbf{v}_2$ is fixed by its vertex values. For degree two {cite:p}`pilkey2002`:

\begin{equation*}
\begin{gathered}
\int_T X\,\mathrm{d}A = \frac{A_T}{3}\sum_k x_k,
\qquad
\int_T X^2\,\mathrm{d}A = \frac{A_T}{6}\Big(\sum_k x_k^2 + \sum_{j<k} x_j x_k\Big),
\\
\int_T XY\,\mathrm{d}A = \frac{A_T}{12}\Big(2\sum_k x_k y_k + \sum_{j\ne k} x_j y_k\Big).
\end{gathered}
\tag{R7.2}
\end{equation*}

Put $\mathbf{v}_0 = \mathbf{0}$, $\mathbf{v}_1 = \mathbf{p}_i$, $\mathbf{v}_2 =
\mathbf{p}_{i+1}$ and $A_T = c_i/2$. Every term with $\mathbf{v}_0$ vanishes, and
summing over $i$ gives the first moments

\begin{equation*}
S_x = \frac16\sum_i (x_i + x_{i+1})\,c_i,
\qquad
S_y = \frac16\sum_i (y_i + y_{i+1})\,c_i ,
\tag{R7.3}
\end{equation*}

and the second moments

\begin{equation*}
\begin{gathered}
S_{xx} = \frac1{12}\sum_i \big(x_i^2 + x_i x_{i+1} + x_{i+1}^2\big)\,c_i,
\\
S_{yy} = \frac1{12}\sum_i \big(y_i^2 + y_i y_{i+1} + y_{i+1}^2\big)\,c_i ,
\end{gathered}
\tag{R7.4}
\end{equation*}

\begin{equation*}
S_{xy} = \frac1{24}\sum_i \big(x_i y_{i+1} + 2x_i y_i + 2x_{i+1}y_{i+1} + x_{i+1}y_i\big)\,c_i .
\tag{R7.5}
\end{equation*}

The integrands are polynomials, and the triangles are exact, so (R7.1)–(R7.5) are
exact for the polygon, with no quadrature error. The centroid is
$\bar X = S_x/A$, $\bar Y = S_y/A$.

### 3.2 Holes and parts

The rules are linear in the region. A polygon with holes $H_j$ integrates as the
outer ring minus each hole, where each hole's sums are taken counter-clockwise. A
multipolygon with parts $P_k$ integrates as the sum of its parts. Both hold for any
raw moment $M$ about a common origin:

\begin{equation*}
M(P \setminus \textstyle\bigcup_j H_j) = M(P) - \sum_j M(H_j),
\qquad
M\big(\textstyle\bigcup_k P_k\big) = \sum_k M(P_k) .
\tag{R7.6}
\end{equation*}

The moments must be summed about the common origin, not about each part's own
centroid. An earlier version shifted them through the combined centroid and lost
each part's offset ([SG-036](../section/bug_log.md#sg-036)).

### 3.3 Reduction to the centroid

The centroidal moments follow from the parallel-axis relations (R1.24) solved for
S′:

\begin{equation*}
I'_{xx} = S_{yy} - A\bar Y^2,
\qquad
I'_{yy} = S_{xx} - A\bar X^2,
\qquad
I'_{xy} = S_{xy} - A\bar X\bar Y .
\tag{R7.7}
\end{equation*}

The radii of gyration are $r_x = \sqrt{I'_{xx}/A}$ and $r_y = \sqrt{I'_{yy}/A}$.
(R7.7) subtracts two nearly equal numbers when the centroid is far from the origin
against the section's size. The relative rounding error then grows by
$(\bar X^2 + \bar Y^2)/r^2$.

### 3.4 A polygon for a curved boundary

The shoelace is exact for the polygon, so its only error is the polygon's. A chord
of length $\ell$ across an arc of curvature $\kappa$ cuts off a segment of area
$\kappa\ell^3/12 + O(\ell^5)$. For a circle of radius $R$ sampled by an inscribed
$N$-gon, (R7.1) gives $\tfrac N2 R^2\sin(2\pi/N)$, so the relative deficit is

\begin{equation*}
1 - \frac{\sin(2\pi/N)}{2\pi/N} = \frac{2\pi^2}{3N^2} + O(N^{-4}) .
\tag{R7.8}
\end{equation*}

The error is second order in the vertex spacing. The aerofoil polygons of R2 are
sampled finely enough that this is the polygonisation error quoted there. The
offset polygons of R6 inherit it.

### 3.5 Interior leaves

A leaf with every corner and its centre inside ($\phi \le 0$) contributes its whole
rectangle $[x_0, x_1]\times[y_0, y_1]$. The integrals separate:

\begin{equation*}
\begin{gathered}
A = (x_1 - x_0)(y_1 - y_0),
\qquad
S_x = A\,\frac{x_0 + x_1}{2},
\\
S_{xx} = A\,\frac{x_0^2 + x_0x_1 + x_1^2}{3},
\\
S_{xy} = \frac{x_1^2 - x_0^2}{2}\cdot\frac{y_1^2 - y_0^2}{2}
= A\,\frac{(x_0 + x_1)(y_0 + y_1)}{4},
\end{gathered}
\tag{R7.9}
\end{equation*}

and likewise in $y$. Before FIXPR the kernel's product moment had a factor
$\tfrac12$ on the mixed terms $x_0y_1 + x_1y_0$
([IN-001](../integration/bug_log.md#in-001), fixed); it now evaluates the factored
form, and a circle of radius 0.2 centred at $(0.2, 0.1)$ gives $I'_{xy}$ at
round-off.

### 3.6 Cut leaves

On a leaf whose corners change sign, the kernel approximates $\phi$ as linear along
each edge. The crossing on the edge from corner $a$ to corner $b$ is then at
(R4.30), $\mathbf{x} = \mathbf{p}_a + \phi_a/(\phi_a - \phi_b)\,(\mathbf{p}_b -
\mathbf{p}_a)$. The inside of the leaf is approximated by the polygon that walks
the corners in order, taking each inside corner and each edge crossing:

\begin{equation*}
Q = \big(\ \dots,\ \mathbf{p}_a\ [\text{if } \phi_a \le 0],\ \mathbf{x}_{ab}\ [\text{if } \operatorname{sign}\phi_a \ne \operatorname{sign}\phi_b],\ \dots\big) .
\tag{R7.10}
\end{equation*}

The leaf contributes the shoelace moments (R7.1)–(R7.5) of $Q$. If $\phi$ is linear
on the leaf, its zero set is a straight line through the two crossings, and $Q$ is
exactly the inside. In a saddle leaf (corners alternating in sign), $Q$ is one
hexagon that joins the two inside corners through the leaf. This is the joined case
of (R4.31), taken without testing the decider.

### 3.7 Leaf classes and refinement

From the corner values $\phi_{1..4}$ and the centre value $\phi_c$, each leaf falls
into one of four classes:

\begin{equation*}
\begin{array}{ll}
\text{exterior:} & \phi_{1..4} > 0,\ \phi_c > 0 \quad\Rightarrow\ \text{nothing},\\
\text{interior:} & \phi_{1..4} \le 0,\ \phi_c \le 0 \quad\Rightarrow\ \text{(R7.9)},\\
\text{cut:} & \text{otherwise, with corners of both signs} \quad\Rightarrow\ \text{(R7.10)},\\
\text{hollow:} & \phi_{1..4} > 0,\ \phi_c \le 0 \quad\Rightarrow\ \text{split up to 2 levels}.
\end{array}
\tag{R7.11}
\end{equation*}

A hollow leaf holds a feature thinner than the leaf, which no corner sees. It is
split into four, and the children are classified again, at most twice.

The tree is refined from square roots. A cell is split, until its side reaches $h$,
when its probes (corners, centre and edge midpoints, or 13 points with
`refine_probes="full"`) either change sign, or come within the band:

\begin{equation*}
\begin{gathered}
\min_{\text{probes}} |\phi| < b(D),
\\
b(D) =
\begin{cases}
\min(b_{\max},\ sD) & \text{filled solids},\\
sD + b_{\max} & \text{thin bands } (\texttt{thin\_band}),
\end{cases}
\end{gathered}
\tag{R7.12}
\end{equation*}

with $s$ = `refine_band_scale` (1.5) and $b_{\max}$ = `refine_band`. The probes form
a $3\times3$ lattice of spacing $w/2$, so every point of the cell is within $D/4$ of
a probe. If the cell holds a boundary point, then by the Lipschitz bound some probe
has $|\phi| \le D/4$. That probe triggers (R7.12) whenever $b(D) > D/4$. This always
holds for thin bands. For filled solids it holds once $D < 4b_{\max}$. Above that
size, the rule relies on a sign change between probes, and a feature narrower than
the probe spacing can be missed until the hollow test of (R7.11) catches it.

### 3.8 Convergence of the quadtree

Every cut leaf has side $h$ by (R7.12), and there are about $L/h$ of them along a
boundary of length $L$. Two errors arise on each:

- Interpolating a 1-Lipschitz field with curvature $\kappa$ moves each crossing by
  $O(\kappa h^2)$.
- The chord between the crossings cuts off $O(\kappa h^3)$ of area, by Section 3.4.

Summed over the cut leaves:

\begin{equation*}
\frac{|A_h - A|}{A} = O\!\left(\frac{\kappa L h^2}{A}\right),
\tag{R7.13}
\end{equation*}

and the same order holds for each moment. Interior leaves are exact by (R7.9). On a circle of radius 0.3, the area error falls from
$5.8\times10^{-4}$ at $h = 0.02$ to $8.8\times10^{-6}$ at $h = 0.0025$, a log-log
slope of 2.02.

### 3.9 The hybrid backend

The kernel builds every wall solid (skin, caps, webs, cores) as a polygon:

- offsets by the polygon buffer of R6, (R6.2);
- the caps' bands and the webs' strips by polygon intersection.

(R7.1)–(R7.7) are exact on those polygons. That leaves only the polygonisation error
(R7.8) of the aerofoil and of the offset joins. Regions defined only by a residual
field, the leading and trailing edge inserts, have no polygon and are integrated on
the quadtree:

\begin{equation*}
P(\text{label}) =
\begin{cases}
\text{(R7.1)–(R7.7), polygon} & \text{walls and cores},\\
\text{(R7.9)–(R7.12), quadtree} & \text{inserts, other labels}.
\end{cases}
\tag{R7.14}
\end{equation*}

At the default leaf size ($h = 0.0015c$), the two backends agree on the walls of a
2 m NACA 2418 section to about $10^{-3}$ in area, which is the quadtree's (R7.13).

## 4. Algorithm

Properties of one section (`compute_section_properties`, backend `"hybrid"`):

1. Split the labels into constructive (walls, cores) and residual (inserts)
   (`classify_property_labels`).
2. Emit a polygon for each constructive label, and drop empty ones.
3. For each polygon, sum (R7.1)–(R7.5) over the outer ring and the holes, (R7.6),
   and reduce by (R7.7).
4. If residual labels exist, build one shared quadtree for them. Roots are squares
   of side about the section's minor extent; refine by (R7.12) to $h$.
5. Sample the centres, classify the leaves by (R7.11), accumulate (R7.9) and
   (R7.10) for each label, and reduce by (R7.7).

## 5. Implementation

| Result | Function (`implicit_geometry_kernel` unless stated) |
|---|---|
| (R7.1)–(R7.5) | `integration.shoelace.ring_sums_S`; `integration.quadtree.moments._polygon_moments`, `_shoelace_moments_batch` |
| (R7.6) | `integration.shoelace.section_properties_S` (holes and `MultiPolygon` parts) |
| (R7.7) | `integration.shoelace.section_properties_Sprime`; `integration.quadtree.quadrature._raw_to_section_properties` |
| (R7.9) | `integration.quadtree.moments._rectangle_moments`, `_rectangle_moments_batch` |
| (R7.10) | `integration.quadtree.moments._inside_polygon_from_cell`, `_boundary_recipe`, `_boundary_moments_batch` |
| (R7.11) | `integration.quadtree.quadrature._accumulate_moments_vectorized`; `integration.quadtree.moments._accumulate_cell` (hollow leaves) |
| (R7.12) | `integration.quadtree.soa._effective_refine_bands`, `_cells_need_refine_batch`; `integration.quadtree.soa_build.build_shared_quadtree`; defaults in `section.properties.backend.quadtree_spec_for_airfoil` |
| (R7.14), algorithm | `section.properties.compute.compute_section_properties`, `_compute_hybrid`; `section.boundary.labels.classify_property_labels`; `section.properties.area_mesh.section_properties_via_quadtree` |

## 6. Limits

1. **Polygonisation** (R7.8) is second order in the vertex spacing. The shoelace is
   exact only for the polygon it is given.
2. **Saddle leaves** (R7.10) are always joined. Two boundaries closer than one leaf
   can be merged in the fill. The error is at most the leaf's area, on the rare
   leaves where it occurs.
3. **Coarse cells.** Above $D = 4b_{\max}$, refinement of a filled solid relies on
   sign changes (R7.12). A sliver thinner than the probe spacing inside such a cell
   is found only by the hollow test, and only two levels deep.
4. **Leaf budget.** `max_cells` stops refinement. If it is reached, cut leaves
   coarser than $h$ remain, and (R7.13) holds with the larger size.
5. **The grid backend** (`SDFGrid.section_properties`) counts cell centres with
   $\phi < 0$. Its error is first order in the grid step, and it is for debugging
   only.
7. **Cancellation** in (R7.7) for sections far from the origin: see Section 3.3.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R7.1)–(R7.5) | A rectangle off the origin against closed forms, and a triangle against (R7.2): $2\times10^{-16}$ | [integration example 02](../../examples/integration/02_property_quadrature/README.md), check 1 |
| (R7.6) | A hole subtracts; two disjoint squares add about the origin | example 02, check 2; `tests/integration/test_area_2d_quadtree.py` |
| (R7.7) | S′ moments equal the moments of the polygon moved to its centroid: $6\times10^{-17}$ | example 02, check 3 |
| (R7.8) | Inscribed $N$-gon area exact; deficit ratio 16.00 per doubling | example 02, check 4; Figure R7.1 |
| (R7.9) | Area, first, second and product moments of a leaf exact; $I'_{xy}$ of an off-axis circle at round-off | example 02, checks 5 and 7; `tests/integration/test_quadtree_quadrature.py` |
| (R7.10) | Straight boundaries with one, two and three inside corners: exact | example 02, check 6 |
| (R7.12), (R7.13) | Every cut leaf at $h$; area and $I_{xx}$ converge with slope 2.02 | example 02, check 7; Figure R7.1; `tests/integration/test_quadtree_quadrature.py` |
| (R7.14) | Hybrid against quadtree on the walls of a 2D-F section: $1.2\times10^{-3}$ | example 02, check 8; [integration example 01](../../examples/integration/01_quadtree_vs_shoelace/README.md) |

```{figure} ../../examples/integration/02_property_quadrature/figures/property_quadrature.png
:alt: Left, quadtree leaves on a quarter circle coloured as interior or cut; centre, area and second moment errors falling as h squared; right, the area deficit of an inscribed N-gon against 2 pi squared over 3 N squared
:width: 100%

**Figure R7.1.** Left: the leaves of a quadtree on a circle at $h = 0.02$. Interior
leaves contribute their rectangle (R7.9), and cut leaves their inside polygon
(R7.10). Centre: the relative errors of area and $I_{xx}$ fall as $h^2$ (R7.13).
Right: the shoelace on an inscribed $N$-gon, against (R7.8). Drawn by
`examples/integration/02_property_quadrature/property_quadrature.py`.
```

These checks move to [V1](../v/index.md) when it is written.

## 8. References

Moments and the parallel-axis relations are as in {cite:t}`megson2017` and
{cite:t}`pilkey2002`. Cut cells follow marching squares {cite:p}`lorensen1987`.
The fields are those of [R4](R4_sdf_csg.md) and the polygons those of
[R6](R6_midsurfaces_cells.md). Full bibliography entries are in
[R0](R0_notation.md#bibliography).
