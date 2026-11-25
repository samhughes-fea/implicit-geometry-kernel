# R9 · Shell mesh generation

## 1. Problem

`mesh.shell.shell_mesh` builds the blade's conforming shell mesh from a
`BladeGeometry`: quadrilaterals on the skin and web mid-surfaces, triangles on the
ribs and the tip closure. Its nodes lie on the mid-surfaces of the cuts
([R6](R6_midsurfaces_cells.md)) at a set of span stations; every element carries
the partition it lies on, with that partition's thickness, layers and cap
eccentricity ([ADR 0004](../decisions/0004-cap-on-skin-mid-surface.md),
[ADR 0005](../decisions/0005-kernel-shell-mesh.md)). This document derives the
choice of span stations, the chordwise element size, the placement of nodes and
mesh lines on each station, the span grading, the element data, the triangle fill
of ribs and the tip, the orientation rule and the hard invariants.

The default mesher since 0.9.6 is the $(t, z)$ partition mesher
([ADR 0015](../decisions/0015-partition-mesher.md)), used whenever the mesh is sized
by $h$ alone and the skin loop has the same partitions at every station. Otherwise
the ruled-strip mesher is used; Section 3.10 gives the rule it adds.

## 2. Definitions and assumptions

1. **Target edge.** $h$ (metres), or $h/c$ as a fraction of the local chord. The
   edge at $z$ is $s(z)$.
2. **Stations.** Span positions $z_0 < z_1 < \dots < z_{m}$ where nodes are
   placed; a bay is $[z_j, z_{j+1}]$.
3. **Skin loop.** The closed mid ring of the shell graph at a station, walked as one
   ring from the partition boundary nearest the trailing edge (TE). Its normalised
   arc length is $u \in [0, 1)$; each partition holds a $u$ interval.
4. **Mesh line.** A curve along the span on which a node sits at every station: a
   web foot, a partition boundary with no web foot on it, or a sub-partition band
   edge ([ADR 0016](../decisions/0016-property-only-sub-partitions.md)).
5. **Sagitta budget.** $\varepsilon = 0.0015\,c$ (`EPS_OVER_C`): the largest
   distance allowed between a straight element edge and the curve it replaces.
6. **Element data.** Each element records `surface` (skin, web, rib, tip closure),
   `partition`, `thickness`, `eccentricity`, `layers`, and the bay it is in.

## 3. Derivation

### 3.1 Target edge and span bays

Given $h$, $h/c$, or both, the edge at $z$ is

\begin{equation*}
s(z) = \min\big(h,\ (h/c)\,c(z)\big),
\tag{R9.1}
\end{equation*}

taking only the term that is given; a call with neither (and no station count) uses
$h/c = 0.05$. The stations start from the **anchors**: the defining stations, the
root and tip, the forced stations of (R8.8), every rib, the start and end of every
partial region, the stations where a web leaves the section (Section 3.2), and the
key crossings. A gap of length $D$ between anchors gets $n$ equal bays with the
step closest to $s$ inside $[s/2, s]$; the step $D/n$ falls as $n$ grows, so the
closest is the smallest admissible $n$:

\begin{equation*}
n = \max\Big(1,\ \Big\lceil \frac{D}{s} \Big\rceil\Big),
\qquad
\frac{s}{2} < \frac{D}{n} \le s \quad (n > 1),
\tag{R9.2}
\end{equation*}

with $s = \min(s(z_a), s(z_b))$ at the gap ends. The lower bound follows from
$n < D/s + 1 \le 2D/s$ for $D > s$.

### 3.2 Where a web leaves the section

A web frozen at the leading edge (LE) or pitch axis is dropped, not slid, once the
inner mould line no longer holds it (R5.3). It carries no stored end, so the end is
found by bisection on whether the web is in the cut's `dropped_web_indices`. With
the bracket $[z_0, z_{\mathrm{tip}}]$, after $k$ halvings

\begin{equation*}
z_{\mathrm{hi}} - z_{\mathrm{lo}} = \frac{z_{\mathrm{tip}} - z_0}{2^k}
\le 10^{-6} \max(z_{\mathrm{tip}} - z_0,\ 1),
\tag{R9.3}
\end{equation*}

which needs at most 20 of the 40 halvings allowed; the station is $z_{\mathrm{hi}}$,
the first point known to be dropped.

### 3.3 Chordwise element size

**Sagitta edge.** A chord of length $\ell$ on a circle of radius $R = 1/\kappa$ is
at most $\delta = R(1 - \cos(\ell/2R)) = \ell^2/(8R) + O(\ell^4/R^3)$ from the arc.
Setting $\delta = \varepsilon$:

\begin{equation*}
\ell_{\mathrm{sag}} = \sqrt{\frac{8\varepsilon}{\kappa}},
\qquad
\ell_0 = \min(h,\ \ell_{\mathrm{sag}}),
\tag{R9.4}
\end{equation*}

where $\kappa \le 10^{-12}$ gives $\ell_{\mathrm{sag}} = \infty$. The neglected term
makes the true height smaller than $\varepsilon$ by $\varepsilon\,\kappa^2\ell^2/48$:
1.2 % at a 15 mm edge on a 20 mm nose radius, 0.15 % on a 0.5 m radius.

**Windowed curvature.** $\kappa$ is the turn of the tangent angle $\theta(s)$ over a
window $w = 0.02\,c$ of arc, divided by the window (clipped at the ends of the run):

\begin{equation*}
\hat\kappa(s) = \frac{\big|\theta(s + w/2) - \theta(s - w/2)\big|}{w},
\tag{R9.5}
\end{equation*}

with $\theta$ linear between segment mid-points. On a circle of radius $R > w$ it
is $1/R$; a single kink of angle $\alpha$ gives $\alpha/w$ instead of $\alpha$ over
one short segment, so a polyline glitch does not shrink the mesh. Each partition
is measured on its own points, so a sharp TE or any other partition end is a mesh
line, not a curvature.

**Grading.** The edge is then limited so that it grows by at most $g = 0.25$ per
metre of arc, by one forward and one backward sweep over the run's points:

\begin{equation*}
\ell_i \leftarrow \min(\ell_i,\ \ell_{i-1} + g\,\Delta s_{i-1}),
\qquad
\ell_i \leftarrow \min(\ell_i,\ \ell_{i+1} + g\,\Delta s_i).
\tag{R9.6}
\end{equation*}

The forward sweep makes $\ell_i - \ell_{i-1} \le g\,\Delta s_{i-1}$ and the
backward sweep, which only lowers values, makes $\ell_i - \ell_{i+1} \le g\,\Delta s_i$
without breaking the first bound, so $\ell$ is $g$-Lipschitz in arc length and never
above $\ell_0$. The refinement at the nose spreads into the flatter skin over a
length $(h - \ell_{\min})/g$.

### 3.4 Edge measure and node placement

The number of edges that fit between the start of a run and arc length $s$ is

\begin{equation*}
W(s) = \int_0^s \frac{\mathrm{d}\sigma}{\ell(\sigma)},
\qquad
W_{i+1} = W_i + \frac{\Delta s_i}{2}\Big(\frac{1}{\ell_i} + \frac{1}{\ell_{i+1}}\Big),
\tag{R9.7}
\end{equation*}

by the trapezoidal rule on the run's points; on the ring, the partitions' measures
are laid end to end. Between two consecutive mesh lines at $u_a < u_b$, segment $k$
gets the largest count any station needs, and its nodes are at equal steps of $W$:

\begin{equation*}
\begin{gathered}
n_k = \max_j \big\lceil W_j(u_b) - W_j(u_a) \big\rceil,
\\
u_{k,i} = W_j^{-1}\Big(W_j(u_a) + \frac{i}{n_k}\big(W_j(u_b) - W_j(u_a)\big)\Big),
\end{gathered}
\tag{R9.8}
\end{equation*}

$i = 0, \dots, n_k$. Each step is at most one unit of $W$, so every edge is at most
the local graded edge (to the trapezoidal rule's accuracy), and the count of each
segment is the same at every station, which makes the patch a structured grid.
Webs get the same treatment with $\ell = h$ (no curvature term). If the whole loop
has fewer than 8 edges, edges are added one at a time to the segment with the
longest edge, so a rib or the tip always has a ring to fill
([ADR 0007](../decisions/0007-rib-cells-as-triangles.md)).

### 3.5 Mesh lines that exist only on part of the span

Every web foot is a mesh line at every station, so the number of lines round the
loop never changes. Where the web is absent, its foot is held at the station $r$
nearest that has the web, at the same fraction of the gap between the two
neighbouring lines that exist at both stations:

\begin{equation*}
u_{\mathrm{here}} = u^{-}_{\mathrm{here}}
+ \frac{u_{\mathrm{ref}} - u^{-}_{\mathrm{ref}}}{u^{+}_{\mathrm{ref}} - u^{-}_{\mathrm{ref}}}
\big(u^{+}_{\mathrm{here}} - u^{-}_{\mathrm{here}}\big),
\tag{R9.9}
\end{equation*}

where $u^{\pm}$ are the nearest shared lines above and below (the ring start and
end if there is none). The held line moves with its neighbours and cannot cross
them. A partition boundary with a web foot on it is that foot's line; one without is
a line of its own. A cap edge that moves from one web foot to another where a web
ends therefore moves between two existing lines. The lines must come in the same
order, and two lines that meet (share a node) must meet, at every station;
otherwise the mesher raises `PartitionMeshError` and `shell_mesh` falls back to the
ruled-strip mesher.

### 3.6 Span grading

Curvature sizing shortens the skin edges near the nose. Without a matching span
step the quads there would be long in span and short across it. With $\bar\ell_j$
the mean loop edge at station $j$ (loop length over $W$ of the whole loop), a bay is
split evenly into

\begin{equation*}
n = \Big\lceil \frac{z_{j+1} - z_j}{1.25\,\min(\bar\ell_j, \bar\ell_{j+1})} \Big\rceil,
\quad\text{reduced while}\quad \frac{z_{j+1} - z_j}{n} < \frac{h}{2},
\tag{R9.10}
\end{equation*}

in one pass. A skin edge shorter than $h/2$ is the geometry (a tight nose), and the
span step does not chase it.

### 3.7 Quads and their data

Between stations $j$ and $j+1$, node $i$ of a row joins node $i$ of the next, so the
skin quad is $(a_i, a_{i+1}, b_{i+1}, b_i)$, and likewise on each web present at both
stations. A quad takes the partition that holds its mid-parameter at the lower
station, so a cap that ends at a station ends on that row. If the partition is the
same at both stations its data are averaged; otherwise the lower station's apply:

\begin{equation*}
t_e = \begin{cases}
\tfrac12\big(t_P(z_j) + t_P(z_{j+1})\big) & P(z_j) = P(z_{j+1}),\\
t_P(z_j) & \text{otherwise},
\end{cases}
\tag{R9.11}
\end{equation*}

where $P$ is the partition at $\tfrac12(u_i + u_{i+1})$, and the same for the
eccentricity $e$. Web feet are not welded: a web's end node
*is* the skin node on its foot line, by index, so the junction is conforming by
construction.

### 3.8 Ribs and the tip closure

At a rib station, and at the tip, the webs split the station's skin ring into $n+1$
cells (R6.17), each bounded by skin and web nodes of that station. Each cell is
filled with triangles whose boundary nodes are those nodes, so only interior nodes
are new. The fill places seeds about $r$ apart (`polygon_size`, default $h$): a
jittered lattice of candidates is kept if it is inside the cell and at least
$0.62\,r$ from its boundary, and then thinned by dart throwing on a background grid
of cell size $r/\sqrt2$, so that no two kept seeds are closer than $r$ (a grid cell's
diagonal is $r$, so it holds at most one seed, and a $5\times5$ block of cells covers
the disk of radius $r$) {cite:p}`bridson2007`. The cloud is triangulated by Delaunay
{cite:p}`bowyer1981,watson1981`, its boundary edges are recovered by flips, and up
to four Laplacian passes move each interior seed to the mean of its neighbours,

\begin{equation*}
\mathbf{p}_i \leftarrow \frac{1}{|N(i)|}\sum_{k \in N(i)} \mathbf{p}_k
\quad\text{if}\quad \mathbf{p}_i^{\mathrm{new}} \in \Omega
\ \text{and}\ d(\mathbf{p}_i^{\mathrm{new}}, \partial\Omega) \ge 0.35\,r,
\tag{R9.12}
\end{equation*}

with a new triangulation after each pass. Triangles whose centroid is outside the
cell, or whose area is below $10^{-16}$, are dropped, and so are slivers on three
collinear ring nodes. Each cell is then a triangulated disk: $V - E + F = 1$, and
the triangle areas sum to the area of its boundary ring.

### 3.9 Orientation

Element normals follow the right-hand rule on the node order. The skin is first made
consistent: walking the sheet across shared edges, a neighbour that runs a shared
edge the same way is reversed. The sign of the whole sheet then follows from the
divergence theorem. For the field $\mathbf{F} = (x - c_x,\ y - c_y,\ 0)$,
$\nabla\cdot\mathbf{F} = 2$; on the end faces of a bay, whose normals are along $z$,
$\mathbf{F}\cdot\mathbf{n} = 0$. So for outward normals the lateral surface of a
bay of volume $V$ gives, with $\mathbf{A}_e$ the vector area of element $e$ and
$\bar{\mathbf{x}}_e$ its vertex mean,

\begin{equation*}
\Phi = \sum_e \mathbf{A}_e\cdot(\bar{\mathbf{x}}_e - \mathbf{c})_{xy}
\approx \oint \mathbf{F}\cdot\mathbf{n}\,\mathrm{d}A = 2V > 0,
\tag{R9.13}
\end{equation*}

whatever the concavity of the section. The constant $\mathbf{c}$ (any point of the
station) cancels, since the in-plane vector areas of a closed band sum to zero. The
sheet is reversed if $\Phi < 0$. A web element is reversed if
$\mathbf{A}_e\cdot(\mathbf{x}_{\mathrm{TE}} - \bar{\mathbf{x}}_e) < 0$, so it faces
the TE; a rib or tip triangle is reversed if $A_{e,z} < 0$, so it faces $+z$.

### 3.10 Hard invariants and the ruled-strip mesher

A finished mesh satisfies, for every mesh the kernel means to build
(`check_shell_mesh`, `shell_mesh(strict=True)`): one entry per element in every
per-element field; at least three distinct nodes per element, all in range; finite
nodes and data with positive thickness; every node used and no two closer than
$10^{-9}$ of the bounding-box diagonal; no duplicate or zero-area element; and, for
each surface kind, every edge $\{a, b\}$ shared by two elements is run in opposite
directions:

\begin{equation*}
(a \to b) \in e_1 \ \Rightarrow\ (b \to a) \in e_2 ,
\tag{R9.14}
\end{equation*}

while an edge used by one element lies on the root station or on a web.

The ruled-strip mesher (`method="ruled"`) is used when the mesh is sized by
`n_span`, $h/c$, `region_h` or `refine`, or when the skin loop's partitions change
along the span. It places stations by (R9.1)–(R9.3) and adds a station at the middle
of any bay whose leading or trailing edge bows past the budget: with $\mathbf{a}$,
$\mathbf{b}$ the edge points at the bay ends and $\mathbf{m}$ the one at mid-bay,

\begin{equation*}
\Big|\mathbf{m} - \mathbf{a} - \frac{(\mathbf{m} - \mathbf{a})\cdot(\mathbf{b} - \mathbf{a})}{|\mathbf{b} - \mathbf{a}|^2}(\mathbf{b} - \mathbf{a})\Big|
> \varepsilon(z_m),
\tag{R9.15}
\end{equation*}

for up to eight rounds, unless the bay is shorter than $s/2$. It then lofts named
region polylines with mapped counts through topology events. Its strip
construction is not derived here; it is kept for one release after the switch
(ADR 0015).

## 4. Algorithm

`shell_mesh(geo, h=h)` with the partition mesher:

1. Stations by (R9.1)–(R9.3), plus every rib; the cut's shell graph at each; the
   span split by (R9.10).
2. Check that the loop has the same partitions at every station, is closed and has
   no branches; otherwise fall back to the ruled-strip mesher.
3. Build the ring of each station from the TE partition boundary, with its graded
   edge measure (R9.4)–(R9.7).
4. Collect the mesh lines: web feet (from the web's skeleton junctions, or projected
   from a section just inside the span where a web ends on a station), partition
   boundaries with no foot, band edges; hold absent feet by (R9.9). Check order
   and meeting at every station.
5. Count each ring segment and web segment by (R9.8); add edges until the loop has
   at least 8.
6. Place nodes; join rows into quads with data (R9.11).
7. Fill rib and tip cells by Section 3.8.
8. Orient by Section 3.9; measure quality; assemble the edge and corner tables.

## 5. Implementation

| Result | Function |
|---|---|
| (R9.1) | `mesh.shell.stations._target_edge`, `_size_fraction`, `TIP_H_OVER_C` |
| (R9.2) | `mesh.shell.stations._n_bays`, `_fill_gaps`, `_anchor_stations`, `mesh_stations` |
| (R9.3) | `mesh.shell.stations._drop_stations` |
| (R9.4)–(R9.7) | `mesh.shell.partition_mesh._window_kappa`, `_edge_measure`; `blade.topology.sizing.EPS_OVER_C` |
| (R9.8) | `mesh.shell.partition_mesh._ts_between`, `_StationLoop`, and the counts in `partition_mesh` |
| (R9.9) | `mesh.shell.partition_mesh._held_u`, `_ring_breaks` |
| (R9.10) | `mesh.shell.partition_mesh._grade_span`, `_mean_loop_edge` |
| (R9.11) | `mesh.shell.partition_mesh.partition_mesh` (element rows) |
| (R9.12) | `mesh.fill.delaunay.triangulate_ring`, `_unstructured_seeds`, `_smooth_interior`, `_recover_boundary`; `partition_mesh._cells` |
| (R9.13) | `mesh.shell.orient._oriented` |
| (R9.14) | `mesh.shell.invariants.check_shell_mesh`, `assert_shell_mesh` |
| (R9.15) | `mesh.shell.spanwise._span_edges`, `_sagitta_pending`, `_freeze_sagitta`; `mesh.shell.shell_mesh._shell_mesh_body` |

## 6. Limits

1. **Same loop everywhere.** The partition mesher needs the skin loop's partitions
   to be the same at every station, a closed loop and no branches (stringers); a
   web may start and stop once. Other blades use the ruled-strip mesher.
2. **Span sagitta.** (R9.15) is applied by the ruled-strip mesher only. The
   partition mesher grades the span by the chordwise edge (R9.10), so a leading or
   trailing edge that curves strongly along the span between anchors is not checked
   against $\varepsilon$.
3. **Partition ends are not refined.** Curvature is measured on each partition
   alone (R9.5), so the mesh is not refined towards a corner. The corner is a mesh
   line, so the geometry is exact there; a consumer that wants a finer mesh near a
   sharp TE passes a smaller $h$.
4. **Edge bound.** "At most the local edge" in (R9.8) holds for the trapezoidal
   measure (R9.7), not exactly for the true integral.
5. **Rib nodes by coordinate.** A rib triangle's corners are matched to the cell's
   nodes by coordinates rounded to $10^{-12}$ m after the fill's projection to the
   cell plane and back. A miss would add a node coincident with a ring node; the
   strict invariants (R9.14) would report it. None has been seen on blades.
6. **Seeds** use a fixed random seed, so the rib fill is reproducible but depends on
   the cell's bounding box.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R9.1), (R9.2) | Default edge $0.05\,c$; bays for five gaps; every anchor kept; largest bay $0.99\,h$ | [blade example 20](../../examples/blade/20_shell_mesh_checks/README.md), check 1; `tests/test_sizing.py` |
| (R9.3) | A web leaving the section found to $10^{-6}$ of the span | example 20, check 2 |
| (R9.4) | Chord-to-arc height over $\varepsilon$: 0.988 to 0.999 | example 20, check 3 |
| (R9.5) | Circle: $\kappa R - 1$ within $3\times10^{-6}$; a kink gives angle over window | example 20, check 4 |
| (R9.6), (R9.7) | Graded edge $\le$ raw edge and 0.25-Lipschitz; measure density equals (R9.7) | example 20, check 5; Figure R9.1 |
| (R9.8), (R9.9) | $W = 1/0.15$ on a 1 m run, $n = 7$, equal steps; a held line at its fraction | example 20, check 6 |
| (R9.10), (R9.11), (R9.14) | 77 stations, all skin and web quads, web feet shared with the skin, no long bay, no invariant broken | example 20, check 7; `tests/test_partition_mesh.py`, `tests/test_mesh_invariants.py` |
| (R9.13) | Skin flux 0.997 of $2V$ (inscribed rings); webs face the TE; ribs face $+z$ | example 20, check 8 |
| (R9.12) | Rib triangles cover their ring exactly; $V - E + F = 1$; seeds $\ge r$ apart and $\ge 0.62\,r$ from the boundary | example 20, check 9 |

```{figure} ../../examples/blade/20_shell_mesh_checks/figures/shell_mesh_checks.png
:alt: Top left, the chordwise edge round the root skin loop, with the sagitta edge and the graded edge dipping at the nose; top right, the root station's skin and web nodes with the shared web feet; bottom left, suction-side skin quads from root to mid-span in planform with the rib station; bottom right, a rib fill of an ellipse with its seeds
:width: 100%

**Figure R9.1.** Top left: the edge round the root skin loop of a 2D-F blade at
$h = 0.08$ m: $\min(h, \ell_{\mathrm{sag}})$ (R9.4) and the graded edge (R9.6);
grey lines are partition ends. Top right: the root station's nodes; the web feet are
skin nodes. Bottom left: suction-side skin quads from the root to mid-span, with the
rib anchor. Bottom right: the unstructured fill (R9.12) of an ellipse at
$r = 0.06$ with its seeds. Drawn by
`examples/blade/20_shell_mesh_checks/shell_mesh_checks.py`.
```

These checks move to V1 when it is written.

## 8. References

The mid-surfaces are those of [R6](R6_midsurfaces_cells.md) and the members those of
[R8](R8_loft_members.md). Delaunay triangulation is by the Bowyer–Watson algorithm
{cite:p}`bowyer1981,watson1981` (or Qhull when SciPy is present); the seed spacing
follows Poisson-disk sampling {cite:p}`bridson2007`. Full bibliography entries are
in [R0](R0_notation.md#bibliography).
