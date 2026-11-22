# R1 · Coordinate frames and twist

## 1. Problem

A blade section is authored in aerofoil coordinates, cut and meshed in a frame
attached to the blade, and integrated about axes through the pitch axis or the
centroid. This document defines those frames and derives every map between them
that the `frames` package evaluates: points, force and moment vectors, second
moments of area, a wall depth, and the spanwise coordinate. It also relates the
`frames` twist to the twist of `BladeGeometry`, which has the opposite sign.

## 2. Definitions and assumptions

1. **Aerofoil coordinates** $(x_c, y_t)$: $x_c$ along the chord from the leading
   edge (LE) towards the trailing edge (TE), $y_t$ normal to the chord towards the
   suction side, both scaled by the chord $c$ {cite:p}`abbott1959`. They are not a
   structural frame.
2. **Station frame S**: origin on the pitch axis, which lies on the chord line at
   $x_c = f c$, with $f$ the pitch-axis fraction (`pivot_fraction`). The unit vector
   $\mathbf{e}_{X_S}$ is normal to the chord towards the suction side (flapwise),
   $\mathbf{e}_{Y_S}$ lies along the chord towards the TE (chordwise), and
   $\mathbf{e}_{Z}$ points along the span from root to tip. S carries no twist: the
   section always appears chord-aligned in it.
3. **Centroidal station frame S′**: S translated so that its origin is the centroid
   of the section; same axes as S.
4. **Blade frame B**: origin on the pitch axis, $\mathbf{e}_{X_B}$ flapwise and
   $\mathbf{e}_{Y_B}$ edgewise, fixed to the blade and independent of the station.
   The section appears twisted in B.
5. **Twist** $\theta$: the angle that takes the axes of B onto the axes of S at a
   station, defined by (R1.3). `BladeGeometry` uses $\tau=-\theta$ (Section 3.5).
6. **Rotor frame R**: the radial coordinate $r$ from the rotor axis and the rotor
   radius $R$. Only the spanwise coordinate is mapped.
7. **Contour frame C**: at a point of a wall mid-line, the arc-length direction
   $\mathbf{e}_s$, the shell normal $\mathbf{e}_n$ and the span direction
   $\mathbf{e}_Z$.
8. **Component order.** In-plane arrays store $(X, Y)$, flapwise first. An array in
   $(Y, X)$ order is *swapped*.
9. **Assumption A1.** Twist is a rotation about $\mathbf{e}_Z$ only: the blade axis
   is straight, with no precone, prebend or sweep. All maps below are therefore
   two-dimensional in the section plane, and $Z_S = Z_B = z$.

Figure R1.1 shows the frames on one station.

```{figure} ../../examples/frames/02_frame_triad/figures/frame_triad.png
:alt: A NACA 2418 section twisted by 20 degrees with the axes of B, S, S′ and C drawn on it
:width: 100%

**Figure R1.1.** Frames B (black), S (blue), S′ (green, dashed, at the centroid) and
C (orange, at a point on the suction side) on a station twisted by
$\theta = 20°$. Horizontal $Y_B$, vertical $X_B$; the horizontal axis is
inverted for display so that the LE is on the right. Drawn by
`examples/frames/02_frame_triad/frame_triad.py`.
```

## 3. Derivation

### 3.1 Aerofoil coordinates to S

The pitch axis is the point $(x_c, y_t) = (f c, 0)$. Translating the origin there
and then listing the flapwise component first gives

\begin{equation*}
\begin{pmatrix} X_S \\ Y_S \end{pmatrix}
= \mathbf{P} \begin{pmatrix} x_c - f c \\ y_t \end{pmatrix},
\qquad\text{that is}\qquad
X_S = y_t,\quad Y_S = x_c - f c,
\tag{R1.1}
\end{equation*}

where $\mathbf{P}$ is the component swap of (R1.8). The LE therefore lies at
$Y_S = -f c < 0$ and the TE at $Y_S = (1-f)c > 0$. Because $\mathbf{P}^2=\mathbf{I}$,
the inverse is

\begin{equation*}
x_c = Y_S + f c,\qquad y_t = X_S .
\tag{R1.2}
\end{equation*}

### 3.2 S to B: the twist rotation

Define the twist by the basis of S written in the basis of B:

\begin{equation*}
\mathbf{e}_{X_S} = \cos\theta\,\mathbf{e}_{X_B} - \sin\theta\,\mathbf{e}_{Y_B},
\qquad
\mathbf{e}_{Y_S} = \sin\theta\,\mathbf{e}_{X_B} + \cos\theta\,\mathbf{e}_{Y_B}.
\tag{R1.3}
\end{equation*}

A point is $\mathbf{x} = X_S\mathbf{e}_{X_S} + Y_S\mathbf{e}_{Y_S}$. Substituting
(R1.3) and collecting the coefficients of $\mathbf{e}_{X_B}$ and $\mathbf{e}_{Y_B}$:

\begin{equation*}
\mathbf{x}_B = \mathbf{Q}(\theta)\,\mathbf{x}_S,
\qquad
\mathbf{Q}(\theta) =
\begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix},
\tag{R1.4}
\end{equation*}

so $X_B = X_S\cos\theta + Y_S\sin\theta$ and
$Y_B = -X_S\sin\theta + Y_S\cos\theta$. A point on the chord, $X_S = 0$, has
$X_B = Y_S\sin\theta$: for $\theta>0$ the LE ($Y_S<0$) moves to negative $X_B$,
towards the pressure side. Direct multiplication gives

\begin{equation*}
\mathbf{Q}(\theta)^{\mathsf{T}}\mathbf{Q}(\theta) = \mathbf{I},
\qquad \det\mathbf{Q}(\theta) = 1,
\tag{R1.5}
\end{equation*}

so $\mathbf{Q}$ is a proper rotation, lengths and areas are preserved, and the
inverse map is the transpose:

\begin{equation*}
\mathbf{x}_S = \mathbf{Q}(\theta)^{\mathsf{T}}\mathbf{x}_B = \mathbf{Q}(-\theta)\,\mathbf{x}_B .
\tag{R1.6}
\end{equation*}

### 3.3 Aerofoil coordinates to B

Composing (R1.1) and (R1.4):

\begin{equation*}
\mathbf{x}_B = \mathbf{Q}(\theta)\,\mathbf{P}
\begin{pmatrix} x_c - f c \\ y_t \end{pmatrix},
\tag{R1.7}
\end{equation*}

and the inverse applies (R1.6) and then (R1.2).

### 3.4 Swapped component order

Let

\begin{equation*}
\mathbf{P} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},
\qquad \mathbf{P}^2 = \mathbf{I}.
\tag{R1.8}
\end{equation*}

A swapped array holds $\mathbf{P}\mathbf{x}$. Rotating it correctly means swapping
back, rotating and swapping again, $\mathbf{P}\mathbf{Q}(\theta)\mathbf{P}$. Multiplying
out,

\begin{equation*}
\mathbf{P}\mathbf{Q}(\theta)\mathbf{P}
= \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix},
\tag{R1.9}
\end{equation*}

which is (R1.4) with the sign of $\theta$ reversed:

\begin{equation*}
\mathbf{P}\mathbf{Q}(\theta)\mathbf{P} = \mathbf{Q}(-\theta).
\tag{R1.10}
\end{equation*}

Applying $\mathbf{Q}(\theta)$ to a swapped array therefore rotates it by the opposite
twist. This is the most common sign error with these frames; it was the cause of
CF-004, fixed in 0.9.6 (Section 6).

### 3.5 The `BladeGeometry` cut frame and its twist

A `SectionCut` stores points as $(x, y, z)$: $x$ along the chord towards the TE,
$y$ towards the suction side, origin on the pitch axis. In the notation above this
is the swapped B order,

\begin{equation*}
\begin{pmatrix} x \\ y \end{pmatrix} = \mathbf{P}\,\mathbf{x}_B,
\qquad \tau = -\theta,
\tag{R1.11}
\end{equation*}

where $\tau$ is the `twist` given to `Station`, positive nose up under the `BladeGeometry` default
`twist_sense="nose_up"` (with `"nose_down"` the sign of the input is reversed
before this point), and the loft passes $-\tau$ to the section geometry. Let $\mathbf{x}_0 = \mathbf{P}\mathbf{x}_S$ be the cut
of the same station at zero twist. Then
$(x, y)^{\mathsf{T}} = \mathbf{P}\mathbf{Q}(\theta)\mathbf{x}_S = \mathbf{P}\mathbf{Q}(\theta)\mathbf{P}\,\mathbf{x}_0$,
and by (R1.10) with $\theta=-\tau$:

\begin{equation*}
\begin{pmatrix} x \\ y \end{pmatrix} = \mathbf{Q}(\tau)\,\mathbf{x}_0 ,
\tag{R1.12}
\end{equation*}

that is $x = x_0\cos\tau + y_0\sin\tau$ and $y = -x_0\sin\tau + y_0\cos\tau$.

For $\tau>0$ the LE ($x_0<0$, $y_0=0$) moves to $y = -x_0\sin\tau > 0$: nose up,
as stated. Users of `BladeGeometry` need only (R1.12); (R1.4) matters when section
properties or loads move between S and B.

### 3.6 Span coordinate and the rotor frame

Under A1 the blade axis is a radial line. With the root at $r = R - L$ and $z$
measured from the root,

\begin{equation*}
r = z + R - L,
\qquad
L = R\left(1 - \frac{r_{\mathrm{root}}}{R}\right),
\tag{R1.13}
\end{equation*}

where $r_{\mathrm{root}}/R$ is the smallest $r/R$ of the supplied radial grids, so that
$r/R = (z + R - L)/R$ at every station.

### 3.7 Force vectors

A force (shear) $\mathbf{V} = V_{X_S}\mathbf{e}_{X_S} + V_{Y_S}\mathbf{e}_{Y_S}$ is expanded
in the basis of B exactly as the point in Section 3.2, so its components obey the
same map:

\begin{equation*}
\mathbf{V}_B = \mathbf{Q}(\theta)\,\mathbf{V}_S .
\tag{R1.14}
\end{equation*}

### 3.8 Bending moment vectors

An axial force $F_z$ acting at the in-plane point $\mathbf{x}$ has the in-plane moment
$\mathbf{x}\times F_z\mathbf{e}_Z$. Its components are linear in $\mathbf{x}$:

\begin{equation*}
\mathbf{M} = F_z\,\mathbf{K}\,\mathbf{x},
\qquad
\mathbf{K} = \pm\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix},
\tag{R1.15}
\end{equation*}

that is $M_X = \pm Y F_z$ and $M_Y = \mp X F_z$, the sign fixed by the handedness of
$(\mathbf{e}_X, \mathbf{e}_Y, \mathbf{e}_Z)$. $\mathbf{K}$ is itself a rotation by a
quarter turn, and plane rotations commute, so $\mathbf{K}\mathbf{Q}=\mathbf{Q}\mathbf{K}$.
Substituting (R1.4) into (R1.15) evaluated in B, and summing over any distribution
of axial force:

\begin{equation*}
\mathbf{M}_B = F_z\mathbf{K}\mathbf{Q}(\theta)\mathbf{x}_S = \mathbf{Q}(\theta)\,F_z\mathbf{K}\mathbf{x}_S
= \mathbf{Q}(\theta)\,\mathbf{M}_S .
\tag{R1.16}
\end{equation*}

Moments, written in $(M_X, M_Y)$ order, rotate exactly like points and forces,
whichever the handedness.

### 3.9 Second moments of area

Define, over the section area $A$,

\begin{equation*}
\begin{gathered}
I_{xx} = \int_A Y^2\,\mathrm{d}A,\qquad
I_{yy} = \int_A X^2\,\mathrm{d}A,\qquad
I_{xy} = \int_A XY\,\mathrm{d}A,
\\
\mathbf{J} = \int_A \mathbf{x}\mathbf{x}^{\mathsf{T}}\,\mathrm{d}A
= \begin{pmatrix} I_{yy} & I_{xy} \\ I_{xy} & I_{xx} \end{pmatrix},
\end{gathered}
\tag{R1.17}
\end{equation*}

about the origin of the frame in which $\mathbf{x}$ is taken (the pitch axis for S and
B). This follows the naming of {cite:t}`megson2017`: $I_{xx}$ is the moment about the
$X$ axis. Substituting $\mathbf{x}_B = \mathbf{Q}\mathbf{x}_S$ and taking the constant
$\mathbf{Q}$ outside the integral, with $\mathrm{d}A$ unchanged by (R1.5):

\begin{equation*}
\mathbf{J}_B = \mathbf{Q}(\theta)\,\mathbf{J}_S\,\mathbf{Q}(\theta)^{\mathsf{T}} .
\tag{R1.18}
\end{equation*}

Writing $c=\cos\theta$, $s=\sin\theta$ and multiplying out the three distinct entries
of (R1.18):

\begin{equation*}
I_{xx,B} = I_{xx,S}\,c^2 + I_{yy,S}\,s^2 - 2 I_{xy,S}\,s c,
\tag{R1.19}
\end{equation*}

\begin{equation*}
I_{yy,B} = I_{xx,S}\,s^2 + I_{yy,S}\,c^2 + 2 I_{xy,S}\,s c,
\tag{R1.20}
\end{equation*}

\begin{equation*}
I_{xy,B} = (I_{xx,S} - I_{yy,S})\,s c + I_{xy,S}\,(c^2 - s^2).
\tag{R1.21}
\end{equation*}

The trace and determinant of $\mathbf{J}$ are unchanged by the similarity (R1.18),
which gives two checks that hold at every twist {cite:p}`pilkey2002`:

\begin{equation*}
\begin{gathered}
I_{xx,B} + I_{yy,B} = I_{xx,S} + I_{yy,S},
\\
I_{xx,B}I_{yy,B} - I_{xy,B}^2 = I_{xx,S}I_{yy,S} - I_{xy,S}^2 .
\end{gathered}
\tag{R1.22}
\end{equation*}

The inverse map is (R1.18) with $\mathbf{Q}(-\theta)$.

### 3.10 S to S′: translation to the centroid

With the centroid $\bar{\mathbf{x}} = (\bar{X}, \bar{Y}) = A^{-1}\int_A \mathbf{x}_S\,\mathrm{d}A$,

\begin{equation*}
\mathbf{x}_{S'} = \mathbf{x}_S - \bar{\mathbf{x}} .
\tag{R1.23}
\end{equation*}

Substituting $\mathbf{x}_S = \mathbf{x}_{S'} + \bar{\mathbf{x}}$ into (R1.17), expanding, and
using $\int_A \mathbf{x}_{S'}\,\mathrm{d}A = \mathbf{0}$ (the definition of the centroid)
to remove the cross terms gives the parallel-axis relations
{cite:p}`megson2017`:

\begin{equation*}
\begin{gathered}
I_{xx,S} = I_{xx,S'} + A\bar{Y}^2,\qquad
I_{yy,S} = I_{yy,S'} + A\bar{X}^2,
\\
I_{xy,S} = I_{xy,S'} + A\bar{X}\bar{Y}.
\end{gathered}
\tag{R1.24}
\end{equation*}

### 3.11 Wall depth between B and S′

A straight segment of length $d$ along $\mathbf{e}_{X_S}$ has, by (R1.4), an extent
along $\mathbf{e}_{X_B}$ of $d_B = d\cos\theta$. Hence

\begin{equation*}
d_{S'} = \frac{d_B}{\cos\theta},\qquad d_B = d_{S'}\cos\theta .
\tag{R1.25}
\end{equation*}

This is a projection, exact only for a segment along $\mathbf{e}_{X_S}$; it is singular
at $\theta = \pm 90°$, where the implementation substitutes 1 for
$|\cos\theta| < 10^{-10}$. No kernel geometry uses it.

### 3.12 Contour frame C

Let $\mathbf{x}_B(s)$ be a wall mid-line parametrised by arc length $s$, in the order in
which its points are stored. Then

\begin{equation*}
\mathbf{e}_s = \frac{\mathrm{d}\mathbf{x}_B}{\mathrm{d}s},\qquad
\mathbf{e}_n \perp \mathbf{e}_s,\qquad
\lVert\mathbf{e}_s\rVert = \lVert\mathbf{e}_n\rVert = 1,
\tag{R1.26}
\end{equation*}

with $\mathbf{e}_n$ pointing into the section on the skin (the sense in which a spar
cap's eccentricity $t_{\mathrm{cap}}/2$ is measured by `LaminateOffsetRegion`). C is a
local frame for shell degrees of freedom; `frames` provides its tag `FRAME_C` but no
transform.

## 4. Algorithm

Every map in Section 3 is a closed-form expression evaluated pointwise on arrays;
there is no iteration. Angles are converted from degrees to radians on entry where
the argument is `*_deg`.

## 5. Implementation

| Result | Function (`implicit_geometry_kernel.frames` unless stated) |
|---|---|
| (R1.1), (R1.2) | `naca_to_S`, `S_to_naca` (long names `naca_to_zero_twist_pitch_axis_station_frame`, …); `transform_points_aerofoil_to_S` |
| (R1.4), (R1.6) | `rotate_S_to_B`, `rotate_B_to_S`; `transform_points_S_to_B`, `transform_points_B_to_S`, `rotate_points_2d` |
| (R1.7) | `naca_to_twisted_pitch_axis_blade_frame` (`chord_to_blade`) and its inverse `blade_to_chord` |
| (R1.10) | `rotate_S_to_B_swapped_order` |
| (R1.12) | `BladeGeometry.section_at` (the loft applies it; no separate function) |
| (R1.13) | `rR_from_zL_scalars`; `CoordinateFrameTransforms(R=…, rR_grids=…)` for grids and named loci |
| (R1.14) | `CoordinateFrameTransforms.transform_forces_S_to_B`, `transform_forces_B_to_S` |
| (R1.16) | `CoordinateFrameTransforms.transform_moments_S_to_B`, `transform_moments_B_to_S` (argument order $(M_Y, M_X)$) |
| (R1.19)–(R1.21) | `CoordinateFrameTransforms.transform_second_moments_S_to_B`, `transform_second_moments_B_to_S` |
| (R1.23) | `transform_points_S_to_Sprime`, `transform_points_Sprime_to_S` |
| (R1.25) | `CoordinateFrameTransforms.transform_thickness_B_to_Sprime`, `transform_thickness_Sprime_to_B` |
| Plot convention | `implicit_geometry_kernel.plotting.frames.section_plot_axes`, `setup_section_axes` |

## 6. Limits

1. **Argument order of the moment transforms.** They take and return
   $(M_Y, M_X)$, not $(M_X, M_Y)$. Before 0.9.6 they applied the $(X, Y)$ matrix to
   that order, which by (R1.10) rotated moments by $-\theta$ (CF-004, fixed; see the
   [frames bug log](../frames/bug_log.md#cf-004)).
2. **Two twist signs.** `frames` and `BladeGeometry` use opposite signs, (R1.11).
   Pass $\theta = -\tau$ when moving a cut's properties with the `frames` functions.
3. **Planar only.** Assumption A1 excludes precone, prebend and sweep.
4. **Depth projection.** (R1.25) holds only for a segment along $\mathbf{e}_{X_S}$.
5. **Out of scope.** `structural_product_moment` closes a product moment from two
   stiffness ratios; it is not a frame transform and is not derived here. The
   laminate frames E and M keep their tags here, but their transforms belong to the
   constitutive package downstream.

## 7. Verification

| Result | Check | Where |
|---|---|---|
| (R1.1)–(R1.7) | Aerofoil → S → B → aerofoil returns every point to 6.7e-16 m for twists from −90° to 90° | [frames example 01](../../examples/frames/01_frames_and_twist/README.md), check 1 |
| (R1.10) | $\mathbf{P}\mathbf{Q}(\theta)\mathbf{P}=\mathbf{Q}(-\theta)$ exactly | [frames example 02](../../examples/frames/02_frame_triad/README.md), check 1 |
| (R1.12) | Twisted cuts equal the rotated untwisted cut to 6.7e-16 m; positive $\tau$ is nose up | example 02, check 2; example 01, check 3 |
| (R1.14) | Round trip and magnitude preservation | `tests/frames/unit/test_coordinate_force_transforms.py` |
| (R1.16) | Moments from rotated points equal the transform to 1e-12 at eight twists; force, moment and point rotate alike (CF-004) | `tests/frames/unit/test_moments_rotate_with_points.py` |
| (R1.19)–(R1.22) | Transform equals integration of the rotated outline to 1.4e-15 of $I_{xx}+I_{yy}$; invariants to 1e-15 | example 01, check 2; example 02, check 4; `tests/frames/unit/test_second_moments_rotate_with_points.py` |
| (R1.24) | Parallel-axis residual 2e-16 of $I_{xx}+I_{yy}$ | example 02, check 3 |

```{figure} ../../examples/frames/01_frames_and_twist/figures/frames_and_twist.png
:alt: Second moments of a NACA 2418 section against twist, transform and direct integration overlaid
:width: 100%

**Figure R1.2.** Top: one NACA 2418 outline in aerofoil coordinates, S and B
($\theta = 20°$). Bottom left: (R1.19)–(R1.21) against second moments integrated
from the rotated outline, twist −90° to 90°. Bottom right: chord-normal webs of
`BladeGeometry` at $\tau=\pm 20°$, showing positive $\tau$ nose up. Drawn by
`examples/frames/01_frames_and_twist/frames_and_twist.py`.
```

These checks move to [V1](../v/index.md) when it is written.

## 8. References

Cited above: {cite:t}`abbott1959` for aerofoil coordinates, {cite:t}`megson2017` for
the naming of second moments and the parallel-axis theorem, and
{cite:t}`pilkey2002` for the transformation and invariants of the second-moment
tensor. Full entries are in the [bibliography](R0_notation.md#bibliography).
