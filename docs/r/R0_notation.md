# R0 · Notation and bibliography

This document fixes the symbols, units and typographic rules used by every R
document, and holds the bibliography for all four manuals. A symbol defined here
keeps its meaning everywhere; a document that needs a new symbol adds it to the
table in the same change.

## Conventions

1. **Units.** SI throughout: lengths in metres, areas in m², second moments of area
   in m⁴, forces in N and moments in N·m. Angles are radians in equations. Public
   functions whose argument ends in `_deg` take degrees and convert on entry.
2. **Vectors and matrices.** Bold lower case for vectors ($\mathbf{x}$), bold
   upper case for matrices ($\mathbf{Q}$). A subscript names the frame in which
   components are taken: $\mathbf{x}_S$ is the column of components of
   $\mathbf{x}$ in frame S.
3. **Component order.** In-plane section arrays are stored as $(X, Y)$: the
   flapwise component first, the chordwise component second (see [R1](R1_frames.md)).
   An array stored the other way round is called *swapped*.
4. **Equation numbers.** Equations are numbered within their document and carry the
   document key: (R1.4) is the fourth numbered equation of R1. Other documents cite
   it by that number.
5. **Figures.** Figures are numbered the same way (Figure R1.1). Each caption states
   what is shown and names the script that draws it.
6. **Citations.** Author and year, for example {cite:t}`megson2017`. Full entries are
   listed at the end of this document.

## Symbols

| Symbol | Meaning | Unit | Defined in |
|---|---|---|---|
| $\mathrm{R}, \mathrm{B}, \mathrm{S}, \mathrm{S}', \mathrm{C}$ | Rotor, blade, station, centroidal station and contour frames | – | [R1](R1_frames.md) |
| $\mathbf{e}_X, \mathbf{e}_Y, \mathbf{e}_Z$ | Unit basis vectors; a frame subscript names their frame | – | R1 |
| $X, Y$ | Flapwise and chordwise in-plane coordinates | m | R1 |
| $z$ | Spanwise coordinate from the blade root | m | R1 |
| $L$ | Blade length, root to tip | m | R1 |
| $r$, $R$ | Radial position from the rotor axis; rotor radius | m | R1 |
| $c$ | Chord | m | R1 |
| $f$ | Pitch-axis position as a fraction of chord from the leading edge | – | R1 |
| $x_c, y_t$ | Aerofoil chordwise and thickness-direction coordinates, leading edge at the origin | m | R1 |
| $\theta$ | Frame twist angle, the rotation from S to B | rad | R1 |
| $\tau$ | `BladeGeometry` twist, positive nose up; $\tau = -\theta$ | rad | R1 |
| $\mathbf{Q}(\theta)$ | In-plane rotation matrix from S components to B components | – | (R1.4) |
| $\mathbf{P}$ | Component swap $(X, Y) \mapsto (Y, X)$ | – | (R1.8) |
| $\bar{X}, \bar{Y}$ | Centroid of a section in frame S | m | R1 |
| $A$ | Section area | m² | R1 |
| $I_{xx}, I_{yy}, I_{xy}$ | Second moments $\int Y^2\,\mathrm{d}A$, $\int X^2\,\mathrm{d}A$ and $\int XY\,\mathrm{d}A$ | m⁴ | (R1.17) |
| $\mathbf{V}$ | In-plane force (shear) vector | N | R1 |
| $\mathbf{M}$ | In-plane bending moment vector | N·m | R1 |
| $s$, $n$ | Arc length along a mid-line; distance along the shell normal | m | R1 |

## Bibliography

```{bibliography}
:all:
```
