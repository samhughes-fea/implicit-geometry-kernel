# Shoelace integrals

↑ [integration](README.md)

`implicit_geometry_kernel.integration.shoelace` evaluates area, first moments,
and second moments of a simple closed polygon in the station plane by Green’s
theorem. Coordinates are \((X,Y)\) in \(\mathrm{S}\) (pitch origin) unless noted.
The implementation is `ring_sums_S`.

## Area

For vertices \((X_i,Y_i)\), \(i=1,\ldots,n\), with \((X_{n+1},Y_{n+1})=(X_1,Y_1)\)
and anti-clockwise orientation,

$$
A = \frac12\sum_{i=1}^{n}\bigl(X_i Y_{i+1}-X_{i+1} Y_i\bigr).
\tag{1}
$$

Holes enter with opposite orientation and subtract.

## First moments and centroid

$$
S_X = \int Y\,\mathrm{d}A,\qquad
S_Y = \int X\,\mathrm{d}A,\qquad
X_c = S_Y/A,\qquad
Y_c = S_X/A.
\tag{2}
$$

On the polyline the discrete forms used by `ring_sums_S` are

$$
S_Y = \frac16\sum_i (X_i+X_{i+1})\,c_i,\qquad
S_X = \frac16\sum_i (Y_i+Y_{i+1})\,c_i,
\tag{3}
$$

where \(c_i=X_i Y_{i+1}-X_{i+1}Y_i\).

## Second moments

About the \(\mathrm{S}\) origin,

$$
I_{xx}=\int Y^2\,\mathrm{d}A,\qquad
I_{yy}=\int X^2\,\mathrm{d}A,\qquad
I_{xy}=\int XY\,\mathrm{d}A.
\tag{4}
$$

Parallel-axis reduction to the centroid \(\mathrm{S}'\) is

$$
I_{xx}' = I_{xx}-A Y_c^2,\qquad
I_{yy}' = I_{yy}-A X_c^2,\qquad
I_{xy}' = I_{xy}-A X_c Y_c.
\tag{5}
$$

The same identities appear after quadtree accumulation
(`section_properties_from_quadtree`). Twist maps \((4)\) into \(\mathrm{B}\) via
[frames](../frames/README.md) equation \((3)\).
