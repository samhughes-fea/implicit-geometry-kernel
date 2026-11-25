# bug_log — integration

**Audience: internal.** Maintainer notes; not part of the library manual.

↑ [`integration`](README.md) · [R7](../r/R7_section_properties.md)

## Role

Package debug log for `implicit_geometry_kernel.integration` (shoelace and quadtree
quadrature). Identifiers are `IN-nnn`. Cross-link
[`section/bug_log.md`](../section/bug_log.md) (`SG-`),
[`blade/bug_log.md`](../blade/bug_log.md) (`BG-`),
[`frames/bug_log.md`](../frames/bug_log.md) (`CF-`) and
[`sdf/bug_log.md`](../sdf/bug_log.md) (`SD-`) when a defect spans packages.

## Index

| ID | Title | Status |
|---|---|---|
| [IN-001](#in-001) | Quadtree interior leaves carry the wrong product moment | Fixed (#124) |

## IN-001

**Symptom.** A circle of radius 0.2 centred at $(0.2, 0.1)$ has $I_{xy} = 0$ about
its centroid; the quadtree backend returns $-5.8\times10^{-4}$, against
$I_{xx} = 1.26\times10^{-3}$. Area, first moments, $I_{xx}$ and $I_{yy}$ converge as $h^2$; $I_{xy}$
does not converge to zero.

**Evidence.** Found while deriving [R7](../r/R7_section_properties.md), (R7.9).
`examples/integration/02_property_quadrature` check 5 evaluates the leaf $[1, 2]^2$: $\int xy\,dA = 1.75$ against the exact 2.25; check 7
prints the circle case. Existing tests use sections centred on the origin or
symmetric about an axis, where the error cancels.

**Root cause.** `_rectangle_moments` and `_rectangle_moments_batch` in
`integration/quadtree/moments.py` evaluate
$A\,(x_0 y_0 + \tfrac12(x_0 y_1 + x_1 y_0) + x_1 y_1)/4$. The integral of $xy$ over a
rectangle factors, $\int_{x_0}^{x_1}x\,dx \int_{y_0}^{y_1}y\,dy$, so the exact value
is $A\,(x_0 + x_1)(y_0 + y_1)/4$ (R7.9): the cross terms carry weight 1, not
$\tfrac12$. The error is $-A\,(x_0 y_1 + x_1 y_0)/8$ per interior leaf, so it grows
with the distance of the section from the origin.

**Scope.** The quadtree backend and the residual labels (LE and TE inserts) that the
hybrid backend integrates with the quadtree. The hybrid wall properties use the
shoelace (R7.1)–(R7.6) and are not affected.

**Fix (#124).** Both functions evaluate $A\,(x_0 + x_1)(y_0 + y_1)/4$. The circle
now gives $I'_{xy} = -1.5\times10^{-16}$; `tests/integration/test_quadtree_quadrature.py`
checks one leaf against the factored integral and an off-axis circle against zero.
