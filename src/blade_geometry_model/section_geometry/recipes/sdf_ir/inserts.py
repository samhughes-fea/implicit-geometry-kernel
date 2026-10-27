"""Leading- and trailing-edge insert subcomponents."""

from __future__ import annotations

from blade_geometry_model.section_geometry.sdf.ir import Circle, HalfPlane, Intersect, Offset, Subtract, subtract_expr

from blade_geometry_model.section_geometry.sections.clip_algebra import (
    airfoil_expr,
)
from blade_geometry_model.section_geometry.sections.labels import (
    ComponentLabel,
)
from .base import SectionComponent


class TEInsert(SectionComponent):
    """Trailing-edge reinforcement zone.

    Occupies the cavity side of the outer skin: ``Offset(airfoil, -t_skin)``,
    flush with the inner-skin / core interface (not the skin mid-surface).
    """

    label = ComponentLabel.TRAILING_EDGE_INSERT

    def __init__(
        self,
        airfoil_sdf,
        skin_thickness,
        x_start=0.75,
        exclusion_sdfs=None,
        *,
        airfoil_expr_node=None,
        exclusion_exprs=None,
    ):
        self.x_start = float(x_start)
        self.skin_thickness = float(skin_thickness)
        self._af = airfoil_sdf

        poly = airfoil_expr_node if airfoil_expr_node is not None else airfoil_expr(airfoil_sdf)
        if poly is None:
            raise ValueError("TEInsert requires an AirfoilSDF with vertices for IR evaluation.")

        expr = Intersect(
            Offset(poly, amount=-float(skin_thickness)),
            HalfPlane(nx=-1.0, ny=0.0, d=float(x_start)),
        )
        cutters = []
        if exclusion_exprs is not None:
            cutters = list(exclusion_exprs)
        elif exclusion_sdfs:
            cutters = [getattr(ex, "expr", getattr(ex, "_expr", None)) for ex in exclusion_sdfs]
            cutters = [c for c in cutters if c is not None]
        if cutters:
            expr = subtract_expr(expr, cutters)
        super().__init__(expr)

    def __repr__(self):
        return f"TEInsert(x_start={self.x_start:.3g})"


class LEInsert(SectionComponent):
    """Leading-edge nose insert (solid glass or carbon nose block).

    Occupies the cavity side of the outer skin: ``Offset(airfoil, -t_skin)``,
    flush with the inner-skin / core interface (not the skin mid-surface).
    """

    label = ComponentLabel.LEADING_EDGE_INSERT

    def __init__(
        self,
        airfoil_sdf,
        skin_thickness,
        x_end=0.10,
        le_x=0.0,
        le_y=0.0,
        radius=None,
        *,
        airfoil_expr_node=None,
    ):
        self.x_end = float(x_end)
        self.skin_thickness = float(skin_thickness)
        self.le_x = float(le_x)
        self.le_y = float(le_y)
        self._af = airfoil_sdf

        if radius is None:
            radius = x_end
        self.radius = float(radius)

        poly = airfoil_expr_node if airfoil_expr_node is not None else airfoil_expr(airfoil_sdf)
        if poly is None:
            raise ValueError("LEInsert requires an AirfoilSDF with vertices for IR evaluation.")

        expr = Intersect(
            Intersect(
                Offset(poly, amount=-float(skin_thickness)),
                HalfPlane(nx=1.0, ny=0.0, d=-float(x_end)),
            ),
            Circle(cx=float(le_x), cy=float(le_y), r=float(radius)),
        )
        super().__init__(expr)

    def __repr__(self):
        return f"LEInsert(x_end={self.x_end:.3g})"
