"""Sandwich core subcomponent."""

from __future__ import annotations

from typing import Literal

from blade_geometry_model.section_geometry.sdf.ir import HalfPlane, Intersect, Offset, subtract_expr

from blade_geometry_model.section_geometry.sections.clip_algebra import (
    airfoil_expr,
    strip_clip_expr,
)
from blade_geometry_model.section_geometry.sections.labels import (
    ComponentLabel,
)
from .base import SectionComponent

BayMode = Literal["chordwise", "flapwise"]
FlapwiseClipMode = Literal["half_lt", "half_gt", "strip"]


class SandwichCore(SectionComponent):
    """Foam/balsa core infill: inside inner skin, outside all wall solids."""

    label = ComponentLabel.SANDWICH_CORE

    def __init__(
        self,
        airfoil_sdf,
        skin_thickness,
        exclusion_sdfs=None,
        x_start=None,
        x_end=None,
        y_min=None,
        y_max=None,
        strip_clip=None,
        *,
        airfoil_expr_node=None,
        exclusion_exprs=None,
        bay_mode: BayMode = "chordwise",
        flapwise_nx: float | None = None,
        flapwise_ny: float | None = None,
        flapwise_c_lo: float | None = None,
        flapwise_c_hi: float | None = None,
        flapwise_clip_mode: FlapwiseClipMode | None = None,
    ):
        poly = airfoil_expr_node if airfoil_expr_node is not None else airfoil_expr(airfoil_sdf)
        if poly is None:
            raise ValueError("SandwichCore requires an AirfoilSDF with vertices for IR evaluation.")

        expr = Offset(poly, amount=-skin_thickness)
        cutters = []
        if exclusion_exprs is not None:
            cutters = list(exclusion_exprs)
        elif exclusion_sdfs:
            cutters = [getattr(ex, "expr", getattr(ex, "_expr", None)) for ex in exclusion_sdfs]
            cutters = [c for c in cutters if c is not None]
        if cutters:
            expr = subtract_expr(expr, cutters)

        if strip_clip is not None:
            expr = Intersect(expr, strip_clip_expr(strip_clip, key=f"core_strip_{id(strip_clip)}"))
            bay_mode = "flapwise"
        else:
            if x_start is not None:
                expr = Intersect(expr, HalfPlane(nx=-1.0, ny=0.0, d=float(x_start)))
            if x_end is not None:
                expr = Intersect(expr, HalfPlane(nx=1.0, ny=0.0, d=-float(x_end)))
        if y_min is not None:
            expr = Intersect(expr, HalfPlane(nx=0.0, ny=-1.0, d=float(y_min)))
        if y_max is not None:
            expr = Intersect(expr, HalfPlane(nx=0.0, ny=1.0, d=-float(y_max)))

        super().__init__(expr)
        self._af = airfoil_sdf
        self.skin_thickness = float(skin_thickness)
        self.x_start = x_start
        self.x_end = x_end
        self.y_min = y_min
        self.y_max = y_max
        self.bay_mode: BayMode = bay_mode
        self.flapwise_nx = flapwise_nx
        self.flapwise_ny = flapwise_ny
        self.flapwise_c_lo = flapwise_c_lo
        self.flapwise_c_hi = flapwise_c_hi
        self.flapwise_clip_mode = flapwise_clip_mode

    def __repr__(self):
        return (
            f"SandwichCore(bay_mode={self.bay_mode!r}, "
            f"x_start={self.x_start}, x_end={self.x_end})"
        )
