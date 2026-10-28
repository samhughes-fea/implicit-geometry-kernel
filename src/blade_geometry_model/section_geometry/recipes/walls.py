"""Frozen wall recipes that delegate polygon emission to constructive emitters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from blade_geometry_model.section_geometry.constructive.emitters import (
    emit_leading_edge_insert_polygon,
    emit_outer_skin_polygon,
    emit_sandwich_core_polygon,
    emit_shear_web_polygon,
    emit_spar_cap_polygon,
    emit_trailing_edge_insert_polygon,
)
from blade_geometry_model.section_geometry.recipes.sdf_ir.outer_skin import (
    OuterSkin,
)
from blade_geometry_model.section_geometry.recipes.sdf_ir.sandwich_core import (
    SandwichCore,
)
from blade_geometry_model.section_geometry.recipes.sdf_ir.shear_web import (
    ShearWeb,
)
from blade_geometry_model.section_geometry.recipes.sdf_ir.spar_cap import (
    SparCap,
)

EmissionMaturity = Literal["constructive", "residual_phi"]
BayMode = Literal["chordwise", "flapwise"]
FlapwiseClipMode = Literal["half_lt", "half_gt", "strip"]


@dataclass(frozen=True)
class SkinWall:
    """Outer-skin annulus recipe (constructive polygon emission)."""

    thickness: float
    skin_side: str | None = None
    emission_maturity: EmissionMaturity = "constructive"

    def _component(self, airfoil, *, camber_xy=None, airfoil_expr_node=None) -> OuterSkin:
        return OuterSkin(
            airfoil,
            float(self.thickness),
            skin_side=self.skin_side,
            camber_xy=camber_xy,
            airfoil_expr_node=airfoil_expr_node,
        )

    def to_sdf_component(self, airfoil, *, camber_xy=None, airfoil_expr_node=None):
        """Return the existing :class:`OuterSkin` instance for this recipe."""
        return self._component(
            airfoil, camber_xy=camber_xy, airfoil_expr_node=airfoil_expr_node
        )

    def to_boundary_polygon(self, airfoil, *, camber_xy=None):
        if self.emission_maturity != "constructive":
            raise NotImplementedError(
                "SkinWall boundary emission currently requires "
                "emission_maturity='constructive'."
            )
        return emit_outer_skin_polygon(self._component(airfoil, camber_xy=camber_xy))

    def to_midline(self, airfoil, *, camber_xy=None):
        return self._component(airfoil, camber_xy=camber_xy).midline_polyline()

    def to_sdf_expr(self, airfoil, *, camber_xy=None, airfoil_expr_node=None):
        return self.to_sdf_component(
            airfoil, camber_xy=camber_xy, airfoil_expr_node=airfoil_expr_node
        ).expr

    def to_sdf(self, airfoil, **kwargs):
        return self.to_sdf_expr(airfoil, **kwargs)



@dataclass(frozen=True)
class WebWall:
    """Shear-web stadium recipe (constructive polygon emission)."""

    skin_thickness: float
    x_top: float
    y_top: float
    x_bot: float
    y_bot: float
    thickness: float
    alignment: str = "chord_normal"
    twist_angle: float = 0.0
    emission_maturity: EmissionMaturity = "constructive"

    def _component(self, airfoil, *, airfoil_expr_node=None) -> ShearWeb:
        return ShearWeb(
            airfoil,
            float(self.skin_thickness),
            float(self.x_top),
            float(self.y_top),
            float(self.x_bot),
            float(self.y_bot),
            float(self.thickness),
            alignment=self.alignment,
            twist_angle=float(self.twist_angle),
            airfoil_expr_node=airfoil_expr_node,
        )

    def to_sdf_component(self, airfoil, *, airfoil_expr_node=None):
        return self._component(airfoil, airfoil_expr_node=airfoil_expr_node)

    def to_boundary_polygon(self, airfoil, *, cavity=None, resolution: int = 16):
        if self.emission_maturity != "constructive":
            raise NotImplementedError(
                "WebWall boundary emission currently requires "
                "emission_maturity='constructive'."
            )
        return emit_shear_web_polygon(
            self._component(airfoil), resolution=resolution, cavity=cavity
        )

    def to_midline(self, airfoil, *, n: int = 20, airfoil_expr_node=None):
        return self._component(
            airfoil, airfoil_expr_node=airfoil_expr_node
        ).midline_polyline(n=n)

    def to_sdf_expr(self, airfoil, *, airfoil_expr_node=None):
        return self.to_sdf_component(
            airfoil, airfoil_expr_node=airfoil_expr_node
        ).expr

    def to_sdf(self, airfoil, **kwargs):
        return self.to_sdf_expr(airfoil, **kwargs)



@dataclass(frozen=True)
class CapWall:
    """Spar-cap slab recipe (constructive polygon emission)."""

    skin_thickness: float
    x_start: float
    x_end: float
    cap_height: float
    surface: str = "upper"
    twist_angle: float = 0.0
    flapwise_nx: float | None = None
    flapwise_ny: float | None = None
    flapwise_c_lo: float | None = None
    flapwise_c_hi: float | None = None
    flapwise_clip_mode: FlapwiseClipMode | None = None
    # Cavity-midline samples from cap_feasibility_clamp (SG-010 / SG-012).
    # Stored on the recipe so constructive emit does not fall back to y=0.
    camber_xy: tuple[tuple[float, float], ...] | None = None
    emission_maturity: EmissionMaturity = "constructive"

    def _resolved_camber_xy(self, camber_xy=None):
        if camber_xy is not None:
            return camber_xy
        return self.camber_xy

    def _component(self, airfoil, *, camber_xy=None, airfoil_expr_node=None) -> SparCap:
        return SparCap(
            airfoil,
            float(self.skin_thickness),
            float(self.x_start),
            float(self.x_end),
            float(self.cap_height),
            surface=self.surface,
            twist_angle=float(self.twist_angle),
            camber_xy=self._resolved_camber_xy(camber_xy),
            airfoil_expr_node=airfoil_expr_node,
            flapwise_nx=self.flapwise_nx,
            flapwise_ny=self.flapwise_ny,
            flapwise_c_lo=self.flapwise_c_lo,
            flapwise_c_hi=self.flapwise_c_hi,
            flapwise_clip_mode=self.flapwise_clip_mode,
        )

    def to_sdf_component(self, airfoil, *, camber_xy=None, airfoil_expr_node=None):
        return self._component(
            airfoil, camber_xy=camber_xy, airfoil_expr_node=airfoil_expr_node
        )

    def to_boundary_polygon(self, airfoil, *, camber_xy=None, **emit_kwargs):
        if self.emission_maturity != "constructive":
            raise NotImplementedError(
                "CapWall boundary emission currently requires "
                "emission_maturity='constructive'."
            )
        return emit_spar_cap_polygon(
            self._component(airfoil, camber_xy=camber_xy), **emit_kwargs
        )

    def to_midline(self, airfoil, *, n: int = 80, camber_xy=None, airfoil_expr_node=None):
        return self._component(
            airfoil, camber_xy=camber_xy, airfoil_expr_node=airfoil_expr_node
        ).midline_polyline(n=n)

    def to_sdf_expr(self, airfoil, *, camber_xy=None, airfoil_expr_node=None):
        return self.to_sdf_component(
            airfoil, camber_xy=camber_xy, airfoil_expr_node=airfoil_expr_node
        ).expr

    def to_sdf(self, airfoil, **kwargs):
        return self.to_sdf_expr(airfoil, **kwargs)



@dataclass(frozen=True)
class InsertWall:
    """LE/TE insert recipe — constructive by default once parity is locked."""

    kind: Literal["leading_edge", "trailing_edge"]
    skin_thickness: float
    x_limit: float
    emission_maturity: EmissionMaturity = "constructive"
    le_x: float = 0.0
    le_y: float = 0.0
    radius: float | None = None

    def to_sdf_component(self, airfoil, *, airfoil_expr_node=None):
        """Build the existing LE/TE insert component (SDF / IR path)."""
        from blade_geometry_model.section_geometry.recipes.sdf_ir.inserts import (
            LEInsert,
            TEInsert,
        )

        if self.kind == "trailing_edge":
            return TEInsert(
                airfoil,
                float(self.skin_thickness),
                x_start=float(self.x_limit),
                airfoil_expr_node=airfoil_expr_node,
            )
        if self.kind == "leading_edge":
            return LEInsert(
                airfoil,
                float(self.skin_thickness),
                x_end=float(self.x_limit),
                le_x=float(self.le_x),
                le_y=float(self.le_y),
                radius=self.radius,
                airfoil_expr_node=airfoil_expr_node,
            )
        raise ValueError(f"Unknown insert kind {self.kind!r}")

    def to_boundary_polygon(self, airfoil, *, airfoil_expr_node=None, **emit_kwargs):
        if self.emission_maturity != "constructive":
            raise NotImplementedError(
                "InsertWall boundary polygon emission requires "
                f"emission_maturity='constructive' "
                f"(kind={self.kind!r}, got {self.emission_maturity!r})."
            )
        comp = self.to_sdf_component(airfoil, airfoil_expr_node=airfoil_expr_node)
        if self.kind == "trailing_edge":
            return emit_trailing_edge_insert_polygon(comp, **emit_kwargs)
        return emit_leading_edge_insert_polygon(comp, **emit_kwargs)

    def to_midline(self, *_args, **_kwargs):
        raise NotImplementedError(
            "InsertWall has no shell midline; inserts are solids, not wall strips."
        )

    def to_sdf_expr(self, airfoil, *, airfoil_expr_node=None):
        return self.to_sdf_component(
            airfoil, airfoil_expr_node=airfoil_expr_node
        ).expr

    def to_sdf(self, airfoil, **kwargs):
        return self.to_sdf_expr(airfoil, **kwargs)



@dataclass(frozen=True)
class CoreWall:
    """Sandwich-core bay recipe (constructive polygon emission)."""

    skin_thickness: float
    bay_mode: BayMode = "chordwise"
    x_start: float | None = None
    x_end: float | None = None
    y_min: float | None = None
    y_max: float | None = None
    flapwise_nx: float | None = None
    flapwise_ny: float | None = None
    flapwise_c_lo: float | None = None
    flapwise_c_hi: float | None = None
    flapwise_clip_mode: FlapwiseClipMode | None = None
    emission_maturity: EmissionMaturity = "constructive"

    def to_sdf_component(
        self,
        airfoil,
        *,
        airfoil_expr_node=None,
        exclusion_exprs=None,
        exclusion_sdfs=None,
        strip_clip=None,
    ) -> SandwichCore:
        return SandwichCore(
            airfoil,
            float(self.skin_thickness),
            exclusion_sdfs=exclusion_sdfs,
            x_start=self.x_start,
            x_end=self.x_end,
            y_min=self.y_min,
            y_max=self.y_max,
            strip_clip=strip_clip,
            airfoil_expr_node=airfoil_expr_node,
            exclusion_exprs=exclusion_exprs,
            bay_mode=self.bay_mode,
            flapwise_nx=self.flapwise_nx,
            flapwise_ny=self.flapwise_ny,
            flapwise_c_lo=self.flapwise_c_lo,
            flapwise_c_hi=self.flapwise_c_hi,
            flapwise_clip_mode=self.flapwise_clip_mode,
        )

    def to_boundary_polygon(
        self,
        airfoil,
        *,
        wall_exclusion_union=None,
        cavity=None,
        airfoil_expr_node=None,
        exclusion_exprs=None,
    ):
        if self.emission_maturity != "constructive":
            raise NotImplementedError(
                "CoreWall boundary emission currently requires "
                "emission_maturity='constructive'."
            )
        comp = self.to_sdf_component(
            airfoil,
            airfoil_expr_node=airfoil_expr_node,
            exclusion_exprs=exclusion_exprs,
        )
        return emit_sandwich_core_polygon(
            comp, wall_exclusion_union, cavity=cavity
        )

    def to_midline(self, *_args, **_kwargs):
        raise NotImplementedError(
            "CoreWall has no shell midline; cores fill bays, they are not wall strips."
        )

    def to_sdf_expr(
        self,
        airfoil,
        *,
        airfoil_expr_node=None,
        exclusion_exprs=None,
        exclusion_sdfs=None,
        strip_clip=None,
    ):
        return self.to_sdf_component(
            airfoil,
            airfoil_expr_node=airfoil_expr_node,
            exclusion_exprs=exclusion_exprs,
            exclusion_sdfs=exclusion_sdfs,
            strip_clip=strip_clip,
        ).expr

    def to_sdf(self, airfoil, **kwargs):
        return self.to_sdf_expr(airfoil, **kwargs)



__all__ = [
    "BayMode",
    "CapWall",
    "CoreWall",
    "EmissionMaturity",
    "FlapwiseClipMode",
    "InsertWall",
    "SkinWall",
    "WebWall",
]
