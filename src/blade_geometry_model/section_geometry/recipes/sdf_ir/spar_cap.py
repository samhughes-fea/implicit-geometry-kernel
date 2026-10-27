"""Spar cap subcomponent."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry.sdf.ir import CallableField, Intersect, Offset, Rotate, Slab

from blade_geometry_model.section_geometry.sections.clip_algebra import (
    ChordBand,
    _airfoil_offset_arc,
    _parallel_web_strip_clip,
    airfoil_expr,
    camber_side_field,
    lower_half_plane_expr,
    strip_clip_expr,
    to_ir,
    upper_half_plane_expr,
)
from blade_geometry_model.section_geometry.sections.labels import (
    ComponentLabel,
)
from .base import SectionComponent


class SparCap(SectionComponent):
    """Constant-thickness spar cap flush against the inner skin surface."""

    label = ComponentLabel.SPAR_CAP

    def __init__(
        self,
        airfoil_sdf,
        skin_thickness,
        x_start,
        x_end,
        cap_height,
        surface="upper",
        strip_clip=None,
        twist_angle=0.0,
        *,
        airfoil_expr_node=None,
        camber_xy=None,
        flapwise_nx=None,
        flapwise_ny=None,
        flapwise_c_lo=None,
        flapwise_c_hi=None,
        flapwise_clip_mode=None,
    ):
        if surface not in ("upper", "lower", "suction_surface", "pressure_surface"):
            raise ValueError(
                "surface must be 'suction_surface'/'pressure_surface' "
                "(aliases 'upper'/'lower')."
            )
        # Camber-split math still uses upper/lower; emit metadata uses long names.
        _camber = {
            "upper": "upper",
            "suction_surface": "upper",
            "lower": "lower",
            "pressure_surface": "lower",
        }[surface]
        self.surface = (
            "suction_surface"
            if _camber == "upper"
            else "pressure_surface"
        )
        self._camber_surface = _camber
        # Flapwise strip-clip parameters, carried explicitly so the Shapely
        # emitter can clip the *solid* directly (``flapwise_bay_clip``) instead
        # of approximating the clip from the midline's chordwise extent. Same
        # contract as ``SandwichCore``; ``twist_component`` copies these.
        self.flapwise_nx = flapwise_nx
        self.flapwise_ny = flapwise_ny
        self.flapwise_c_lo = flapwise_c_lo
        self.flapwise_c_hi = flapwise_c_hi
        self.flapwise_clip_mode = flapwise_clip_mode
        self.x_start = float(x_start)
        self.x_end = float(x_end)
        self.cap_height = float(cap_height)
        self.skin_thickness = float(skin_thickness)
        self.twist_angle = float(twist_angle)
        self._af = airfoil_sdf
        # CapWall carries flapwise_* without an explicit callable strip; derive
        # the same midplane-to-midplane clip cores already honor (SG-011).
        if strip_clip is None and self._flapwise_strip_ready():
            strip_clip = _parallel_web_strip_clip(
                float(flapwise_nx),
                float(flapwise_ny),
                float(flapwise_c_lo),
                float(flapwise_c_hi),
            )
        self._strip_clip = strip_clip
        self._camber_xy = None if camber_xy is None else np.asarray(camber_xy, dtype=float)

        poly = airfoil_expr_node if airfoil_expr_node is not None else airfoil_expr(airfoil_sdf)
        if poly is None:
            raise ValueError("SparCap requires an AirfoilSDF with vertices for IR evaluation.")

        inner_skin = Offset(poly, amount=-skin_thickness)
        cap_lam = Slab(inner_skin, depth=cap_height)
        if self._camber_xy is not None and self._camber_xy.shape[0] >= 2:
            side_expr = CallableField(
                key=f"cap_side_{self._camber_surface}_{id(self)}",
                fn=camber_side_field(self._camber_xy, self._camber_surface),
            )
        else:
            side_expr = (
                upper_half_plane_expr()
                if self._camber_surface == "upper"
                else lower_half_plane_expr()
            )
        if self._strip_clip is not None:
            clip_expr = strip_clip_expr(
                self._strip_clip, key=f"spar_strip_{id(self._strip_clip)}"
            )
        else:
            clip_expr = to_ir(ChordBand(x_start=float(x_start), x_end=float(x_end)))
        expr_raw = Intersect(Intersect(cap_lam, side_expr), clip_expr)
        if abs(twist_angle) > 1e-10:
            cx = 0.5 * (self.x_start + self.x_end)
            expr = Rotate(expr_raw, angle=float(twist_angle), cx=float(cx), cy=0.0)
        else:
            expr = expr_raw
        super().__init__(expr)

    def _flapwise_strip_ready(self) -> bool:
        mode = str(self.flapwise_clip_mode or "strip")
        return (
            self.flapwise_nx is not None
            and self.flapwise_ny is not None
            and self.flapwise_c_lo is not None
            and self.flapwise_c_hi is not None
            and mode == "strip"
        )

    def midline_polyline(self, n=80):
        arc = _airfoil_offset_arc(
            self._af,
            skin_thickness=self.skin_thickness,
            cap_height=self.cap_height,
            x_start=self.x_start,
            x_end=self.x_end,
            surface=self._camber_surface,
            n=n,
            strip_clip=self._strip_clip,
        )
        if abs(self.twist_angle) > 1e-10:
            cx = 0.5 * (self.x_start + self.x_end)
            c = float(np.cos(self.twist_angle))
            s = float(np.sin(self.twist_angle))
            rel = arc - np.array([cx, 0.0])
            arc = np.column_stack(
                [c * rel[:, 0] - s * rel[:, 1], s * rel[:, 0] + c * rel[:, 1]]
            ) + np.array([cx, 0.0])
        return arc

    def __repr__(self):
        return (
            f"SparCap(surface={self.surface!r}, "
            f"x=[{self.x_start:.3g},{self.x_end:.3g}], "
            f"cap_height={self.cap_height:.3g})"
        )
