"""Frozen geometry input for one station: layout, airfoil, twist, gauges, webs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from implicit_geometry_kernel.section.spec.kwarg_aliases import (
    resolve_shear_web_thickness,
    resolve_skin_gauges,
    resolve_spar_cap_thickness,
    resolve_twist_rad,
)


@dataclass(frozen=True)
class Gauges:
    """Wall thicknesses in metres: skin suction/pressure, spar cap, shear web."""

    skin_ss: float
    skin_ps: float
    cap: float
    web: float | tuple[float, ...]

    @property
    def skin(self) -> float | tuple[float, float]:
        if self.skin_ss == self.skin_ps:
            return self.skin_ss
        return (self.skin_ss, self.skin_ps)

    @property
    def skin_max(self) -> float:
        return max(self.skin_ss, self.skin_ps)

    @property
    def web_arg(self) -> float | list[float]:
        if isinstance(self.web, tuple):
            return list(self.web)
        return self.web


@dataclass(frozen=True)
class SectionSpec:
    """One name per geometry quantity. Validates in ``__post_init__``."""

    layout: Any
    airfoil: Any
    twist_rad: float
    gauges: Gauges
    webs: tuple[float, ...] | None = None
    web_y_anchors: tuple[tuple[float, float], ...] | None = None
    spar_cap_ply_pitch_m: float | None = None
    zone_gauges_m: tuple[tuple[str, float], ...] | None = None

    def __post_init__(self) -> None:
        if self.layout is None:
            raise ValueError("layout is required")
        if self.airfoil is None:
            raise ValueError("airfoil is required")
        object.__setattr__(self, "twist_rad", float(self.twist_rad))
        if self.webs is not None:
            object.__setattr__(self, "webs", tuple(float(x) for x in self.webs))
        if self.web_y_anchors is not None:
            object.__setattr__(
                self,
                "web_y_anchors",
                tuple((float(a), float(b)) for a, b in self.web_y_anchors),
            )
        if self.spar_cap_ply_pitch_m is not None:
            object.__setattr__(self, "spar_cap_ply_pitch_m", float(self.spar_cap_ply_pitch_m))
        g = self.gauges
        if min(g.skin_ss, g.skin_ps, g.cap) < 0.0:
            raise ValueError("gauges must be non-negative")
        web = g.web
        if isinstance(web, tuple):
            if any(float(v) < 0.0 for v in web):
                raise ValueError("gauges must be non-negative")
        elif float(web) < 0.0:
            raise ValueError("gauges must be non-negative")

    @classmethod
    def from_kwargs(
        cls,
        layout: Any,
        airfoil: Any,
        *,
        skin_thickness: Any = None,
        spar_cap_thickness: float | None = None,
        shear_web_thickness: Any = None,
        web_thickness: Any = None,
        cap_height: float | None = None,
        twist_rad: float | None = None,
        twist_angle_rad: float | None = None,
        twist_angle: float | None = None,
        web_x_positions_m: Sequence[float] | None = None,
        web_offsets: Sequence[float] | None = None,
        web_y_coords: Sequence[tuple[float, float]] | None = None,
        spar_cap_ply_pitch_m: float | None = None,
        zone_gauges_m: Mapping[str, float] | None = None,
        **_ignored: Any,
    ) -> SectionSpec:
        """Build a spec from the historical keyword form. Validation is in ``__post_init__``."""
        zone_items = None
        if zone_gauges_m:
            from implicit_geometry_kernel.section.spec.prep import (
                _compat_kwargs_from_zone_gauges,
            )

            if any(
                v is not None
                for v in (skin_thickness, spar_cap_thickness, shear_web_thickness)
            ):
                raise ValueError(
                    "geometry received both zone_gauges_m and component thickness kwargs"
                )
            derived = _compat_kwargs_from_zone_gauges(zone_gauges_m)
            skin_thickness = derived["skin_thickness"]
            spar_cap_thickness = derived["spar_cap_thickness"]
            shear_web_thickness = derived["shear_web_thickness"]
            zone_items = tuple(sorted((str(k), float(v)) for k, v in zone_gauges_m.items()))
        ss, ps = resolve_skin_gauges(skin_thickness=skin_thickness, required=True)
        cap = resolve_spar_cap_thickness(
            spar_cap_thickness=spar_cap_thickness, cap_height=cap_height, required=True
        )
        web = resolve_shear_web_thickness(
            shear_web_thickness=shear_web_thickness,
            web_thickness=web_thickness,
            required=True,
        )
        twist = resolve_twist_rad(
            twist_rad=twist_rad,
            twist_angle_rad=twist_angle_rad,
            twist_angle=twist_angle,
        )
        if web_offsets is not None and web_x_positions_m is not None:
            raise ValueError("pass web_offsets or web_x_positions_m, not both")
        xs = web_offsets if web_offsets is not None else web_x_positions_m
        web_t: float | tuple[float, ...]
        if isinstance(web, list):
            web_t = tuple(float(v) for v in web)
        else:
            web_t = float(web)
        return cls(
            layout=layout,
            airfoil=airfoil,
            twist_rad=twist,
            gauges=Gauges(skin_ss=ss, skin_ps=ps, cap=float(cap), web=web_t),
            webs=None if xs is None else tuple(float(x) for x in xs),
            web_y_anchors=None if web_y_coords is None else tuple((float(a), float(b)) for a, b in web_y_coords),
            spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
            zone_gauges_m=zone_items,
        )
