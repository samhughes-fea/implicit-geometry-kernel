"""Frozen geometry input for one station: layout, airfoil, twist, gauges, webs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, TypeAlias

SkinGaugeArg: TypeAlias = float | Sequence[float] | None
WebThicknessArg: TypeAlias = float | Sequence[float] | None


def _pick_exclusive(
    preferred: float | None,
    aliases: Sequence[float | None],
    *,
    names: str,
    default: float | None = None,
    required: bool = False,
) -> float:
    """Return the single provided value among preferred + aliases.

    Conflict (two unequal non-None values) → ``ValueError``.
    """
    present: list[float] = []
    if preferred is not None:
        present.append(float(preferred))
    for a in aliases:
        if a is not None:
            present.append(float(a))
    if not present:
        if required:
            raise TypeError(f"Missing required thickness kwargs ({names})")
        if default is None:
            raise TypeError(f"Missing thickness kwargs ({names})")
        return float(default)
    first = present[0]
    for v in present[1:]:
        if v != first:
            raise ValueError(f"Conflicting thickness kwargs for {names}")
    return first


def resolve_skin_gauges(
    *,
    skin_thickness: SkinGaugeArg = None,
    default: float | None = None,
    required: bool = False,
) -> tuple[float, float]:
    """Return ``(suction_surface, pressure_surface)`` gauges [m].

    A scalar broadcasts to both sides. A length-2 sequence is ``(SS, PS)``.
    """
    if skin_thickness is None:
        if required and default is None:
            raise TypeError("Missing required thickness kwargs (skin_thickness)")
        if default is None:
            raise TypeError("Missing thickness kwargs (skin_thickness)")
        t = float(default)
        return t, t
    if isinstance(skin_thickness, (int, float)):
        t = float(skin_thickness)
        return t, t
    vals = [float(x) for x in skin_thickness]
    if len(vals) != 2:
        raise ValueError(
            "skin_thickness sequence must be length 2 (suction, pressure); "
            f"got {len(vals)}."
        )
    return vals[0], vals[1]


def resolve_skin_thickness(
    *,
    skin_thickness: SkinGaugeArg = None,
    default: float | None = None,
    required: bool = False,
) -> float:
    """Scalar skin thickness for cavity / anchor probes (``max(SS, PS)``)."""
    ss, ps = resolve_skin_gauges(
        skin_thickness=skin_thickness, default=default, required=required
    )
    return max(ss, ps)


def resolve_spar_cap_thickness(
    *,
    spar_cap_thickness: float | tuple[float, float] | None = None,
    cap_height: float | None = None,
    default: float | None = None,
    required: bool = False,
) -> float | tuple[float, float]:
    """Scalar, or ``(suction, pressure)`` when the two caps differ."""
    raw = spar_cap_thickness if spar_cap_thickness is not None else cap_height
    if isinstance(raw, (tuple, list)) and not isinstance(raw, str):
        vals = [float(v) for v in raw]
        if len(vals) != 2:
            raise ValueError(
                "spar_cap_thickness sequence must be (suction, pressure); "
                f"got {len(vals)}"
            )
        return (vals[0], vals[1])
    scalar = spar_cap_thickness if isinstance(spar_cap_thickness, (int, float)) else None
    return _pick_exclusive(
        scalar,
        (cap_height,),
        names="spar_cap_thickness / cap_height",
        default=default,
        required=required,
    )


def resolve_shear_web_thickness(
    *,
    shear_web_thickness: WebThicknessArg = None,
    web_thickness: WebThicknessArg = None,
    default: float | None = None,
    required: bool = False,
) -> float | list[float]:
    """Scalar or per-web list of shear-web thicknesses [m]."""
    preferred = shear_web_thickness if shear_web_thickness is not None else web_thickness
    if preferred is None:
        if required and default is None:
            raise TypeError(
                "Missing required thickness kwargs "
                "(shear_web_thickness / web_thickness)"
            )
        if default is None:
            raise TypeError(
                "Missing thickness kwargs (shear_web_thickness / web_thickness)"
            )
        return float(default)
    if shear_web_thickness is not None and web_thickness is not None:
        a = (
            float(shear_web_thickness)
            if isinstance(shear_web_thickness, (int, float))
            else [float(x) for x in shear_web_thickness]
        )
        b = (
            float(web_thickness)
            if isinstance(web_thickness, (int, float))
            else [float(x) for x in web_thickness]
        )
        if a != b:
            raise ValueError(
                "Conflicting thickness kwargs for shear_web_thickness / web_thickness"
            )
    if isinstance(preferred, (int, float)):
        return float(preferred)
    return [float(x) for x in preferred]


def resolve_twist_rad(
    *,
    twist_rad: float | None = None,
    twist_angle_rad: float | None = None,
    twist_angle: float | None = None,
    default: float = 0.0,
) -> float:
    """Canonical twist [rad]. ``twist_angle`` / ``twist_angle_rad`` are aliases (also radians)."""
    return _pick_exclusive(
        twist_rad,
        (twist_angle_rad, twist_angle),
        names="twist_rad / twist_angle_rad / twist_angle",
        default=default,
    )


def resolve_geometry_mode(token: str | None) -> str:
    """Map legacy ``airfoil_sdf_only`` → ``zero_web_catalogue``."""
    if token is None:
        return "multicell"
    t = str(token)
    if t == "airfoil_sdf_only":
        return "zero_web_catalogue"
    return t


@dataclass(frozen=True)
class Gauges:
    """Wall thicknesses in metres: skin suction/pressure, spar cap, shear web."""

    skin_ss: float
    skin_ps: float
    cap: float | tuple[float, float]
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
    web_indices: tuple[int, ...] | None = None
    web_y_anchors: tuple[tuple[float, float], ...] | None = None
    spar_cap_ply_pitch_m: float | None = None
    zone_gauges_m: tuple[tuple[str, float], ...] | None = None
    pivot_fraction: float | None = None

    def __post_init__(self) -> None:
        if self.layout is None:
            raise ValueError("layout is required")
        if self.airfoil is None:
            raise ValueError("airfoil is required")
        object.__setattr__(self, "twist_rad", float(self.twist_rad))
        if self.webs is not None:
            object.__setattr__(self, "webs", tuple(float(x) for x in self.webs))
        if self.web_indices is not None:
            object.__setattr__(self, "web_indices", tuple(int(i) for i in self.web_indices))
        if self.web_y_anchors is not None:
            object.__setattr__(
                self,
                "web_y_anchors",
                tuple((float(a), float(b)) for a, b in self.web_y_anchors),
            )
        if self.spar_cap_ply_pitch_m is not None:
            object.__setattr__(self, "spar_cap_ply_pitch_m", float(self.spar_cap_ply_pitch_m))
        if self.pivot_fraction is not None:
            f = float(self.pivot_fraction)
            if not 0.0 < f < 1.0:
                raise ValueError("pivot_fraction must be a chord fraction in (0, 1)")
            object.__setattr__(self, "pivot_fraction", f)
        g = self.gauges
        caps = g.cap if isinstance(g.cap, tuple) else (g.cap,)
        if min(g.skin_ss, g.skin_ps, *caps) < 0.0:
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
        web_offsets: Sequence[float] | None = None,
        web_indices: Sequence[int] | None = None,
        web_y_coords: Sequence[tuple[float, float]] | None = None,
        spar_cap_ply_pitch_m: float | None = None,
        zone_gauges_m: Mapping[str, float] | None = None,
        pivot_fraction: float | None = None,
        **_ignored: Any,
    ) -> SectionSpec:
        """Build a spec from the historical keyword form.

        Alias conflicts (``cap_height`` vs ``spar_cap_thickness``, twist
        names, web offsets) are resolved here; canonical fields are then
        checked in ``__post_init__``.
        """
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
        cap_field: float | tuple[float, float]
        if isinstance(cap, tuple):
            cap_field = cap
        else:
            cap_field = float(cap)
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
        if "web_x_positions_m" in _ignored:
            raise TypeError("web_x_positions_m was removed in 0.6.0; use web_offsets")
        xs = web_offsets
        web_t: float | tuple[float, ...]
        if isinstance(web, list):
            web_t = tuple(float(v) for v in web)
        else:
            web_t = float(web)
        return cls(
            layout=layout,
            airfoil=airfoil,
            twist_rad=twist,
            gauges=Gauges(skin_ss=ss, skin_ps=ps, cap=cap_field, web=web_t),
            webs=None if xs is None else tuple(float(x) for x in xs),
            web_indices=None if web_indices is None else tuple(int(i) for i in web_indices),
            web_y_anchors=None if web_y_coords is None else tuple((float(a), float(b)) for a, b in web_y_coords),
            spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
            zone_gauges_m=zone_items,
            pivot_fraction=pivot_fraction,
        )
