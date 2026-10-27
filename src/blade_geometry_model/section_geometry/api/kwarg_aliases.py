"""Resolve preferred long thickness / twist kwargs (optional secondary aliases)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TypeAlias

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
    spar_cap_thickness: float | None = None,
    cap_height: float | None = None,
    default: float | None = None,
    required: bool = False,
) -> float:
    return _pick_exclusive(
        spar_cap_thickness,
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
    """Canonical twist [rad]. Aliases are also radians (not degrees)."""
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


__all__ = [
    "SkinGaugeArg",
    "WebThicknessArg",
    "resolve_skin_gauges",
    "resolve_skin_thickness",
    "resolve_spar_cap_thickness",
    "resolve_shear_web_thickness",
    "resolve_twist_rad",
    "resolve_geometry_mode",
]
