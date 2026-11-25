"""Spar-cap run-out: terminate caps when they cannot be built at gauge.

Web cascade and cap run-out are separate events. This module is the single
pre-assemble gate for the latter: global cavity room *and* placement
buildability (same feasibility ``build_spar_caps`` uses). Callers must apply
it before assembling a capped layout.
"""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import Any

import numpy as np

from implicit_geometry_kernel.catalogue.system_layout import (
    layout_without_spar_caps,
)
from implicit_geometry_kernel.section.assembly.cap_feasibility import (
    MIN_CAP_VOID_ABS_M,
    MIN_CAP_VOID_FRAC,
    cap_feasibility_clamp,
)

logger = logging.getLogger(__name__)


def max_spar_from_cavity(
    h_local: float,
    skin_thickness: float,
    *,
    gap_min: float = MIN_CAP_VOID_ABS_M,
) -> float:
    """Upper bound on each spar cap so upper/lower never meet [m].

    ``h_local`` is the depth between the two skin mid-surfaces, so the cavity depth is
    ``d = h_local - skin_thickness``. The caps keep the void of the cap clamp,
    ``max(gap_min, MIN_CAP_VOID_FRAC * d)`` (R10.7), and share the rest equally
    (R10.14). With the default ``gap_min`` the bound equals half the clamp's allowance,
    so the gate and the build agree (SG-039).
    """
    if not np.isfinite(h_local) or h_local <= 0.0:
        return float("inf")
    depth = float(h_local) - float(skin_thickness)
    void = max(float(gap_min), MIN_CAP_VOID_FRAC * depth)
    room = 0.5 * (depth - void)
    if not np.isfinite(room):
        return float("inf")
    return float(max(0.0, room))


def _cap_rule(layout: Any) -> Any:
    rule = getattr(layout, "cap_placement", None)
    if rule is not None:
        return rule
    from implicit_geometry_kernel.section.layout.cap_placement import (
        CapPlacementRule,
    )

    fam = getattr(layout, "structural_family", None) or "A"
    return CapPlacementRule.from_structural_family(str(fam))


def _fixed_bands(rule: Any, airfoil: Any, chord: float, pitch_fraction: float | None):
    """The chordwise band(s) ``_caps_fixed_band`` clamps, and its extent anchor."""
    from implicit_geometry_kernel.section.layout.section_axes import (
        max_thickness_chord_x,
        pitch_axis_x_from_le,
    )

    side_windows = getattr(rule, "side_windows", None)
    if side_windows is not None:
        return [(min(w) * chord, max(w) * chord) for w in side_windows], "deepest"
    window = getattr(rule, "chord_window", None)
    if window is not None:
        return [(min(window) * chord, max(window) * chord)], "deepest"
    half = float(getattr(rule, "half_width", None) or (0.05 * chord))
    if (getattr(rule, "anchor", None) or "max_thickness") == "pitching":
        if pitch_fraction is None:
            import warnings

            warnings.warn(
                "apply_spar_cap_runout: a fixed band anchored on the pitch axis needs "
                "pitch_fraction to test the band the build clamps; testing it at 0.5c as "
                "before. pitch_fraction will be required for this anchor in a later "
                "release (SG-041).",
                DeprecationWarning,
                stacklevel=4,
            )
            x_c = 0.5 * chord
        else:
            x_c = pitch_axis_x_from_le(airfoil, float(pitch_fraction))
    else:
        x_c = max_thickness_chord_x(airfoil)
    return [(x_c - half, x_c + half)], "band_center"


def _placement_buildable_at_gauge(
    resolution: Any,
    *,
    airfoil: Any,
    skin_thickness: float,
    spar_min_gauge: float,
    pitch_fraction: float | None = None,
    pad_half_width: float | None = None,
) -> bool:
    """True when min-gauge caps survive the same feasibility build_spar_caps uses.

    The bands are those the build clamps (SG-041): per-web pads of
    ``pad_half_width`` (default ``0.04 * chord``, the builder's default), and a fixed
    band at its published window(s) or centred on the rule's anchor.
    """
    layout = resolution.layout
    if not bool(getattr(layout, "has_spar_caps", False)):
        return True
    gauge = float(spar_min_gauge)
    if not (gauge > 0.0):
        return True

    rule = _cap_rule(layout)
    kind = str(getattr(rule, "kind", "none") or "none")
    if kind == "none":
        return True

    chord = float(getattr(airfoil, "chord", 1.0) or 1.0)
    xs = tuple(float(x) for x in (resolution.web_x_positions_m or ()))
    skin = float(skin_thickness)

    if kind == "per_web_pads":
        if not xs:
            return False
        dhw = float(pad_half_width) if pad_half_width is not None else 0.04 * chord
        live = 0
        for x in xs:
            _x_lo, _x_hi, h_u, h_l, info = cap_feasibility_clamp(
                airfoil,
                skin,
                float(x) - dhw,
                float(x) + dhw,
                gauge,
                gauge,
                y_search=None,
                phi_inner_ref=None,
                aft_half_depth_rule=False,
                extent_anchor="band_center",
                center_x=float(x),
            )
            if (
                float(h_u) >= gauge - 1e-12
                and float(h_l) >= gauge - 1e-12
                and (float(_x_hi) - float(_x_lo)) > 1e-12
                and not bool(info.get("pad_run_out"))
            ):
                live += 1
        return live > 0

    if kind == "continuous_box":
        if len(xs) < 2:
            return False
        _x_lo, _x_hi, h_u, h_l, _info = cap_feasibility_clamp(
            airfoil,
            skin,
            float(xs[0]),
            float(xs[-1]),
            gauge,
            gauge,
            y_search=None,
            phi_inner_ref=None,
            aft_half_depth_rule=True,
        )
        return float(h_u) >= gauge - 1e-12 and float(h_l) >= gauge - 1e-12

    if kind == "fixed_band":
        bands, extent_anchor = _fixed_bands(rule, airfoil, chord, pitch_fraction)
        for x_lo, x_hi in bands:
            x_c = 0.5 * (x_lo + x_hi)
            _x_lo, _x_hi, h_u, h_l, info = cap_feasibility_clamp(
                airfoil,
                skin,
                x_lo,
                x_hi,
                gauge,
                gauge,
                y_search=None,
                phi_inner_ref=None,
                aft_half_depth_rule=False,
                extent_anchor=extent_anchor,
                center_x=x_c,
            )
            if bool(info.get("pad_run_out")):
                return False
            if not (float(h_u) >= gauge - 1e-12 and float(h_l) >= gauge - 1e-12):
                return False
        return True

    return True


def apply_spar_cap_runout(
    resolution: Any,
    *,
    h_local: float | None,
    gap_min: float,
    skin_thickness_floor: float,
    spar_min_gauge: float,
    station_idx: int = -1,
    airfoil: Any | None = None,
    pitch_fraction: float | None = None,
    pad_half_width: float | None = None,
) -> Any:
    """Return resolution, possibly stepped to ``layout_without_spar_caps``.

    Fires when:
    1. Global cavity room at the skin floor cannot hold min gauge, or
    2. ``airfoil`` is provided and placement feasibility at min gauge fails
       (pad run-out / continuous-box clamp / fixed-band collapse).

    ``h_local`` is the depth between the skin mid-surfaces (see
    :func:`max_spar_from_cavity`). ``pitch_fraction`` places a fixed band anchored on the pitch
    axis (without it that band is tested at ``0.5 * chord``, with a
    ``DeprecationWarning``); ``pad_half_width`` is the per-web pad half-width the
    section is built with (default ``0.04 * chord``).
    """
    if resolution is None:
        return resolution
    layout = resolution.layout
    if not bool(getattr(layout, "has_spar_caps", False)):
        return resolution
    if not (float(spar_min_gauge) > 0.0):
        return resolution

    reasons: list[str] = []
    if h_local is not None and float(h_local) > 0.0:
        cap_ceiling = max_spar_from_cavity(
            float(h_local), float(skin_thickness_floor), gap_min=float(gap_min)
        )
        if np.isfinite(cap_ceiling) and float(cap_ceiling) < float(spar_min_gauge) - 1e-12:
            reasons.append(
                f"cavity allows {float(cap_ceiling) * 1e3:.3f} mm "
                f"< {float(spar_min_gauge) * 1e3:.3f} mm gauge"
            )

    if airfoil is not None and not reasons:
        if not _placement_buildable_at_gauge(
            resolution,
            airfoil=airfoil,
            skin_thickness=float(skin_thickness_floor),
            spar_min_gauge=float(spar_min_gauge),
            pitch_fraction=pitch_fraction,
            pad_half_width=pad_half_width,
        ):
            reasons.append("placement cannot host min-gauge caps (pad/box feasibility)")

    if not reasons:
        return resolution

    capless = layout_without_spar_caps(layout)
    logger.warning(
        "[Station %d] spar cap runs out: %s. Terminating the cap and building %s "
        "instead of %s -- every web is kept. See docs/section/bug_log.md.",
        int(station_idx),
        "; ".join(reasons),
        getattr(capless, "key", capless),
        getattr(layout, "key", layout),
    )
    return replace(resolution, layout=capless)


__all__ = [
    "apply_spar_cap_runout",
    "max_spar_from_cavity",
]
