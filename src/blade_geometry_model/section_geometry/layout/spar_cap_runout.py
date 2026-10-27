"""BL-006 spar-cap run-out: terminate caps when they cannot be built at gauge.

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

from blade_geometry_model.section_geometry.layout.system_layout import (
    layout_without_spar_caps,
)
from blade_geometry_model.section_geometry.sections.cap_feasibility import (
    cap_feasibility_clamp,
)

logger = logging.getLogger(__name__)


def max_spar_from_cavity(
    h_local: float,
    skin_thickness: float,
    *,
    gap_min: float = 0.002,
) -> float:
    """Upper bound on each spar cap so upper/lower never meet [m]."""
    if not np.isfinite(h_local) or h_local <= 0.0:
        return float("inf")
    room = 0.5 * (float(h_local) - float(skin_thickness) - float(gap_min))
    if not np.isfinite(room):
        return float("inf")
    return float(max(0.0, room))


def _cap_rule(layout: Any) -> Any:
    rule = getattr(layout, "cap_placement", None)
    if rule is not None:
        return rule
    from blade_geometry_model.section_geometry.layout.cap_placement import (
        CapPlacementRule,
    )

    fam = getattr(layout, "structural_family", None) or "A"
    return CapPlacementRule.from_structural_family(str(fam))


def _placement_buildable_at_gauge(
    resolution: Any,
    *,
    airfoil: Any,
    skin_thickness: float,
    spar_min_gauge: float,
) -> bool:
    """True when min-gauge caps survive the same feasibility build_spar_caps uses."""
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
        dhw = 0.04 * chord
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
        half = float(getattr(rule, "half_width", None) or (0.05 * chord))
        x_c = 0.5 * chord
        _x_lo, _x_hi, h_u, h_l, info = cap_feasibility_clamp(
            airfoil,
            skin,
            x_c - half,
            x_c + half,
            gauge,
            gauge,
            y_search=None,
            phi_inner_ref=None,
            aft_half_depth_rule=False,
            extent_anchor="band_center",
            center_x=x_c,
        )
        if bool(info.get("pad_run_out")):
            return False
        return float(h_u) >= gauge - 1e-12 and float(h_l) >= gauge - 1e-12

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
) -> Any:
    """Return resolution, possibly stepped to ``layout_without_spar_caps``.

    Fires when:
    1. Global cavity room at the skin floor cannot hold min gauge, or
    2. ``airfoil`` is provided and placement feasibility at min gauge fails
       (pad run-out / continuous-box clamp / fixed-band collapse).
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
        ):
            reasons.append("placement cannot host min-gauge caps (pad/box feasibility)")

    if not reasons:
        return resolution

    capless = layout_without_spar_caps(layout)
    logger.warning(
        "[Station %d] spar cap runs out: %s. Terminating the cap and building %s "
        "instead of %s -- every web is kept. See docs/bug_log.md BL-006.",
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
