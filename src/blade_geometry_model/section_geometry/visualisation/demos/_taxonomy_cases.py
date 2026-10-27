"""Shared twist / airfoil helpers for taxonomy-driven section_geometry examples."""

from __future__ import annotations

from typing import Any

import numpy as np

from blade_geometry_model.section_geometry.layout.system_layout import SystemLayoutSpec


def twist_deg_for_layout(
    layout: SystemLayoutSpec,
    *,
    twist_deg: float = 15.0,
) -> float:
    """Return section twist [deg] for a registry layout (uniform across taxonomy)."""
    del layout  # reserved for future per-layout overrides
    return float(twist_deg)


def twist_rad_for_layout(
    layout: SystemLayoutSpec,
    *,
    twist_deg: float = 15.0,
) -> float:
    return float(np.deg2rad(twist_deg_for_layout(layout, twist_deg=twist_deg)))


def airfoil_for_twist(airfoil: Any, twist_rad: float) -> Any:
    """Chord-frame airfoil rotated into B-frame when twist is non-zero."""
    if abs(float(twist_rad)) > 1e-12:
        return airfoil.rotate(twist_rad)
    return airfoil


def system_type_display_name(key: str) -> str:
    return f"SystemType{key}"
