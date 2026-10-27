"""Constructive spar-cap solids must split on the cavity midline (SG-012).

``cap_feasibility_clamp`` already samples ``camber_xy`` and the SDF SparCap path
honours it. CapWall recipes used by ``emit_wall_solid_polygons`` previously
omitted those samples, so constructive emission fell back to the chord line
``y=0`` and suction stole pressure-side material on cambered sections.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SectionGeometryRunner,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.constructive.export import (
    emit_wall_solid_polygons,
)


@pytest.mark.parametrize("twist_deg", [0.0, 2.7])
def test_constructive_caps_follow_cavity_midline_not_y0(twist_deg: float) -> None:
    """Equal-height caps on NACA 2412 must not meet on y=0 with a 2.7x area skew."""
    chord = 0.1228
    result = (
        SectionGeometryRunner(
            layout=resolve_system_type("2D-F"),
            airfoil=AirfoilSDF.from_naca("2412", chord=chord),
            skin_thickness=0.003,
            spar_cap_thickness=0.006125,
            shear_web_thickness=0.004,
            twist_rad=float(np.deg2rad(twist_deg)),
            outputs=("properties",),
        )
        .execute()
        .get_results()
    )
    recipes = result.section.spar_cap_wall_recipes
    assert recipes["spar_cap_suction_surface"].camber_xy is not None
    assert recipes["spar_cap_pressure_surface"].camber_xy is not None

    polys = emit_wall_solid_polygons(result.section)
    su = polys["spar_cap_suction_surface"]
    pr = polys["spar_cap_pressure_surface"]
    assert not su.is_empty and not pr.is_empty
    assert su.intersection(pr).area < 1e-10

    # Pre-fix: ratio ~ 2.7 with the shared boundary pinned at y=0.
    ratio = float(su.area / pr.area)
    assert 0.75 <= ratio <= 1.35, f"camber split still skewed: area ratio={ratio:.3f}"

    if twist_deg == 0.0:
        cam = np.asarray(recipes["spar_cap_suction_surface"].camber_xy, dtype=float)
        cam_y = float(np.mean(cam[:, 1]))
        assert cam_y > 1e-4
        # Untwisted solids: suction sits above the midline, pressure below it.
        assert su.bounds[1] > 0.5 * cam_y
        assert pr.bounds[3] > 0.5 * cam_y
        assert abs(0.5 * (su.bounds[1] + pr.bounds[3]) - cam_y) < 1e-3
