"""Family-D box caps must stay flush against both outer webs under twist.

At twist=0, ChordBand(x_start, x_end) coincides with flapwise web midplanes, so
a CapWall that ignored flapwise_* still looked flush. Under twist the webs sit
on slanted midplanes (c_station); caps clipped chordwise left a bay gap of
several millimetres on 2D-F.
"""

from __future__ import annotations

import numpy as np
import pytest

from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF
from implicit_geometry_kernel.section import (
    SectionGeometry,
    SectionOptions,
    SectionSpec,
    resolve_system_type,
)
from implicit_geometry_kernel.section.boundary.polygons import OPS
from implicit_geometry_kernel.section.boundary.walls import emit_wall_solid_polygons

_MAX_SOLID_GAP_M = 1e-6


@pytest.mark.parametrize("twist_deg", [0.0, 15.0])
def test_family_d_cap_web_solid_flush_under_twist(twist_deg):
    twist_rad = float(np.deg2rad(twist_deg))
    result = (
        SectionGeometry(
            SectionSpec.from_kwargs(
                resolve_system_type("2D-F-LC-0R"),
                AirfoilSDF.from_naca("2412", chord=1.0),
                skin_thickness=0.003,
                spar_cap_thickness=0.012,
                shear_web_thickness=0.004,
                twist_rad=twist_rad,
            ),
            SectionOptions(outputs=("properties",)),
        )
        .execute()
        .get_results()
    )
    polys = emit_wall_solid_polygons(result.section)
    webs = {
        k: polys[k]
        for k in sorted(polys)
        if str(k).startswith("shear_web_")
    }
    assert len(webs) == 2, f"expected two webs, got {list(webs)}"
    caps = {
        k: polys[k]
        for k in ("spar_cap_suction_surface", "spar_cap_pressure_surface")
    }
    assert all(not c.is_empty for c in caps.values())

    for cap_key, cap in caps.items():
        for web_key, web in webs.items():
            gap = float(OPS.distance(cap, web))
            assert gap <= _MAX_SOLID_GAP_M, (
                f"2D-F @ {twist_deg} deg: solid {cap_key} to {web_key} gap "
                f"{gap * 1e3:.3f} mm (limit {_MAX_SOLID_GAP_M * 1e3:.3f} mm). "
                "CapWall flapwise strip must match web midplanes under twist."
            )

