"""Family-D box caps must stay flush against both outer webs under twist (SG-011).

At twist=0, ChordBand(x_start, x_end) coincides with flapwise web midplanes, so
a CapWall that ignored flapwise_* still looked flush. Under twist the webs sit
on slanted midplanes (c_station); caps clipped chordwise left a bay gap of
several millimetres on 2D-F.
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

_MAX_SOLID_GAP_M = 1e-6
_MAX_MIDLINE_GAP_M = 1e-6


def _point_to_polyline_min(point: np.ndarray, poly: np.ndarray) -> float:
    """Min Euclidean distance from point to an open polyline (vertices + segments)."""
    p = np.asarray(point, dtype=float).reshape(2)
    xy = np.asarray(poly, dtype=float)
    if xy.shape[0] < 2:
        return float(np.linalg.norm(xy[0] - p)) if xy.shape[0] else float("nan")
    d_min = float(np.linalg.norm(xy - p, axis=1).min())
    for a, b in zip(xy[:-1], xy[1:]):
        ab = b - a
        lab2 = float(np.dot(ab, ab))
        if lab2 < 1e-30:
            continue
        t = float(np.clip(np.dot(p - a, ab) / lab2, 0.0, 1.0))
        d_min = min(d_min, float(np.linalg.norm(a + t * ab - p)))
    return d_min


@pytest.mark.parametrize("twist_deg", [0.0, 15.0])
def test_family_d_cap_web_solid_and_midline_flush_under_twist(twist_deg):
    twist_rad = float(np.deg2rad(twist_deg))
    result = (
        SectionGeometryRunner(
            layout=resolve_system_type("2D-F"),
            airfoil=AirfoilSDF.from_naca("2412", chord=1.0),
            skin_thickness=0.003,
            spar_cap_thickness=0.012,
            shear_web_thickness=0.004,
            twist_rad=twist_rad,
            outputs=("properties", "midlines"),
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
            gap = float(cap.distance(web))
            assert gap <= _MAX_SOLID_GAP_M, (
                f"2D-F @ {twist_deg} deg: solid {cap_key} to {web_key} gap "
                f"{gap * 1e3:.3f} mm (limit {_MAX_SOLID_GAP_M * 1e3:.3f} mm). "
                "CapWall flapwise strip must match web midplanes under twist (SG-011)."
            )

    web_mid = [
        np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for s in (result.midline_strips or ())
        if s.kind == "web"
    ]
    assert len(web_mid) == 2
    for strip in result.midline_strips or ():
        if strip.kind != "cap":
            continue
        pts = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for end_i, end in enumerate((pts[0], pts[-1])):
            d_web = min(_point_to_polyline_min(end, w) for w in web_mid)
            assert d_web <= _MAX_MIDLINE_GAP_M, (
                f"2D-F @ {twist_deg} deg: midline {strip.label} end{end_i} is "
                f"{d_web * 1e3:.3f} mm from nearest web mid-surface "
                f"(limit {_MAX_MIDLINE_GAP_M * 1e3:.3f} mm). SG-011."
            )
