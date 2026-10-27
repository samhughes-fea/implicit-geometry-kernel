"""Axis × skin-mid offset ring lands near production web-foot endpoints."""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.midlines.junction_primitives import (
    solve_web_foot_on_skin_mid,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)


def test_web_foot_oracle_near_production_midline_ends():
    skin = 0.006
    section = build_section_view(
        AirfoilSDF.from_naca("2412", n_points=200, chord=1.0),
        resolve_system_type("2D-F"),
        twist_rad=0.0,
        skin_thickness=skin,
        spar_cap_thickness=0.014,
        shear_web_thickness=0.008,
    )
    host = getattr(section, "host", section)
    comps = dict(getattr(host, "components", {}) or {})
    webs = [(k, v) for k, v in comps.items() if str(k).startswith("shear_web")]
    assert webs, "expected shear webs on 2D-F"

    chord = float(getattr(host._af if hasattr(host, "_af") else section, "chord", 1.0) or 1.0)
    # Production feet prefer the same offset-ring oracle; allow micron-scale
    # residual from polyline discretisation of the ring.
    tol = 5e-5 * max(chord, 1.0)

    for label, web in webs:
        mid = np.asarray(web.midline_polyline(n=8), dtype=float)
        assert mid.shape[0] >= 2
        # Axis from raw IML anchors (pre-extend) through extended ends.
        p0 = np.array([web.x_top, web.y_top], dtype=float)
        p1 = np.array([web.x_bot, web.y_bot], dtype=float)
        # Extend axis beyond IML so ring hits are interior to the segment.
        direction = p0 - p1
        nrm = float(np.linalg.norm(direction))
        if nrm < 1e-15:
            continue
        direction = direction / nrm
        axis_p0 = p0 + 2.0 * skin * direction
        axis_p1 = p1 - 2.0 * skin * direction
        af = getattr(web, "_af", None)
        hits = solve_web_foot_on_skin_mid(axis_p0, axis_p1, af, skin)
        assert len(hits) >= 1, f"{label}: no axis×ring hits"
        for end in (mid[0], mid[-1]):
            gap = min(float(np.linalg.norm(np.asarray(h) - end)) for h in hits)
            assert gap <= tol, (
                f"{label}: production foot {end} is {gap*1e3:.3f} mm from "
                f"nearest oracle hit (tol={tol*1e3:.3f} mm)"
            )
