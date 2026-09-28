"""Frame-invariance harness: planar rotation in B."""

from __future__ import annotations

import numpy as np
import pytest

from blade_structure_utilities.geometry_utilities.coordinate_frames.testing import (
    FRAME_B,
    assert_planar_pierce_rotates,
    assert_scalar_invariant_under_rotation,
    rotate_points_in_b,
)


def test_rotate_points_in_b_round_trip():
    pts = np.array([[1.0, 0.0], [0.0, 2.0], [-1.5, 0.5]])
    back = rotate_points_in_b(rotate_points_in_b(pts, 0.4), -0.4)
    assert np.allclose(back, pts)


def test_shear_centre_and_gj_are_frame_invariant_on_two_cell_box():
    from blade_structure_utilities.geometry_utilities.section_properties.reference_axes import (
        shear_centre_from_midlines,
    )
    from blade_structure_utilities.geometry_utilities.section_properties.torsion import (
        closed_section_torsional_stiffness,
    )

    t = 0.01
    # Asymmetric box: thicker nose wall shifts the shear centre off the origin
    # so planar rotation is observable (a symmetric box would stay at 0).
    midlines = [
        {"label": "skin_top:left", "kind": "skin", "xy": [(-1.0, 0.5), (0.0, 0.5)], "thickness_m": t},
        {"label": "skin_top:right", "kind": "skin", "xy": [(0.0, 0.5), (1.0, 0.5)], "thickness_m": t},
        {"label": "skin_bottom:left", "kind": "skin", "xy": [(-1.0, -0.5), (0.0, -0.5)], "thickness_m": t},
        {"label": "skin_bottom:right", "kind": "skin", "xy": [(0.0, -0.5), (1.0, -0.5)], "thickness_m": t},
        {"label": "skin_nose", "kind": "skin", "xy": [(-1.0, -0.5), (-1.0, 0.5)], "thickness_m": 4.0 * t},
        {"label": "skin_tail", "kind": "skin", "xy": [(1.0, -0.5), (1.0, 0.5)], "thickness_m": t},
        {"label": "shear_web_0", "kind": "web", "xy": [(0.0, -0.5), (0.0, 0.5)], "thickness_m": t},
    ]
    loops = [
        {
            "xy": [(-1.0, -0.5), (0.0, -0.5), (0.0, 0.5), (-1.0, 0.5), (-1.0, -0.5)],
            "area": 1.0,
            "wall_labels": ("skin_bottom:left", "shear_web_0", "skin_top:left", "skin_nose"),
        },
        {
            "xy": [(0.0, -0.5), (1.0, -0.5), (1.0, 0.5), (0.0, 0.5), (0.0, -0.5)],
            "area": 1.0,
            "wall_labels": ("skin_bottom:right", "skin_tail", "skin_top:right", "shear_web_0"),
        },
    ]

    def pierce(mids):
        res = shear_centre_from_midlines(mids, cell_wall_loops=loops)
        assert res.meta is not None and res.meta.get("frame") == FRAME_B
        return (float(res.xs), float(res.ys))

    def gj(mids):
        return float(
            closed_section_torsional_stiffness(
                mids,
                cell_wall_loops=loops,
                shear_modulus_by_component={"skin": 1.0, "web": 1.0},
            ).GJ
        )

    x0, y0 = pierce(midlines)
    assert abs(x0) > 1e-3 or abs(y0) > 1e-3
    assert_planar_pierce_rotates(midlines, pierce=pierce)
    assert_scalar_invariant_under_rotation(midlines, scalar=gj)


def test_require_frame_tag_on_harness():
    with pytest.raises(ValueError):
        rotate_points_in_b(np.zeros((2, 2)), 0.1, frame="not-a-frame")
