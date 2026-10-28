"""Smoke tests for long preferred kwargs with secondary aliases."""

from __future__ import annotations

import pytest

from blade_geometry_model.coordinate_frames import (
    TWISTED_PITCH_AXIS_BLADE_FRAME,
    ZERO_TWIST_CENTROID_AXIS_STATION_FRAME,
    ZERO_TWIST_PITCH_AXIS_STATION_FRAME,
    naca_to_S,
    naca_to_twisted_pitch_axis_blade_frame,
    naca_to_zero_twist_pitch_axis_station_frame,
)
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    build_section_with_thicknesses,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.api.kwarg_aliases import (
    resolve_shear_web_thickness,
    resolve_spar_cap_thickness,
    resolve_twist_rad,
)


def test_frame_id_constants():
    assert ZERO_TWIST_PITCH_AXIS_STATION_FRAME == "zero_twist_pitch_axis_station_frame"
    assert ZERO_TWIST_CENTROID_AXIS_STATION_FRAME == "zero_twist_centroid_axis_station_frame"
    assert TWISTED_PITCH_AXIS_BLADE_FRAME == "twisted_pitch_axis_blade_frame"


def test_geometry_transform_aliases_match():
    import numpy as np

    verts = np.array([[0.0, 0.0], [0.5, 0.05], [1.0, 0.0]], dtype=float)
    a = naca_to_zero_twist_pitch_axis_station_frame(verts, 1.0)
    b = naca_to_S(verts, 1.0)
    np.testing.assert_allclose(a, b)
    c = naca_to_twisted_pitch_axis_blade_frame(verts, 0.1, 1.0)
    assert c.shape == (3, 2)


def test_thickness_alias_conflict():
    with pytest.raises(ValueError, match="Conflicting"):
        resolve_spar_cap_thickness(spar_cap_thickness=0.01, cap_height=0.02)


def test_build_prefers_long_thickness_names():
    af = AirfoilSDF.from_naca("2412", n_points=60, chord=1.0)
    layout = resolve_system_type("1B-CN")
    section = build_section_with_thicknesses(
        af,
        layout,
        twist_rad=0.0,
        skin_thickness=0.005,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.006,
    )
    assert "skin_suction_surface" in section.labels


def test_t_star_aliases_raise_type_error():
    af = AirfoilSDF.from_naca("2412", n_points=60, chord=1.0)
    layout = resolve_system_type("1B-CN")
    with pytest.raises(TypeError):
        build_section_with_thicknesses(
            af, layout, twist_rad=0.0, t_skin=0.005, t_spar=0.008, t_web=0.006
        )


def test_twist_alias_resolution():
    assert resolve_twist_rad(twist_angle_rad=0.2) == pytest.approx(0.2)
    assert resolve_spar_cap_thickness(cap_height=0.01) == pytest.approx(0.01)
    assert resolve_shear_web_thickness(web_thickness=0.004) == pytest.approx(0.004)
