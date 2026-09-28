"""section_contract speaks coordinate_frames B/S polarity and product channels."""

from __future__ import annotations

import numpy as np
import pytest

from blade_structure_utilities.geometry_utilities.coordinate_frames import (
    blade_resultants_to_shell_channels,
    rotate_S_to_B,
)
from blade_structure_utilities.geometry_utilities.section_contract import (
    BEAM_LOAD_AXIS_XY,
    BEAM_LOAD_AXIS_YZ,
    FRAME_B,
    FRAME_S,
    FRAME_SECTION_YZ,
    BeamResultantVector,
)


def test_frame_aliases_match_canon_tags() -> None:
    assert FRAME_B == "B"
    assert FRAME_S == "S"
    assert FRAME_SECTION_YZ == FRAME_B
    assert BEAM_LOAD_AXIS_YZ == BEAM_LOAD_AXIS_XY == (0.0, 0.0)


def test_from_flap_edge_matches_shell_channels() -> None:
    ch = blade_resultants_to_shell_channels(
        M_flapwise=1000.0,
        M_edgewise=200.0,
        V_flapwise=300.0,
        V_edgewise=40.0,
        N=1.0,
        T=2.0,
    )
    vec = BeamResultantVector.from_flap_edge(
        M_flapwise=1000.0,
        M_edgewise=200.0,
        V_flapwise=300.0,
        V_edgewise=40.0,
        N=1.0,
        T=2.0,
    )
    assert vec.frame == FRAME_B
    assert vec.Mz == ch.Mz == 1000.0
    assert vec.My == ch.My == 200.0
    assert vec.Vz == ch.Vz == 300.0
    assert vec.Vy == ch.Vy == 40.0
    assert vec.about_xy == BEAM_LOAD_AXIS_XY


def test_rotated_S_to_B_about_matches_rotate_S_to_B() -> None:
    twist = 0.35
    about_S = (0.12, -0.25)
    vec_S = BeamResultantVector.from_named(
        My=10.0,
        Mz=20.0,
        Vy=3.0,
        Vz=4.0,
        frame=FRAME_S,
        about_xy=about_S,
    )
    vec_B = vec_S.rotated_S_to_B(twist)
    assert vec_B.frame == FRAME_B
    expect = rotate_S_to_B(np.asarray([about_S], dtype=float), twist)[0]
    np.testing.assert_allclose(vec_B.about_xy, expect, rtol=0.0, atol=1e-12)


def test_rotated_S_to_B_roundtrip() -> None:
    twist = 0.21
    vec_S = BeamResultantVector.from_named(
        N=1.0,
        My=11.0,
        Mz=-7.0,
        Vy=2.5,
        Vz=-3.5,
        T=0.5,
        frame=FRAME_S,
        about_xy=(0.05, 0.08),
    )
    back = vec_S.rotated_S_to_B(twist).rotated_B_to_S(twist)
    assert back.frame == FRAME_S
    np.testing.assert_allclose(back.R, vec_S.R, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(back.about_xy, vec_S.about_xy, rtol=0.0, atol=1e-12)


def test_default_from_named_tags_B() -> None:
    vec = BeamResultantVector.from_named(Mz=1.0)
    assert vec.frame == FRAME_B
    assert vec.about_xy == BEAM_LOAD_AXIS_XY


def test_transported_to_uses_canonical_xy_slots() -> None:
    vec = BeamResultantVector.from_named(N=100.0, about_xy=(0.0, 0.0))
    out = vec.transported_to(0.1, -0.2)
    assert out.about_xy == pytest.approx((0.1, -0.2))
    # dX = 0 - 0.1 = -0.1 → Mz' = Mz - dX N = 10
    # dY = 0 - (-0.2) = 0.2 → My' = My + dY N = 20
    assert out.Mz == pytest.approx(10.0)
    assert out.My == pytest.approx(20.0)


def test_transported_to_torque_uses_product_shear_map() -> None:
    """T' = T + ΔX Vy - ΔY Vz with Vz=F_XB (flap), Vy=F_YB (edge)."""
    # Pure flap shear through origin → about (0.1, 0.2): T' = -ΔY Vz = 0.2 * Vz
    flap = BeamResultantVector.from_named(Vz=1000.0, about_xy=(0.0, 0.0))
    out = flap.transported_to(0.1, 0.2)
    assert out.T == pytest.approx(0.2 * 1000.0)
    np.testing.assert_allclose(out.R[:5], [0.0, 0.0, 1000.0, 0.0, 0.0])
    # Pure edge shear → T' = ΔX Vy = -0.1 * Vy
    edge = BeamResultantVector.from_named(Vy=500.0, about_xy=(0.0, 0.0))
    out_e = edge.transported_to(0.1, 0.2)
    assert out_e.T == pytest.approx(-0.1 * 500.0)
    assert out_e.Mz == pytest.approx(0.0)
    assert out_e.My == pytest.approx(0.0)
    assert out_e.N == pytest.approx(0.0)


def test_with_about_xy_retags_only() -> None:
    vec = BeamResultantVector.from_named(N=100.0, about_xy=(0.0, 0.0))
    out = vec.with_about_xy(0.3, 0.4)
    assert out.about_xy == (0.3, 0.4)
    np.testing.assert_allclose(out.R, vec.R)


def test_about_yz_deprecated_alias_warns() -> None:
    vec = BeamResultantVector.from_named(N=1.0, about_xy=(0.1, 0.2))
    with pytest.warns(DeprecationWarning, match="about_yz is deprecated"):
        assert vec.about_yz == (0.1, 0.2)
    with pytest.warns(DeprecationWarning, match="with_about_yz is deprecated"):
        out = vec.with_about_yz(0.3, 0.4)
    assert out.about_xy == (0.3, 0.4)
