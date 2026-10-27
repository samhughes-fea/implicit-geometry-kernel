"""
Unit tests for force/shear vector B-to-S and S-to-B coordinate transforms.

Validates that transform_forces_B_to_S / transform_forces_S_to_B apply the
same 2D rotation as transform_moments_B_to_S / transform_moments_S_to_B,
and that the inverse round-trips correctly.
"""

import numpy as np
import pytest
import sys
from pathlib import Path

_workspace_root = Path(__file__).parent.parent.parent.parent.parent
if str(_workspace_root) not in sys.path:
    sys.path.insert(0, str(_workspace_root))

from blade_geometry_model.coordinate_frames import CoordinateFrameTransforms


@pytest.fixture
def transform():
    return CoordinateFrameTransforms()


class TestTransformForcesBtoS:
    """Tests for transform_forces_B_to_S."""

    def test_identity_at_zero_twist(self, transform):
        V_XS, V_YS = transform.transform_forces_B_to_S(100.0, 50.0, twist_deg=0.0)
        assert np.isclose(V_XS, 100.0)
        assert np.isclose(V_YS, 50.0)

    def test_90_deg_twist(self, transform):
        V_XS, V_YS = transform.transform_forces_B_to_S(100.0, 0.0, twist_deg=90.0)
        assert np.isclose(V_XS, 0.0, atol=1e-10)
        assert np.isclose(V_YS, 100.0, atol=1e-10)

    def test_matches_moment_transform(self, transform):
        """Force transform must produce the same rotation as moment transform."""
        V_XB, V_YB, twist = 1234.5, -678.9, 17.3
        M_YS, M_XS = transform.transform_moments_B_to_S(
            M_YB=V_XB, M_XB=V_YB, twist_deg=twist,
        )
        V_XS, V_YS = transform.transform_forces_B_to_S(
            V_XB=V_XB, V_YB=V_YB, twist_deg=twist,
        )
        assert np.isclose(V_XS, M_YS)
        assert np.isclose(V_YS, M_XS)

    def test_array_input(self, transform):
        V_XB = np.array([100.0, 200.0, 300.0])
        V_YB = np.array([10.0, 20.0, 30.0])
        twist = np.array([0.0, 45.0, 90.0])
        V_XS, V_YS = transform.transform_forces_B_to_S(V_XB, V_YB, twist)
        assert V_XS.shape == (3,)
        assert V_YS.shape == (3,)
        assert np.isclose(V_XS[0], 100.0)
        assert np.isclose(V_YS[0], 10.0)
        assert np.isclose(V_XS[2], -30.0, atol=1e-10)
        assert np.isclose(V_YS[2], 300.0, atol=1e-10)

    def test_preserves_magnitude(self, transform):
        V_XB, V_YB = 300.0, 400.0
        mag_in = np.hypot(V_XB, V_YB)
        for twist in [0.0, 15.0, 45.0, 90.0, 135.0, -30.0]:
            V_XS, V_YS = transform.transform_forces_B_to_S(V_XB, V_YB, twist)
            mag_out = np.hypot(V_XS, V_YS)
            assert np.isclose(mag_in, mag_out, rtol=1e-12), (
                f"Magnitude not preserved at twist={twist}"
            )


class TestTransformForcesStoB:
    """Tests for transform_forces_S_to_B (inverse)."""

    def test_round_trip_scalar(self, transform):
        V_XB, V_YB, twist = 500.0, -200.0, 23.7
        V_XS, V_YS = transform.transform_forces_B_to_S(V_XB, V_YB, twist)
        V_XB_rt, V_YB_rt = transform.transform_forces_S_to_B(V_XS, V_YS, twist)
        assert np.isclose(V_XB_rt, V_XB, rtol=1e-12)
        assert np.isclose(V_YB_rt, V_YB, rtol=1e-12)

    def test_round_trip_array(self, transform):
        rng = np.random.default_rng(42)
        V_XB = rng.normal(0, 1000, size=20)
        V_YB = rng.normal(0, 500, size=20)
        twist = rng.uniform(-30, 30, size=20)
        V_XS, V_YS = transform.transform_forces_B_to_S(V_XB, V_YB, twist)
        V_XB_rt, V_YB_rt = transform.transform_forces_S_to_B(V_XS, V_YS, twist)
        np.testing.assert_allclose(V_XB_rt, V_XB, rtol=1e-12)
        np.testing.assert_allclose(V_YB_rt, V_YB, rtol=1e-12)

    def test_identity_at_zero_twist(self, transform):
        V_XS, V_YS = 100.0, 50.0
        V_XB, V_YB = transform.transform_forces_S_to_B(V_XS, V_YS, twist_deg=0.0)
        assert np.isclose(V_XB, V_XS)
        assert np.isclose(V_YB, V_YS)
