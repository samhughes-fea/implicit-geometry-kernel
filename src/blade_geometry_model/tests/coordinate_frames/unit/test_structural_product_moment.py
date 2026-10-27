"""The structural ``Ixy`` closure must stay admissible at any edgewise stiffness.

The stress model is handed two scalars per station -- ``EI_flap_B`` and
``EI_edge_B`` -- and needs a full inertia tensor, so the product moment has to be
closed from the outer geometry. The closure it used to apply scaled ``Ixy`` by the
**flapwise** ratio alone, which is correct only while the section is described by a
single scale parameter.

Once ``EI_edge`` is selected in its own right that assumption is false, and it fails
in a way that does not look like a modelling error: the implied correlation
``rho = Ixy / sqrt(Ixx * Iyy)`` grows without bound as ``EI_edge`` falls, so the
B->S transform rejects the section as physically impossible. Measured on a real
3.2 m blade, that put the softest "valid" edgewise stiffness at up to 2.6x the
flapwise one -- at a station whose *geometric* ``I_edge/I_flap`` is 4.4. The
rejection reflected the closure, not the aerofoil, and it took the whole ULS floor
solve down with a message about section geometry.

These tests pin both halves: that the closure is admissible at any positive pair,
and that it still agrees with the old flapwise-only form wherever that form was
right.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.coordinate_frames import (
    structural_product_moment,
)
from blade_geometry_model.coordinate_frames.section_frames.tensors import (
    transform_second_moments_B_to_S,
)

#: Spans "far softer than flapwise" to "far stiffer", covering the whole range the
#: edgewise ULS bracket searches over.
_EDGE_OVER_FLAP = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4)


def _geometric_tensor(twist_deg: float) -> tuple[float, float, float]:
    """A plausible asymmetric B-frame section: ``Ixx > Iyy`` with real ``Ixy``."""
    Ixx, Iyy = 4.0e-5, 8.0e-6
    # Correlation of 0.6 -- firmly asymmetric, and safely positive-definite.
    Ixy = 0.6 * np.sqrt(Ixx * Iyy)
    return Ixx, Iyy, Ixy


class TestAgreesWithTheFormItReplaces:
    def test_identical_ratios_cost_no_ulp(self):
        """When both ratios land on the *same float*, the answer is exactly that float.

        This is why the closure is written ``flap * sqrt(edge / flap)`` and not
        ``sqrt(flap * edge)``: the quotient is exactly 1.0, ``sqrt(1.0)`` is exactly
        1.0, and nothing is lost. The sqrt-of-a-square form would return
        ``sqrt(flap**2)`` and shed a ULP for no reason.
        """
        rng = np.random.default_rng(20260820)
        # Equal geometric and equal structural moments force both divisions to
        # produce bit-identical results, which is the condition being pinned.
        I_geom = rng.uniform(1e-7, 1e-4, 512)
        I_str = rng.uniform(1e-9, 1e-3, 512)
        Ixy_geom = rng.uniform(-1e-5, 1e-5, 512)

        got = structural_product_moment(
            Ixy_B_geom=Ixy_geom,
            Ixx_B_geom=I_geom,
            Iyy_B_geom=I_geom,
            Ixx_B_structural=I_str,
            Iyy_B_structural=I_str,
        )
        np.testing.assert_array_equal(got, Ixy_geom * (I_str / I_geom))

    def test_the_ratio_fallback_agrees_to_within_a_few_ulp(self):
        """The realistic version of the case above, which is *not* bit-exact.

        The fallback reaches the edgewise ratio as
        ``EI_flap * (Ixx_geom / Iyy_geom) / Ixx_geom`` where the flapwise one is
        ``EI_flap / Iyy_geom`` -- three roundings against one, so the two ratios
        agree only to within a ULP or so and the geometric mean inherits that. The
        agreement is close enough to be irrelevant physically and is asserted at
        that strength rather than claimed to be exact.
        """
        rng = np.random.default_rng(20260820)
        Ixx_geom = rng.uniform(1e-6, 1e-4, 512)
        Iyy_geom = rng.uniform(1e-7, 1e-5, 512)
        Ixy_geom = rng.uniform(-1e-5, 1e-5, 512)
        Iyy_str = rng.uniform(1e-9, 1e-3, 512)
        # Exactly how `_compute_EI_edge_B_from_ratio` builds it.
        Ixx_str = Iyy_str * (Ixx_geom / Iyy_geom)

        got = structural_product_moment(
            Ixy_B_geom=Ixy_geom,
            Ixx_B_geom=Ixx_geom,
            Iyy_B_geom=Iyy_geom,
            Ixx_B_structural=Ixx_str,
            Iyy_B_structural=Iyy_str,
        )
        old_flapwise_only = Ixy_geom * (Iyy_str / Iyy_geom)
        np.testing.assert_allclose(got, old_flapwise_only, rtol=1e-14)

    def test_it_is_the_geometric_mean_for_independent_stiffnesses(self):
        rng = np.random.default_rng(7)
        Ixx_geom = rng.uniform(1e-6, 1e-4, 256)
        Iyy_geom = rng.uniform(1e-7, 1e-5, 256)
        Ixy_geom = rng.uniform(-1e-5, 1e-5, 256)
        Ixx_str = rng.uniform(1e-9, 1e-3, 256)
        Iyy_str = rng.uniform(1e-9, 1e-3, 256)

        got = structural_product_moment(
            Ixy_B_geom=Ixy_geom,
            Ixx_B_geom=Ixx_geom,
            Iyy_B_geom=Iyy_geom,
            Ixx_B_structural=Ixx_str,
            Iyy_B_structural=Iyy_str,
        )
        want = Ixy_geom * np.sqrt((Ixx_str / Ixx_geom) * (Iyy_str / Iyy_geom))
        np.testing.assert_allclose(got, want, rtol=1e-14)


class TestItIsAdmissibleByConstruction:
    def test_the_correlation_coefficient_is_preserved(self):
        """``rho`` is what positive-definiteness depends on, so it is the invariant."""
        rng = np.random.default_rng(11)
        Ixx_geom = rng.uniform(1e-6, 1e-4, 256)
        Iyy_geom = rng.uniform(1e-7, 1e-5, 256)
        rho_geom = rng.uniform(-0.95, 0.95, 256)
        Ixy_geom = rho_geom * np.sqrt(Ixx_geom * Iyy_geom)
        Ixx_str = rng.uniform(1e-9, 1e-3, 256)
        Iyy_str = rng.uniform(1e-9, 1e-3, 256)

        Ixy_str = structural_product_moment(
            Ixy_B_geom=Ixy_geom,
            Ixx_B_geom=Ixx_geom,
            Iyy_B_geom=Iyy_geom,
            Ixx_B_structural=Ixx_str,
            Iyy_B_structural=Iyy_str,
        )
        np.testing.assert_allclose(
            Ixy_str / np.sqrt(Ixx_str * Iyy_str), rho_geom, rtol=1e-13
        )

    @pytest.mark.parametrize("twist_deg", [0.0, 5.0, 15.0, 34.1, -20.0])
    @pytest.mark.parametrize("edge_over_flap", _EDGE_OVER_FLAP)
    def test_the_B_to_S_transform_accepts_any_edgewise_stiffness(
        self, twist_deg, edge_over_flap
    ):
        """The regression guard for the failure this closure was written to fix.

        34.1 deg is the measured root twist on ``new_3.2m_rotor``, the station where
        the old closure rejected everything below ``2.6 x EI_flap``.
        """
        Ixx_geom, Iyy_geom, Ixy_geom = _geometric_tensor(twist_deg)
        Iyy_str = 1.234e-5
        Ixx_str = Iyy_str * edge_over_flap

        Ixy_str = structural_product_moment(
            Ixy_B_geom=Ixy_geom,
            Ixx_B_geom=Ixx_geom,
            Iyy_B_geom=Iyy_geom,
            Ixx_B_structural=Ixx_str,
            Iyy_B_structural=Iyy_str,
        )
        # Raises rather than returning a sentinel when the tensor is inadmissible.
        Ixx_S, Iyy_S, _ = transform_second_moments_B_to_S(
            Ixx_B=Ixx_str, Iyy_B=Iyy_str, Ixy_B=Ixy_str, twist_deg=twist_deg
        )
        assert Ixx_S > 0.0 and Iyy_S > 0.0

    @pytest.mark.parametrize("edge_over_flap", _EDGE_OVER_FLAP)
    def test_the_old_flapwise_only_closure_is_what_used_to_fail(self, edge_over_flap):
        """Pins *why* the fix was needed, so a revert cannot pass quietly.

        Scaling ``Ixy`` by the flapwise ratio alone leaves the product moment
        untouched as ``EI_edge`` falls, so a soft enough edgewise stiffness stops
        describing a section at all. If this ever stops failing, the premise of the
        closure has changed and its docstring is stale.
        """
        twist_deg = 34.1
        Ixx_geom, Iyy_geom, Ixy_geom = _geometric_tensor(twist_deg)
        Iyy_str = 1.234e-5
        Ixx_str = Iyy_str * edge_over_flap
        Ixy_str_old = Ixy_geom * (Iyy_str / Iyy_geom)  # the form that was replaced

        implied_rho = Ixy_str_old / np.sqrt(Ixx_str * Iyy_str)
        if edge_over_flap < 1e-2:
            assert abs(implied_rho) > 1.0, (
                "the old closure should imply an impossible correlation at a soft "
                "edgewise stiffness; if it no longer does, this test's premise moved"
            )
            with pytest.raises(ValueError):
                transform_second_moments_B_to_S(
                    Ixx_B=Ixx_str,
                    Iyy_B=Iyy_str,
                    Ixy_B=Ixy_str_old,
                    twist_deg=twist_deg,
                )


class TestDegenerateGeometry:
    def test_a_vanishing_geometric_Ixx_falls_back_to_the_flapwise_form(self):
        """Those stations kept the historical answer; nothing new is asserted there."""
        got = structural_product_moment(
            Ixy_B_geom=np.array([1.0e-6]),
            Ixx_B_geom=np.array([0.0]),
            Iyy_B_geom=np.array([2.0e-6]),
            Ixx_B_structural=np.array([5.0e-6]),
            Iyy_B_structural=np.array([4.0e-6]),
        )
        np.testing.assert_array_equal(got, np.array([1.0e-6 * (4.0e-6 / 2.0e-6)]))

    def test_a_zero_flapwise_stiffness_gives_a_zero_product_moment(self):
        """``EI_flap = 0`` is the free tip; the geometric mean there is zero, not NaN."""
        got = structural_product_moment(
            Ixy_B_geom=np.array([1.0e-6]),
            Ixx_B_geom=np.array([4.0e-5]),
            Iyy_B_geom=np.array([8.0e-6]),
            Ixx_B_structural=np.array([1.0e-5]),
            Iyy_B_structural=np.array([0.0]),
        )
        assert np.all(np.isfinite(got))
        np.testing.assert_array_equal(got, np.array([0.0]))
