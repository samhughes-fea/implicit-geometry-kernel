"""E (laminate/element) ↔ M (ply material) transform contract.

The governing identity is ``sigma_M == Q_material @ eps_M``: stress and strain
must be rotated with matrices that stay mutually consistent with the ply's own
material stiffness. Rotating both with the same matrix satisfies neither, which
is the defect these tests pin.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.coordinate_frames import (
    ply_strain_transform_E_to_M,
    ply_strain_transform_M_to_E,
    ply_stress_transform_3d_E_to_M,
    ply_stress_transform_E_to_M,
    ply_stress_transform_M_to_E,
    rotate_ply_strain_E_to_M,
    rotate_ply_stress_E_to_M,
    rotate_ply_stress_M_to_E,
    rotate_ply_vectors,
)

E1, E2, G12, NU12 = 120e9, 9e9, 4.42e9, 0.3
ANGLES = (0.0, 30.0, 45.0, -45.0, 60.0, 90.0, -73.5)


def _Q_material() -> np.ndarray:
    """Plane-stress stiffness in the ply's own material axes."""
    nu21 = NU12 * E2 / E1
    d = 1.0 - NU12 * nu21
    return np.array(
        [[E1 / d, NU12 * E2 / d, 0.0], [NU12 * E2 / d, E2 / d, 0.0], [0.0, 0.0, G12]]
    )


def _single_ply_laminate(theta_deg: float):
    from blade_structure_utilities.material_utilities import (
        LaminateDefinition,
        OrthotropicPly,
    )

    ply = OrthotropicPly(
        name="probe", E1=E1, E2=E2, G12=G12, nu12=NU12, rho=1530.0, t_ply=0.001,
        Xt=1e9, Xc=1e9, Yt=5e7, Yc=2e8, S12=5e7, Zt=0.0, S13=0.0, S23=0.0,
    )
    return LaminateDefinition(plies=[(ply, theta_deg)])


class TestTransformShape:
    def test_theta_zero_is_identity_for_both(self):
        assert np.allclose(ply_stress_transform_E_to_M(0.0), np.eye(3))
        assert np.allclose(ply_strain_transform_E_to_M(0.0), np.eye(3))

    def test_stress_and_strain_transforms_differ_off_axis(self):
        """The whole point: they are not interchangeable away from 0/90."""
        for th in (30.0, 45.0, -60.0):
            assert not np.allclose(
                ply_stress_transform_E_to_M(th), ply_strain_transform_E_to_M(th)
            ), th

    def test_shear_terms_carry_the_factor_two_differently(self):
        th = 45.0
        ts, te = ply_stress_transform_E_to_M(th), ply_strain_transform_E_to_M(th)
        # T_sigma: 2sc on the normal rows, sc on the shear row. T_eps: reversed.
        assert ts[0, 2] == pytest.approx(2.0 * te[0, 2])
        assert te[2, 0] == pytest.approx(2.0 * ts[2, 0])

    def test_scalar_theta_gives_3x3_array_theta_gives_stacked(self):
        assert ply_stress_transform_E_to_M(12.0).shape == (3, 3)
        assert ply_strain_transform_E_to_M(np.zeros(7)).shape == (7, 3, 3)


class TestInverses:
    def test_M_to_E_is_the_negative_angle_build(self):
        for th in ANGLES:
            assert np.allclose(
                ply_stress_transform_M_to_E(th), ply_stress_transform_E_to_M(-th)
            )
            assert np.allclose(
                ply_strain_transform_M_to_E(th), ply_strain_transform_E_to_M(-th)
            )

    def test_matrices_are_true_inverses(self):
        for th in ANGLES:
            for fwd, inv in (
                (ply_stress_transform_E_to_M, ply_stress_transform_M_to_E),
                (ply_strain_transform_E_to_M, ply_strain_transform_M_to_E),
            ):
                assert np.allclose(fwd(th) @ inv(th), np.eye(3), atol=1e-12), th

    def test_vector_round_trip(self):
        v = np.array([100e6, 20e6, 40e6])
        for th in ANGLES:
            assert np.allclose(rotate_ply_stress_M_to_E(rotate_ply_stress_E_to_M(v, th), th), v)


class TestPhysics:
    def test_45_degree_ply_under_uniaxial_stress(self):
        """Textbook: pure sigma_x at 45 deg splits evenly with equal shear."""
        sigma_E = np.array([100e6, 0.0, 0.0])
        got = rotate_ply_stress_E_to_M(sigma_E, 45.0)
        assert got == pytest.approx(np.array([50e6, 50e6, -50e6]))

    def test_90_degree_ply_swaps_normal_components(self):
        got = rotate_ply_stress_E_to_M(np.array([100e6, 20e6, 0.0]), 90.0)
        assert got == pytest.approx(np.array([20e6, 100e6, 0.0]), abs=1.0)

    @pytest.mark.parametrize("theta_deg", [30.0, 45.0, -60.0, 73.5])
    def test_stress_strain_consistency_against_material_stiffness(self, theta_deg):
        """sigma_M == Q_material @ eps_M — the identity that exposed the defect.

        Uses a laminate stress state with non-zero shear. After the coordinated
        ``Q_bar`` / ``T_ply`` sign fix, consistency holds under ``+theta``.
        """
        from blade_structure_utilities.material_utilities.laminate import (
            clpt_recovery_from_laminate,
        )

        lam = _single_ply_laminate(theta_deg)
        h = lam.total_thickness()
        _, _, sig_lam, eps_lam = clpt_recovery_from_laminate(
            lam, np.array([100e6 * h, 20e6 * h, 40e6 * h]), np.zeros(3)
        )
        sig_E = np.asarray(sig_lam)[0]
        eps_E = np.asarray(eps_lam)[0]

        sig_M = rotate_ply_stress_E_to_M(sig_E, theta_deg)
        eps_M = rotate_ply_strain_E_to_M(eps_E, theta_deg)

        assert sig_M == pytest.approx(_Q_material() @ eps_M, rel=1e-9)

    @pytest.mark.parametrize("theta_deg", [30.0, 45.0, -60.0])
    def test_rotating_strain_with_the_stress_matrix_is_inconsistent(self, theta_deg):
        """Guard against the defect returning: the wrong matrix must not pass."""
        from blade_structure_utilities.material_utilities.laminate import (
            clpt_recovery_from_laminate,
        )

        lam = _single_ply_laminate(theta_deg)
        h = lam.total_thickness()
        _, _, sig_lam, eps_lam = clpt_recovery_from_laminate(
            lam, np.array([100e6 * h, 20e6 * h, 40e6 * h]), np.zeros(3)
        )
        sig_M = rotate_ply_stress_E_to_M(np.asarray(sig_lam)[0], theta_deg)
        eps_M_wrong = rotate_ply_stress_E_to_M(np.asarray(eps_lam)[0], theta_deg)

        assert not np.allclose(sig_M, _Q_material() @ eps_M_wrong, rtol=1e-6)


class TestBatching:
    def test_per_ply_angles_match_a_loop(self):
        angles = np.array([0.0, 45.0, -45.0, 90.0, 30.0])
        vecs = np.arange(15, dtype=float).reshape(5, 3) * 1e6
        batched = rotate_ply_stress_E_to_M(vecs, angles)
        loop = np.stack(
            [ply_stress_transform_E_to_M(a) @ vecs[i] for i, a in enumerate(angles)]
        )
        assert np.allclose(batched, loop)

    def test_face_axis_broadcasts_with_column_angles(self):
        angles = np.array([0.0, 45.0, -45.0])
        vecs = np.arange(27, dtype=float).reshape(3, 3, 3) * 1e6  # (n_ply, n_face, 3)
        batched = rotate_ply_strain_E_to_M(vecs, angles[:, None])
        assert batched.shape == (3, 3, 3)
        for i, a in enumerate(angles):
            expected = (ply_strain_transform_E_to_M(a) @ vecs[i].T).T
            assert np.allclose(batched[i], expected)

    def test_rotate_ply_vectors_handles_both_documented_layouts(self):
        rng = np.random.default_rng(0)
        T = ply_stress_transform_E_to_M(np.array([[0.0, 45.0]]))  # (1, 2, 3, 3)

        flat = rng.normal(size=(1, 2, 3))          # (n_comp, n_ply, 3)
        assert rotate_ply_vectors(flat, T).shape == (1, 2, 3)

        faced = rng.normal(size=(1, 2, 3, 3))      # (n_comp, n_ply, n_face, 3)
        assert rotate_ply_vectors(faced, T).shape == (1, 2, 3, 3)


class TestThreeDimensionalTransforms:
    def test_in_plane_block_matches_3x3(self):
        th = 37.0
        t3 = ply_stress_transform_3d_E_to_M(th)
        t2 = ply_stress_transform_E_to_M(th)
        idx = np.ix_([0, 1, 5], [0, 1, 5])
        assert np.allclose(t3[idx], t2)

    def test_transverse_shear_round_trip(self):
        from blade_geometry_model.coordinate_frames import (
            ply_stress_transform_3d_E_to_M,
            rotate_transverse_shear_E_to_M,
        )

        tau13, tau23 = 12e6, -3e6
        m13, m23 = rotate_transverse_shear_E_to_M(tau13, tau23, 30.0)
        vec = np.array([0.0, 0.0, 0.0, tau23, tau13, 0.0])
        back = ply_stress_transform_3d_E_to_M(-30.0) @ (
            ply_stress_transform_3d_E_to_M(30.0) @ vec
        )
        assert back == pytest.approx(vec)
        assert m13 != pytest.approx(tau13) or m23 != pytest.approx(tau23)

    def test_strain_work_invariance(self):
        from blade_geometry_model.coordinate_frames import (
            ply_strain_transform_3d_E_to_M,
            ply_stress_transform_3d_E_to_M,
        )

        th = 25.0
        sigma_E = np.array([10.0, 4.0, 1.0, 0.5, -0.2, 3.0])
        # engineering strain dummy
        eps_E = np.array([0.001, 0.0004, 0.0001, 0.0002, -0.0001, 0.0006])
        ts = ply_stress_transform_3d_E_to_M(th)
        te = ply_strain_transform_3d_E_to_M(th)
        # Voigt work: σ·ε with engineering shear is invariant if T_ε^T T_σ^{-T} = I
        # equivalently σ_M · ε_M == σ_E · ε_E
        sigma_M = ts @ sigma_E
        eps_M = te @ eps_E
        assert float(sigma_M @ eps_M) == pytest.approx(float(sigma_E @ eps_E), rel=1e-12)
    def test_build_T_ply_is_the_stress_transform(self):
        lam = _single_ply_laminate(37.0)
        assert np.allclose(lam.build_T_ply()[0], ply_stress_transform_E_to_M(37.0))

    def test_build_T_ply_eps_is_the_strain_transform(self):
        lam = _single_ply_laminate(37.0)
        assert np.allclose(lam.build_T_ply_eps()[0], ply_strain_transform_E_to_M(37.0))

    def test_stiffness_bundle_carries_both(self):
        bundle = _single_ply_laminate(-22.0).stiffness_bundle()
        assert bundle.T_ply_eps is not None
        assert not np.allclose(bundle.T_ply, bundle.T_ply_eps)
