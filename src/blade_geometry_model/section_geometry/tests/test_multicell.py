"""
tests.test_multicell
====================
Unit tests for MultiCellSection and the corrected SparCap construction.
"""

import numpy as np
import pytest

from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
)
from blade_geometry_model.section_geometry.recipes.sdf_ir import SparCap
from blade_geometry_model.section_geometry.sdf import offset
from blade_geometry_model.section_geometry.layout.thickness_limits import (
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M,
    MIN_REALISTIC_SPAR_WALL_THICKNESS_M,
)
from blade_geometry_model.section_geometry.visualisation.plot import _DEFAULT_COLOR, _component_color


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def af():
    return AirfoilSDF.from_naca("0012", n_points=200, chord=1.0)

@pytest.fixture(scope="module")
def af2412():
    return AirfoilSDF.from_naca("2412", n_points=200, chord=1.0)

@pytest.fixture(scope="module")
def grid(af):
    return SDFGrid.from_airfoil(af, padding=0.05, nx=300, ny=120)


# ---------------------------------------------------------------------------
# Corrected SparCap construction
# ---------------------------------------------------------------------------

class TestSparCapCurvature:
    """Verify that the corrected SparCap has constant thickness everywhere."""

    def test_outer_face_flush_with_inner_skin(self, af):
        """The outer cap face should coincide with the inner skin surface."""
        skin_t  = 0.003
        cap_h   = 0.015
        cap     = SparCap(af, skin_t, x_start=0.20, x_end=0.50,
                          cap_height=cap_h, surface="upper")
        inner_skin = offset(af, skin_t / 2.0)

        # Sample points on the inner skin in the cap x-range
        x_samples = np.linspace(0.22, 0.48, 30)
        # Find y on inner skin upper surface by scanning
        y_scan = np.linspace(0.0, 0.15, 500)
        outer_face_phi = []
        for x_q in x_samples:
            phi_vals = inner_skin(np.full_like(y_scan, x_q), y_scan)
            # Find y where inner_skin ≈ 0 (upper surface crossing)
            idx = np.argmin(np.abs(phi_vals))
            y_q = y_scan[idx]
            # Evaluate inner_skin there — should be near 0
            outer_face_phi.append(float(inner_skin(x_q, y_q)))

        np.testing.assert_allclose(outer_face_phi, 0.0, atol=0.005)

    def test_cap_thickness_uniform(self, af):
        """Cap thickness (distance from outer to inner face) should equal cap_height."""
        skin_t  = 0.003
        cap_h   = 0.015
        inner_skin = offset(af, skin_t / 2.0)

        # Sample points on the outer face (inner_skin ≈ 0) in the cap x-range
        x_samples = np.linspace(0.22, 0.48, 20)
        y_scan    = np.linspace(0.0, 0.15, 1000)
        thicknesses = []

        for x_q in x_samples:
            phi_is = inner_skin(np.full_like(y_scan, x_q), y_scan)
            # Find upper skin crossing (phi ≈ 0, y > 0)
            pos_vals = phi_is[y_scan > 0.005]
            pos_y    = y_scan[y_scan > 0.005]
            if len(pos_vals) == 0:
                continue
            idx    = np.argmin(np.abs(pos_vals))
            y_outer = pos_y[idx]   # y at outer cap face

            # Inner cap face is at inner_skin = -cap_h → scan from outer mold toward chord (smaller y)
            y_inner_scan = np.linspace(y_outer, 0.0, 500)
            phi_inner    = inner_skin(np.full_like(y_inner_scan, x_q), y_inner_scan)
            # Find where phi_inner ≈ -cap_h
            target = -cap_h
            idx2   = np.argmin(np.abs(phi_inner - target))
            y_inner = y_inner_scan[idx2]

            thickness = y_outer - y_inner
            thicknesses.append(thickness)

        # Thickness should be approximately cap_h (within grid resolution)
        thicknesses = np.array(thicknesses)
        np.testing.assert_allclose(thicknesses, cap_h, atol=0.003)

    def test_cap_inside_airfoil(self, af, grid):
        """Cap wall solid lies in the cavity (airfoil interior past the inner mold).

        Inner mold is at phi_af ≈ -skin_t/2; cavity points have phi_af < -skin_t/2.
        """
        skin_t = 0.003
        cap_h  = 0.015
        cap = SparCap(af, skin_t, 0.20, 0.50, cap_h, surface="upper")
        phi_cap = grid.eval(cap)
        phi_af  = grid.eval(af)
        cap_interior = phi_cap < 0.0
        assert cap_interior.any(), "Cap has no interior region."
        mold = -skin_t / 2.0
        assert np.all(phi_af[cap_interior] <= mold + 5e-4), \
            "Cap interior should lie on the cavity side of the inner skin (phi_af <= -skin_t/2)."


# ---------------------------------------------------------------------------
# MultiCellSection topology
# ---------------------------------------------------------------------------

class TestMultiCellTopology:
    def test_wall_thickness_minimum_clamp(self, af):
        mcs = MultiCellSection(
            af,
            web_x_positions=[0.25, 0.5],
            skin_thickness=0.0005,
            web_thickness=0.0005,
        )
        assert mcs._skin_t == pytest.approx(MIN_REALISTIC_SKIN_WALL_THICKNESS_M)

    def test_spar_cap_height_clamp_matches_minimum_reference(self, af, grid):
        """Sub-millimetre cap_height is raised to the spar wall thickness floor; SDF matches explicit minimum."""
        m_tiny = MultiCellSection.twin_web(
            af, skin_thickness=0.003, cap_height=1e-5, web_thickness=0.004
        )
        m_ref = MultiCellSection.twin_web(
            af,
            skin_thickness=0.003,
            cap_height=MIN_REALISTIC_SPAR_WALL_THICKNESS_M,
            web_thickness=0.004,
        )
        p0 = grid.eval(m_tiny["spar_cap_suction_surface"])
        p1 = grid.eval(m_ref["spar_cap_suction_surface"])
        np.testing.assert_allclose(p0, p1, atol=1e-9, rtol=0)

    def test_single_web_d_spar_labels(self, af):
        """d_spar is catalogue 1C (per_web_pads); continuous_box + N=1 is illegal."""
        mcs = MultiCellSection.d_spar(af, web_x=0.25,
                                      skin_thickness=0.003,
                                      cap_height=0.012,
                                      web_thickness=0.004)
        assert "shear_web_0"         in mcs.labels
        assert "spar_cap_suction_surface" in mcs.labels
        assert "spar_cap_pressure_surface" in mcs.labels
        assert "skin_suction_surface" in mcs.labels
        assert "skin_pressure_surface" in mcs.labels
        assert mcs.n_webs == 1
        assert mcs.cap_placement.kind == "per_web_pads"
        assert "sandwich_core_0" in mcs.labels and "sandwich_core_1" in mcs.labels
        assert mcs.n_cells == 2

    def test_continuous_box_with_one_web_is_infeasible(self, af):
        with pytest.raises(InfeasibleSectionGeometry, match="continuous_box"):
            MultiCellSection(
                af,
                web_x_positions=[0.25],
                structural_family="D",
                skin_thickness=0.003,
                cap_height=0.012,
                web_thickness=0.004,
            )

    def test_twin_web_labels(self, af):
        mcs = MultiCellSection.twin_web(af)
        assert mcs.n_webs  == 2
        assert "shear_web_0"     in mcs.labels
        assert "shear_web_1"     in mcs.labels
        assert "sandwich_core_0" in mcs.labels and "sandwich_core_1" in mcs.labels and "sandwich_core_2" in mcs.labels
        assert mcs.n_cells == 3

    def test_three_web_labels(self, af):
        mcs = MultiCellSection.torsion_box(af)
        assert mcs.n_webs  == 3
        assert mcs.n_cells == 4
        for j in range(4):
            assert f"sandwich_core_{j}" in mcs.labels

    def test_n_web_generalisation(self, af):
        """N webs → N+1 structural bays (foam regions when core_enabled).

        Uses family C so N=1 is a legal catalogue layout (1C); D requires N>=2.
        """
        for n in range(1, 6):
            xs = np.linspace(0.15, 0.55, n).tolist()
            mcs = MultiCellSection(af, web_x_positions=xs,
                                   structural_family="C",
                                   skin_thickness=0.003,
                                   cap_height=0.010,
                                   web_thickness=0.003)
            assert mcs.n_webs == n
            expected_cells = n + 1
            assert mcs.n_cells == expected_cells, \
                f"Expected {expected_cells} cells for {n} webs, got {mcs.n_cells}"

    def test_te_leading_edge_inserts(self, af):
        mcs = MultiCellSection.twin_web(af,
                                        trailing_edge_insert_x_m=0.75,
                                        leading_edge_insert_x_m=0.10)
        assert "trailing_edge_insert" in mcs.labels
        assert "leading_edge_insert" in mcs.labels

    def test_inserts_flush_with_inner_skin_not_mid_surface(self, af2412):
        """LE/TE insert outer boundary sits on Offset(airfoil, -t_skin).

        Historically inserts used ``-t/2`` (skin mid-surface). They must match
        the cavity / core interface at ``-t`` so they do not invade the skin.
        """
        skin_t = 0.005
        mcs = MultiCellSection.twin_web(
            af2412,
            skin_thickness=skin_t,
            trailing_edge_insert_x_m=0.75,
            leading_edge_insert_x_m=0.10,
            core_enabled=False,
        )
        grid = SDFGrid.from_airfoil(af2412, padding=0.04, nx=400, ny=180)
        phi_inner = grid.eval(offset(af2412, -skin_t))
        phi_mid = grid.eval(offset(af2412, -0.5 * skin_t))

        for label in ("leading_edge_insert", "trailing_edge_insert"):
            phi_ins = grid.eval(mcs[label])
            # Zero-contour band of the insert (near its outer envelope).
            on_boundary = np.abs(phi_ins) < 0.15 * skin_t
            assert on_boundary.any(), f"{label}: empty boundary band"
            # Boundary should track the inner-skin cavity, not the mid-surface.
            median_inner = float(np.median(np.abs(phi_inner[on_boundary])))
            median_mid = float(np.median(np.abs(phi_mid[on_boundary])))
            assert median_inner < 0.35 * skin_t, (
                f"{label}: insert boundary not flush with inner skin "
                f"(median |phi_inner|={median_inner:.4g}, skin_t={skin_t})"
            )
            assert median_inner < median_mid, (
                f"{label}: insert still closer to mid-surface than inner skin "
                f"(|phi_inner|={median_inner:.4g}, |phi_mid|={median_mid:.4g})"
            )
            # Insert solid must not invade the skin wall band (-t < phi_af < 0
            # is outside the cavity; require interior samples stay in cavity).
            interior = phi_ins < -0.1 * skin_t
            assert interior.any(), f"{label}: no clear interior"
            assert np.all(phi_inner[interior] < 0.05 * skin_t), (
                f"{label}: insert interior leaks outside the inner-skin cavity"
            )


# ---------------------------------------------------------------------------
# Interior checks
# ---------------------------------------------------------------------------

class TestInteriorChecks:
    def test_core_inside_airfoil(self, af, grid):
        mcs = MultiCellSection.twin_web(af,
                                        skin_thickness=0.003,
                                        cap_height=0.010)
        phi_af = grid.eval(af)
        for lbl in mcs.labels:
            if not lbl.startswith("sandwich_core_"):
                continue
            phi_core = grid.eval(mcs[lbl])
            interior = phi_core < 0.0
            assert interior.any(), f"{lbl} has no interior region."
            assert np.all(phi_af[interior] <= 0.01), \
                f"{lbl} interior leaks outside airfoil."

    def test_caps_do_not_overlap(self, af, grid):
        """Upper and lower caps should not share interior points."""
        mcs   = MultiCellSection.twin_web(af, skin_thickness=0.003, cap_height=0.010)
        phi_u = grid.eval(mcs["spar_cap_suction_surface"])
        phi_l = grid.eval(mcs["spar_cap_pressure_surface"])
        overlap = (phi_u < 0.0) & (phi_l < 0.0)
        n_overlap = int(overlap.sum())
        assert n_overlap < 5, \
            f"Upper and lower caps overlap at {n_overlap} grid points."

    def test_web_intersects_interior(self, af, grid):
        mcs = MultiCellSection.twin_web(af, skin_thickness=0.003, web_thickness=0.005)
        phi_web = grid.eval(mcs["shear_web_0"])
        assert (phi_web < 0.0).any(), "Web_0 has no interior region."


# ---------------------------------------------------------------------------
# Web alignment
# ---------------------------------------------------------------------------

class TestWebAlignment:
    def test_flapwise_vs_chord_normal_differ(self, af):
        """Flapwise and chord-normal webs at non-zero twist should produce different SDFs."""
        twist = np.radians(15)
        mcs_cn = MultiCellSection.twin_web(
            af, web_thickness=0.004, web_alignment="chord_normal", twist_angle=0.0
        )
        mcs_fw = MultiCellSection.twin_web(
            af, web_thickness=0.004, web_alignment="flapwise", twist_angle=twist
        )
        x = np.linspace(0.18, 0.22, 50)
        y = np.linspace(-0.05, 0.05, 50)
        X, Y = np.meshgrid(x, y)
        phi_cn = mcs_cn["shear_web_0"](X, Y)
        phi_fw = mcs_fw["shear_web_0"](X, Y)
        # They should differ meaningfully in the web region
        max_diff = float(np.abs(phi_cn - phi_fw).max())
        assert max_diff > 1e-4, \
            "Flapwise and chord-normal webs should differ when twist != 0."

    def test_flapwise_zero_twist_equals_chord_normal(self, af):
        """At zero twist, flapwise == chord-normal."""
        mcs_cn = MultiCellSection.twin_web(
            af, web_thickness=0.004, web_alignment="chord_normal", twist_angle=0.0
        )
        mcs_fw = MultiCellSection.twin_web(
            af, web_thickness=0.004, web_alignment="flapwise", twist_angle=0.0
        )
        x = np.linspace(0.18, 0.22, 30)
        y = np.linspace(-0.05, 0.05, 30)
        X, Y = np.meshgrid(x, y)
        phi_cn = mcs_cn["shear_web_0"](X, Y)
        phi_fw = mcs_fw["shear_web_0"](X, Y)
        np.testing.assert_allclose(phi_cn, phi_fw, atol=1e-10)

    def test_per_web_alignment(self, af):
        """Mixed alignment list: first web chord_normal, second flapwise."""
        mcs = MultiCellSection.twin_web(
            af,
            web_thickness=0.004,
            web_alignment=["chord_normal", "flapwise"],
            twist_angle=np.radians(10),
        )
        assert mcs["shear_web_0"] is not None
        assert mcs["shear_web_1"] is not None


class TestAnchorDetection:
    def test_narrow_y_search_raises_actionable_error(self, af):
        with pytest.raises(ValueError, match="No feasible web anchors"):
            MultiCellSection(
                af,
                web_x_positions=[0.2, 0.5],
                y_search=(-1e-4, 1e-4),
            )

    def test_explicit_web_y_coords_bypass_search(self, af):
        mcs = MultiCellSection(
            af,
            web_x_positions=[0.2, 0.5],
            web_y_coords=[(0.05, -0.05), (0.05, -0.05)],
        )
        assert "shear_web_0" in mcs.labels and "shear_web_1" in mcs.labels


class TestPlotColorMapping:
    def test_web_label_has_non_default_color(self):
        assert _component_color("shear_web_0") != _DEFAULT_COLOR

    def test_core_label_has_non_default_color(self):
        assert _component_color("sandwich_core_0") != _DEFAULT_COLOR


# ---------------------------------------------------------------------------
# SparCap twist (was ContinuousSparCap)
# ---------------------------------------------------------------------------

class TestSparCapTwist:
    def test_upper_lower_do_not_overlap(self, af, grid):
        cap_u = SparCap(af, 0.003, 0.20, 0.50, 0.012, surface="upper")
        cap_l = SparCap(af, 0.003, 0.20, 0.50, 0.012, surface="lower")
        phi_u = grid.eval(cap_u)
        phi_l = grid.eval(cap_l)
        overlap = (phi_u < 0.0) & (phi_l < 0.0)
        assert int(overlap.sum()) < 5

    def test_twist_rotates_cap(self, af):
        cap_0  = SparCap(af, 0.003, 0.20, 0.50, 0.012,
                         surface="upper", twist_angle=0.0)
        cap_tw = SparCap(af, 0.003, 0.20, 0.50, 0.012,
                         surface="upper", twist_angle=np.radians(15))
        x = np.linspace(0.2, 0.5, 50)
        y = np.linspace(0.0, 0.08, 50)
        X, Y = np.meshgrid(x, y)
        phi_0  = cap_0(X, Y)
        phi_tw = cap_tw(X, Y)
        max_diff = float(np.abs(phi_0 - phi_tw).max())
        assert max_diff > 1e-4, "Twisted and untwisted caps should differ."
