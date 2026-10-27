"""Hybrid vs quadtree parity across taxonomy; overlap policy; insert retention."""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.constructive import (
    classify_property_labels,
    is_constructive_label,
    is_core_label,
    is_wall_solid_label,
    is_residual_label,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    props_backend_from_env,
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.sections.multicell import (
    MultiCellSection,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF


def _assert_props_close(hy, qt, *, wall_solid: bool, le_te: bool = False):
    a_qt = float(qt["area"])
    if a_qt < 1e-5:
        return
    if le_te:
        rtol, atol = 0.15, 1e-3
    elif wall_solid:
        rtol, atol = 1e-3, 5e-4
    else:
        rtol, atol = 5e-3, 5e-4
    np.testing.assert_allclose(hy["area"], a_qt, rtol=rtol, atol=atol)
    if wall_solid:
        np.testing.assert_allclose(hy["cx"], qt["cx"], atol=0.02)
        np.testing.assert_allclose(hy["cy"], qt["cy"], atol=0.02)
        if np.isfinite(qt["Ixx"]):
            # Small second moments (thin webs/caps) need a larger absolute floor.
            np.testing.assert_allclose(hy["Ixx"], qt["Ixx"], rtol=0.05, atol=2e-5)
            np.testing.assert_allclose(hy["Iyy"], qt["Iyy"], rtol=0.05, atol=2e-5)


@pytest.mark.parametrize(
    "layout_key,twist_deg",
    [
        ("2D-CN", 0.0),
        ("2D-F", 0.0),
        ("2C-CN", 0.0),
        ("1B-CN", 0.0),
        ("2D-CN", 15.0),
        ("0A", 0.0),
    ],
)
def test_hybrid_vs_quadtree_taxonomy(layout_key: str, twist_deg: float):
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    twist = float(np.deg2rad(twist_deg))
    section = build_section_view(
        airfoil, resolve_system_type(layout_key), twist_angle_rad=twist
    )
    af_spec = airfoil.rotate(twist) if abs(twist) > 1e-12 else airfoil
    spec = quadtree_spec_for_airfoil(af_spec, chord=1.0, max_depth=7, max_cells=120000)
    hy, hy_meta = compute_section_properties(
        section, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=af_spec
    )
    qt, _ = compute_section_properties(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=af_spec
    )
    assert hy_meta["backend"] == "hybrid"
    # Cores and inserts are constructive; residual list empty for standard builds
    for lab in section.labels:
        key = str(getattr(lab, "value", lab))
        assert key in hy
        assert key in qt
        if key in {"leading_edge_insert", "trailing_edge_insert"}:
            assert key in hy_meta["constructive_labels"]
            assert key not in hy_meta["residual_labels"]
            _assert_props_close(hy[key], qt[key], wall_solid=True)
        if is_core_label(key):
            assert key in hy_meta["constructive_labels"]
            assert key not in hy_meta["residual_labels"]
            _assert_props_close(hy[key], qt[key], wall_solid=True)
        elif is_wall_solid_label(key):
            _assert_props_close(hy[key], qt[key], wall_solid=True)
        elif is_residual_label(key):
            _assert_props_close(
                hy[key],
                qt[key],
                wall_solid=False,
                le_te=key in {"leading_edge_insert", "trailing_edge_insert"},
            )


def test_overlap_policy_independent_label_moments():
    """Hybrid must not unique-partition skin∩cap; per-label areas stay independent."""
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=120000)
    hy, _ = compute_section_properties(
        section, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    qt, _ = compute_section_properties(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    # Sum of wall-solid areas need not equal solid union; both backends agree per label.
    for lab in ("skin_suction_surface", "skin_pressure_surface", "spar_cap_suction_surface", "spar_cap_pressure_surface"):
        if lab in hy and lab in qt and qt[lab]["area"] > 1e-5:
            np.testing.assert_allclose(hy[lab]["area"], qt[lab]["area"], rtol=1e-3, atol=5e-4)


def test_inserts_are_constructive_classifier():
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    mcs = MultiCellSection.twin_web(
        airfoil,
        skin_thickness=0.01,
        web_thickness=0.008,
        cap_height=0.02,
        trailing_edge_insert_x_m=0.85,
        leading_edge_insert_x_m=0.05,
        core_enabled=True,
    )
    constructive, res = classify_property_labels(mcs.labels)
    assert "leading_edge_insert" in constructive and "trailing_edge_insert" in constructive
    assert "leading_edge_insert" not in res and "trailing_edge_insert" not in res
    assert any(k.startswith("sandwich_core_") for k in constructive)
    assert not any(k.startswith("sandwich_core_") for k in res)
    assert is_constructive_label("sandwich_core_0")
    assert is_constructive_label("leading_edge_insert")
    assert not is_residual_label("leading_edge_insert")


def test_default_env_backend_is_hybrid(monkeypatch):
    monkeypatch.delenv("SECTION_GEOMETRY_PROPS_BACKEND", raising=False)
    assert props_backend_from_env() == "hybrid"
    monkeypatch.setenv("SECTION_GEOMETRY_PROPS_BACKEND", "quadtree")
    assert props_backend_from_env() == "quadtree"
