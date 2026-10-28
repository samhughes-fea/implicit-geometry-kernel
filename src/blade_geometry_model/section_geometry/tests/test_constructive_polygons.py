"""Unit tests for constructive wall-solid polygons and shoelace properties."""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import Polygon

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.constructive import (
    classify_property_labels,
    emit_constructive_polygons,
    emit_wall_solid_polygons,
    is_constructive_label,
    is_core_label,
    is_wall_solid_label,
    is_residual_label,
    section_properties_from_polygons,
)
from blade_geometry_model.section_geometry.constructive.emitters import (
    emit_outer_skin_polygon,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.shoelace_integration import section_properties_S


def test_label_classifier_cores_and_inserts_constructive():
    constructive, res = classify_property_labels(
        [
            "skin_suction_surface",
            "skin_pressure_surface",
            "shear_web_0",
            "spar_cap_suction_surface",
            "sandwich_core_0",
            "leading_edge_insert",
            "trailing_edge_insert",
        ]
    )
    assert "skin_suction_surface" in constructive and "skin_pressure_surface" in constructive
    assert "shear_web_0" in constructive
    assert "sandwich_core_0" in constructive and is_core_label("sandwich_core_0")
    assert is_constructive_label("sandwich_core_1")
    assert "leading_edge_insert" in constructive and "trailing_edge_insert" in constructive
    assert not is_residual_label("leading_edge_insert")
    assert not is_residual_label("sandwich_core_0")
    assert is_wall_solid_label("spar_cap_pressure_surface")
    assert is_wall_solid_label("skin_suction_surface")
    assert is_wall_solid_label("leading_edge_insert")
    assert not is_wall_solid_label("sandwich_core_1")
    assert res == []


def test_section_properties_from_polygons_matches_shoelace_rectangle():
    rect = Polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
    props = section_properties_from_polygons({"rect": rect})["rect"]
    A, (cx, cy), Ixx, Iyy, Ixy = section_properties_S(rect)
    # Centroidal moments
    ixx_c = Ixx - A * cy**2
    iyy_c = Iyy - A * cx**2
    np.testing.assert_allclose(props["area"], A, rtol=1e-12)
    np.testing.assert_allclose(props["cx"], cx, rtol=1e-12)
    np.testing.assert_allclose(props["cy"], cy, rtol=1e-12)
    np.testing.assert_allclose(props["Ixx"], ixx_c, rtol=1e-12)
    np.testing.assert_allclose(props["Iyy"], iyy_c, rtol=1e-12)
    # Analytic: A=2, cx=1, cy=0.5, Ixx=A*h^2/12=2/12, Iyy=A*w^2/12=8/12
    np.testing.assert_allclose(props["area"], 2.0, rtol=1e-12)
    np.testing.assert_allclose(props["Ixx"], 2.0 / 12.0, rtol=1e-10)
    np.testing.assert_allclose(props["Iyy"], 8.0 / 12.0, rtol=1e-10)


def test_outer_skin_annulus_positive_area():
    airfoil = AirfoilSDF.from_naca("0012", n_points=80, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    a_ss = emit_outer_skin_polygon(section["skin_suction_surface"]).area
    a_ps = emit_outer_skin_polygon(section["skin_pressure_surface"]).area
    assert a_ss > 0.0 and a_ps > 0.0


def test_emit_wall_solid_polygons_2d_cn():
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    polys = emit_wall_solid_polygons(section)
    assert "skin_suction_surface" in polys and "skin_pressure_surface" in polys
    assert any(k.startswith("shear_web_") for k in polys)
    assert "spar_cap_suction_surface" in polys and "spar_cap_pressure_surface" in polys
    for p in polys.values():
        assert p.area > 0.0


def test_twist_preserves_wall_solid_area():
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    layout = resolve_system_type("2D-CN")
    s0 = build_section_view(airfoil, layout, twist_angle_rad=0.0)
    s15 = build_section_view(airfoil, layout, twist_angle_rad=float(np.deg2rad(15.0)))
    p0 = emit_wall_solid_polygons(s0)
    p15 = emit_wall_solid_polygons(s15)
    for key in p0:
        np.testing.assert_allclose(p15[key].area, p0[key].area, rtol=1e-3, atol=1e-5)


@pytest.mark.parametrize("layout_key", ["2D-CN", "2D-F"])
def test_hybrid_constructive_area_vs_quadtree(layout_key: str):
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type(layout_key), twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=150000)
    hy, meta = compute_section_properties(
        section, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    qt, _ = compute_section_properties(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    # No inserts → residual φ skipped
    assert meta["residual_labels"] == []
    assert meta.get("residual_quadtree_s", 0.0) < 1e-3
    for label in section.labels:
        key = str(getattr(label, "value", label))
        if not is_constructive_label(key):
            continue
        if qt[key]["area"] < 1e-5:
            continue
        np.testing.assert_allclose(
            hy[key]["area"], qt[key]["area"], rtol=1e-3, atol=5e-4
        )


def test_emit_core_polygons_bay_partition():
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    polys = emit_constructive_polygons(section)
    core_keys = sorted(
        k for k in polys if k.startswith("sandwich_core_")
    )
    assert len(core_keys) >= 2
    core_area = sum(float(polys[k].area) for k in core_keys)
    assert core_area > 0.0


def test_section_properties_from_polygons_handles_degenerate():
    """Zero-area / bow-tie leftovers must not raise (e2e FoS-thick crash guard)."""
    empty = Polygon()
    # Self-touching bow-tie: Shapely may report tiny area; shoelace A can be 0.
    bow = Polygon([(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)])
    props = section_properties_from_polygons({"empty": empty, "bow": bow})
    assert props["empty"]["area"] == 0.0
    assert props["bow"]["area"] == 0.0 or np.isfinite(props["bow"]["area"])


def test_hybrid_5d_f_thick_skin_does_not_crash():
    """Heavy wall thicknesses (FoS-thickened geometry) must stay on hybrid path."""
    from blade_geometry_model.section_geometry.sections import (
        MultiCellSection,
    )

    # Chord must leave cavity depth for 12 mm skin + five webs; tip-scale
    # chords raise InfeasibleSectionGeometry before hybrid props can run.
    airfoil = AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)
    layout = resolve_system_type("5D-F")
    xs = [float(f) * float(airfoil.chord) for f in layout.web_chord_fracs]
    section = MultiCellSection(
        airfoil,
        web_x_positions=xs,
        web_thickness=0.008,
        web_alignment="flapwise",
        cap_height=0.015,
        skin_thickness=0.012,
        twist_angle=0.0,
        core_enabled=True,
        structural_family=layout.structural_family,
    )
    polys = emit_constructive_polygons(section)
    props = section_properties_from_polygons(polys)
    assert "skin_suction_surface" in props and "skin_pressure_surface" in props
    assert props["skin_suction_surface"]["area"] > 0.0
    assert props["skin_pressure_surface"]["area"] > 0.0
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=80000)
    hy, meta = compute_section_properties(
        section, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    assert meta["backend"] == "hybrid"
    assert hy["skin_suction_surface"]["area"] > 0.0
    assert hy["skin_pressure_surface"]["area"] > 0.0
    total = sum(float(p["area"]) for p in hy.values() if np.isfinite(p["area"]))
    assert total > 0.0


def test_hybrid_retains_core_and_optional_inserts():
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    # MultiCellSection via build_section_view; enable inserts through factory if available.
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    # Direct MultiCellSection with inserts for constructive coverage
    from blade_geometry_model.section_geometry.sections.multicell import (
        MultiCellSection,
    )

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
    assert any(k.startswith("sandwich_core_") for k in constructive)
    assert not any(k.startswith("sandwich_core_") for k in res)
    assert "leading_edge_insert" in constructive and "trailing_edge_insert" in constructive
    assert "leading_edge_insert" not in res and "trailing_edge_insert" not in res
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=80000)
    props, meta = compute_section_properties(
        mcs, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    assert meta["backend"] == "hybrid"
    assert "leading_edge_insert" not in meta["residual_labels"]
    assert "trailing_edge_insert" not in meta["residual_labels"]
    assert "leading_edge_insert" in meta["constructive_labels"]
    assert "trailing_edge_insert" in meta["constructive_labels"]
    assert any(
        k.startswith("sandwich_core_")
        for k in meta["constructive_labels"]
    )
    for key in ("leading_edge_insert", "trailing_edge_insert"):
        assert key in props
        assert np.isfinite(props[key]["area"]) or props[key]["area"] == 0.0
    assert any(k.startswith("sandwich_core_") for k in props)
