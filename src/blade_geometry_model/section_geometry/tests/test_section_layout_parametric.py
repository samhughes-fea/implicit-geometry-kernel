"""Parametric SectionLayout + CapPlacementRule catalogue conversion."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry.api.runner import (
    SectionGeometryRunner,
)
from blade_geometry_model.section_geometry.api.worker import (
    build_section_geometry,
    resolve_station_layout,
)
from blade_geometry_model.section_geometry.layout.cap_placement import (
    CapPlacementRule,
)
from blade_geometry_model.section_geometry.layout.section_layout import (
    SectionLayout,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)
from blade_geometry_model.section_geometry.sections.from_layout import (
    multicell_from_layout,
    multicell_from_section_layout,
)


def test_from_catalogue_key_2d_f_has_two_web_fracs():
    layout = SectionLayout.from_catalogue_key("2D-F")
    assert len(layout.web_chord_fracs) == 2
    assert layout.orientation == "flapwise"
    assert layout.catalogue_key == "2D-F"
    assert layout.cap_placement.kind == "continuous_box"
    spec = layout.to_system_layout_spec()
    assert spec.key == "2D-F"
    assert spec.n_webs == 2


def test_cap_placement_family_d_is_continuous_box():
    rule = CapPlacementRule.from_structural_family("D")
    assert rule.kind == "continuous_box"
    assert CapPlacementRule.from_structural_family("A").kind == "none"
    assert CapPlacementRule.from_structural_family("C").kind == "per_web_pads"
    assert CapPlacementRule.from_structural_family("B").kind == "fixed_band"


def test_assemble_from_section_layout_without_catalogue_key():
    airfoil = AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)
    layout = SectionLayout(
        web_chord_fracs=(0.25, 0.45),
        orientation="flapwise",
        cap_placement=CapPlacementRule.continuous_box(),
        catalogue_key=None,
        structural_family="D",
    )
    host = multicell_from_section_layout(
        airfoil,
        layout,
        skin_thickness=0.006,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    labels = [str(getattr(k, "value", k)) for k in host.labels]
    assert any(k.startswith("shear_web_") for k in labels)
    assert "spar_cap_suction_surface" in labels
    assert host.n_webs == 2


def test_multicell_from_layout_accepts_section_layout():
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = SectionLayout.from_catalogue_key("2D-F")
    a = multicell_from_layout(
        airfoil,
        layout,
        skin_thickness=0.006,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    b = multicell_from_section_layout(
        airfoil,
        layout,
        skin_thickness=0.006,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    assert a.n_webs == b.n_webs == 2
    # Areas of a skin label should match across the two assemble paths.
    key = "skin_suction_surface"
    xa = np.linspace(0.1, 0.6, 20)
    ya = np.zeros_like(xa)
    np.testing.assert_allclose(a.eval(key, xa, ya), b.eval(key, xa, ya), atol=1e-12)


def test_resolve_and_build_accept_keyless_section_layout():
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = SectionLayout(
        web_chord_fracs=(0.25, 0.45),
        orientation="flapwise",
        cap_placement=CapPlacementRule.continuous_box(),
        catalogue_key=None,
        structural_family="D",
    )
    res = resolve_station_layout(layout, airfoil, skin_thickness=0.006)
    assert isinstance(res.layout, SectionLayout)
    assert res.layout.n_webs == 2
    assert res.stepped_down is False

    result = build_section_geometry(
        res.layout,
        airfoil,
        skin_thickness=0.006,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
        outputs=("midlines",),
        web_x_positions_m=res.web_x_positions_m,
        web_y_coords=res.web_y_coords,
    )
    n_web_strips = sum(
        1 for s in (result.midline_strips or []) if str(s.kind) == "web"
    )
    assert n_web_strips == 2


def test_keyless_section_layout_steps_down_without_catalogue_key():
    """Drop cascade is N survivors + CapPlacementRule, not a SYSTEM_TYPE key."""
    # Same station geometry that steps catalogue 5D-F → 1C-F.
    root_chord = 0.2409
    airfoil = AirfoilSDF.from_naca("0012", n_points=200, chord=0.0878)
    layout = SectionLayout.from_catalogue_key("5D-F")
    layout = SectionLayout(
        web_chord_fracs=layout.web_chord_fracs,
        orientation=layout.orientation,
        cap_placement=layout.cap_placement,
        catalogue_key=None,
        structural_family=layout.structural_family,
    )
    anchored = [float(f) * root_chord for f in layout.web_chord_fracs]
    res = resolve_station_layout(
        layout,
        airfoil,
        skin_thickness=0.001,
        web_x_positions_m=anchored,
    )
    assert isinstance(res.layout, SectionLayout)
    assert res.layout.catalogue_key is None
    assert res.stepped_down is True
    assert res.layout.n_webs == len(res.web_x_positions_m) == 1
    assert res.dropped
    assert res.layout.cap_placement.kind == "per_web_pads"
    assert res.layout.structural_family == "C"


def test_runner_accepts_section_layout():
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    layout = SectionLayout.from_catalogue_key("2D-F")
    runner = SectionGeometryRunner(
        layout=layout,
        airfoil=airfoil,
        skin_thickness=0.006,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
        outputs=("midlines",),
    ).execute()
    res = runner.get_layout_resolution()
    result = runner.get_results()
    assert len(res.web_x_positions_m) >= 1
    assert result.midline_strips is not None


def test_keyless_section_layout_matches_catalogue_2d_f_recipes():
    """Keyless CapPlacementRule layout equals catalogue preset for 2D-F."""
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    catalog = SectionLayout.from_catalogue_key("2D-F")
    keyless = SectionLayout(
        web_chord_fracs=catalog.web_chord_fracs,
        orientation=catalog.orientation,
        cap_placement=catalog.cap_placement,
        catalogue_key=None,
        structural_family=catalog.structural_family,
    )
    kwargs = dict(
        skin_thickness=0.006,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    a = multicell_from_section_layout(airfoil, catalog, **kwargs)
    b = multicell_from_section_layout(airfoil, keyless, **kwargs)
    assert a.cap_placement.kind == b.cap_placement.kind == "continuous_box"
    assert set(a.wall_recipes) == set(b.wall_recipes)
    assert any(k.startswith("spar_cap") for k in a.wall_recipes)
    # Cap recipes authored before / alongside SDF components.
    assert a.spar_cap_wall_recipes
    xa = np.linspace(0.15, 0.55, 16)
    ya = np.zeros_like(xa)
    np.testing.assert_allclose(
        a.eval("skin_suction_surface", xa, ya),
        b.eval("skin_suction_surface", xa, ya),
        atol=1e-12,
    )
