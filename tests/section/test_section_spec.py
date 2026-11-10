"""SectionSpec / SectionOptions are the geometry input; kwargs still work."""

from __future__ import annotations

import pytest

from implicit_geometry_kernel import Airfoil
from implicit_geometry_kernel.section import Gauges, SectionGeometry, SectionOptions, SectionSpec, resolve_system_type
from implicit_geometry_kernel.section.spec.worker import build_section_geometry


def _airfoil():
    return Airfoil.naca4("2412").to_sdf(chord=1.0)


def test_section_spec_validates_and_normalises():
    spec = SectionSpec(
        layout=resolve_system_type("2D-F-LC-0R"),
        airfoil=_airfoil(),
        twist_rad=0.1,
        gauges=Gauges(skin_ss=0.004, skin_ps=0.004, cap=0.012, web=0.005),
        webs=(0.15, 0.45),
    )
    assert spec.webs == (0.15, 0.45)
    assert spec.gauges.skin == 0.004


def test_section_spec_rejects_missing_layout():
    with pytest.raises(ValueError, match="layout"):
        SectionSpec(
            layout=None,
            airfoil=_airfoil(),
            twist_rad=0.0,
            gauges=Gauges(0.004, 0.004, 0.012, 0.005),
        )


def test_build_from_spec_matches_runner():
    layout = resolve_system_type("2D-F-LC-0R")
    airfoil = _airfoil()
    spec = SectionSpec.from_kwargs(
        layout,
        airfoil,
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
    )
    opts = SectionOptions(outputs=("properties",))
    runner = SectionGeometry(spec, opts).execute().get_results()
    via_spec = build_section_geometry(spec, options=opts)
    assert via_spec.properties_map.keys() == runner.properties_map.keys()
    for label, props in runner.properties_map.items():
        assert props["area"] == pytest.approx(via_spec.properties_map[label]["area"], rel=1e-9)


def test_circle_sdf_warns():
    from implicit_geometry_kernel.sdf.nodes.tagged import CircleSDF

    with pytest.warns(DeprecationWarning, match="CircleSDF"):
        CircleSDF(0.0, 0.0, 1.0)


def test_kwargs_constructor_and_from_kwargs_warn():
    layout = resolve_system_type("2D-F-LC-0R")
    airfoil = _airfoil()
    kwargs = dict(
        layout=layout,
        airfoil=airfoil,
        skin_thickness=0.004,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.005,
        outputs=("properties",),
    )
    with pytest.warns(DeprecationWarning, match="layout="):
        SectionGeometry(**kwargs)
    with pytest.warns(DeprecationWarning, match="from_kwargs"):
        SectionGeometry.from_kwargs(**kwargs)


def test_from_kwargs_rejects_conflicting_aliases():
    layout = resolve_system_type("2D-F-LC-0R")
    airfoil = _airfoil()
    with pytest.raises(ValueError, match="Conflicting"):
        SectionSpec.from_kwargs(
            layout,
            airfoil,
            skin_thickness=0.004,
            spar_cap_thickness=0.012,
            cap_height=0.02,
            shear_web_thickness=0.005,
        )
    with pytest.raises(ValueError, match="Conflicting"):
        SectionSpec.from_kwargs(
            layout,
            airfoil,
            skin_thickness=0.004,
            spar_cap_thickness=0.012,
            shear_web_thickness=0.005,
            twist_rad=0.1,
            twist_angle=0.2,
        )


def test_keyword_form_removed():
    with pytest.raises(TypeError, match="SectionSpec"):
        build_section_geometry(resolve_system_type("2D-F-LC-0R"))
