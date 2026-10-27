"""Independent SS/PS skin gauges and per-web thicknesses."""

from __future__ import annotations

import pytest

from blade_geometry_model.section_geometry.api.kwarg_aliases import (
    resolve_shear_web_thickness,
    resolve_skin_gauges,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)
from blade_geometry_model.section_geometry.sections.multicell import (
    MultiCellSection,
)

def test_resolve_skin_gauges_scalar_broadcasts():
    assert resolve_skin_gauges(skin_thickness=0.005) == (0.005, 0.005)

def test_resolve_skin_gauges_pair():
    assert resolve_skin_gauges(skin_thickness=(0.004, 0.006)) == (0.004, 0.006)

def test_resolve_shear_web_thickness_list():
    assert resolve_shear_web_thickness(shear_web_thickness=[0.005, 0.007]) == [
        0.005,
        0.007,
    ]

def test_multicell_asymmetric_skin_recipes():
    af = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    mcs = MultiCellSection(
        af,
        web_x_positions=[0.25, 0.5],
        skin_thickness=(0.004, 0.008),
        web_thickness=[0.005, 0.007],
        cap_height=0.01,
        structural_family="D",
        core_enabled=False,
    )
    ss = mcs.wall_recipes["skin_suction_surface"]
    ps = mcs.wall_recipes["skin_pressure_surface"]
    assert float(ss.thickness) == 0.004
    assert float(ps.thickness) == 0.008
    w0 = mcs.wall_recipes["shear_web_0"]
    w1 = mcs.wall_recipes["shear_web_1"]
    assert float(w0.thickness) == 0.005
    assert float(w1.thickness) == 0.007


def test_zone_gauges_m_preserve_independent_webs(airfoil=None):
    from blade_geometry_model.section_geometry.api.runner import (
        SectionGeometryRunner,
    )
    from blade_geometry_model.section_geometry.layout.system_layout import (
        resolve_system_type,
    )
    from blade_geometry_model.section_geometry.sdf.airfoil import (
        AirfoilSDF,
    )

    af = airfoil or AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    result = (
        SectionGeometryRunner(
            layout=resolve_system_type("2D-F"),
            airfoil=af,
            zone_gauges_m={
                "skin:suction": 0.004,
                "skin:pressure": 0.004,
                "cap:suction": 0.010,
                "cap:pressure": 0.010,
                "web:0": 0.005,
                "web:1": 0.008,
            },
            outputs=("midlines",),
        )
        .execute()
        .get_results()
    )
    webs = {
        str(getattr(s, "label", "")): float(s.thickness_m)
        for s in result.midline_strips
        if str(getattr(s, "kind", "")) == "web"
    }
    assert len(webs) >= 2
    thicknesses = sorted(webs.values())
    assert thicknesses[0] == pytest.approx(0.005, abs=1e-9)
    assert thicknesses[-1] == pytest.approx(0.008, abs=1e-9)


def test_prep_key_distinguishes_unequal_zone_maps():
    from blade_geometry_model.section_geometry.api.worker import (
        section_geometry_prep_key,
    )
    from blade_geometry_model.section_geometry.layout.system_layout import (
        resolve_system_type,
    )
    from blade_geometry_model.section_geometry.sdf.airfoil import (
        AirfoilSDF,
    )

    af = AirfoilSDF.from_naca("2412", n_points=40, chord=1.0)
    layout = resolve_system_type("2D-F")
    common = {
        "skin:suction": 0.004,
        "skin:pressure": 0.004,
        "cap:suction": 0.010,
        "cap:pressure": 0.010,
    }
    a = section_geometry_prep_key(
        layout,
        af,
        zone_gauges_m={**common, "web:0": 0.005, "web:1": 0.005},
        outputs=("properties",),
    )
    b = section_geometry_prep_key(
        layout,
        af,
        zone_gauges_m={**common, "web:0": 0.005, "web:1": 0.008},
        outputs=("properties",),
    )
    assert a != b
