"""Non-parallel (slanted) web axes via SectionLayout.web_axes_m."""

from __future__ import annotations

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
    multicell_from_section_layout,
)

def test_slanted_web_axis_on_section_layout():
    airfoil = AirfoilSDF.from_naca("0012", n_points=100, chord=1.0)
    layout = SectionLayout(
        web_chord_fracs=(0.35,),
        orientation="chord_normal",
        cap_placement=CapPlacementRule.per_web_pads(),
        catalogue_key=None,
        structural_family="C",
        web_axes_m=(((0.30, 0.04), (0.40, -0.04)),),
    )
    host = multicell_from_section_layout(
        airfoil,
        layout,
        skin_thickness=0.005,
        spar_cap_thickness=0.01,
        shear_web_thickness=0.004,
    )
    wall = host.wall_recipes["shear_web_0"]
    assert abs(float(wall.x_top) - 0.30) < 1e-12
    assert abs(float(wall.x_bot) - 0.40) < 1e-12
    assert host.n_webs == 1
