"""Assemble fail-closed: requested webs and cap rule must match.

``resolve_station_layout`` may step ``2D-F`` → ``1C-F`` when webs drop.
``MultiCellSection`` must *not* silently shrink the web set while keeping a
``continuous_box`` rule (illegal 1D → ``None`` caps → twist ``AttributeError``).
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    build_section_with_thicknesses,
    resolve_station_layout,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.layout.cap_placement import (
    CapPlacementRule,
)
from blade_geometry_model.section_geometry.layout.system_layout import (
    layout_without_spar_caps,
)
from blade_geometry_model.section_geometry.sections.build_spar_caps import (
    build_spar_caps,
)


ROOT_CHORD = 0.2409
T_SKIN_FLOOR = 0.0025


def _naca_chord(chord: float) -> AirfoilSDF:
    return AirfoilSDF.from_naca("0012", n_points=200, chord=float(chord))


def _anchored_2d_f(root_chord: float = ROOT_CHORD) -> list[float]:
    layout = resolve_system_type("2D-F")
    return [float(f) * root_chord for f in layout.web_chord_fracs]


def test_partial_web_drop_at_assemble_raises_not_attribute_error():
    """2D-F + 2 xs + thick skin that loses a web → InfeasibleSectionGeometry.

    Must never reach twist_component(None) / ``NoneType.expr``.
    """
    chord = 0.1353
    af = _naca_chord(chord)
    xs = _anchored_2d_f()
    assert len(xs) == 2

    with pytest.raises(InfeasibleSectionGeometry, match="cannot all anchor"):
        MultiCellSection(
            af,
            web_x_positions=xs,
            structural_family="D",
            skin_thickness=0.004,
            cap_height=0.005,
            web_thickness=0.004,
            twist_angle=float(np.radians(4.0)),
        )


def test_resolve_then_assemble_1c_succeeds_with_non_none_components():
    """After resolve steps to 1C-F, assemble with frozen anchors is legal."""
    chord = 0.1353
    af = _naca_chord(chord)
    xs = _anchored_2d_f()
    res = resolve_station_layout(
        resolve_system_type("2D-F"),
        af,
        skin_thickness=0.004,
        web_x_positions_m=xs,
    )
    assert res.layout.key == "1C-F"
    assert len(res.web_x_positions_m) == 1
    assert res.layout.cap_placement.kind == "per_web_pads"

    section = build_section_with_thicknesses(
        af,
        res.layout,
        skin_thickness=T_SKIN_FLOOR,
        spar_cap_thickness=0.005,
        shear_web_thickness=0.004,
        twist_rad=float(np.radians(4.0)),
        web_x_positions_m=res.web_x_positions_m,
        web_y_coords=res.web_y_coords,
    )
    comps = section.components
    assert comps
    assert all(c is not None for c in comps.values())
    assert sum(
        1 for k in comps if str(getattr(k, "value", k)).startswith("shear_web_")
    ) == 1


def test_frozen_web_y_coords_do_not_silently_drop_a_web():
    """Thicker skin with frozen y-coords keeps N webs (no mid-assemble shrink)."""
    af = _naca_chord(ROOT_CHORD)
    xs = _anchored_2d_f()
    res = resolve_station_layout(
        resolve_system_type("2D-F"),
        af,
        skin_thickness=T_SKIN_FLOOR,
        web_x_positions_m=xs,
    )
    assert res.layout.key == "2D-F"
    assert len(res.web_x_positions_m) == 2

    section = build_section_with_thicknesses(
        af,
        res.layout,
        skin_thickness=0.008,
        spar_cap_thickness=0.006,
        shear_web_thickness=0.004,
        twist_rad=float(np.radians(3.0)),
        web_x_positions_m=res.web_x_positions_m,
        web_y_coords=res.web_y_coords,
    )
    n_webs = sum(
        1
        for k in section.components
        if str(getattr(k, "value", k)).startswith("shear_web_")
    )
    assert n_webs == 2
    assert all(c is not None for c in section.components.values())


def test_build_spar_caps_rejects_continuous_box_with_one_web():
    af = _naca_chord(1.0)
    with pytest.raises(InfeasibleSectionGeometry, match="continuous_box"):
        build_spar_caps(
            af,
            0.003,
            [0.25],
            [(0.05, -0.05)],
            0.01,
            0.01,
            [],
            cap_placement=CapPlacementRule.continuous_box(),
            airfoil_expr_node=None,
            all_flapwise=False,
            c_station=None,
            nx=None,
            ny=None,
            N=1,
            pitch_fr=1.0 / 3.0,
            fhw=0.05,
            dhw=0.04,
            y_search=None,
            phi_inner_ref=None,
        )


def test_cap_runout_via_layout_without_spar_caps_omits_cap_keys():
    af = _naca_chord(ROOT_CHORD)
    layout = layout_without_spar_caps(resolve_system_type("1C-F"))
    assert layout.key == "1A-F"
    assert not layout.has_spar_caps
    section = build_section_with_thicknesses(
        af,
        layout,
        skin_thickness=0.003,
        spar_cap_thickness=0.005,
        shear_web_thickness=0.004,
        twist_rad=0.0,
        web_x_positions_m=[0.25 * ROOT_CHORD],
    )
    labels = {str(getattr(k, "value", k)) for k in section.components}
    assert "spar_cap_suction_surface" not in labels
    assert "spar_cap_pressure_surface" not in labels
    assert all(c is not None for c in section.components.values())

