"""``0B`` spar caps must respond to ``spar_cap_thickness`` and fit the section they sit in.

``0B`` is the terminal rung of the spanwise step-down: once a section can no
longer carry its family's minimum web count it drops here — skin plus a
continuous, web-independent cap band on :class:`MultiCellSection` (family B,
zero webs, max-thickness anchor). That only works if the cap is actually sized.

Historically ``OuterSkinOnlySection`` built caps from a hardcoded height and
ignored ``spar_cap_thickness`` (docs/bug_log.md BL-008). The sized worker now always builds
``MultiCellSection`` and passes ``spar_cap_thickness`` as ``cap_height``.

Above the height at which the caps stop fitting the local cavity, the clamp
restricts the band and scales the height, so cap area plateaus. These tests pin
responsiveness below the plateau and containment above it.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    build_section_geometry,
    build_section_with_thicknesses,
    resolve_system_type,
)

ROOT_CHORD = 0.2409
TIP_CHORD = 0.0878
T_SKIN = 0.001


def _cap_area(chord: float, spar_cap_thickness: float) -> float:
    layout = resolve_system_type("0B")
    airfoil = AirfoilSDF.from_naca("0012", n_points=200, chord=chord)
    result = build_section_geometry(
        layout,
        airfoil,
        skin_thickness=T_SKIN,
        spar_cap_thickness=spar_cap_thickness,
        shear_web_thickness=0.0015,
        outputs=("properties",),
    )
    props = result.properties_map or {}
    return sum(
        float(v.get("area", 0.0)) for k, v in props.items() if "spar_cap" in k
    )


def test_cap_area_responds_to_t_spar():
    """spar_cap_thickness must size the 0B band (above the MultiCell spar wall thickness floor)."""
    # MultiCell clamps spar wall thickness to >= 5 mm; sample above that floor.
    areas = [_cap_area(ROOT_CHORD, t) for t in (0.0055, 0.007, 0.0085)]

    assert all(b > a for a, b in zip(areas, areas[1:])), areas
    assert areas[-1] == pytest.approx(areas[0] * (0.0085 / 0.0055), rel=0.2)


def test_cap_is_contained_by_a_shallow_section():
    """A requested cap larger than the cavity must be clamped, not built."""
    deep = _cap_area(ROOT_CHORD, 0.050)
    shallow = _cap_area(TIP_CHORD, 0.050)

    assert shallow > 0.0, "clamping must not delete the cap outright"
    assert shallow < deep
    airfoil = AirfoilSDF.from_naca("0012", n_points=200, chord=TIP_CHORD)
    _xc, thickness = airfoil.thickness_distribution(n_points=200)
    max_depth = float(np.max(thickness))
    assert shallow < max_depth * TIP_CHORD


def test_zero_web_layout_still_has_caps():
    """0B keeps its spar caps — it is 0 *webs*, not 0 structure."""
    layout = resolve_system_type("0B")
    assert layout.n_webs == 0
    assert layout.has_spar_caps is True

    assert _cap_area(ROOT_CHORD, 0.004) > 0.0


def test_0B_sized_build_is_multicell():
    layout = resolve_system_type("0B")
    airfoil = AirfoilSDF.from_naca("0012", n_points=200, chord=ROOT_CHORD)
    section = build_section_with_thicknesses(
        airfoil,
        layout,
        twist_rad=0.0,
        skin_thickness=T_SKIN,
        spar_cap_thickness=0.004,
        shear_web_thickness=0.0,
    )
    assert isinstance(section, MultiCellSection)
    assert section.n_webs == 0
    assert "spar_cap_suction_surface" in section.labels
