"""What a station actually builds, once webs that no longer fit are dropped.

``resolve_station_layout`` is the single entry point composing the two halves
that already existed separately: ``resolve_web_anchors`` decides which anchored
webs the local cavity can hold, and ``effective_layout_for_web_count`` maps that
count onto a catalogue rung.

It is deliberately cheap — inner-skin intercept probes, no geometry assembly —
because the effective layout has to be known *before* the build, so it can drive
both the section and the design variables the sizing loop allocates.
"""

from __future__ import annotations

import pytest

from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.api.worker import (
    resolve_station_layout,
)

ROOT_CHORD = 0.2409
T_SKIN = 0.001


def _airfoil(chord: float):
    return AirfoilSDF.from_naca("0012", n_points=200, chord=chord)


def _anchored(layout, root_chord: float = ROOT_CHORD) -> list[float]:
    return [float(f) * root_chord for f in layout.web_chord_fracs]


def test_root_station_keeps_every_web():
    layout = resolve_system_type("5D-F")

    res = resolve_station_layout(
        layout, _airfoil(ROOT_CHORD), skin_thickness=T_SKIN,
        web_x_positions_m=_anchored(layout),
    )

    assert res.layout.key == "5D-F"
    assert len(res.web_x_positions_m) == 5
    assert len(res.web_y_coords) == 5
    for yt, yb in res.web_y_coords:
        assert yt > yb
    assert res.dropped == ()
    assert res.stepped_down is False


@pytest.mark.parametrize(
    "chord, expected_key, expected_webs",
    [
        (0.2409, "5D-F", 5),
        (0.2026, "4D-F", 4),
        (0.1353, "3D-F", 3),
        (0.1121, "2D-F", 2),
        (0.0878, "1C-F", 1),
    ],
)
def test_layout_steps_down_as_the_section_tapers(chord, expected_key, expected_webs):
    """The measured ladder from A1, now resolved in one call."""
    layout = resolve_system_type("5D-F")

    res = resolve_station_layout(
        layout, _airfoil(chord), skin_thickness=T_SKIN,
        web_x_positions_m=_anchored(layout),
    )

    assert res.layout.key == expected_key
    assert len(res.web_x_positions_m) == expected_webs
    assert len(res.dropped) == 5 - expected_webs


def test_survivors_stay_on_their_anchored_offsets():
    """No migration: a surviving web is exactly where the root put it."""
    layout = resolve_system_type("5D-F")
    anchored = _anchored(layout)

    res = resolve_station_layout(
        layout, _airfoil(0.1121), skin_thickness=T_SKIN, web_x_positions_m=anchored
    )

    for x in res.web_x_positions_m:
        assert min(abs(x - a) for a in anchored) == pytest.approx(0.0, abs=1e-12)


def test_dropped_webs_carry_a_reason():
    layout = resolve_system_type("5D-F")

    res = resolve_station_layout(
        layout, _airfoil(0.0878), skin_thickness=T_SKIN,
        web_x_positions_m=_anchored(layout),
    )

    assert res.dropped
    for entry in res.dropped:
        assert isinstance(entry.x_m, float)
        assert entry.reason in {
            "no_inner_skin_intercept", "merged_into_previous_web"
        }


def test_omitting_anchored_positions_keeps_the_legacy_behaviour():
    """Callers that pass nothing must see no drops and no step-down.

    This is what keeps every existing call site unaffected until it opts in.
    """
    layout = resolve_system_type("5D-F")

    res = resolve_station_layout(layout, _airfoil(0.0878), skin_thickness=T_SKIN)

    assert res.layout.key == "5D-F"
    assert len(res.web_x_positions_m) == 5
    assert res.dropped == ()
    assert res.stepped_down is False


def test_zero_web_layout_resolves_without_probing():
    """0A / 0B are selected, never stepped into, and have no webs to check."""
    for key in ("0A", "0B"):
        layout = resolve_system_type(key)
        res = resolve_station_layout(layout, _airfoil(0.0878), skin_thickness=T_SKIN)

        assert res.layout.key == key
        assert res.web_x_positions_m == ()
        assert res.stepped_down is False


def test_a_section_that_cannot_hold_any_web_is_infeasible():
    """Losing the last web is over-constrained, not a rung.

    Anchoring a 5-web layout on a chord far larger than the station's puts every
    web outside the buildable cavity.
    """
    layout = resolve_system_type("5D-F")
    far_anchored = _anchored(layout, root_chord=2.0)

    with pytest.raises(InfeasibleSectionGeometry):
        resolve_station_layout(
            layout, _airfoil(0.0878), skin_thickness=T_SKIN,
            web_x_positions_m=far_anchored,
        )


def test_nominal_layout_is_retained_for_comparison():
    """Callers need to report what was asked for as well as what was built."""
    layout = resolve_system_type("5D-F")

    res = resolve_station_layout(
        layout, _airfoil(0.0878), skin_thickness=T_SKIN,
        web_x_positions_m=_anchored(layout),
    )

    assert res.nominal_layout.key == "5D-F"
    assert res.layout.key == "1C-F"
    assert res.stepped_down is True


@pytest.mark.parametrize("start", ["5C-CN", "5B-F", "5A-F"])
def test_other_families_step_down_too(start):
    """The entry point is family-general, not 5D-F-specific."""
    layout = resolve_system_type(start)

    res = resolve_station_layout(
        layout, _airfoil(0.1353), skin_thickness=T_SKIN,
        web_x_positions_m=_anchored(layout),
    )

    assert res.layout.n_webs < layout.n_webs
    assert res.layout.structural_family == layout.structural_family
    assert bool(res.layout.has_spar_caps) is bool(layout.has_spar_caps)
