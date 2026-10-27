"""Stepping a layout down as webs drop, generalised across families A/B/C/D.

Web positions are anchored spanwise, so a tapering section eventually stops
being able to carry every web. ``effective_layout_for_web_count`` maps the
surviving count back onto a catalogue layout.

The rule is family-general: keep the parent family and orientation where the
catalogue defines that web count, otherwise fall to the next less-elaborate
family that does, and never cross the cap-bearing boundary. Family ``D`` is the
only one that has to cross, because ``1D`` does not exist -- a continuous box
needs a fore and an aft web.

The cascade **floors at one web**. ``0A`` / ``0B`` are selectable layouts but
not degradation targets: a section with no web has no shear path at all, so
losing the last one is a different kind of event from losing one of several.
A station that cannot anchor a single web raises instead.

Expected cascades::

    A: 5A -> 4A -> 3A -> 2A -> 1A
    B: 5B -> 4B -> 3B -> 2B -> 1B
    C: 5C -> 4C -> 3C -> 2C -> 1C
    D: 5D -> 4D -> 3D -> 2D -> 1C
"""

from __future__ import annotations

import pytest

from blade_geometry_model.section_geometry import (
    SYSTEM_TYPE_KEYS,
    resolve_system_type,
)
from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)
from blade_geometry_model.section_geometry.layout.system_layout import (
    effective_layout_for_web_count,
    family_min_webs,
    step_down_cascade,
)


def test_family_floors_match_the_catalogue():
    """Catalogue floors, which are not the cascade floor.

    ``family_min_webs`` reports the fewest webs the *registry* defines per
    family; the cascade separately floors at one web for every family. The two
    coincide only for C.
    """
    assert family_min_webs("D") == 2
    assert family_min_webs("C") == 1
    assert family_min_webs("B") == 0
    assert family_min_webs("A") == 0


@pytest.mark.parametrize(
    "start, expected",
    [
        ("5A-F", ("5A-F", "4A-F", "3A-F", "2A-F", "1A-F")),
        ("5B-F", ("5B-F", "4B-F", "3B-F", "2B-F", "1B-F")),
        ("5C-F", ("5C-F", "4C-F", "3C-F", "2C-F", "1C-F")),
        ("5D-F", ("5D-F", "4D-F", "3D-F", "2D-F", "1C-F")),
        ("5D-CN", ("5D-CN", "4D-CN", "3D-CN", "2D-CN", "1C-CN")),
        ("3D-F", ("3D-F", "2D-F", "1C-F")),
        ("2D-CN", ("2D-CN", "1C-CN")),
        ("1C-F", ("1C-F",)),
        ("0B", ("0B",)),
        ("0A", ("0A",)),
    ],
)
def test_cascades_are_as_specified(start, expected):
    assert step_down_cascade(resolve_system_type(start)) == expected


def test_family_d_crosses_to_c_at_one_web():
    """1D does not exist, so D must hand over rather than fail.

    The catalogue omits it deliberately: a continuous box between the first and
    last web degenerates to zero-span caps with only one web.
    """
    assert "1D-F" not in SYSTEM_TYPE_KEYS
    assert "1D-CN" not in SYSTEM_TYPE_KEYS

    stepped = effective_layout_for_web_count(resolve_system_type("5D-F"), 1)

    assert stepped.key == "1C-F"
    assert stepped.n_webs == 1
    assert stepped.has_spar_caps is True


@pytest.mark.parametrize("key", [k for k in SYSTEM_TYPE_KEYS])
def test_every_catalogue_layout_has_a_valid_cascade(key):
    """Exhaustive: no layout in the catalogue can fail to step down."""
    layout = resolve_system_type(key)
    cascade = step_down_cascade(layout)

    assert cascade[0] == key
    # Floors at one web, so an n-web layout yields n rungs; a 0-web layout is
    # its own terminal rung.
    assert len(cascade) == max(layout.n_webs, 1)
    for step in cascade:
        assert step in SYSTEM_TYPE_KEYS


@pytest.mark.parametrize("key", [k for k in SYSTEM_TYPE_KEYS])
def test_cascade_web_count_is_strictly_decreasing(key):
    layout = resolve_system_type(key)
    counts = [resolve_system_type(k).n_webs for k in step_down_cascade(layout)]

    assert counts == sorted(counts, reverse=True)
    assert len(set(counts)) == len(counts), "a rung repeated a web count"
    assert counts[-1] == min(layout.n_webs, 1), "cascade must floor at one web"


@pytest.mark.parametrize("key", [k for k in SYSTEM_TYPE_KEYS])
def test_cap_bearing_is_never_silently_lost(key):
    """Stepping down must not delete the spar caps.

    Family A has no caps and must never acquire them; B/C/D have caps and must
    keep them all the way to the terminal rung. With the cascade floored at one
    web this is what stops D handing over to 1A instead of 1C.
    """
    layout = resolve_system_type(key)
    want = bool(layout.has_spar_caps)

    for step in step_down_cascade(layout):
        assert bool(resolve_system_type(step).has_spar_caps) is want


@pytest.mark.parametrize("key", [k for k in SYSTEM_TYPE_KEYS if "-" in k])
def test_orientation_is_preserved_until_there_are_no_webs(key):
    """A flapwise blade must not become chord-normal partway down the span."""
    layout = resolve_system_type(key)

    for step in step_down_cascade(layout):
        spec = resolve_system_type(step)
        if spec.n_webs == 0:
            continue  # a selected 0-web layout has no web to orient
        assert spec.web_orientation == layout.web_orientation


def test_asking_for_the_same_count_returns_the_same_layout():
    layout = resolve_system_type("3C-CN")
    assert effective_layout_for_web_count(layout, 3) is layout


def test_stepping_up_is_rejected():
    """The cascade is one-way; a section cannot regain a web outboard."""
    layout = resolve_system_type("2D-F")

    with pytest.raises(ValueError, match="step-down only"):
        effective_layout_for_web_count(layout, 3)


def test_negative_web_count_is_rejected():
    with pytest.raises(ValueError):
        effective_layout_for_web_count(resolve_system_type("2D-F"), -1)


@pytest.mark.parametrize("key", ["5D-F", "5C-CN", "5B-F", "5A-F", "1C-F", "1A-CN"])
def test_losing_the_last_web_is_infeasible_not_a_rung(key):
    """0A/0B are selectable layouts, not degradation targets.

    A section with no web has no shear path at all, so losing the last one is a
    different kind of event from losing one of several. It must surface as an
    over-constrained station rather than quietly becoming a web-less section.
    """
    with pytest.raises(InfeasibleSectionGeometry):
        effective_layout_for_web_count(resolve_system_type(key), 0)


@pytest.mark.parametrize("key", ["0A", "0B"])
def test_an_explicitly_selected_zero_web_layout_is_its_own_terminal_rung(key):
    """Selecting 0A/0B is a design choice and must still work."""
    layout = resolve_system_type(key)

    assert step_down_cascade(layout) == (key,)
    assert effective_layout_for_web_count(layout, 0) is layout
