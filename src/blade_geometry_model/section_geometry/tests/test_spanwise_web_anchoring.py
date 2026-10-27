"""Webs anchor spanwise and drop when they no longer fit; they never migrate.

Web chordwise positions used to be re-derived per station as
``frac x local_chord``, so every web slid inboard as the chord tapered and the
web centroid moved station to station. ``resolve_web_anchors`` compounded that
by walking any web that could not find inner-skin intercepts *further* LE-ward,
up to 12 steps, before dropping it.

The model these tests pin instead: the root station's web positions define the
centroidal axis for each web, held constant in the section plane along the span.
A web that no longer fits the local cavity is **dropped**, not moved.

Measured on the reference rotor (see
``examples/output/spanwise_web_anchors/README.md``): under the legacy model the
LE-ward walk never fired at all across 64 stations, so removing it is not what
changes behaviour here -- anchoring the positions is.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.api.worker import (
    resolve_web_x_positions,
    section_geometry_prep_key,
)
from blade_geometry_model.section_geometry.sections.web_layout import (
    resolve_web_anchors,
)

T_SKIN = 0.001
ROOT_CHORD = 0.2409
TIP_CHORD = 0.0878


def _airfoil(chord: float):
    return AirfoilSDF.from_naca("0012", n_points=200, chord=chord)


def _root_anchored(layout, root_chord: float = ROOT_CHORD) -> list[float]:
    """The root station's absolute web offsets, reused at every station."""
    return [float(f) * root_chord for f in layout.web_chord_fracs]


def test_anchored_positions_do_not_scale_with_local_chord():
    """The whole point: the same offsets, whatever the local chord."""
    layout = resolve_system_type("5D-F")
    anchored = _root_anchored(layout)

    root = resolve_web_x_positions(layout, _airfoil(ROOT_CHORD), anchored)
    tip = resolve_web_x_positions(layout, _airfoil(TIP_CHORD), anchored)

    assert root == pytest.approx(tip)
    assert root == pytest.approx(anchored)


def test_omitting_anchored_positions_keeps_the_legacy_per_station_form():
    """Callers that pass nothing must be unaffected."""
    layout = resolve_system_type("5D-F")

    tip = resolve_web_x_positions(layout, _airfoil(TIP_CHORD), None)

    assert tip == pytest.approx([f * TIP_CHORD for f in layout.web_chord_fracs])
    # ... and that is precisely the migration the anchored form removes.
    root = resolve_web_x_positions(layout, _airfoil(ROOT_CHORD), None)
    assert root != pytest.approx(tip)


def test_a_web_that_does_not_fit_is_dropped_not_moved():
    """Default behaviour: no LE-ward walk.

    Anchored positions from a 241 mm root chord include webs at 160.6 mm and
    200.7 mm, which on an 87.8 mm tip section are outside the buildable cavity.
    They must disappear, and the survivors must stay exactly where they were.
    """
    layout = resolve_system_type("5D-F")
    anchored = _root_anchored(layout)
    airfoil = _airfoil(TIP_CHORD)

    resolved, dropped = resolve_web_anchors(
        airfoil, T_SKIN, anchored, return_dropped=True
    )

    xs = [a[0] for a in resolved]
    assert dropped, "expected the outboard anchored webs to drop at the tip"
    assert len(xs) < len(anchored)
    # Every survivor sits exactly on an anchored position -- none was walked.
    for x in xs:
        assert min(abs(x - a) for a in anchored) == pytest.approx(0.0, abs=1e-12)
    for d in dropped:
        assert d.reason in {"no_inner_skin_intercept", "merged_into_previous_web"}


def test_legacy_migration_is_still_reachable_and_moves_webs():
    """The opt-in flag must actually restore the old behaviour.

    Guards against the flag becoming a no-op: under migration a web that cannot
    anchor is relocated rather than lost, so more webs survive and at least one
    sits away from its requested position.
    """
    layout = resolve_system_type("5D-F")
    anchored = _root_anchored(layout)
    airfoil = _airfoil(TIP_CHORD)

    dropped_only = resolve_web_anchors(airfoil, T_SKIN, anchored)
    migrated = resolve_web_anchors(
        airfoil, T_SKIN, anchored, allow_le_ward_migration=True
    )

    assert len(migrated) >= len(dropped_only)
    moved = [
        x for x, *_ in migrated
        if min(abs(x - a) for a in anchored) > 1e-9
    ]
    assert moved, "migration enabled but nothing moved; the flag is inert"


def test_anchored_webs_survive_where_the_section_is_deep():
    """No spurious drops: at the root every anchored web must be buildable."""
    layout = resolve_system_type("5D-F")
    anchored = _root_anchored(layout)

    resolved, dropped = resolve_web_anchors(
        _airfoil(ROOT_CHORD), T_SKIN, anchored, return_dropped=True
    )

    assert len(resolved) == len(anchored)
    assert dropped == []


def test_survivor_count_is_non_increasing_outboard():
    """The layout may only step down along the span, never back up."""
    layout = resolve_system_type("5D-F")
    anchored = _root_anchored(layout)

    counts = [
        len(resolve_web_anchors(_airfoil(c), T_SKIN, anchored))
        for c in np.linspace(ROOT_CHORD, TIP_CHORD, 12)
    ]

    assert all(b <= a for a, b in zip(counts, counts[1:])), counts
    assert counts[0] > counts[-1], "fixture no longer exercises any drop"


def test_prep_key_separates_sections_that_differ_only_in_web_positions():
    """Anchored positions change the section without changing layout or DVs.

    Leaving them out of the cache key would let two different sections collide
    on one entry, which is a silently wrong result rather than a slow one.
    """
    layout = resolve_system_type("5D-F")
    airfoil = _airfoil(TIP_CHORD)
    common = dict(skin_thickness=T_SKIN, spar_cap_thickness=0.004, shear_web_thickness=0.0015, twist_rad=0.0)

    legacy = section_geometry_prep_key(layout, airfoil, **common)
    anchored = section_geometry_prep_key(
        layout, airfoil, web_x_positions_m=_root_anchored(layout), **common
    )
    shifted = section_geometry_prep_key(
        layout, airfoil,
        web_x_positions_m=[x + 1e-3 for x in _root_anchored(layout)],
        **common,
    )

    assert legacy != anchored
    assert anchored != shifted
    # Same positions must still hit the same entry.
    assert anchored == section_geometry_prep_key(
        layout, airfoil, web_x_positions_m=_root_anchored(layout), **common
    )
