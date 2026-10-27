"""A web must sit somewhere its family can actually build a cap.

Anchored web positions are absolute offsets from the root, so as the chord
shrinks outboard a surviving web migrates *aft* in ``x/c`` into ever-shallower
cavity. The only drop criteria were topological — no inner-skin intercept, or
merged into its neighbour — and having an intercept is not the same as being
worth anchoring. On the reference rotor a web sits at 0.85c in a 1.0 mm cavity
and still probes two crossings, so it survives; family D then keys its cap
height on it (``h_allow = 0.5 x depth``) and reports a 0.5 mm allowance against
a 2 mm minimum gauge.

The visible symptom was a sawtooth in the cap allowance — it decays as the web
migrates aft, collapses below gauge, then jumps when that web finally drops and
a forward web inherits the aft end (1.080 mm at station 38, 5.742 mm at 42).
The cause is the section carrying a web it cannot build a cap against.

These tests pin the criterion and the invariant it buys: every surviving web
leaves an allowance at or above minimum gauge.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.sections.cap_feasibility import (
    inner_y_at_xs,
)
from blade_geometry_model.section_geometry.sections.web_layout import (
    resolve_web_anchors,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)

T_SKIN = 0.001
GAUGE = 0.002
# Family D: h_allow = 0.5 x depth, so a minimum-gauge cap needs 2 x gauge of it.
THRESHOLD = 2.0 * GAUGE

ROOT_CHORD = 0.2409
# 5D-F uniform fracs, frozen at the root into absolute offsets.
ANCHORS = [ROOT_CHORD * (i / 6.0) for i in range(1, 6)]

# Reference-rotor sizing stations and chords [mm].
STATIONS = (0, 4, 8, 13, 17, 21, 25, 29, 34, 38, 42, 46, 50, 55, 59, 63)
CHORD_MM = (240.9, 287.6, 292.0, 268.8, 241.9, 214.6, 191.7, 172.8,
            153.6, 141.0, 130.1, 120.5, 112.1, 102.6, 94.9, 87.8)
# The three whose aft web sits in cavity too shallow for a gauge cap.
SUB_GAUGE_STATIONS = {38, 55, 59}


def _airfoil(chord_mm: float) -> AirfoilSDF:
    return AirfoilSDF.from_naca_series(
        6, 63.0, 4.0, 15.0, n_points=200, chord=chord_mm / 1e3, closed_te=True
    )


def _cap_allowance(af: AirfoilSDF, resolved) -> float:
    """``0.5 x`` cavity depth at the aft web — family D's cap height rule."""
    if not resolved:
        return float("nan")
    probe = inner_y_at_xs(
        af, T_SKIN, [resolved[-1][0]], y_search=None, phi_inner=None,
        bisection_iters=0, raise_on_miss=False,
    )[0]
    if probe is None:
        return float("nan")
    return 0.5 * (probe[0] - probe[1])


def _resolve(chord_mm: float, threshold: float):
    af = _airfoil(chord_mm)
    return af, resolve_web_anchors(af, T_SKIN, ANCHORS, min_cavity_depth_m=threshold)


class TestCriterion:
    def test_zero_threshold_changes_nothing(self):
        """The opt-in contract: callers that do not ask keep exact behaviour."""
        for chord in CHORD_MM:
            af = _airfoil(chord)
            assert resolve_web_anchors(af, T_SKIN, ANCHORS) == resolve_web_anchors(
                af, T_SKIN, ANCHORS, min_cavity_depth_m=0.0
            )

    def test_drops_exactly_the_stations_that_cannot_carry_a_gauge_cap(self):
        """Not a blunt instrument: only where the rule is actually violated.

        A criterion that dropped more would be trading structure for tidiness;
        one that dropped fewer would leave a sub-gauge cap in the export.
        """
        for station, chord in zip(STATIONS, CHORD_MM):
            _af, before = _resolve(chord, 0.0)
            _af2, after = _resolve(chord, THRESHOLD)
            lost = len(before) - len(after)
            if station in SUB_GAUGE_STATIONS:
                assert lost == 1, (
                    f"station {station} should lose its aft web; lost {lost}"
                )
            else:
                assert lost == 0, (
                    f"station {station} lost {lost} web(s) it could have carried"
                )

    def test_reports_why_the_web_went(self):
        """Callers stepping the layout down need the reason, not just the count.

        ``cavity_too_shallow`` has to be distinguishable from the topological
        drops: it means the section could hold the web but not a cap against it,
        which is a different engineering statement.
        """
        af, _ = _resolve(94.9, THRESHOLD)  # station 59
        _resolved, dropped = resolve_web_anchors(
            af, T_SKIN, ANCHORS, min_cavity_depth_m=THRESHOLD, return_dropped=True
        )
        reasons = [d.reason for d in dropped]
        assert "cavity_too_shallow" in reasons, reasons
        shallow = next(d for d in dropped if d.reason == "cavity_too_shallow")

        # The aft-most web that could anchor at all -- not the aft-most
        # *requested*. By station 59 the outer three anchors lie beyond the
        # 94.9 mm chord and have already gone for no_inner_skin_intercept, so
        # the shallow one is the last survivor of what remained.
        anchorable = resolve_web_anchors(af, T_SKIN, ANCHORS)
        assert shallow.x_m == pytest.approx(anchorable[-1][0], abs=1e-9), (
            "the aft-most anchorable web is the one that should go"
        )
        assert shallow.x_m < max(ANCHORS), (
            "fixture drifted: the aft-most requested web now anchors, so this "
            "no longer tests the survivor case"
        )

    def test_a_web_that_still_has_an_intercept_can_still_be_dropped(self):
        """The whole point: the old criterion would have kept this web.

        If the shallow station ever stops resolving an intercept on its own,
        this fixture has stopped exercising the gap between 'anchorable' and
        'usable' and the criterion is untested.
        """
        af, _ = _resolve(94.9, 0.0)
        kept = resolve_web_anchors(af, T_SKIN, ANCHORS)  # no threshold
        probe = inner_y_at_xs(
            af, T_SKIN, [kept[-1][0]], y_search=None, phi_inner=None,
            bisection_iters=0, raise_on_miss=False,
        )[0]
        assert probe is not None, "aft web has no intercept; fixture no longer bites"
        depth = probe[0] - probe[1]
        assert 0.0 < depth < THRESHOLD, (
            f"aft web cavity is {depth * 1e3:.3f} mm — needs to be anchorable but "
            f"below the {THRESHOLD * 1e3:.3f} mm threshold for this test to mean anything"
        )


class TestInvariant:
    def test_every_surviving_web_leaves_a_buildable_cap(self):
        """The property the criterion exists to establish, across the span."""
        worst = float("inf")
        for station, chord in zip(STATIONS, CHORD_MM):
            af, resolved = _resolve(chord, THRESHOLD)
            assert resolved, f"station {station} lost every web"
            allow = _cap_allowance(af, resolved)
            worst = min(worst, allow)
            assert allow >= GAUGE - 1e-12, (
                f"station {station} allows only {allow * 1e3:.3f} mm of cap "
                f"against a {GAUGE * 1e3:.1f} mm gauge"
            )
        assert np.isfinite(worst)

    def test_the_invariant_did_not_hold_before(self):
        """Guard on the guard: without the threshold the export is sub-gauge."""
        violations = {
            station
            for station, chord in zip(STATIONS, CHORD_MM)
            if (lambda af_r: _cap_allowance(*af_r) < GAUGE)(_resolve(chord, 0.0))
        }
        assert violations == SUB_GAUGE_STATIONS, (
            f"expected {sorted(SUB_GAUGE_STATIONS)} to be sub-gauge without the "
            f"threshold; got {sorted(violations)}"
        )
