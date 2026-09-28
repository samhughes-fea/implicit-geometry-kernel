"""Topology preparation must not flatten the spar-cap midline.

``resample_cap_interiors_from_skin_rays`` replaces a cap's interior vertices with
points found by casting inward rays from the skin vertices. Before casting, it
rejects skin vertices outside the cap's bounding box plus a margin — and that
margin used to be the fixed constant ``_SKIN_RAY_CAP_BBOX_MARGIN_M = 5 mm``.

The skin midline sits ``(skin_thickness + cap_thickness) / 2`` from the cap
midline. Whenever that separation exceeded the margin (e.g. 10 mm for a 6 mm
skin and a 14 mm cap) the filter discarded nearly every skin vertex above the
cap: on ``2D-F`` only 6 of 24 survived, all of them clustered where the skin
curves toward the web feet, so the prepared cap became a straight chord across
its middle 80 %. Sagitta fell 10.14 mm -> 6.12 mm and the mesh built on it put
nodes 5.97 mm off the true cap midline — 43 % of the cap thickness, and not
recoverable by refining, since refinement subdivides the wrong chord.

Nothing pinned that, which is why it went unnoticed. These tests pin it.

Families B and C were originally excluded here, because
``split_skin_at_junctions`` snapped their free, mid-bay cap ends onto the skin
midline and their caps were displaced ~10 mm regardless of this margin. That
relocation has since been removed — a cap end is either already at a web
junction or is free, and neither belongs on the skin — so B and C are covered
too and the suite no longer has a known-bad hole in it.

The bounding-box pre-filter itself is also gone: cap resampling now selects skin
vertices by the strips' declared suction/pressure side, so there is no margin
left to get wrong. These tests pin the *outcome* (the prepared cap tracks the
raw cap) rather than the mechanism, so they survived that change unaltered.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_structure_utilities.geometry_utilities.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.airfoil import AirfoilSDF
from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import (
    prepare_shell_midline_junctions,
)

# D terminates its caps on webs; B (single band) and C (discrete) terminate them
# mid-bay with free ends. Both shapes are covered.
BOX_CAP_LAYOUTS = ["2D-F", "2D-CN", "3D-F", "5D-CN", "2B-F", "2C-F"]

# (skin, cap) gauges in metres. The first pair is inside the old 5 mm margin and
# always passed; the rest exceed it and are the regression cases.
GAUGES = [(0.003, 0.006), (0.006, 0.014), (0.010, 0.030)]

_CHORD_M = 1.0
_MAX_DEVIATION_FRAC_OF_CAP = 0.02


def _prepared_and_raw_caps(layout_key: str, skin_thickness: float, cap_thickness: float):
    airfoil = AirfoilSDF.from_naca("2412", n_points=200, chord=_CHORD_M)
    section = build_section_view(
        airfoil,
        resolve_system_type(layout_key),
        twist_rad=0.0,
        skin_thickness=skin_thickness,
        spar_cap_thickness=cap_thickness,
        shear_web_thickness=0.008,
    )
    raw = build_shell_midline_strips(section, twist_rad=0.0, chord_m=_CHORD_M)
    prepared = prepare_shell_midline_junctions(raw)

    def caps(strips):
        return {
            str(s.label): np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
            for s in strips
            if str(s.kind) == "cap"
        }

    return caps(raw), caps(prepared)


def _max_distance_to_polyline(points: np.ndarray, polyline: np.ndarray) -> float:
    """Worst distance from any of ``points`` to the ``polyline`` as a curve.

    Sampling the prepared vertices against the raw curve instead would pass
    trivially: every resampled vertex is a ray hit and so lies on the raw cap by
    construction. The failure is the raw curve departing from the chords drawn
    between those vertices, which is what this direction measures.
    """
    start, end = polyline[:-1], polyline[1:]
    segment = end - start
    worst = 0.0
    for point in points:
        t = np.clip(
            ((point - start) * segment).sum(1) / ((segment * segment).sum(1) + 1e-30),
            0.0,
            1.0,
        )
        projected = start + t[:, None] * segment
        worst = max(worst, float(np.linalg.norm(projected - point, axis=1).min()))
    return worst


@pytest.mark.parametrize("skin_thickness,cap_thickness", GAUGES)
@pytest.mark.parametrize("layout_key", BOX_CAP_LAYOUTS)
def test_prepared_cap_tracks_the_raw_cap(layout_key, skin_thickness, cap_thickness):
    raw_caps, prepared_caps = _prepared_and_raw_caps(
        layout_key, skin_thickness, cap_thickness
    )
    assert raw_caps, f"{layout_key} produced no cap strips"
    tolerance = _MAX_DEVIATION_FRAC_OF_CAP * cap_thickness

    for label, raw_points in raw_caps.items():
        prepared_points = prepared_caps.get(label)
        assert prepared_points is not None, f"cap {label} vanished during preparation"
        deviation = _max_distance_to_polyline(raw_points, prepared_points)
        assert deviation <= tolerance, (
            f"{layout_key} {label}: prepared cap departs from the raw cap by "
            f"{deviation * 1e3:.3f} mm (tolerance {tolerance * 1e3:.3f} mm). "
            "The skin-ray bbox margin is likely narrower than the "
            "skin-to-cap separation again."
        )


@pytest.mark.parametrize("skin_thickness,cap_thickness", GAUGES)
@pytest.mark.parametrize("layout_key", BOX_CAP_LAYOUTS)
def test_prepared_cap_keeps_its_curvature(layout_key, skin_thickness, cap_thickness):
    """Sagitta is the quantity that collapsed; pin it directly."""

    def sagitta(points: np.ndarray) -> float:
        chord = points[-1] - points[0]
        length = float(np.linalg.norm(chord))
        if length < 1e-12:
            return 0.0
        normal = np.array([-chord[1], chord[0]]) / length
        return float(np.abs((points - points[0]) @ normal).max())

    raw_caps, prepared_caps = _prepared_and_raw_caps(
        layout_key, skin_thickness, cap_thickness
    )
    for label, raw_points in raw_caps.items():
        raw_sagitta = sagitta(raw_points)
        if raw_sagitta < 1e-4:
            continue  # a nearly straight cap has no curvature to lose
        prepared_sagitta = sagitta(prepared_caps[label])
        assert prepared_sagitta == pytest.approx(raw_sagitta, rel=0.02), (
            f"{layout_key} {label}: cap sagitta {raw_sagitta * 1e3:.3f} mm -> "
            f"{prepared_sagitta * 1e3:.3f} mm"
        )


@pytest.mark.parametrize("skin_thickness,cap_thickness", GAUGES)
@pytest.mark.parametrize("layout_key", BOX_CAP_LAYOUTS)
def test_topology_never_relocates_a_cap_endpoint(
    layout_key, skin_thickness, cap_thickness
):
    """Cap ends are at their true extent on emit; preparation must not move them.

    This is the invariant, not a tolerance: a cap end is either already on a web
    mid-surface (measured 0.000 mm for family D) or a free end mid-bay (40-149 mm
    from any web on B and C). Relocating either was the ~10 mm displacement bug.
    """
    raw_caps, prepared_caps = _prepared_and_raw_caps(
        layout_key, skin_thickness, cap_thickness
    )
    for label, raw_points in raw_caps.items():
        prepared_points = prepared_caps[label]
        for name, raw_end, prepared_end in (
            ("start", raw_points[0], prepared_points[0]),
            ("end", raw_points[-1], prepared_points[-1]),
        ):
            moved = float(np.linalg.norm(prepared_end - raw_end))
            assert moved <= 1e-12, (
                f"{layout_key} {label} {name}: topology moved the cap endpoint "
                f"{moved * 1e3:.4f} mm"
            )
