"""A shear web must stay straight when its feet are snapped to the skin.

The web midline is a capsule axis — dead straight on emit. Junction preparation
has to extend each foot until it meets the skin midline, and the only direction
that preserves straightness is the web's **own axis**.

``split_skin_at_junctions`` used to snap the foot with
``point_to_open_polyline_nearest``, i.e. the nearest point on the *skin*. That is
a projection perpendicular to the skin, which coincides with the web axis only
where the two meet at exactly 90 degrees — never, on a cambered or twisted
section. The foot slid sideways and the last segment of the web picked up a
kink: 1.0-3.4 deg at zero twist, **8.45 deg at 15 deg twist**, putting the
junction node 1.19 mm off the web axis on an 8 mm web.

Nothing pinned that, which is why it went unnoticed. These tests pin it.

The second test guards the coupling that makes the fix correct rather than
merely straight: the same block derives ``best_s``, the arc length at which the
skin is **split**, from the foot point. Move the foot without moving the split
and the kink is traded for a gap between the web foot and the skin node.
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

LAYOUTS = ["2D-F", "2D-CN", "3D-F", "3D-CN", "5D-F"]
# Twist is the aggravating factor: it inclines the skin relative to a flapwise
# web, so the perpendicular and axial directions diverge further.
TWIST_DEG = [0.0, 15.0]
AIRFOILS = ["2412", "0012", "4415"]

_CHORD_M = 1.0
_MAX_KINK_DEG = 0.05
_MAX_FOOT_GAP_M = 1e-9


def _prepared_strips(layout_key: str, twist_deg: float, naca: str):
    twist_rad = float(np.deg2rad(twist_deg))
    section = build_section_view(
        AirfoilSDF.from_naca(naca, n_points=200, chord=_CHORD_M),
        resolve_system_type(layout_key),
        twist_rad=twist_rad,
        skin_thickness=0.006,
        spar_cap_thickness=0.014,
        shear_web_thickness=0.008,
    )
    raw = build_shell_midline_strips(section, twist_rad=twist_rad, chord_m=_CHORD_M)
    return raw, prepare_shell_midline_junctions(raw)


def _webs(strips):
    return {
        str(s.label): np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for s in strips
        if str(s.kind) == "web"
    }


@pytest.mark.parametrize("naca", AIRFOILS)
@pytest.mark.parametrize("twist_deg", TWIST_DEG)
@pytest.mark.parametrize("layout_key", LAYOUTS)
def test_prepared_web_stays_straight(layout_key, twist_deg, naca):
    raw, prepared = _prepared_strips(layout_key, twist_deg, naca)
    raw_webs, prepared_webs = _webs(raw), _webs(prepared)
    assert raw_webs, f"{layout_key} produced no web strips"

    for label, prepared_points in prepared_webs.items():
        assert prepared_points.shape[0] >= 2
        axis = prepared_points[-1] - prepared_points[0]
        axis = axis / np.linalg.norm(axis)

        segments = np.diff(prepared_points, axis=0)
        lengths = np.linalg.norm(segments, axis=1)
        keep = lengths > 1e-12
        directions = segments[keep] / lengths[keep, None]
        # Signed area against the axis gives sin(angle) directly.
        cross = np.abs(directions[:, 0] * axis[1] - directions[:, 1] * axis[0])
        worst_deg = float(np.degrees(np.arcsin(np.clip(cross, 0.0, 1.0))).max())

        assert worst_deg <= _MAX_KINK_DEG, (
            f"{layout_key} {label} @ {twist_deg} deg twist, NACA {naca}: a segment "
            f"leans {worst_deg:.3f} deg off the web axis. The foot is probably "
            "being snapped perpendicular to the skin again instead of along the "
            "web axis."
        )


@pytest.mark.parametrize("twist_deg", TWIST_DEG)
@pytest.mark.parametrize("layout_key", LAYOUTS)
def test_web_foot_lands_exactly_on_a_skin_node(layout_key, twist_deg):
    """Straightening the foot must not decouple it from the skin split."""
    _raw, prepared = _prepared_strips(layout_key, twist_deg, "2412")
    skin_vertices = np.concatenate(
        [
            np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
            for s in prepared
            if str(s.kind) == "skin"
        ]
    )
    for label, points in _webs(prepared).items():
        for name, endpoint in (("start", points[0]), ("end", points[-1])):
            gap = float(np.linalg.norm(skin_vertices - endpoint, axis=1).min())
            assert gap <= _MAX_FOOT_GAP_M, (
                f"{layout_key} {label} {name} @ {twist_deg} deg: foot sits "
                f"{gap * 1e3:.4f} mm from the nearest skin node — the skin split "
                "arc length no longer matches the foot point."
            )


def test_tip_absolute_2a_f_loops_build_under_twist():
    """Tip-scale 2A-F with absolute LE web + tip twist must emit SG-002-clean feet.

    Captured from WP6 confirmatory tip mono DV (chord ≈ 88 mm, skin 3.536 mm,
    web 4.42 mm, twist ≈ 0.41 deg). Fractional layouts already passed; this pins
    the absolute-web tip path that BL-006 steps into (TE web dropped → 1A-F).
    """
    from blade_structure_utilities.geometry_utilities.section_geometry.api.worker import (
        build_section_geometry,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.cell_wall_loops.cell_wall_loop_builders import (
        build_cell_wall_loops,
    )

    chord = 0.08775168
    twist_rad = float(np.deg2rad(0.407709))
    af = AirfoilSDF.from_naca("2412", n_points=200, chord=chord)
    # Tip plan webs: LE kept, TE past tip chord dropped → 1A-F survivors.
    result = build_section_geometry(
        layout=resolve_system_type("1A-F"),
        airfoil=af,
        skin_thickness=0.003536,
        shear_web_thickness=0.00442,
        spar_cap_thickness=0.0,
        twist_rad=twist_rad,
        web_x_positions_m=[0.036130035],
        outputs=("midlines", "cell_wall_loops"),
    )
    strips = result.midline_strips or ()
    webs = [s for s in strips if s.kind == "web"]
    assert webs, "tip must keep the LE web"
    assert result.cell_wall_loops and len(result.cell_wall_loops) >= 2
    build_cell_wall_loops(strips)


def test_tip_fractional_2a_f_loops_build_under_twist():
    """Tip-scale catalogue 2A-F (fractional webs) still builds loops under tip twist."""
    from blade_structure_utilities.geometry_utilities.section_geometry.api.worker import (
        build_section_geometry,
    )

    chord = 0.08775168
    twist_rad = float(np.deg2rad(0.407709))
    af = AirfoilSDF.from_naca("2412", n_points=200, chord=chord)
    result = build_section_geometry(
        layout=resolve_system_type("2A-F"),
        airfoil=af,
        skin_thickness=0.003536,
        shear_web_thickness=0.00442,
        spar_cap_thickness=0.0,
        twist_rad=twist_rad,
        outputs=("midlines", "cell_wall_loops"),
    )
    webs = [s for s in (result.midline_strips or ()) if s.kind == "web"]
    assert len(webs) == 2
    assert result.cell_wall_loops and len(result.cell_wall_loops) == 3


def test_stale_frozen_web_y_under_twist_raises_infeasible_not_sg002_gap():
    """Frozen thin-skin web_y + thicker skin + twist must not emit a blind foot.

    Pre-fix: flapwise_web_sframe_inputs rotated stale IML anchors and
    ShearWeb._extend_to_skin_mid_surface returned origin+t/2, leaving the foot
    ~3.5 mm off the skin mid (SG-002). Now the assemble path raises
    InfeasibleSectionGeometry so mono can keep prior gauges.
    """
    from blade_structure_utilities.geometry_utilities.section_contract.errors import (
        InfeasibleSectionGeometry,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.api.worker import (
        build_section_geometry,
        resolve_station_layout,
    )

    chord = 0.08775168
    twist_rad = float(np.deg2rad(0.407709))
    af = AirfoilSDF.from_naca("2412", n_points=200, chord=chord)
    skin_cold = 0.003536
    skin_mono = 0.007072
    web_x = [0.036130035]
    res = resolve_station_layout(
        resolve_system_type("1A-F"),
        af,
        skin_thickness=skin_cold,
        web_x_positions_m=web_x,
    )
    assert list(res.web_x_positions_m), "cold tip must host the LE web"
    # Matched gauges still build.
    ok = build_section_geometry(
        layout=resolve_system_type("1A-F"),
        airfoil=af,
        skin_thickness=skin_cold,
        shear_web_thickness=0.00442,
        spar_cap_thickness=0.0,
        twist_rad=twist_rad,
        web_x_positions_m=list(res.web_x_positions_m),
        web_y_coords=list(res.web_y_coords),
        outputs=("midlines", "cell_wall_loops"),
    )
    assert ok.cell_wall_loops
    with pytest.raises(InfeasibleSectionGeometry, match="frozen web_y|re-anchor|inner skin"):
        build_section_geometry(
            layout=resolve_system_type("1A-F"),
            airfoil=af,
            skin_thickness=skin_mono,
            shear_web_thickness=0.00442,
            spar_cap_thickness=0.0,
            twist_rad=twist_rad,
            web_x_positions_m=list(res.web_x_positions_m),
            web_y_coords=list(res.web_y_coords),
            outputs=("midlines", "cell_wall_loops"),
        )
