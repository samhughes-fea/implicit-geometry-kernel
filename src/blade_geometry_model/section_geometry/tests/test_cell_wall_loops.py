"""Cell wall loops are the faces of the labelled wall graph (SG-014).

Fixtures are in canonical B, as production midlines are: ``X_B`` is flapwise
(thickness) and ``Y_B`` chordwise. The previous fixtures were chord-first --
webs at constant ``x`` -- which is the one convention under which the old
builder happened to work, and why nothing caught it pairing the wrong web ends
on every B-frame section with flapwise webs.
"""

from __future__ import annotations

import itertools
import warnings

import numpy as np
import pytest

from blade_structure_utilities.geometry_utilities.section_contract import CellWallLoop, ShellMidlineStrip
from blade_structure_utilities.geometry_utilities.section_geometry.cell_wall_loops.cell_wall_geom import (
    shoelace_area,
    stitch_skin_ring,
)
from blade_structure_utilities.geometry_utilities.section_geometry.cell_wall_loops.cell_wall_loops import (
    build_cell_wall_loops,
    cell_wall_loops_case_payload,
)

HALF_DEPTH = 0.5
"""Flapwise half-depth of the box section, along ``X_B``."""
HALF_CHORD = 1.0
"""Chordwise half-length, along ``Y_B``; the LE is at ``Y_B = -HALF_CHORD``."""


def _box_skin(*, samples_per_side: int = 9) -> ShellMidlineStrip:
    """Closed B-frame box: suction at ``X_B = +HALF_DEPTH``, pressure at ``-HALF_DEPTH``."""
    chord = np.linspace(-HALF_CHORD, HALF_CHORD, samples_per_side)
    depth = np.linspace(HALF_DEPTH, -HALF_DEPTH, samples_per_side)
    ring = np.vstack(
        [
            np.column_stack([np.full_like(chord, HALF_DEPTH), chord[::-1]]),  # suction, TE -> LE
            np.column_stack([depth[1:], np.full(samples_per_side - 1, -HALF_CHORD)]),  # nose
            np.column_stack([np.full_like(chord[1:], -HALF_DEPTH), chord[1:]]),  # pressure, LE -> TE
            np.column_stack([depth[::-1][1:], np.full(samples_per_side - 1, HALF_CHORD)]),  # tail
        ]
    )
    return ShellMidlineStrip(
        label="outer_skin",
        kind="skin",
        midline_twisted_pitch_axis_blade_frame=ring,
        thickness_m=0.01,
        closed=True,
    )


def _flapwise_web(chord_position: float, label: str, *, reverse: bool = False) -> ShellMidlineStrip:
    """A web at constant ``Y_B`` spanning the depth -- skin to skin, as emitted (SG-002)."""
    points = np.column_stack(
        [np.linspace(-HALF_DEPTH, HALF_DEPTH, 5), np.full(5, chord_position)]
    )
    return ShellMidlineStrip(
        label=label,
        kind="web",
        midline_twisted_pitch_axis_blade_frame=points[::-1].copy() if reverse else points,
        thickness_m=0.008,
        closed=False,
        alignment="flapwise",
    )


def _box_area() -> float:
    return (2.0 * HALF_DEPTH) * (2.0 * HALF_CHORD)


def _membership(loops: tuple[CellWallLoop, ...]) -> dict[str, list[int]]:
    owners: dict[str, list[int]] = {}
    for loop in loops:
        for label in loop.wall_labels:
            owners.setdefault(label, []).append(loop.id)
    return owners


def _assert_planar_subdivision(strips, loops) -> None:
    """Every web bounds two cells, every skin arc one, and the cells tile the ring."""
    ring = stitch_skin_ring([s for s in strips if s.kind == "skin"])
    assert sum(loop.area_m2 for loop in loops) == pytest.approx(abs(shoelace_area(ring)), rel=1e-9)
    owners = _membership(loops)
    for label, cells in owners.items():
        kind = next(s.kind for s in strips if s.label == label) if any(s.label == label for s in strips) else "skin"
        assert len(cells) == (2 if kind == "web" else 1), (label, cells)
    assert all(loop.is_closed() for loop in loops)


def test_no_web_is_one_cell_bounded_by_the_skin():
    skin = _box_skin()
    loops = build_cell_wall_loops([skin])
    assert len(loops) == 1
    assert loops[0].web_labels == ()
    assert loops[0].wall_labels == ("outer_skin",)
    assert loops[0].area_m2 == pytest.approx(_box_area(), rel=1e-12)


@pytest.mark.parametrize("n_webs", [1, 2, 3, 4, 5])
def test_flapwise_webs_give_contiguous_le_to_te_cells(n_webs):
    """The SG-014 failure: flapwise webs have both ends at the same ``Y_B``."""
    positions = np.linspace(-0.7, 0.7, n_webs)
    webs = [
        _flapwise_web(float(y), f"shear_web_{i}", reverse=bool(i % 2))
        for i, y in enumerate(positions)
    ]
    strips = [_box_skin(), *webs]
    loops = build_cell_wall_loops(strips)
    assert len(loops) == n_webs + 1
    edges = [-HALF_CHORD, *positions, HALF_CHORD]
    for loop, (lo, hi) in zip(loops, zip(edges[:-1], edges[1:])):
        chord = loop.polyline_twisted_pitch_axis_blade_frame[:, 1]
        assert chord.min() == pytest.approx(lo, abs=1e-12)
        assert chord.max() == pytest.approx(hi, abs=1e-12)
        assert loop.area_m2 == pytest.approx((hi - lo) * 2.0 * HALF_DEPTH, rel=1e-12)
    for index, loop in enumerate(loops[1:-1], start=1):
        assert loop.web_labels == (f"shear_web_{index - 1}", f"shear_web_{index}")
    _assert_planar_subdivision(strips, loops)


def test_web_list_order_and_direction_do_not_matter():
    webs = [_flapwise_web(y, f"shear_web_{i}") for i, y in enumerate((-0.5, 0.0, 0.5))]
    reference = build_cell_wall_loops([_box_skin(), *webs])
    shuffled = build_cell_wall_loops(
        [webs[2], _box_skin(), _flapwise_web(-0.5, "shear_web_0", reverse=True), webs[1]]
    )
    assert [l.area_m2 for l in shuffled] == pytest.approx([l.area_m2 for l in reference], rel=1e-12)
    assert [set(l.wall_labels) for l in shuffled] == [set(l.wall_labels) for l in reference]


def test_membership_does_not_depend_on_the_frame():
    """Rotating the whole section must not change which walls bound which cell.

    Which end is the *leading* edge is defined in B, so under an arbitrary
    rotation only the LE-to-TE numbering may change -- never the cells.
    """
    webs = [_flapwise_web(y, f"shear_web_{i}") for i, y in enumerate((-0.6, -0.1, 0.3, 0.8))]
    strips = [_box_skin(), *webs]
    reference = build_cell_wall_loops(strips)
    for angle in (0.3, np.pi / 2, 2.0, -1.1):
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        rotated = [
            ShellMidlineStrip(
                label=s.label,
                kind=s.kind,
                midline_twisted_pitch_axis_blade_frame=np.asarray(s.midline_twisted_pitch_axis_blade_frame) @ rotation.T,
                thickness_m=s.thickness_m,
                closed=s.closed,
            )
            for s in strips
        ]
        loops = build_cell_wall_loops(rotated)
        cells = {frozenset(l.wall_labels) - _skin_arc_labels(l): round(l.area_m2, 12) for l in loops}
        expected = {frozenset(l.wall_labels) - _skin_arc_labels(l): round(l.area_m2, 12) for l in reference}
        assert cells == expected, angle


def _skin_arc_labels(loop: CellWallLoop) -> frozenset[str]:
    """Skin-arc names carry the suction/pressure and LE/TE naming, which is B-relative."""
    return frozenset(label for label in loop.wall_labels if label.startswith("skin_"))


def test_a_closed_b_frame_skin_is_split_into_suction_and_pressure_by_x_b():
    """SG-015: the closed-ring split reads B, where suction is +X_B and the LE is -Y_B."""
    from blade_structure_utilities.geometry_utilities.section_geometry.midlines.skin_split import (
        split_skin_at_junctions,
    )

    prepared = split_skin_at_junctions((_box_skin(), _flapwise_web(0.0, "shear_web_0")))
    for strip in (s for s in prepared if s.kind == "skin"):
        depth = np.asarray(strip.midline_twisted_pitch_axis_blade_frame)[:, 0]
        if strip.skin_side == "suction_surface":
            assert depth.mean() > 0.0, strip.label
        else:
            assert depth.mean() < 0.0, strip.label
    loops = build_cell_wall_loops(list(prepared))
    assert loops[0].polyline_twisted_pitch_axis_blade_frame[:, 1].min() == pytest.approx(-HALF_CHORD)


def test_caps_bound_no_cell():
    cap = ShellMidlineStrip(
        label="spar_cap_suction_surface",
        kind="cap",
        midline_twisted_pitch_axis_blade_frame=np.array([[0.45, -0.3], [0.45, 0.3]]),
        thickness_m=0.02,
        closed=False,
        surface="suction_surface",
    )
    loops = build_cell_wall_loops([_box_skin(), _flapwise_web(0.0, "shear_web_0"), cap])
    assert len(loops) == 2
    assert all("spar_cap_suction_surface" not in loop.wall_labels for loop in loops)


def test_a_web_that_stops_short_of_the_skin_is_refused():
    """Extrapolating it is relocation, which SG-002 forbids; the builder reports it."""
    short = ShellMidlineStrip(
        label="shear_web_0",
        kind="web",
        midline_twisted_pitch_axis_blade_frame=np.array([[-0.4, 0.0], [0.4, 0.0]]),
        thickness_m=0.008,
        closed=False,
    )
    with pytest.raises(ValueError, match="skin junction"):
        build_cell_wall_loops([_box_skin(), short])


def test_crossing_webs_are_refused():
    crossing = ShellMidlineStrip(
        label="shear_web_1",
        kind="web",
        midline_twisted_pitch_axis_blade_frame=np.array([[-HALF_DEPTH, -0.5], [HALF_DEPTH, 0.5]]),
        thickness_m=0.008,
        closed=False,
    )
    with pytest.raises(ValueError, match="cross"):
        build_cell_wall_loops([_box_skin(), _flapwise_web(0.0, "shear_web_0"), crossing])


def test_two_webs_on_one_junction_are_refused():
    with pytest.raises(ValueError, match="same junction"):
        build_cell_wall_loops(
            [_box_skin(), _flapwise_web(0.0, "shear_web_0"), _flapwise_web(0.0, "shear_web_1")]
        )


def test_hierarchical_skin_arcs_use_ordered_ring_path():
    """Long skin_side/skin_arc tokens stitch via the hierarchical ring order."""
    from blade_structure_utilities.geometry_utilities.section_geometry.midlines.skin_split import (
        split_skin_at_junctions,
    )

    prepared = list(split_skin_at_junctions((_box_skin(), _flapwise_web(0.0, "shear_web_0"))))
    skins = [s for s in prepared if s.kind == "skin"]
    assert {s.skin_arc for s in skins} == {"leading_edge", "trailing_edge"}
    ring = stitch_skin_ring(list(reversed(skins)))
    assert abs(shoelace_area(ring)) == pytest.approx(_box_area(), rel=1e-12)
    loops = build_cell_wall_loops(list(reversed(prepared)))
    assert len(loops) == 2
    assert {"skin_suction_surface:leading_edge", "skin_pressure_surface:leading_edge"} <= set(loops[0].wall_labels)
    _assert_planar_subdivision(prepared, loops)


@pytest.mark.parametrize(
    "key, twist_rad, chord",
    list(itertools.product(["1C-F", "2D-F", "2D-CN", "3D-F", "3C-F", "4D-F", "4D-CN", "2B-F"], (0.0, 0.2, 0.6, 1.2), (1.0, 0.2))),
)
def test_worker_loops_tile_every_layout_twist_and_chord(key, twist_rad, chord):
    """The SG-014 stress test: 36 of these 64 emitted overlapping cells before the fix."""
    from blade_structure_utilities.load_utilities.section_loads.visualisation.demos._common import (
        build_midlines_for_key,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        strips, loops, _ = build_midlines_for_key(
            key,
            chord=chord,
            twist_rad=twist_rad,
            skin_thickness=0.004 * chord,
            spar_cap_thickness=0.01 * chord,
            shear_web_thickness=0.004 * chord,
        )
    n_webs = sum(1 for s in strips if s.kind == "web")
    assert len(loops) == n_webs + 1
    _assert_planar_subdivision(strips, loops)


@pytest.mark.parametrize("n_webs", [0, 1, 2])
def test_cell_wall_loops_case_payload_shape(n_webs):
    strips = [_box_skin(), *(_flapwise_web(y, f"shear_web_{i}") for i, y in enumerate(np.linspace(-0.4, 0.4, n_webs)))]
    payload = cell_wall_loops_case_payload(strips, name=f"n_webs_{n_webs}")
    assert payload["name"] == f"n_webs_{n_webs}"
    assert payload["n_webs"] == n_webs
    assert payload["n_loops"] == n_webs + 1
    assert len(payload["loops"]) == n_webs + 1
    for loop in payload["loops"]:
        assert isinstance(loop["id"], int)
        assert isinstance(loop["area_m2"], float)
        assert isinstance(loop["web_labels"], list)
        assert isinstance(loop["wall_labels"], list)
        assert isinstance(loop["polyline_twisted_pitch_axis_blade_frame"], list)
        assert all(isinstance(v, float) for pt in loop["polyline_twisted_pitch_axis_blade_frame"] for v in pt)
