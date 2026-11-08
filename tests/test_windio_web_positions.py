"""WindIO nd_arc anchors → web chord fractions (G8 median + span tables)."""

import numpy as np
import pytest

from implicit_geometry_kernel.blade.windio import (
    cap_chord_windows_from_windio_structure,
    cap_window_table_from_windio_structure,
    nd_arc_to_chord_frac,
    web_position_table_from_windio_structure,
    web_positions_fn_from_windio_structure,
    web_positions_from_windio_structure,
)


def _web(name: str, start: float, end: float) -> dict:
    return {
        "name": name,
        "start_nd_arc": {"grid": [0.0, 1.0], "values": [start, start]},
        "end_nd_arc": {"grid": [0.0, 1.0], "values": [end, end]},
    }


def test_nd_arc_le_te():
    assert nd_arc_to_chord_frac(0.5) == 0.0
    assert nd_arc_to_chord_frac(0.0) == 1.0
    assert nd_arc_to_chord_frac(1.0) == 1.0
    assert nd_arc_to_chord_frac(0.39) == nd_arc_to_chord_frac(0.61)


def test_web_positions_from_nd_arc_match_old_hardcoded_pair():
    structure = {
        "webs": [{"name": "fore_web"}, {"name": "rear_web"}],
        "anchors": [
            _web("fore_web", 0.39, 0.61),
            _web("rear_web", 0.25, 0.75),
        ],
    }
    assert web_positions_from_windio_structure(structure) == (0.22, 0.50)


def test_three_webs_sorted_le_to_te():
    structure = {
        "webs": [{"name": "a"}, {"name": "b"}, {"name": "c"}],
        "anchors": [
            _web("a", 0.35, 0.65),
            _web("c", 0.15, 0.85),
            _web("b", 0.40, 0.60),
        ],
    }
    pos = web_positions_from_windio_structure(structure)
    assert pos == (0.20, 0.30, 0.70)


def test_missing_anchor_returns_none():
    structure = {"webs": [{"name": "fore_web"}], "anchors": []}
    assert web_positions_from_windio_structure(structure) is None


def test_cap_window_from_layer_nd_arc():
    structure = {
        "layers": [
            {
                "name": "spar_cap_ss",
                "start_nd_arc": {
                    "anchor": {"name": "spar_cap_ss", "handle": "start_nd_arc"},
                },
                "end_nd_arc": {
                    "anchor": {"name": "spar_cap_ss", "handle": "end_nd_arc"},
                },
            }
        ],
        "anchors": [
            {
                "name": "spar_cap_ss",
                "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.40, 0.40]},
                "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.30, 0.30]},
            }
        ],
    }
    windows = cap_chord_windows_from_windio_structure(structure)
    assert windows == ((0.20, 0.40),)


def test_web_position_table_varies_with_nd_arc():
    structure = {
        "webs": [{"name": "fore_web"}, {"name": "rear_web"}],
        "anchors": [
            {
                "name": "fore_web",
                "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.39, 0.42]},
                "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.61, 0.58]},
            },
            {
                "name": "rear_web",
                "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.25, 0.25]},
                "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.75, 0.75]},
            },
        ],
    }
    z_of_s = np.column_stack([[0.0, 1.0], [0.0, 10.0]])
    table = web_position_table_from_windio_structure(structure, z_of_s)
    assert table is not None
    z, vals = table
    assert vals.shape[1] == 2
    root = tuple(float(np.interp(0.0, z, vals[:, j])) for j in range(2))
    tip = tuple(float(np.interp(10.0, z, vals[:, j])) for j in range(2))
    assert root[0] == pytest.approx(0.22)
    assert root[1] == pytest.approx(0.50)
    assert tip[0] == pytest.approx(0.16)
    assert tip[1] == pytest.approx(0.50)
    assert web_positions_from_windio_structure(structure) == (0.19, 0.50)
    fn = web_positions_fn_from_windio_structure(structure, z_of_s)
    assert fn is not None
    assert fn(0.0) == root


def test_cap_window_table_varies_with_nd_arc():
    structure = {
        "layers": [
            {
                "name": "spar_cap_ss",
                "start_nd_arc": {"grid": [0.0, 1.0], "values": [0.40, 0.38]},
                "end_nd_arc": {"grid": [0.0, 1.0], "values": [0.30, 0.28]},
            }
        ],
    }
    z_of_s = np.column_stack([[0.0, 1.0], [0.0, 20.0]])
    table = cap_window_table_from_windio_structure(structure, z_of_s)
    assert table is not None
    z, vals = table
    root = tuple(float(np.interp(0.0, z, vals[:, j])) for j in range(2))
    tip = tuple(float(np.interp(20.0, z, vals[:, j])) for j in range(2))
    assert root[0] == pytest.approx(0.20)
    assert root[1] == pytest.approx(0.40)
    assert tip[0] == pytest.approx(0.24)
    assert tip[1] == pytest.approx(0.44)

