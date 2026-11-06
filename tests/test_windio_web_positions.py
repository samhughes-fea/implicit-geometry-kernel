"""WindIO nd_arc anchors → span-constant web chord fractions (G8)."""

from implicit_geometry_kernel.blade.windio import (
    cap_chord_windows_from_windio_structure,
    nd_arc_to_chord_frac,
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
