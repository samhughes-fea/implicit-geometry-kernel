"""Shell nodes sit on the mid-surface, or on the outer mould line."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station


def _blade(reference="mid"):
    foil = Airfoil.ellipse(1.0)
    stations = [
        Station(0.0, 2.0, foil, pitch_axis=0.5),
        Station(4.0, 2.0, foil, pitch_axis=0.5),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.02,
        spar_cap_thickness=0.01,
        shear_web_thickness=0.008,
        web_positions=(0.3, 0.7),
        pivot_fraction=0.5,
        wall_reference=reference,
    )


def _span(graph):
    pts = np.vstack([run.points[:, :2] for run in graph.loop if run.kind == "skin"])
    return float(pts[:, 0].max() - pts[:, 0].min()), float(pts[:, 1].max() - pts[:, 1].min())


def test_mid_is_the_default():
    assert _blade().wall_reference == "mid"


def test_oml_reference_moves_the_skin_out_by_half_the_thickness():
    mid_x, mid_y = _span(_blade("mid").section_at(1.0).shell_graph())
    oml_x, oml_y = _span(_blade("oml").section_at(1.0).shell_graph())
    # both sides of a circle step out by half the 0.02 m skin
    assert oml_x - mid_x == pytest.approx(0.02, abs=2e-3)
    assert oml_y - mid_y == pytest.approx(0.02, abs=2e-3)


def test_unknown_wall_reference_is_rejected():
    with pytest.raises(ValueError):
        _blade("inner")


def _web_line_offset(web, end):
    """Distance of a web end from the line through the web's two nearest points."""
    pts = web.points[:, :2]
    a, b = (pts[1], pts[2]) if end == 0 else (pts[-2], pts[-3])
    d = (a - b) / np.linalg.norm(a - b)
    rel = pts[end] - a
    return abs(float(rel[0] * d[1] - rel[1] * d[0]))


@pytest.mark.parametrize("reference", ["mid", "oml"])
def test_skin_partitions_and_web_feet_share_one_point(reference):
    """SG-035: OML offset per run used to split every junction (19 mm here)."""
    g = _blade(reference).section_at(1.0).shell_graph()
    assert g.closed_loop_chain(rel_tol=1e-12) is not None
    ends = np.array([p for r in g.loop for p in (r.points[0, :2], r.points[-1, :2])])
    for web in g.webs:
        for end in (0, -1):
            foot = web.points[end, :2]
            assert float(np.min(np.linalg.norm(ends - foot, axis=1))) == 0.0


def test_oml_web_feet_follow_the_web_line():
    """The web reaches the outer skin along itself; it is not slid sideways."""
    g = _blade("oml").section_at(1.0).shell_graph()
    for web in g.webs:
        if len(web.points) < 3:
            continue
        for end in (0, -1):
            assert _web_line_offset(web, end) < 1e-9
