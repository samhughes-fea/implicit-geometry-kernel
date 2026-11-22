"""Tagged shell skeleton, normalised run t and the web-to-skin angle (decision 0015)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade import blade_skeleton, run_t, section_skeleton
from implicit_geometry_kernel.blade.skeleton import phi_issues, run_point_at


def _blade(**kw):
    foil = Airfoil.naca4("0018")
    stations = [
        Station(0.0, 1.0, foil, twist=5.0, pitch_axis=0.3),
        Station(6.0, 0.6, foil, twist=0.0, pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3, **kw,
    )


def test_run_t_is_normalised_arc_length():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 3.0]])
    assert run_t(pts).tolist() == [0.0, 0.25, 1.0]
    assert np.allclose(run_point_at(pts, 0.5), [1.0, 1.0])
    assert run_t(np.zeros((1, 2))).tolist() == [0.0]


def test_every_loop_partition_edge_and_web_foot_is_a_node():
    graph = _blade().section_at(2.0).shell_graph()
    sk = graph.skeleton()
    boundaries = sk.of_class("partition_boundary")
    assert len({n.key for n in boundaries if n.key[0] == "partition_boundary"}) == len(graph.partition_edges())
    feet = [n for n in sk.nodes if n.key[0] == "junction"]
    assert {(n.key[1], n.key[2]) for n in feet} == {(w.name, e) for w in graph.webs for e in (0, 1)}
    for node in feet:
        assert node.key[1] in node.runs
        assert node.runs[node.key[1]] == (0.0 if node.key[2] == 0 else 1.0)
        assert len(node.runs) >= 2


def test_web_feet_at_cap_edges_carry_both_classes():
    sk = section_skeleton(_blade().section_at(2.0).shell_graph())
    foot = sk.by_key()[("junction", "shear_web_0", 0)]
    assert {"junction", "partition_boundary"} <= foot.classes


def test_phi_is_near_normal_and_follows_web_tilt():
    def phis(tilt):
        sk = _blade(web_tilt_deg=tilt).section_at(2.0).shell_graph().skeleton()
        return {n.key[1:]: n.phi_deg for n in sk.nodes if n.key[0] == "junction"}

    base = phis(0.0)
    assert all(60.0 < v < 120.0 for v in base.values())
    tilted = phis(20.0)
    for key, value in base.items():
        assert tilted[key] == pytest.approx(value - 20.0, abs=6.0)


def test_phi_outside_the_range_is_reported():
    sk = _blade().section_at(2.0).shell_graph().skeleton()
    assert phi_issues(sk) == []
    narrow = phi_issues(sk, (89.0, 91.0))
    assert narrow and "web-to-skin angle" in narrow[0]
    assert _blade().validate() == []


def test_blade_skeleton_follows_nodes_along_the_span():
    curves = blade_skeleton(_blade(), [0.0, 2.0, 4.0, 6.0])
    feet = [c for c in curves if c.key[0] == "junction"]
    assert len(feet) == 4
    for curve in feet:
        assert curve.z.tolist() == [0.0, 2.0, 4.0, 6.0]
        assert curve.points.shape == (4, 3)
        assert np.allclose(curve.points[:, 2], curve.z)
        assert curve.starts == curve.ends == "blade_end"
        assert "termination" not in curve.classes
        assert np.all(np.isfinite(curve.phi_deg))


def test_a_partial_span_stringer_starts_and_ends_at_terminations():
    from implicit_geometry_kernel import Stringer

    geo = _blade(stringers=[Stringer("stiff", 0.6, 0.004, 0.02, extent=(1.5, 4.5))])
    curves = {c.key: c for c in blade_skeleton(geo, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])}
    stringer = [c for k, c in curves.items() if k[1] == "stiff"]
    assert stringer
    for curve in stringer:
        assert curve.starts == "termination" and curve.ends == "termination"
        assert "termination" in curve.classes
        assert curve.z.tolist() == [2.0, 3.0, 4.0]
    classes = {cls for c in stringer for cls in c.classes}
    assert {"junction", "free_edge"} <= classes
