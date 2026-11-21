"""StackThickness can take a layup or a windIO ply stack."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology.fields import StackThickness, _polyline_distance


def _blade():
    foil = Airfoil.naca4("0024")
    stations = [
        Station(0.0, 1.0, foil, pitch_axis=0.3),
        Station(4.0, 1.0, foil, pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
        web_positions=(0.2, 0.55),
        pivot_fraction=0.3,
    )


def _point_on(cut, label):
    for run in cut.shell_graph().runs():
        if run.name == label:
            point = np.asarray(run.points, float).mean(axis=0)
            if point.size == 2:
                point = np.array([point[0], point[1], cut.z])
            return point[:3]
    raise KeyError(label)


def _polyline_distance_python(xy: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Reference unsigned distance (pure NumPy / Python loops)."""
    xy = np.asarray(xy, float)
    poly = np.asarray(poly, float)
    if len(poly) < 2 or len(xy) == 0:
        return np.full(len(xy), np.inf)
    out = np.full(len(xy), np.inf)
    for i in range(len(poly) - 1):
        a = poly[i]
        b = poly[i + 1]
        ab = b - a
        denom = float(np.dot(ab, ab)) + 1e-30
        for j, p in enumerate(xy):
            t = float(np.clip(np.dot(p - a, ab) / denom, 0.0, 1.0))
            d = float(np.linalg.norm(p - (a + t * ab)))
            if d < out[j]:
                out[j] = d
    return out


def test_polyline_distance_matches_python_reference():
    rng = np.random.default_rng(0)
    poly = np.column_stack([np.linspace(0.0, 1.0, 17), 0.1 * np.sin(np.linspace(0.0, np.pi, 17))])
    xy = rng.normal(size=(64, 2)) * 0.3 + np.array([0.5, 0.0])
    got = _polyline_distance(xy, poly)
    ref = _polyline_distance_python(xy, poly)
    np.testing.assert_allclose(got, ref, rtol=0.0, atol=1e-12)


def test_stack_thickness_samples_match_nearest_run_assignment():
    geo = _blade()
    cut = geo.section_at(1.0)
    graph = cut.shell_graph()
    field = StackThickness(geo)
    samples = []
    for run in list(graph.loop) + list(graph.webs):
        mid = np.asarray(run.points, float)[len(run.points) // 2]
        if mid.size == 2:
            mid = np.array([mid[0], mid[1], cut.z])
        samples.append(mid[:3])
    pts = np.asarray(samples, float)
    T = field(pts)
    assert T.shape == (len(pts),)
    assert np.all(np.isfinite(T))
    assert float(np.min(T)) >= 0.003 - 1e-12
    assert float(np.max(T)) <= 0.003 + 0.008 + 1e-12


def test_a_layup_replaces_the_cap_stack_and_leaves_the_skin():
    geo = _blade()
    cut = geo.section_at(1.0)
    field = StackThickness.from_layup(
        geo, lambda z, name: 0.021 if "spar_cap_suction" in name else None,
    )
    cap = _point_on(cut, "spar_cap_suction_surface")
    skin_pt = _point_on(cut, "skin_suction_surface:leading_edge")
    assert field(np.asarray([cap])) == pytest.approx([0.021])
    geometric = StackThickness(geo)
    assert field(np.asarray([skin_pt])) == pytest.approx(geometric(np.asarray([skin_pt])))


def test_cap_edge_vertices_share_one_thickness():
    geo = _blade()
    cut = geo.section_at(1.0)
    cap = next(r for r in cut.shell_graph().loop if r.kind == "cap")
    field = StackThickness(geo)
    ends = field(cap.points[[0, -1]])
    mid = float(field(cap.points[len(cap.points) // 2].reshape(1, 3))[0])
    assert ends[0] == pytest.approx(ends[1], abs=1e-9)
    assert float(ends[0]) < mid
    assert mid == pytest.approx(0.003 + 0.008)


def test_a_windio_stack_sums_the_plies_of_each_wall():
    geo = _blade()
    cut = geo.section_at(1.0)
    structure = {
        "layers": [
            {"name": "shell", "plies": [{"thickness": 0.001}, {"thickness": 0.002}]},
            {"name": "spar_cap_ss", "thickness": 0.008},
            {"name": "spar_cap_ps", "plies": [{"thickness": 0.006}, {"thickness": 0.005}]},
            {"name": "fore", "web": "0", "thickness": {"grid": [0.0, 4.0], "values": [0.004, 0.006]}},
        ],
    }
    field = StackThickness.from_windio(geo, structure)
    suction = _point_on(cut, "spar_cap_suction_surface")
    pressure = _point_on(cut, "spar_cap_pressure_surface")
    web = _point_on(cut, "shear_web_0")
    assert field(np.asarray([suction])) == pytest.approx([0.003 + 0.008])
    assert field(np.asarray([pressure])) == pytest.approx([0.003 + 0.011])
    assert field(np.asarray([web])) == pytest.approx([0.0045])
