"""windIO arcs on the blade's own sections; two-foot webs; per-side cap windows."""

from collections import Counter

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.cut import _twist_offset
from implicit_geometry_kernel.blade.topology.region_counts import MIN_CELL_NODES
from implicit_geometry_kernel.blade.windio import nd_arc_to_chord_frac
from implicit_geometry_kernel.blade.windio_arc import (
    arc_point,
    cap_side_windows_fn_from_windio_structure,
    section_arc_loop,
    web_feet_fn_from_windio_structure,
    web_feet_series_from_windio_structure,
)
from implicit_geometry_kernel.mesh import shell_mesh


def _stations(a="2418", b="2415"):
    return [Station(0.0, 1.0, Airfoil.naca4(a), twist=10.0),
            Station(3.0, 0.7, Airfoil.naca4(b), twist=2.0)]


def _oml(**kw):
    return BladeGeometry(_stations(), "0A-LC-0R", 0.003, pivot_fraction=0.3, **kw)


def _table(v0, v1=None):
    return {"grid": [0.0, 1.0], "values": [v0, v0 if v1 is None else v1]}


def _section(geo, z):
    *_, (x, upper, lower) = geo.interpolate(z)
    return np.asarray(x, float), np.asarray(upper, float), np.asarray(lower, float)


def _to_local(geo, z, pts):
    """Blade-frame section points back to chord fractions and unit-chord y."""
    c, t, p, o, _ = geo.interpolate(z)
    tt = np.radians(t)
    rot = np.array([[np.cos(tt), np.sin(tt)], [-np.sin(tt), np.cos(tt)]])
    loc = np.linalg.solve(rot, (np.asarray(pts, float)[:, :2] - np.asarray(o, float)).T).T
    return np.column_stack([loc[:, 0] / c + p, loc[:, 1] / c])


def test_arc_loop_runs_te_suction_to_le_to_te_pressure():
    x, upper, lower = _section(_oml(), 0.0)
    pts, arc, i_le = section_arc_loop(x, upper, lower)
    assert arc[0] == 0.0 and arc[-1] == 1.0 and np.all(np.diff(arc) >= 0.0)
    assert pts[i_le, 0] == pytest.approx(x.min())
    assert pts[0, 1] >= pts[-1, 1]  # starts on the suction side


def test_arc_point_is_at_its_arc_length_and_differs_from_the_fixed_formula():
    x, upper, lower = _section(_oml(), 0.0)
    pts, arc, i_le = section_arc_loop(x, upper, lower)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    total = seg.sum()
    worst = 0.0
    for s in np.linspace(0.05, 0.95, 19):
        px, py, side = arc_point(x, upper, lower, s)
        assert side == ("suction" if s <= arc[i_le] else "pressure")
        # walk the polyline to the point: its arc length is s of the whole
        k = int(np.searchsorted(arc, s, side="right")) - 1
        walked = seg[:k].sum() + np.hypot(px - pts[k, 0], py - pts[k, 1])
        assert walked / total == pytest.approx(s, abs=1e-12)
        worst = max(worst, abs(px - nd_arc_to_chord_frac(s)))
    # an 18 % section: arc length is not chord position (the old reader's mapping)
    assert worst > 0.02


def test_two_foot_webs_lie_on_the_line_through_their_oml_feet():
    feet = lambda z: ((0.20 + 0.01 * z, 0.26), (0.55, 0.50))  # noqa: E731
    geo = BladeGeometry(_stations(), "2D-CN-LC-0R", 0.003, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, pivot_fraction=0.3, web_positions=feet)
    for z in (0.0, 1.5, 3.0):
        c, t, p, o, (x, upper, lower) = geo.interpolate(z)
        webs = geo.section_at(z).shell_graph().webs
        assert len(webs) == 2
        for web, (fs, fp) in zip(webs, feet(z)):
            a = _twist_offset(np.array([[(fs - p) * c, np.interp(fs, x, upper) * c]]), t, o)[0]
            b = _twist_offset(np.array([[(fp - p) * c, np.interp(fp, x, lower) * c]]), t, o)[0]
            d = b - a
            n = np.array([-d[1], d[0]]) / np.linalg.norm(d)
            P = np.asarray(web.points, float)[:, :2]
            assert np.abs((P - a) @ n).max() < 1e-9, (z, web.name)
    assert "web" in set(shell_mesh(geo, None, h=0.1).surface)


def test_each_side_keeps_its_own_cap_window():
    caps = lambda z: {"suction": (0.20, 0.45), "pressure": (0.25, 0.52)}  # noqa: E731
    geo = BladeGeometry(_stations(), "2B-CN-LC-0R", 0.003, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, pivot_fraction=0.3,
                        web_positions=(0.15, 0.6), cap_windows=caps)
    for z in (0.0, 1.5):
        got = {}
        for run in geo.section_at(z).shell_graph().loop:
            if "cap" in run.name:
                f = _to_local(geo, z, run.points)[:, 0]
                side = "suction" if "suction" in run.name else "pressure"
                lo, hi = got.get(side, (np.inf, -np.inf))
                got[side] = (min(lo, f.min()), max(hi, f.max()))
        assert set(got) == {"suction", "pressure"}
        for side, (lo, hi) in caps(z).items():
            # cap runs lie on the skin mid-surface and their edges are cut on the cap
            # mid-plane, so the chord fractions sit within a wall thickness (0.1 % of chord here)
            assert got[side][0] == pytest.approx(lo, abs=2e-3), (z, side)
            assert got[side][1] == pytest.approx(hi, abs=2e-3), (z, side)
    parts = Counter(shell_mesh(geo, None, h=0.1).partition)
    assert any("cap" in p for p in parts)


def _structure():
    return {
        "webs": [
            {"name": "fore", "start_nd_arc": _table(0.36), "end_nd_arc": _table(0.66)},
            {"name": "aft", "start_nd_arc": {"grid": [0.0, 0.3, 0.6], "values": [0.2, 0.2, 0.2]},
             "end_nd_arc": {"grid": [0.0, 0.3, 0.6], "values": [0.83, 0.83, 0.83]},
             "start_nd_span": 0.0, "end_nd_span": 0.6},
        ],
        "layers": [
            {"name": "spar_cap_ss", "start_nd_arc": _table(0.22), "end_nd_arc": _table(0.36)},
            {"name": "spar_cap_ps", "start_nd_arc": _table(0.64), "end_nd_arc": _table(0.80)},
            {"name": "le_reinf", "start_nd_arc": _table(0.45), "end_nd_arc": _table(0.55)},
        ],
    }


def test_windio_web_feet_are_the_oml_points_at_their_arcs():
    geo = _oml()
    z_of_s = np.array([[0.0, 0.0], [1.0, 3.0]])
    series = web_feet_series_from_windio_structure(_structure(), geo, z_of_s)
    assert [s[0] for s in series] == ["fore", "aft"]
    for name, zs, feet, _window in series:
        web = next(w for w in _structure()["webs"] if w["name"] == name)
        for z, (fs, fp) in zip(zs, feet):
            sec = _section(geo, float(z))
            assert fs == pytest.approx(arc_point(*sec, web["start_nd_arc"]["values"][0])[0])
            assert fp == pytest.approx(arc_point(*sec, web["end_nd_arc"]["values"][0])[0])
    fn = web_feet_fn_from_windio_structure(_structure(), geo, z_of_s)
    # between the grid points too: the arc is placed on the section at that z
    for z in (0.4, 1.1, 2.7):
        sec = _section(geo, z)
        fore = min(fn(z), key=lambda p: p[0])
        assert fore == pytest.approx((arc_point(*sec, 0.36)[0], arc_point(*sec, 0.66)[0]), abs=1e-12)
    assert len(fn(1.0)) == 2  # both webs live inboard
    assert len(fn(2.5)) == 1  # the aft web ends at 60 % span
    pairs = fn(1.0)
    assert pairs[0][0] + pairs[0][1] < pairs[1][0] + pairs[1][1]  # LE to TE


def test_windio_cap_windows_per_side_leave_out_a_band_round_the_le():
    geo = _oml()
    fn = cap_side_windows_fn_from_windio_structure(_structure(), geo, np.array([[0.0, 0.0], [1.0, 3.0]]))
    for z in (0.0, 3.0):  # grid points: between them the windows are linear in z
        w = fn(z)
        assert set(w) == {"suction", "pressure"}
        sec = _section(geo, z)
        assert w["suction"] == pytest.approx(
            tuple(sorted(arc_point(*sec, s)[0] for s in (0.22, 0.36))))
        assert w["pressure"] == pytest.approx(
            tuple(sorted(arc_point(*sec, s)[0] for s in (0.64, 0.80))))
    assert w["suction"] != pytest.approx(w["pressure"])


def test_a_single_cell_tip_still_gets_its_tip_closure():
    # no webs at the tip and a coarse h: the loop alone must close on enough nodes
    geo = _oml()
    mesh = shell_mesh(geo, None, h=5.0)
    tips = [e for e, s in zip(mesh.elements, mesh.surface) if s == "tip_closure"]
    assert tips
    ring = {k for e in tips for k in e}
    assert len(ring) >= MIN_CELL_NODES


def test_station_airfoil_passes_through_its_points_at_a_blunt_te_corner():
    # a surface that flares into a blunt trailing edge over the last 0.005 % of chord:
    # sharper than the spacing of the airfoil's own spline samples at that end
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 241)))
    y = 0.6 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    y[-1] -= 4e-4  # the corner
    pts = np.vstack([np.column_stack([x[::-1], y[::-1]]), np.column_stack([x[1:], -y[1:]])])
    foil = Airfoil(pts, "flared")
    # the spline's own most forward point sets x = 0 (2e-11 of chord ahead of the given
    # nose point), so the test is aft of the nose and to 1e-9; the cut corner was 4e-5
    aft = x > 0.5
    assert np.abs(foil.surface("upper", x) - y)[aft].max() < 1e-9
    assert np.abs(foil.surface("lower", x) + y)[aft].max() < 1e-9


def test_round_trailing_end_ends_at_its_end_point():
    # a circle given on a cosine grid: both surfaces arrive vertically at (1, 0)
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 241)))
    y = np.sqrt(x * (1.0 - x))
    pts = np.vstack([np.column_stack([x[::-1], y[::-1]]), np.column_stack([x[1:], -y[1:]])])
    foil = Airfoil(pts, "circle")
    for side, sign in (("upper", 1.0), ("lower", -1.0)):
        got = foil.surface(side, x)
        aft = x > 0.5
        assert np.abs(got - sign * y)[aft].max() < 1e-9  # through the points, to the end
        dense = np.linspace(0.999, 1.0, 2001)  # single-valued and falling to the end
        assert np.all(np.diff(sign * foil.surface(side, dense)) < 0.0)


def test_blunt_base_listed_as_points_keeps_its_corner():
    # a flatback whose surfaces flare, then run down the base to its mid-point: the last
    # segment is the base, not part of the smooth surface
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 241)))
    y = 0.6 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    y += 0.01 * x**8  # half the base height, reached with a flare
    y[-1] = 0.0
    pts = np.vstack([np.column_stack([x[::-1], y[::-1]]), np.column_stack([x[1:], -y[1:]])])
    foil = Airfoil(pts, "flatback")
    aft = x > 0.5
    assert np.abs(foil.surface("upper", x) - y)[aft].max() < 1e-9
    assert np.abs(foil.surface("lower", x) + y)[aft].max() < 1e-9
    # straight down the base between the corner and the end point
    xm = 0.5 * (x[-2] + x[-1])
    assert foil.surface("upper", xm) == pytest.approx(0.5 * y[-2], rel=1e-9)


def test_arc_runs_corner_to_corner_on_a_section_that_lists_its_base():
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 241)))
    y = 0.6 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    y += 0.01 * x**8
    open_upper, open_lower = y.copy(), -y.copy()  # the surfaces alone: an open base
    listed_upper, listed_lower = y.copy(), -y.copy()
    listed_upper[-1] = listed_lower[-1] = 0.0  # each runs on down the base to its mid-point
    for s in (0.1, 0.26, 0.5, 0.68, 0.9):
        px, py, side = arc_point(x, listed_upper, listed_lower, s)
        # the same point as on the surfaces up to their corners (x[-2]), arc for arc
        qx, qy, q_side = arc_point(x[:-1], open_upper[:-1], open_lower[:-1], s)
        assert side == q_side
        assert (px, py) == pytest.approx((qx, qy), abs=1e-12)
