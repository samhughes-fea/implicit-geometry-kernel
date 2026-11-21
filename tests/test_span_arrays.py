"""Span arrays: one placed airfoil stack, sagitta before any section, skin and web splits."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.loft.span_place import offset_rings, split_at_arc, web_lines
from implicit_geometry_kernel.blade.topology.sizing import EPS_OVER_C
from implicit_geometry_kernel.mesh.shell import _freeze_sagitta, _span_edges, _stacked_shell


def _blade():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=8.0),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3, twist=0.0),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )


def _bowed():
    """Twist that bows the leading edge well past the sagitta budget."""
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=0.0),
        Station(4.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=40.0),
    ]
    return BladeGeometry(stations, "0A-LC-0R", 0.004, pivot_fraction=0.25)


def test_stacked_placement_matches_oml_polyline():
    geo = _blade()
    zs = np.linspace(0.0, 10.0, 5)
    rings = geo.oml_polylines(zs)
    for i, z in enumerate(zs):
        one = geo.oml_polyline(float(z))
        assert rings[i].shape == one.shape
        assert np.allclose(rings[i], one)


def test_span_sagitta_matches_the_scalar_ruling():
    geo = _bowed()
    zs = np.linspace(0.0, 4.0, 3)
    _le, _te, le, te, budget = _span_edges(geo, zs)
    for i in range(len(zs) - 1):
        za, zb = float(zs[i]), float(zs[i + 1])
        mid = 0.5 * (za + zb)
        expect_budget = EPS_OVER_C * geo._span_scalar(geo._chord, [s.chord for s in geo.stations], mid)
        assert le[i] == pytest.approx(float(_scalar_sag(geo, za, zb, "le")))
        assert te[i] == pytest.approx(float(_scalar_sag(geo, za, zb, "te")))
        assert budget[i] == pytest.approx(expect_budget)


def test_sagitta_freeze_matches_the_per_bay_loop():
    geo = _bowed()
    coarse = (0.0, 4.0)
    frozen, _inserted = _freeze_sagitta(geo, coarse, None, 1.0)
    expect = _scalar_freeze(geo, coarse, None, 1.0)
    assert frozen == pytest.approx(expect)
    assert len(frozen) > 2


def test_vertex_normal_offset_is_the_hand_built_midline():
    square = np.array([[
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0],
    ]])
    moved = offset_rings(square, 0.1)
    step = 0.1 / np.sqrt(2.0)
    expect = np.array([
        [step, step, 1.0],
        [1.0 - step, step, 1.0],
        [1.0 - step, 1.0 - step, 1.0],
        [step, 1.0 - step, 1.0],
        [step, step, 1.0],
    ])
    assert np.allclose(moved[0], expect)


def test_web_chord_and_arc_split_match_a_hand_built_ring():
    x = np.array([0.0, 0.5, 1.0])
    upper = np.array([[0.0, 0.2, 0.0], [0.0, 0.2, 0.0]])
    lower = -upper
    chord = np.array([2.0, 2.0])
    pitch = np.zeros(2)
    twist = np.zeros(2)
    offset = np.zeros((2, 2))
    zs = np.array([3.0, 4.0])
    webs = web_lines(x, upper, lower, chord, pitch, twist, offset, zs, np.array([0.5]))
    assert webs.shape == (2, 1, 2, 3)
    assert np.allclose(webs[0, 0, 0], [1.0, 0.4, 3.0])
    assert np.allclose(webs[0, 0, 1], [1.0, -0.4, 3.0])
    assert np.allclose(webs[1, 0, 0], [1.0, 0.4, 4.0])

    square = np.array([[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ]] * 2)
    cuts = split_at_arc(square, np.array([0.0, 0.25, 0.5]))
    assert np.allclose(cuts[:, 0], [0.0, 0.0, 0.0])
    assert np.allclose(cuts[:, 1], [1.0, 0.0, 0.0])
    assert np.allclose(cuts[:, 2], [1.0, 1.0, 0.0])


def test_stacked_shell_splits_every_station_without_a_section():
    geo = _blade()
    zs = (0.0, 5.0, 10.0)
    skin, webs, splits, le, te, le_sag, te_sag, budget = _stacked_shell(geo, zs)
    assert skin.shape[0] == 3 and skin.shape[-1] == 3
    assert np.allclose(skin[:, 0], skin[:, -1])
    assert webs.shape[0] == 3 and webs.shape[2] == 2
    assert splits.shape[0] == 3 and splits.shape[1] == webs.shape[1]
    assert le.shape == (3, 3) and te.shape == (3, 3)
    assert le_sag.shape == te_sag.shape == budget.shape == (2,)
    for i, z in enumerate(zs):
        assert np.allclose(skin[i], _hand_midline(geo, z))


def test_stacked_skin_against_the_section_midline_when_geos_is_present():
    """The mesh keeps the section sampler until this offset matches the midline."""
    geo = _blade()
    try:
        cut = geo.section_at(5.0)
    except OSError:
        pytest.skip("libgeos_c is not available")
    from implicit_geometry_kernel.blade.topology.map_walk import _skin_mid_loop

    mid = np.asarray(_skin_mid_loop(cut.shell_graph(), 5.0), float)
    skin, *_rest = _stacked_shell(geo, (5.0,))
    assert skin.shape[-1] == 3 and len(mid) >= 3
    # Recorded before any switch of the conforming sampler. A vertex-normal
    # offset is not a mitre buffer, so this is a comparison, not a switch.
    _distance = float(np.min(np.linalg.norm(skin[0][:, None, :] - mid[None, :, :], axis=2)))
    assert np.isfinite(_distance)


def _scalar_sag(geo, za: float, zb: float, which: str) -> float:
    from implicit_geometry_kernel.blade.cut import _twist_offset

    mid = 0.5 * (za + zb)

    def point(z: float) -> np.ndarray:
        chord, twist, pitch, offset, (x, upper, lower) = geo.interpolate(z)
        i = int(np.argmin(x) if which == "le" else np.argmax(x))
        y = 0.5 * (float(upper[i]) + float(lower[i]))
        local = np.array([[float(x[i]) * chord - float(pitch) * chord, y * chord]])
        xy = _twist_offset(local, geo._placement_twist(twist), offset)
        return np.array([float(xy[0, 0]), float(xy[0, 1]), z])

    a, m, b = point(za), point(mid), point(zb)
    ab = b - a
    length = float(np.linalg.norm(ab))
    if length < 1e-12:
        return 0.0
    t = float(np.dot(m - a, ab) / (length * length))
    return float(np.linalg.norm(m - (a + t * ab)))


def _scalar_freeze(geo, zs, h, fraction) -> tuple[float, ...]:
    from implicit_geometry_kernel.mesh.shell import _target_edge, _unique_z

    chosen = [float(z) for z in zs]
    for _round in range(8):
        present = {round(z, 9) for z in chosen}
        fresh = []
        for i in range(len(chosen) - 1):
            za, zb = chosen[i], chosen[i + 1]
            le = _scalar_sag(geo, za, zb, "le")
            te = _scalar_sag(geo, za, zb, "te")
            mid_z = 0.5 * (za + zb)
            budget = EPS_OVER_C * geo._span_scalar(
                geo._chord, [s.chord for s in geo.stations], mid_z)
            if max(le, te) <= budget:
                continue
            if zb - za < 0.5 * _target_edge(geo, zb, h, fraction):
                continue
            key = round(mid_z, 9)
            if key in present or any(abs(mid_z - s) <= 1e-6 for s in chosen):
                continue
            fresh.append(mid_z)
            present.add(key)
        if not fresh:
            break
        chosen = list(_unique_z(chosen + fresh))
    return tuple(chosen)


def _hand_midline(geo, z: float) -> np.ndarray:
    """Outer mould at ``z``, moved inward along vertex normals by half the skin."""
    ring = geo.oml_polyline(float(z))[None, :, :]
    return offset_rings(ring, 0.5 * float(geo._skin(float(z))))[0]


def test_chord_normal_placement_matches_the_scalar_map():
    from implicit_geometry_kernel.blade.cut import _twist_offset
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal, to_chord_normal

    unit = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 0.2], [0.0, 0.1]])
    chord = np.array([2.0, 1.0])
    twist = np.array([0.0, 25.0])
    pitch = np.array([0.3, 0.4])
    offset = np.array([[0.1, -0.2], [0.0, 0.05]])
    zs = np.array([0.0, 4.0])
    placed = place_chord_normal(unit, chord, twist, pitch, offset, zs)
    for i in range(2):
        local = np.column_stack([
            unit[:, 0] * chord[i] - pitch[i] * chord[i],
            unit[:, 1] * chord[i],
        ])
        xy = _twist_offset(local, float(twist[i]), offset[i])
        expect = np.column_stack([xy, np.full(len(xy), zs[i])])
        assert np.allclose(placed[i], expect)
        back = to_chord_normal(placed[i, :, :2], chord[i], twist[i], pitch[i], offset[i])
        assert np.allclose(back, unit)


def test_shape_key_follows_airfoil_and_thickness():
    from implicit_geometry_kernel.mesh.shell import _shape_keys

    same = BladeGeometry(
        [
            Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=0.0),
            Station(10.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=20.0),
        ],
        "0A-LC-0R", 0.004, pivot_fraction=0.25,
    )
    keys = _shape_keys(same, (0.0, 5.0, 10.0))
    assert keys[0] == keys[1] == keys[2]

    stepped = BladeGeometry(
        [
            Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25),
            Station(10.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25),
        ],
        "0A-LC-0R", lambda z: 0.004 if float(z) < 5.0 else 0.012, pivot_fraction=0.25,
    )
    thick = _shape_keys(stepped, (0.0, 4.0, 8.0))
    assert thick[0] == thick[1]
    assert thick[0] != thick[2]

    morph = _shape_keys(_blade(), (0.0, 5.0, 10.0))
    assert len(set(morph)) == 3

    similar = BladeGeometry(
        [
            Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=0.0),
            Station(4.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=10.0),
        ],
        "0A-LC-0R", lambda z: 0.01 * (2.0 - 0.25 * float(z)), pivot_fraction=0.25,
    )
    shared = _shape_keys(similar, (0.0, 4.0))
    assert shared[0] == shared[1]


def test_open_resample_matches_one_row_at_a_time():
    from implicit_geometry_kernel.blade.loft.locked_loft import arc_resample_open
    from implicit_geometry_kernel.blade.loft.span_place import resample_open_many

    rng = np.random.default_rng(0)
    stack = rng.normal(size=(4, 7, 3))
    together = resample_open_many(stack, 5)
    for i in range(4):
        assert np.allclose(together[i], arc_resample_open(stack[i], 5))
    bunched = resample_open_many(stack, 6, bunch=True)
    reference = []
    for row in stack:
        seg = np.linalg.norm(np.diff(row, axis=0), axis=1)
        s = np.concatenate([[0.0], np.cumsum(seg)])
        fraction = np.sin(np.linspace(0.0, 0.5 * np.pi, 6)) ** 2
        target = fraction * s[-1]
        out = np.empty((6, 3))
        for j, tj in enumerate(target):
            k = int(np.searchsorted(s, tj, side="right") - 1)
            k = min(max(k, 0), len(seg) - 1)
            span = seg[k] if seg[k] > 1e-15 else 1.0
            w = (tj - s[k]) / span
            out[j] = row[k] + w * (row[k + 1] - row[k])
        reference.append(out)
    assert np.allclose(bunched, np.stack(reference))


def test_grouped_nearest_station_matches_the_scan():
    from implicit_geometry_kernel.mesh.shell import _nearest_at, _nearest_runs, _nearest_stations

    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.2, 0.4, 1.0],
        [0.8, 0.1, 1.0],
    ])
    origins = ["station", "station", "station", "web"]
    queries = np.array([[0.1, 0.0, 0.0], [0.9, 0.2, 1.0]])
    found = _nearest_stations(nodes, origins, queries[:, 2], queries)
    for i, z in enumerate(queries[:, 2]):
        assert found[i] == _nearest_at(nodes, origins, float(z), queries[i])

    class _Run:
        def __init__(self, points):
            self.points = np.asarray(points, float)

    runs = (_Run([[0.0, 0.0], [1.0, 0.0]]), _Run([[0.0, 1.0], [1.0, 1.0]]))
    pts = np.array([[0.2, 0.1], [0.2, 0.9]])
    chosen = _nearest_runs(runs, pts)
    assert chosen[0] is runs[0]
    assert chosen[1] is runs[1]


def test_sibling_section_matches_section_at_when_geos_is_present():
    from implicit_geometry_kernel.mesh.shell import _placed_cuts

    geo = BladeGeometry(
        [
            Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=0.0),
            Station(4.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.25, twist=10.0),
        ],
        "0A-LC-0R", lambda z: 0.01 * (2.0 - 0.25 * float(z)), pivot_fraction=0.25,
    )
    try:
        cuts, _keys = _placed_cuts(geo, (0.0, 4.0))
    except OSError:
        pytest.skip("libgeos_c is not available")
    direct = geo.section_at(4.0)
    placed = cuts[1]
    assert len(placed.walls) == len(direct.walls)
    for wall, other in zip(placed.walls, direct.walls):
        assert wall.label == other.label
        assert np.allclose(wall.points, other.points, atol=1e-6, rtol=1e-5)


def test_cap_passing_a_web_is_a_station_without_a_section():
    from implicit_geometry_kernel.blade.topology.crossings import key_crossing_stations

    geo = BladeGeometry(
        [Station(0.0, 1.0, Airfoil.naca4("2418"), twist=10.0),
         Station(3.0, 0.7, Airfoil.naca4("2415"), twist=2.0)],
        "2B-CN-LC-0R", 0.003, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.3, web_positions=(0.26, 0.7),
        cap_windows=lambda z: (0.20 + 0.04 * z, 0.45 + 0.04 * z),
    )
    extra = key_crossing_stations(geo, [1.0, 2.0])
    assert len(extra) == 1 and extra[0] == pytest.approx(1.5, abs=1e-6)


def test_section_sample_keeps_a_rib_and_the_ends_of_a_web():
    from dataclasses import replace

    from implicit_geometry_kernel.blade.stations import TransverseMember
    from implicit_geometry_kernel.blade.topology.graph import WebMember
    from implicit_geometry_kernel.mesh.shell import _section_sites

    geo = BladeGeometry(
        [Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25),
         Station(10.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25)],
        "0A-LC-0R", 0.004, pivot_fraction=0.25,
    )

    def foot(z: float) -> float:
        span = (float(z) - 2.0) / 6.0
        return 0.25 + 0.2 * span * span

    web = WebMember(0, foot, z0=2.0, z1=8.0, name="web_0")
    geo.region_graph = replace(geo.region_graph, webs=(web,))
    geo.transverse = (TransverseMember("rib_1", 4.0, 0.012),)
    sites = _section_sites(geo)
    for z in (0.0, 2.0, 4.0, 8.0, 10.0):
        assert any(abs(s - z) < 1e-9 for s in sites)
    assert not any(0.0 < s < 2.0 or 8.0 < s < 10.0 for s in sites)
    assert any(2.0 < s < 8.0 and abs(s - 4.0) > 1e-6 for s in sites)


def test_a_straight_section_keeps_its_defining_stations():
    from implicit_geometry_kernel.mesh.shell import _refine_shape_sites

    geo = BladeGeometry(
        [Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25),
         Station(10.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25)],
        "0A-LC-0R", 0.004, pivot_fraction=0.25,
    )
    assert _refine_shape_sites(geo, (0.0, 10.0)) == (0.0, 10.0)


def test_a_bowed_cap_law_adds_section_sites():
    from implicit_geometry_kernel.mesh.shell import _loft_gap, _refine_shape_sites

    geo = BladeGeometry(
        [Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25),
         Station(10.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.25)],
        "0A-LC-0R", 0.004, pivot_fraction=0.25,
        cap_windows=lambda z: (0.2 + 0.15 * (float(z) / 10.0) ** 2, 0.7),
    )
    sites = _refine_shape_sites(geo, (0.0, 10.0))
    assert sites[0] == 0.0 and sites[-1] == 10.0 and len(sites) > 2
    for za, zb in zip(sites, sites[1:]):
        assert _loft_gap(geo, za, zb) <= 0.0


def test_lofted_polyline_is_the_blend_of_the_bay_ends():
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal
    from implicit_geometry_kernel.mesh.shell import _lerp_polylines

    class End:
        def __init__(self, z, chord):
            self.z, self.chord, self.twist, self.pitch_axis, self.offset = z, chord, 0.0, 0.0, (0.0, 0.0)

    left, right, mid = End(0.0, 2.0), End(10.0, 1.0), End(5.0, 1.5)
    unit_a = np.array([[0.0, 0.0], [1.0, 0.2]])
    unit_b = np.array([[0.0, 0.0], [1.0, -0.2]])
    samples_a = {"skin": place_chord_normal(unit_a, left.chord, 0.0, 0.0, (0.0, 0.0), left.z)[0]}
    samples_b = {"skin": place_chord_normal(unit_b, right.chord, 0.0, 0.0, (0.0, 0.0), right.z)[0]}
    blended = _lerp_polylines(samples_a, samples_b, left, right, mid, 0.5)
    expect = place_chord_normal(0.5 * (unit_a + unit_b), mid.chord, 0.0, 0.0, (0.0, 0.0), mid.z)[0]
    assert np.allclose(blended["skin"], expect)


def test_region_loft_keeps_branch_order_when_the_near_end_is_the_other_foot():
    """A short region that slides keeps start-to-start pairing.

    The first point of one station can lie nearer the other station's end
    than its start. The branch still runs start to end on both, and reversing
    it folds the strip.
    """
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal
    from implicit_geometry_kernel.mesh.shell import _lerp_polylines

    class End:
        def __init__(self, z):
            self.z, self.chord, self.twist, self.pitch_axis, self.offset = z, 1.0, 0.0, 0.0, (0.0, 0.0)

    left, right, mid = End(0.0), End(10.0), End(5.0)
    unit_a = np.array([[0.0, 1.0], [0.2, 0.0]])
    unit_b = np.array([[0.15, 0.0], [0.25, 1.0]])
    samples_a = {"cap_pr": place_chord_normal(unit_a, 1.0, 0.0, 0.0, (0.0, 0.0), left.z)[0]}
    samples_b = {"cap_pr": place_chord_normal(unit_b, 1.0, 0.0, 0.0, (0.0, 0.0), right.z)[0]}
    blended = _lerp_polylines(samples_a, samples_b, left, right, mid, 0.5)
    expect = place_chord_normal(0.5 * (unit_a + unit_b), 1.0, 0.0, 0.0, (0.0, 0.0), mid.z)[0]
    assert np.allclose(blended["cap_pr"], expect)
