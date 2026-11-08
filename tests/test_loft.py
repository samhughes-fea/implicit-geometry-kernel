"""BladeGeometry: the implicit section geometry lofted between stations and cut at any z."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, section_at


AIRFOIL = Airfoil.naca4("0024")


def _prismatic(twist=0.0, offset=(0.0, 0.0)):
    return [Station(0.0, 1.0, AIRFOIL, twist=twist, pitch_axis=0.3, offset=offset),
            Station(4.0, 1.0, AIRFOIL, twist=twist, pitch_axis=0.3, offset=offset)]


def _cut(stations=None, z=1.0, layout="2D-F-LC-0R", **kw):
    options = dict(
        pivot_fraction=0.3,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48),
    )
    options.update(kw)
    return BladeGeometry(stations or _prismatic(), layout, 0.003, **options).section_at(z)


@pytest.fixture(scope="module")
def cut():
    return _cut()


def test_walls_labels_and_thicknesses(cut):
    kinds = {w.kind for w in cut.walls}
    assert kinds == {"skin", "web", "cap"}
    assert {w.label for w in cut.walls_of("web")} == {"shear_web_0", "shear_web_1"}
    assert {w.label for w in cut.walls_of("cap")} == {"spar_cap_suction_surface", "spar_cap_pressure_surface"}
    # thin walls are below section_geometry's manufacturing floors and still come back as given
    t = {w.kind: w.thickness for w in cut.walls}
    assert t == pytest.approx({"skin": 0.003, "web": 0.002, "cap": 0.008})
    # webs at their chord fractions, measured from the pitch axis
    xs = sorted(float(np.mean(w.xy[:, 0])) for w in cut.walls_of("web"))
    assert xs == pytest.approx([0.18 - 0.3, 0.48 - 0.3], abs=1e-9)
    # caps bonded inside the skin between the webs
    for cap in cut.walls_of("cap"):
        assert cap.bonds and all(host.startswith("skin_") for host in cap.bonds)
        assert cap.xy[:, 0].min() == pytest.approx(-0.12) and cap.xy[:, 0].max() == pytest.approx(0.18)
    assert cut.layout == "2D-F-LC-0R" and cut.dropped_webs == ()


def test_cell_areas(cut):
    assert len(cut.cells) == 3
    # the wall midline encloses a little less than the NACA 0024 outline (0.685 t c^2)
    total = sum(c.area for c in cut.cells)
    assert 0.95 * 0.685 * 0.24 < total < 0.685 * 0.24
    for c in cut.cells:
        x, y = c.points[:, 0], c.points[:, 1]
        assert abs(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)) == pytest.approx(c.area, rel=1e-6)


def test_twist_and_offset_place_the_untwisted_section():
    cut = _cut(layout="2D-CN-LC-0R")
    placed = _cut(_prismatic(twist=15.0, offset=(0.1, -0.05)), layout="2D-CN-LC-0R")
    for wall in cut.walls:
        assert np.allclose(placed.to_blade(wall.xy[:, 0], wall.xy[:, 1]), placed.wall(wall.label).xy, atol=1e-9)
    # nose up: the leading edge moves to +y
    le = placed.to_blade([-0.3], [0.0])[0]
    assert le[1] > -0.05


@pytest.mark.parametrize("layout, tilt", [("2D-F-LC-0R", 0.0), ("2D-CN-LC-0R", 12.0)])
def test_flapwise_webs_stay_normal_to_the_rotor_plane(layout, tilt):
    twisted = _cut(_prismatic(twist=12.0), layout=layout)
    for web in twisted.walls_of("web"):
        d = web.xy[-1] - web.xy[0]
        assert abs(np.degrees(np.arctan2(abs(d[0]), abs(d[1])))) == pytest.approx(tilt, abs=1e-6)
    assert len(twisted.cells) == 3


def test_ellipse_to_naca_morph_follows_smoothstep_blend():
    ellipse, naca = Airfoil.ellipse(1.0), Airfoil.naca4("0012")
    stations = [
        Station(0.0, 2.0, ellipse, pitch_axis=0.5),
        Station(10.0, 1.0, naca, pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, geo.n_airfoil_points)))
    *_, (_x, upper, lower) = geo.interpolate(2.5)
    w = 0.25
    s = w * w * (3.0 - 2.0 * w)
    blended = Airfoil.blend(ellipse, naca, s)
    lerp_u = 0.75 * ellipse.surface("upper", x) + 0.25 * naca.surface("upper", x)
    d_morph = float(np.max(np.abs(upper - blended.surface("upper", x))))
    d_lerp = float(np.max(np.abs(upper - lerp_u)))
    assert d_morph < d_lerp
    assert d_morph < 0.05


def test_naca_to_naca_hits_station_thickness():
    a, b = Airfoil.naca4("0024"), Airfoil.naca4("0012")
    stations = [Station(0.0, 2.0, a, pitch_axis=0.3), Station(10.0, 1.0, b, pitch_axis=0.3)]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    *_, (_x, u0, l0) = geo.interpolate(0.0)
    *_, (_x, u1, l1) = geo.interpolate(10.0)
    assert float(np.max(u0 - l0)) == pytest.approx(a.max_thickness, rel=0.05)
    assert float(np.max(u1 - l1)) == pytest.approx(b.max_thickness, rel=0.05)


def test_interpolate_keeps_sharp_te_a_point():
    stations = [
        Station(0.0, 2.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
        Station(10.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    *_, (_x, upper, lower) = geo.interpolate(5.0)
    assert abs(float(upper[-1]) - float(lower[-1])) < 1e-6


def test_interpolate_preserves_blunt_te_gap():
    blunt = Airfoil.naca4("0018").with_te_thickness(0.04)
    stations = [
        Station(0.0, 2.0, blunt, pitch_axis=0.3),
        Station(10.0, 1.0, blunt, pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    *_, (_x, upper, lower) = geo.interpolate(5.0)
    assert abs(float(upper[-1]) - float(lower[-1])) == pytest.approx(0.04, rel=0.15)


def test_taper_interpolates_chord_twist_and_offset():
    stations = [Station(0.0, 2.0, AIRFOIL, pitch_axis=0.3),
                Station(10.0, 1.0, AIRFOIL, twist=10.0, pitch_axis=0.3, offset=(0.0, 0.2))]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    mid = geo.section_at(5.0)
    assert (mid.chord, mid.twist, *mid.offset) == pytest.approx((1.5, 5.0, 0.0, 0.1))
    root, tip = geo.section_at(0.0), geo.section_at(10.0)
    total = [sum(c.area for c in s.cells) for s in (root, mid, tip)]
    assert total[0] > total[1] > total[2]
    assert total[0] / total[2] == pytest.approx(4.0, rel=0.05)  # area scales with chord^2
    with pytest.raises(ValueError):
        geo.section_at(11.0)


def test_fixed_web_offsets_drop_a_web_that_no_longer_fits():
    stations = [Station(0.0, 2.0, AIRFOIL, pitch_axis=0.3), Station(10.0, 1.0, AIRFOIL, pitch_axis=0.3)]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
                        web_offsets=(0.4, 1.2), pivot_fraction=0.3)
    root, tip = geo.section_at(0.0), geo.section_at(10.0)
    assert len(root.walls_of("web")) == 2 and root.dropped_webs == ()
    assert len(tip.walls_of("web")) == 1 and len(tip.dropped_webs) == 1
    assert "1.2 m" in tip.dropped_webs[0]
    xs = [float(np.mean(w.xy[:, 0])) for w in root.walls_of("web")]
    assert sorted(xs) == pytest.approx([0.4 - 0.6, 1.2 - 0.6], abs=1e-9)  # pitch axis 0.3 x 2 m


def test_thickness_as_a_function_of_span():
    stations = _prismatic()
    cut = section_at(stations, 3.0, "2D-F-LC-0R", lambda z: 0.002 + 0.001 * z, pivot_fraction=0.3,
                     spar_cap_thickness=0.008, shear_web_thickness=0.004)
    assert [w.thickness for w in cut.walls_of("skin")] == pytest.approx([0.005] * 6)


def test_bad_options():
    with pytest.raises(ValueError):
        BladeGeometry(_prismatic(), "2D-F-LC-0R", 0.003, web_positions=(0.2, 0.5), web_offsets=(0.2, 0.5), pivot_fraction=0.3)
    with pytest.raises(ValueError):
        BladeGeometry(_prismatic(), "2D-F-LC-0R", 0.003, web_offsets=(0.2,), pivot_fraction=0.3)


def test_mid_surface_graph(cut):
    mid = cut.mid_surface()
    names = {s.name for s in mid.segments}
    assert {"spar_cap_suction_surface", "shear_web_0", "skin_pressure_surface:trailing_edge"} <= names
    capped = [s for s in mid.segments if s.name == "spar_cap_suction_surface"]
    # the cap's sub-millimetre reach past the webs onto the next skins is left to the skins
    assert {s.layers for s in capped} == {("spar_cap_suction_surface", "skin_suction_surface:bay_0")}
    # oriented runs: every point is left by at most one segment of each wall
    starts = [(s.start, s.name) for s in mid.segments]
    assert len(starts) == len(set(starts))


def test_oml_polyline_closed_and_contains_skin():
    from implicit_geometry_kernel.section.boundary import LineString, OPS, Polygon

    geo = BladeGeometry(_prismatic(), "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008, shear_web_thickness=0.002,
                        web_positions=(0.18, 0.48), pivot_fraction=0.3)
    z = 1.0
    oml = geo.oml_polyline(z)
    assert oml.shape[1] == 3
    assert np.allclose(oml[:, 2], z)
    assert np.allclose(oml[0], oml[-1])
    cut = geo.section_at(z)
    assert (cut.chord, cut.twist, cut.pitch_axis) == pytest.approx(geo.interpolate(z)[:3])
    hull = OPS.buffer(Polygon(oml[:, :2]), 1e-6, join_style=1)
    for wall in cut.walls_of("skin"):
        assert hull.contains(LineString(wall.xy))


def test_span_tables_override_station_scalars_not_airfoils():
    stations = [Station(0.0, 2.0, AIRFOIL, pitch_axis=0.3),
                Station(10.0, 1.0, AIRFOIL, twist=10.0, pitch_axis=0.3, offset=(0.0, 0.2))]
    z = np.array([0.0, 5.0, 10.0])
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        chord=(z, np.array([2.0, 2.0, 1.0])),
        twist=(z, np.array([0.0, 0.0, 10.0])),
        offset=(z, np.array([[0.0, 0.0], [0.4, 0.0], [0.0, 0.2]])), pivot_fraction=0.3)
    chord, twist, pitch, offset, _ = geo.interpolate(5.0)
    assert (chord, twist, *offset) == pytest.approx((2.0, 0.0, 0.4, 0.0))
    assert pitch == pytest.approx(0.3)
    plain = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    assert plain.interpolate(5.0)[0] == pytest.approx(1.5)
    cut = geo.section_at(5.0)
    assert cut.chord == pytest.approx(2.0)
    assert np.allclose(geo.oml_polyline(5.0)[:, 2], 5.0)


def test_bad_span_table():
    with pytest.raises(ValueError, match="chord"):
        BladeGeometry(_prismatic(), "2D-F-LC-0R", 0.003, chord=([0.0], [1.0]), pivot_fraction=0.3)


def test_pchip_span_hits_knots_and_is_smoother_than_linear():
    stations = [
        Station(0.0, 2.0, AIRFOIL, pitch_axis=0.3),
        Station(5.0, 1.2, AIRFOIL, pitch_axis=0.3),
        Station(10.0, 1.0, AIRFOIL, pitch_axis=0.3),
    ]
    linear = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)
    pchip = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
                          span_kind="pchip", pivot_fraction=0.3)
    for z, c in ((0.0, 2.0), (5.0, 1.2), (10.0, 1.0)):
        assert pchip.interpolate(z)[0] == pytest.approx(c)
    chords = [pchip.interpolate(z)[0] for z in np.linspace(0.0, 10.0, 21)]
    assert all(chords[i] >= chords[i + 1] - 1e-12 for i in range(len(chords) - 1))
    h = 1e-3
    d_lin = (linear.interpolate(5.0 + h)[0] - linear.interpolate(5.0 - h)[0]) / (2 * h)
    d_pchip = (pchip.interpolate(5.0 + h)[0] - pchip.interpolate(5.0 - h)[0]) / (2 * h)
    left_lin = (linear.interpolate(5.0)[0] - linear.interpolate(5.0 - h)[0]) / h
    right_lin = (linear.interpolate(5.0 + h)[0] - linear.interpolate(5.0)[0]) / h
    left_p = (pchip.interpolate(5.0)[0] - pchip.interpolate(5.0 - h)[0]) / h
    right_p = (pchip.interpolate(5.0 + h)[0] - pchip.interpolate(5.0)[0]) / h
    assert abs(right_lin - left_lin) > abs(right_p - left_p)
    assert d_pchip != pytest.approx(d_lin, abs=1e-6)


def test_pivot_fraction_is_s_frame_origin():
    """BladeGeometry.pivot_fraction is the XY stacking point and the S origin."""
    stations = _prismatic(twist=0.0)
    geo = BladeGeometry(stations, "0A-LC-0R", 0.003, pivot_fraction=0.40)
    assert geo.interpolate(1.0)[2] == pytest.approx(0.40)
    oml = geo.oml_polyline(1.0)
    assert float(oml[:, 0].min()) == pytest.approx(-0.40, abs=2e-3)
    skin_x = np.concatenate([w.xy[:, 0] for w in geo.section_at(1.0).walls_of("skin")])
    assert float(skin_x.min()) == pytest.approx(float(oml[:, 0].min()), abs=2e-3)


def test_pivot_fraction_moves_the_twist_axis():
    stations = _prismatic(twist=20.0)
    kw = dict(spar_cap_thickness=0.008, shear_web_thickness=0.002)
    a = BladeGeometry(stations, "2D-F-LC-0R", 0.003, pivot_fraction=0.25, **kw)
    b = BladeGeometry(stations, "2D-F-LC-0R", 0.003, pivot_fraction=0.333, **kw)
    assert not np.allclose(a.oml_polyline(2.0), b.oml_polyline(2.0), atol=1e-6)


def test_pivot_fraction_conflicts_with_pitch_axis_table():
    with pytest.raises(ValueError, match="pitch_axis or pivot_fraction"):
        BladeGeometry(
            _prismatic(),
            "0A-LC-0R",
            0.003,
            pitch_axis=lambda z: 0.3,
            pivot_fraction=0.25,
        )


def test_missing_pivot_fraction_raises():
    with pytest.raises(TypeError, match="pivot_fraction is required"):
        BladeGeometry(_prismatic(), "0A-LC-0R", 0.003)
