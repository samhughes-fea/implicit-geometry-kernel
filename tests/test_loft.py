"""BladeGeometry: the implicit section geometry lofted between stations and cut at any z."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, section_at


AIRFOIL = Airfoil.naca4("0024")


def _prismatic(twist=0.0, offset=(0.0, 0.0)):
    return [Station(0.0, 1.0, AIRFOIL, twist=twist, pitch_axis=0.3, offset=offset),
            Station(4.0, 1.0, AIRFOIL, twist=twist, pitch_axis=0.3, offset=offset)]


def _cut(stations=None, z=1.0, layout="2D-F", **kw):
    options = dict(spar_cap_thickness=0.008, shear_web_thickness=0.002, web_positions=(0.18, 0.48))
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
    assert cut.layout == "2D-F" and cut.dropped_webs == ()


def test_cell_areas(cut):
    assert len(cut.cells) == 3
    # the wall midline encloses a little less than the NACA 0024 outline (0.685 t c^2)
    total = sum(c.area for c in cut.cells)
    assert 0.95 * 0.685 * 0.24 < total < 0.685 * 0.24
    for c in cut.cells:
        x, y = c.points[:, 0], c.points[:, 1]
        assert abs(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)) == pytest.approx(c.area, rel=1e-6)


def test_twist_and_offset_place_the_untwisted_section():
    cut = _cut(layout="2D-CN")
    placed = _cut(_prismatic(twist=15.0, offset=(0.1, -0.05)), layout="2D-CN")
    for wall in cut.walls:
        assert np.allclose(placed.to_blade(wall.xy[:, 0], wall.xy[:, 1]), placed.wall(wall.label).xy, atol=1e-9)
    # nose up: the leading edge moves to +y
    le = placed.to_blade([-0.3], [0.0])[0]
    assert le[1] > -0.05


@pytest.mark.parametrize("layout, tilt", [("2D-F", 0.0), ("2D-CN", 12.0)])
def test_flapwise_webs_stay_normal_to_the_rotor_plane(layout, tilt):
    twisted = _cut(_prismatic(twist=12.0), layout=layout)
    for web in twisted.walls_of("web"):
        d = web.xy[-1] - web.xy[0]
        assert abs(np.degrees(np.arctan2(abs(d[0]), abs(d[1])))) == pytest.approx(tilt, abs=1e-6)
    assert len(twisted.cells) == 3


def test_taper_interpolates_chord_twist_and_offset():
    stations = [Station(0.0, 2.0, AIRFOIL, pitch_axis=0.3),
                Station(10.0, 1.0, AIRFOIL, twist=10.0, pitch_axis=0.3, offset=(0.0, 0.2))]
    geo = BladeGeometry(stations, "2D-F", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004)
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
    geo = BladeGeometry(stations, "2D-F", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
                        web_offsets=(0.4, 1.2))
    root, tip = geo.section_at(0.0), geo.section_at(10.0)
    assert len(root.walls_of("web")) == 2 and root.dropped_webs == ()
    assert len(tip.walls_of("web")) == 1 and len(tip.dropped_webs) == 1
    assert "1.2 m" in tip.dropped_webs[0]
    xs = [float(np.mean(w.xy[:, 0])) for w in root.walls_of("web")]
    assert sorted(xs) == pytest.approx([0.4 - 0.6, 1.2 - 0.6], abs=1e-9)  # pitch axis 0.3 x 2 m


def test_thickness_as_a_function_of_span():
    stations = _prismatic()
    cut = section_at(stations, 3.0, "2D-F", lambda z: 0.002 + 0.001 * z, spar_cap_thickness=0.008,
                     shear_web_thickness=0.004)
    assert [w.thickness for w in cut.walls_of("skin")] == pytest.approx([0.005] * 6)


def test_bad_options():
    with pytest.raises(ValueError):
        BladeGeometry(_prismatic(), "2D-F", 0.003, web_positions=(0.2, 0.5), web_offsets=(0.2, 0.5))
    with pytest.raises(ValueError):
        BladeGeometry(_prismatic(), "2D-F", 0.003, web_offsets=(0.2,))


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
