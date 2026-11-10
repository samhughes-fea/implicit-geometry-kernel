"""In-plane size field ``h`` for CSG conformal skin/web nodes."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology import (
    _template_from_cut,
    place_station_keys,
    sample_region_polylines,
)
from implicit_geometry_kernel.blade.topology.sizing import (
    H_MAX_OVER_C,
    H_MIN_OVER_T,
    resample_open_h,
    size_field_h,
    station_live_thickness,
)


def test_size_field_h_zero_kappa_is_h_max_not_wall_gauge():
    c, t = 2.0, 0.008
    expect = max(t * H_MIN_OVER_T, H_MAX_OVER_C * c)
    assert float(size_field_h(0.0, c, t)) == pytest.approx(expect)
    assert float(size_field_h(np.zeros(4), c, t)[0]) == pytest.approx(expect)
    assert expect == pytest.approx(H_MAX_OVER_C * c)


def test_resample_open_h_flat_count_follows_h_max():
    c, t = 1.0, 0.004
    pts = np.column_stack([np.linspace(0.0, c, 41), np.zeros(41), np.zeros(41)])
    out = resample_open_h(pts, chord=c, t=t)
    n_expect = int(np.ceil(c / (H_MAX_OVER_C * c) - 1e-12)) + 1
    assert len(out) == n_expect
    assert len(out) < 40


def test_resample_open_h_keeps_endpoints():
    pts = np.column_stack([np.linspace(0.0, 1.0, 21), np.zeros(21), np.zeros(21)])
    out = resample_open_h(pts, chord=1.0, t=0.004)
    assert np.allclose(out[0], pts[0])
    assert np.allclose(out[-1], pts[-1])
    assert len(out) >= 2


def test_size_field_sampler_pins_web_feet_and_varies_n():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.004,
        pivot_fraction=0.3, spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )
    cut = geo.section_at(5.0)
    sg = cut.shell_graph()
    tmpl = _template_from_cut(geo.region_graph, cut, 5.0, 8)
    keys = place_station_keys(cut, sg, geo.region_graph, 5.0)
    fixed = sample_region_polylines(cut, sg, tmpl, geo.region_graph)
    sized = sample_region_polylines(cut, sg, tmpl, geo.region_graph, size_field=True)
    assert sized
    assert station_live_thickness(cut) > 0.0
    for web in tmpl.webs:
        assert web.name in sized
        pts = np.asarray(sized[web.name], float)
        su = np.asarray(keys[(f"web_{web.index}", "suction")], float)
        pr = np.asarray(keys[(f"web_{web.index}", "pressure")], float)
        d0 = min(float(np.linalg.norm(pts[0] - su)), float(np.linalg.norm(pts[0] - pr)))
        d1 = min(float(np.linalg.norm(pts[-1] - su)), float(np.linalg.norm(pts[-1] - pr)))
        assert d0 == pytest.approx(0.0, abs=1e-12)
        assert d1 == pytest.approx(0.0, abs=1e-12)
    sized_n = {name: len(pts) for name, pts in sized.items()}
    fixed_n = {name: len(pts) for name, pts in fixed.items()}
    assert len(set(sized_n.values())) > 1
    assert sized_n != fixed_n
