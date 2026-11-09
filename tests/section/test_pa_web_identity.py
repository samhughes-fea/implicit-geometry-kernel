"""PA freeze identity: omit/seed by intercept, keep authored web indices."""

from __future__ import annotations

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.concepts.web_anchoring import pa_web_x_from_pitch
from implicit_geometry_kernel.section import SectionGeometry, resolve_system_type
from implicit_geometry_kernel.section.spec.layout_resolution import resolve_station_layout
from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF


def _stations_grow():
    foil = Airfoil.naca4("0012")
    return [
        Station(0.0, 1.0, foil, twist=0.0, pitch_axis=0.5),
        Station(5.0, 2.0, foil, twist=0.0, pitch_axis=0.5),
        Station(10.0, 0.6, foil, twist=0.0, pitch_axis=0.5),
    ]


def _web_names(cut) -> list[str]:
    return [str(w.name) for w in cut.shell_graph().webs]


def test_pa_four_web_freeze_matches_helper():
    st = _stations_grow()
    geo = BladeGeometry(
        st, "4D-F-PA-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.5,
    )
    expected = pa_web_x_from_pitch(st, lambda z: 0.5, 4)
    assert geo._web_x_from_pitch == pytest.approx(expected)
    assert expected[-1] > expected[0]


def test_pa_five_web_freeze_matches_helper():
    st = _stations_grow()
    geo = BladeGeometry(
        st, "5D-F-PA-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.5,
    )
    expected = pa_web_x_from_pitch(st, lambda z: 0.5, 5)
    assert geo._web_x_from_pitch == pytest.approx(expected)
    assert len(expected) == 5
    dx = expected[1] - expected[0]
    assert expected[4] == pytest.approx(expected[2] + 2 * dx)


def test_lc_uniform_fracs_unchanged():
    spec = resolve_system_type("4D-F-LC-0R")
    assert spec.web_chord_fracs == pytest.approx((0.2, 0.4, 0.6, 0.8))


def test_pa_three_web_freeze_stays_in_root_band():
    st = _stations_grow()
    geo = BladeGeometry(
        st, "3D-F-PA-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.5,
    )
    expected = pa_web_x_from_pitch(st, lambda z: 0.5, 3)
    assert geo._web_x_from_pitch == pytest.approx(expected)
    assert all(abs(x) <= 0.5 + 1e-9 for x in expected)


def test_two_web_pa_straddles_pitch_axis():
    st = _stations_grow()
    geo = BladeGeometry(
        st, "2D-F-PA-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.5,
    )
    xs = geo._web_x_from_pitch
    assert xs == pytest.approx(pa_web_x_from_pitch(st, lambda z: 0.5, 2))
    assert xs[0] < 0.0 < xs[1]


def test_pa_onset_keeps_shear_web_3():
    geo = BladeGeometry(
        _stations_grow(), "4D-F-PA-0R", 0.004, spar_cap_thickness=0.008,
        shear_web_thickness=0.004, pivot_fraction=0.5,
    )
    root = _web_names(geo.section_at(0.0))
    mid = _web_names(geo.section_at(5.0))
    tip = _web_names(geo.section_at(10.0))
    assert "shear_web_3" not in root
    assert "shear_web_3" in mid
    assert "shear_web_3" not in tip
    for name in ("shear_web_0", "shear_web_1", "shear_web_2"):
        if name in root:
            assert name in mid


def test_pa_drop_keeps_higher_indices():
    foil = AirfoilSDF.from_naca("0012", n_points=121, chord=1.0)
    xs = (0.15, 0.45, 1.20)
    res = resolve_station_layout(
        resolve_system_type("3D-F-PA-0R"), foil, skin_thickness=0.004, web_x_positions_m=xs,
    )
    assert 2 in {d.index for d in res.dropped}
    assert tuple(res.web_indices) == (0, 1)
    runner = SectionGeometry(
        layout=resolve_system_type("3D-F-PA-0R"),
        airfoil=foil,
        skin_thickness=0.004,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.004,
        web_offsets=xs,
        outputs=("midlines",),
    ).execute()
    labels = [s.label for s in (runner.get_results().midline_strips or ()) if "web" in str(s.label)]
    assert "shear_web_0" in labels and "shear_web_1" in labels
    assert "shear_web_2" not in labels


def test_le_te_drop_keeps_leading_indices():
    st = [
        Station(0.0, 2.0, Airfoil.naca4("0012"), twist=0.0, pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0012"), twist=0.0, pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        st, "3D-F-LE-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.3,
    )
    tip = _web_names(geo.section_at(10.0))
    assert "shear_web_0" in tip
    assert all(n.startswith("shear_web_") for n in tip)
    nums = [int(n.rsplit("_", 1)[1]) for n in tip]
    assert nums == sorted(nums)
    assert nums[0] == 0
    assert "shear_web_2" not in tip


def _mean_blade_x(cut, name: str) -> float:
    run = next(w for w in cut.shell_graph().webs if str(w.name) == name)
    return float(np.mean(np.asarray(run.points)[:, 0]))


def test_pa_named_web_keeps_blade_x():
    geo = BladeGeometry(
        _stations_grow(), "4D-F-PA-0R", 0.004, spar_cap_thickness=0.008,
        shear_web_thickness=0.004, pivot_fraction=0.5,
    )
    root = geo.section_at(0.0)
    mid = geo.section_at(5.0)
    for name in ("shear_web_0", "shear_web_1", "shear_web_2"):
        if name not in _web_names(root):
            continue
        assert name in _web_names(mid)
        assert _mean_blade_x(mid, name) == pytest.approx(_mean_blade_x(root, name), abs=0.08)


def test_flapwise_retry_drops_failing_index_not_te():
    from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
    from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
    from implicit_geometry_kernel.section.spec.layout_resolution import StationLayoutResolution

    calls = []

    class _Host:
        def _assemble(self, res, live_y=True):
            calls.append(tuple(res.web_indices))
            return "ok"

    nominal = SectionLayout.from_catalogue_key("3D-F-PA-0R")
    res = StationLayoutResolution(
        nominal, (0.2, 0.5, 0.8), (), nominal, web_indices=(0, 1, 2),
    )
    exc = InfeasibleSectionGeometry("re-anchor", web_index=0, web_x=0.2)
    out = SectionGeometry._assemble_dropping_webs(_Host(), res, exc)
    assert out == "ok"
    assert calls[0] == (1, 2)
