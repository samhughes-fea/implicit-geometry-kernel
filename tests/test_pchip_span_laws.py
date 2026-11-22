"""PCHIP span laws and the station-density check (decision 0014)."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.structure import (
    PchipLaw,
    TableLaw,
    UnderSampledSpanLaw,
    check_pchip_sampling,
)

Z = np.array([0.0, 1.0, 2.5, 4.0, 6.0, 9.0])
V = np.array([1.0, 1.3, 1.2, 0.9, 0.5, 0.3])


def _chord(z):
    return 1.0 + 0.3 * np.sin(np.asarray(z) / 3.0)


def test_pchip_matches_the_fritsch_carlson_reference():
    interp = pytest.importorskip("scipy.interpolate")
    law = PchipLaw(Z, V, sampling="flag")
    ref = interp.PchipInterpolator(Z, V)
    zs = np.linspace(Z[0], Z[-1], 401)
    assert np.allclose(law(zs), ref(zs), atol=1e-14)
    assert np.allclose(law.derivative(zs), ref(zs, 1), atol=1e-13)
    assert np.allclose(law.derivative(zs, 2), ref(zs, 2), atol=1e-12)


def test_pchip_passes_through_knots_and_holds_the_ends():
    law = PchipLaw(Z, V, sampling="flag")
    assert np.allclose(law(Z), V)
    assert law.at(-1.0) == pytest.approx(V[0])
    assert law.at(10.0) == pytest.approx(V[-1])
    assert law.derivative([-1.0, 10.0]).tolist() == [0.0, 0.0]


def test_pchip_never_overshoots_and_keeps_plateaus_flat():
    z = [0.0, 1.0, 2.0, 3.0, 4.0]
    v = [0.0, 0.0, 1.0, 1.0, 1.0]
    law = PchipLaw(z, v, sampling="flag")
    zs = np.linspace(0.0, 4.0, 801)
    out = law(zs)
    assert out.min() >= 0.0 and out.max() <= 1.0
    assert np.all(np.diff(out) >= -1e-15)
    assert np.all(out[zs <= 1.0] == 0.0)
    assert np.allclose(out[zs >= 2.0], 1.0)


def test_pchip_is_c1_at_interior_knots():
    zz = np.linspace(0.0, 9.0, 41)
    law = PchipLaw(zz, _chord(zz))
    eps = 1e-9
    for k in zz[1:-1]:
        left = law.derivative([k - eps])[0]
        right = law.derivative([k + eps])[0]
        assert left == pytest.approx(right, abs=1e-6)


def test_dense_smooth_knots_pass_the_leave_one_out_check():
    zz = np.linspace(0.0, 9.0, 41)
    report = check_pchip_sampling(zz, _chord(zz), name="chord")
    assert report.ok
    assert report.loo_error[0] == report.loo_error[-1] == 0.0
    assert report.bay_gap.shape == (40,)


def test_sparse_knots_are_refused_with_the_quantity_and_station_named():
    zz = np.linspace(0.0, 9.0, 4)
    with pytest.raises(UnderSampledSpanLaw) as info:
        PchipLaw(zz, _chord(zz), name="chord")
    assert "chord" in str(info.value)
    assert not info.value.report.ok
    worst = info.value.report.worst
    assert f"z = {zz[worst]:g}" in str(info.value)


def test_flag_mode_keeps_the_law_and_the_failed_report():
    zz = np.linspace(0.0, 9.0, 4)
    law = PchipLaw(zz, _chord(zz), name="chord", sampling="flag")
    assert not law.sampling.ok
    assert np.allclose(law(zz), _chord(zz))


def test_linear_laws_are_not_checked():
    zz = np.linspace(0.0, 9.0, 4)
    assert np.allclose(TableLaw(zz, _chord(zz))(zz), _chord(zz))


def test_two_and_three_knot_laws():
    two = PchipLaw([0.0, 2.0], [1.0, 3.0])
    assert two.at(1.0) == pytest.approx(2.0)
    assert check_pchip_sampling([0.0, 1.0, 2.0], [0.0, 1.0, 2.0]).ok


def _stations(n):
    foil = Airfoil.naca4("0018")
    zs = np.linspace(0.0, 9.0, n)
    return [
        Station(float(z), float(_chord(z)), foil, twist=float(10.0 - z), pitch_axis=0.3)
        for z in zs
    ]


def _blade(n, **kw):
    return BladeGeometry(
        _stations(n), "2D-F-LC-0R", 0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
        pivot_fraction=0.3, **kw,
    )


def test_linear_stays_the_default():
    geo = _blade(4)
    assert geo.span_interp == "linear"
    assert geo.span_sampling_reports == []
    assert geo.validate() == []


def test_pchip_blade_matches_stations_and_curves_between_them():
    lin = _blade(31)
    smooth = _blade(31, span_interp="pchip")
    for st in smooth.stations:
        chord, twist, *_ = smooth.interpolate(st.z)
        assert chord == pytest.approx(st.chord)
        assert twist == pytest.approx(st.twist)
    mid = 0.5 * (smooth.stations[3].z + smooth.stations[4].z)
    exact = float(_chord(mid))
    err_pchip = abs(smooth.interpolate(mid)[0] - exact)
    err_linear = abs(lin.interpolate(mid)[0] - exact)
    assert err_pchip < err_linear
    assert all(report.ok for report in smooth.span_sampling_reports)
    names = {report.name for report in smooth.span_sampling_reports}
    assert {"chord", "twist"} <= names


def test_pchip_blade_with_sparse_stations_is_refused():
    with pytest.raises(UnderSampledSpanLaw, match="chord"):
        _blade(4, span_interp="pchip")


def test_flagged_pchip_blade_lists_the_failure_in_validate():
    geo = _blade(4, span_interp="pchip", span_sampling="flag")
    issues = [i for i in geo.validate() if "PCHIP" in i]
    assert issues and "chord" in issues[0]


def test_pchip_span_table_and_options_are_checked():
    zt = np.linspace(0.0, 9.0, 31)
    geo = _blade(2, span_interp="pchip", chord=(zt, _chord(zt)))
    assert geo.interpolate(4.5)[0] == pytest.approx(float(_chord(4.5)), abs=1e-4)
    with pytest.raises(ValueError, match="span_interp"):
        _blade(2, span_interp="cubic")
    with pytest.raises(ValueError, match="span_sampling"):
        _blade(2, span_interp="pchip", span_sampling="warn")
