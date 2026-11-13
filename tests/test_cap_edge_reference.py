"""A cap window is a chord fraction on the mid-plane, or on the outer mould line."""

import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station


def _blade(reference="midplane"):
    foil = Airfoil.ellipse(1.0)
    stations = [
        Station(0.0, 2.0, foil, pitch_axis=0.5),
        Station(4.0, 2.0, foil, pitch_axis=0.5),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.04,
        spar_cap_thickness=0.06,
        shear_web_thickness=0.01,
        web_positions=(0.2, 0.8),
        cap_windows=((0.35, 0.65),),
        pivot_fraction=0.5,
        cap_edge_reference=reference,
    )


def _suction_fracs(cut):
    cap = next(w for w in cut.walls_of("cap") if "suction" in w.label)
    chord, pivot = cut.chord, cut.pitch_axis
    return sorted(float(x) / chord + pivot for x in (cap.xy[:, 0].min(), cap.xy[:, 0].max()))


def test_midplane_is_the_default_and_keeps_the_window():
    geo = _blade()
    assert geo.cap_edge_reference == "midplane"
    lo, hi = _suction_fracs(geo.section_at(1.0))
    assert lo == pytest.approx(0.35, abs=2e-3)
    assert hi == pytest.approx(0.65, abs=2e-3)


def test_outer_reference_moves_the_midplane_edge_onto_the_normal():
    mid = _suction_fracs(_blade("midplane").section_at(1.0))
    geo = _blade("outer")
    landed = geo._outer_fraction_to_midplane(1.0, 0.35, "suction")
    lo, hi = _suction_fracs(geo.section_at(1.0))
    assert abs(landed - 0.35) > 1e-3
    assert lo == pytest.approx(landed, abs=2e-3)
    assert abs(lo - mid[0]) > 1e-3


def test_unknown_cap_edge_reference_is_rejected():
    with pytest.raises(ValueError):
        _blade("node")
