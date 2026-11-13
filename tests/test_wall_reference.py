"""Shell nodes sit on the mid-surface, or on the outer mould line."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station


def _blade(reference="mid"):
    foil = Airfoil.ellipse(1.0)
    stations = [
        Station(0.0, 2.0, foil, pitch_axis=0.5),
        Station(4.0, 2.0, foil, pitch_axis=0.5),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.02,
        spar_cap_thickness=0.01,
        shear_web_thickness=0.008,
        web_positions=(0.3, 0.7),
        pivot_fraction=0.5,
        wall_reference=reference,
    )


def _span(graph):
    pts = np.vstack([run.points[:, :2] for run in graph.loop if run.kind == "skin"])
    return float(pts[:, 0].max() - pts[:, 0].min()), float(pts[:, 1].max() - pts[:, 1].min())


def test_mid_is_the_default():
    assert _blade().wall_reference == "mid"


def test_oml_reference_moves_the_skin_out_by_half_the_thickness():
    mid_x, mid_y = _span(_blade("mid").section_at(1.0).shell_graph())
    oml_x, oml_y = _span(_blade("oml").section_at(1.0).shell_graph())
    # both sides of a circle step out by half the 0.02 m skin
    assert oml_x - mid_x == pytest.approx(0.02, abs=2e-3)
    assert oml_y - mid_y == pytest.approx(0.02, abs=2e-3)


def test_unknown_wall_reference_is_rejected():
    with pytest.raises(ValueError):
        _blade("inner")
