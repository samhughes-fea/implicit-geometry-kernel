"""StackThickness can take a layup or a windIO ply stack."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology.fields import StackThickness


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
