"""Optional Gmsh wall-strip mesh."""

import importlib.util

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.extract.curves import locate_shell_vertex
from implicit_geometry_kernel.mesh.gmsh_shell import junction_web_errors, mesh_shell


def _geo():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3)


def test_locate_shell_vertex_tracks_te():
    geo = _geo()
    zs = np.linspace(0.0, 10.0, 6)
    pts = locate_shell_vertex(geo, ("TE", "suction"), zs)
    assert len(pts) == 6
    jumps = np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)
    assert float(np.max(jumps)) < 0.5


@pytest.mark.skipif(importlib.util.find_spec("gmsh") is None, reason="gmsh not installed")
def test_mesh_shell_is_all_quad_on_prismatic_2df():
    geo = _geo()
    stats = mesh_shell(geo, np.linspace(0.0, 10.0, 5), n_chord=4)
    assert stats["n_elements"] > 0
    assert stats["n_nodes"] > 0
    assert stats["faces"]
    assert all(len(f) == 4 for f in stats["faces"])
    assert junction_web_errors(geo, np.linspace(0.0, 10.0, 5), n_chord=4) < 0.05
