"""Mesh digest and provenance."""

import json
from dataclasses import replace

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, __version__
from implicit_geometry_kernel.mesh import shell_mesh


@pytest.fixture(scope="module")
def mesh():
    foil = Airfoil.naca4("0018")
    stations = [Station(0.0, 1.0, foil, pitch_axis=0.3), Station(4.0, 0.7, foil, pitch_axis=0.3)]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, pivot_fraction=0.3)
    return shell_mesh(geo, h=0.25)


def test_digest_is_stable_and_content_sensitive(mesh):
    d = mesh.digest()
    assert len(d) == 64 and d == mesh.digest()
    assert replace(mesh, nodes=np.round(np.array(mesh.nodes), 12)).digest() == d  # hashed after rounding
    moved = np.array(mesh.nodes)
    moved[0, 0] += 1e-6
    assert replace(mesh, nodes=moved).digest() != d
    th = np.array(mesh.thickness)
    th[0] *= 2.0
    assert replace(mesh, thickness=th).digest() != d


def test_provenance_is_json_ready(mesh):
    prov = mesh.provenance()
    assert prov["kernel_version"] == __version__
    assert prov["mesher"] == "partition"  # the default since 0.9.6
    assert prov["h"] == pytest.approx(0.25)
    assert prov["n_elements"] == len(mesh.elements)
    assert prov["digest"] == mesh.digest()
    json.dumps(prov)
