"""Strict invariant gate for shell meshes, and a property fuzzer over generic blades."""
from __future__ import annotations

import dataclasses
import os

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.mesh import shell_mesh
from implicit_geometry_kernel.mesh.shell import (
    PartitionMeshError,
    ShellMeshInvariantError,
    assert_shell_mesh,
    check_shell_mesh,
)


def _blade(code="2D-F-LC-2R", *, digits="0018", root=1.0, taper=0.7, length=4.0, twist=5.0, pivot=0.3):
    foil = Airfoil.naca4(digits)
    stations = [
        Station(0.0, root, foil, twist=twist, pitch_axis=pivot),
        Station(length, root * taper, foil, twist=0.0, pitch_axis=pivot),
    ]
    return BladeGeometry(stations, code, 0.004, spar_cap_thickness=0.008,
                         shear_web_thickness=0.004, pivot_fraction=pivot)


@pytest.fixture(scope="module", params=["ruled", "partition"])
def mesh(request):
    return shell_mesh(_blade(), h=0.25, method=request.param)


def test_kernel_meshes_pass_the_gate(mesh):
    assert check_shell_mesh(mesh) == []
    assert assert_shell_mesh(mesh) is mesh


def test_strict_shell_mesh_returns_the_same_mesh():
    geo = _blade()
    assert shell_mesh(geo, h=0.25, strict=True).digest() == shell_mesh(geo, h=0.25).digest()


def _problems(mesh, **changes):
    return check_shell_mesh(dataclasses.replace(mesh, **changes))


def test_gate_catches_field_length_mismatch(mesh):
    out = _problems(mesh, thickness=np.asarray(mesh.thickness)[:-1])
    assert any("thickness has" in p for p in out)


def test_gate_catches_bad_values(mesh):
    th = np.asarray(mesh.thickness, float).copy()
    th[3] = 0.0
    nodes = np.asarray(mesh.nodes, float).copy()
    nodes[0, 1] = np.nan
    out = _problems(mesh, thickness=th, nodes=nodes)
    assert any("non-positive thickness: 3" in p for p in out)
    assert any("non-finite nodes: 0" in p for p in out)
    surface = list(mesh.surface)
    surface[0] = "wing"
    assert any("unknown surface: 0" in p for p in _problems(mesh, surface=tuple(surface)))


def test_gate_catches_bad_connectivity(mesh):
    elements = list(mesh.elements)
    elements[0] = (elements[0][0], elements[0][0], elements[0][1], elements[0][2])
    assert any("repeated node" in p for p in _problems(mesh, elements=tuple(elements)))
    elements = list(mesh.elements)
    elements[1] = (len(mesh.nodes),) + tuple(elements[1][1:])
    assert any("out of range" in p for p in _problems(mesh, elements=tuple(elements)))


def test_gate_catches_flipped_element(mesh):
    i = next(k for k, s in enumerate(mesh.surface) if s == "skin" and np.mean(mesh.nodes[list(mesh.elements[k]), 2]) > 1.0)
    elements = list(mesh.elements)
    elements[i] = tuple(reversed(elements[i]))
    out = _problems(mesh, elements=tuple(elements))
    assert any("inconsistent normals" in p for p in out)


def test_gate_catches_unshared_nodes(mesh):
    """A duplicated node that one element uses instead of the original: coincident pair and a free seam."""
    i = next(k for k, s in enumerate(mesh.surface) if s == "skin" and np.mean(mesh.nodes[list(mesh.elements[k]), 2]) > 1.0)
    k = mesh.elements[i][0]
    nodes = np.vstack([mesh.nodes, mesh.nodes[k]])
    elements = list(mesh.elements)
    elements[i] = (len(mesh.nodes),) + tuple(elements[i][1:])
    out = _problems(mesh, nodes=nodes, elements=tuple(elements))
    assert any("coincident node pairs" in p for p in out)
    assert any("free edges away from the root" in p for p in out)


def test_gate_catches_duplicate_and_unused():
    geo = _blade()
    m = shell_mesh(geo, h=0.25)
    n = len(m.elements)

    def extend(values):
        return tuple(values) + (tuple(values)[0],)

    dup = dataclasses.replace(
        m, elements=extend(m.elements), name=extend(m.name), surface=extend(m.surface),
        partition=extend(m.partition), thickness=np.append(m.thickness, m.thickness[0]),
        eccentricity=np.concatenate([np.asarray(m.eccentricity), np.asarray(m.eccentricity)[:1]]),
        layers=extend(m.layers), layer_thicknesses=extend(m.layer_thicknesses),
        constructor=extend(m.constructor) if len(m.constructor) else ())
    assert any(f"duplicate elements: (0, {n})" in p for p in check_shell_mesh(dup))
    extra = dataclasses.replace(m, nodes=np.vstack([m.nodes, [[9.0, 9.0, 9.0]]]))
    assert any(f"unused nodes: {len(m.nodes)}" in p for p in check_shell_mesh(extra))


def test_invariant_error_lists_problems(mesh):
    bad = dataclasses.replace(mesh, thickness=-np.asarray(mesh.thickness))
    with pytest.raises(ShellMeshInvariantError) as err:
        assert_shell_mesh(bad)
    assert err.value.problems and "non-positive thickness" in str(err.value)


# --- property fuzzer ---------------------------------------------------------------

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import HealthCheck, given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

_blades = st.fixed_dictionaries({
    "code": st.sampled_from(["2D-F-LC-0R", "2D-F-LC-2R", "2D-CN-LC-0R", "3D-F-PA-10R", "5D-F-LE-0R", "3D-F-LE-10R"]),
    "digits": st.sampled_from(["0012", "0018", "0024", "2418", "4421"]),
    "root": st.floats(0.6, 1.5),
    "taper": st.floats(0.5, 1.0),
    "length": st.floats(2.0, 6.0),
    "twist": st.floats(-10.0, 10.0),
    "pivot": st.floats(0.25, 0.4),
})


def _fuzz(params, h_over_root, method):
    geo = _blade(**params)
    try:
        mesh = shell_mesh(geo, h=h_over_root * params["root"], method=method)
    except PartitionMeshError:
        return  # the partition mesher may refuse a blade; refusing is not a defect
    assert check_shell_mesh(mesh) == []


@settings(max_examples=6, deadline=None, derandomize=True, suppress_health_check=[HealthCheck.too_slow])
@given(params=_blades, h_over_root=st.floats(0.18, 0.35), method=st.sampled_from(["ruled", "auto"]))
def test_fuzz_generic_blades_pass_the_gate(params, h_over_root, method):
    _fuzz(params, h_over_root, method)


@pytest.mark.skipif(not os.environ.get("IGK_FUZZ"), reason="long fuzz run: set IGK_FUZZ=1 (about 4 minutes)")
@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(params=_blades, h_over_root=st.floats(0.08, 0.35), method=st.sampled_from(["ruled", "partition", "auto"]))
def test_fuzz_generic_blades_pass_the_gate_long(params, h_over_root, method):
    _fuzz(params, h_over_root, method)
