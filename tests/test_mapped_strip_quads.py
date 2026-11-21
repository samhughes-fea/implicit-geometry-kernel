"""Ruled and split loft quads match the per-edge Python loops."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology.map_loft import (
    _append_face,
    _quad_bowtie,
    _ruled_quads,
    _split_quads,
)
from implicit_geometry_kernel.mesh.shell import TIP_H_OVER_C, _sample_stations, mesh_stations


def _ruled_loops(name, pa, pb):
    out = []
    for i in range(len(pa) - 1):
        quad = np.stack([pa[i], pa[i + 1], pb[i + 1], pb[i]])
        if _quad_bowtie(quad):
            a, b, c, d = quad
            _append_face(out, name, np.stack([a, b, c]), True)
            _append_face(out, name, np.stack([a, c, d]), True)
        else:
            _append_face(out, name, quad, False)
    return out


def _split_loops(name, pa, pb):
    mids = [0.5 * (pa[i] + pb[i]) for i in range(len(pa))]
    out = []
    for i in range(len(pa) - 1):
        _append_face(out, name, np.stack([pa[i], pa[i + 1], mids[i + 1], mids[i]]), False)
        _append_face(out, name, np.stack([mids[i], mids[i + 1], pb[i + 1], pb[i]]), False)
    return out


def _same_faces(got, expected) -> None:
    assert len(got) == len(expected)
    for (name, face, collapsing), (ename, eface, ecollapsing) in zip(got, expected):
        assert name == ename
        assert collapsing is ecollapsing
        assert face.shape == eface.shape
        assert np.allclose(face, eface)


def _pair(n: int, z0: float, z1: float, bow: bool) -> tuple[np.ndarray, np.ndarray]:
    s = np.linspace(0.0, 1.0, n)
    pa = np.column_stack([s, np.zeros(n), np.full(n, z0)])
    if bow:
        pb = np.column_stack([s[::-1], np.linspace(0.2, -0.2, n), np.full(n, z1)])
    else:
        pb = np.column_stack([s, np.full(n, 0.3), np.full(n, z1)])
    return pa, pb


def test_a_ruled_strip_does_not_ask_for_a_section():
    from implicit_geometry_kernel.blade.topology.map_loft import loft_mapped_strips
    from implicit_geometry_kernel.blade.topology.template import LoopBranch, ShellTopologyTemplate

    def station(z):
        loop = (LoopBranch("panel", ("TE",), ("anchor", "le"), 4, "suction"),)
        return ShellTopologyTemplate(
            z=z, signature=(), loop=loop, webs=(), vertices=(("TE",), ("anchor", "le")),
        )

    def line(z):
        return np.column_stack([np.linspace(0.0, 1.0, 4), np.zeros(4), np.full(4, z)])

    def on_section(*_args, **_kwargs):
        raise AssertionError("ruled bay projected a seam")

    faces = loft_mapped_strips(
        {"panel": line(0.0)}, {"panel": line(1.0)},
        station(0.0), station(1.0), on_section=on_section,
    )
    assert faces
    assert {face.constructor for face in faces} == {"ruled"}
    assert all(len(face.face) == 4 for face in faces)


def test_array_quads_match_python_loops():
    for n, bow in ((2, False), (5, False), (8, True)):
        pa, pb = _pair(n, 0.0, 1.0, bow)
        _same_faces(_ruled_quads("skin", pa, pb), _ruled_loops("skin", pa, pb))
        _same_faces(_split_quads("web", pa, pb), _split_loops("web", pa, pb))


def _parametric():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=8.0),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )


def test_parametric_blade_quads_match_python_loops():
    try:
        geo = _parametric()
        zs = mesh_stations(geo, h_over_c=TIP_H_OVER_C)
        _cuts, _tmpls, samples = _sample_stations(geo, zs, 8, None, TIP_H_OVER_C, None)
    except OSError as exc:
        pytest.skip(f"GEOS is required for the parametric blade: {exc}")
    compared = 0
    for i in range(len(samples) - 1):
        for name, pa in samples[i].items():
            pb = samples[i + 1].get(name)
            if pb is None or len(pa) != len(pb) or len(pa) < 2:
                continue
            pa = np.asarray(pa, float)
            pb = np.asarray(pb, float)
            _same_faces(_ruled_quads(name, pa, pb), _ruled_loops(name, pa, pb))
            _same_faces(_split_quads(name, pa, pb), _split_loops(name, pa, pb))
            compared += 1
    assert compared > 0
