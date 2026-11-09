"""Wall-clock bake-off: per-rib Gmsh vs one-session Gmsh vs native Delaunay.

Does not assert timings (CI-flaky). Prints a table when run with ``-s``.
"""

from __future__ import annotations

import importlib.util
import time

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeConcept, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology import (
    ShellTopologyTemplate,
    _cut_tc,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
    section_cut_near,
)
from implicit_geometry_kernel.blade.topology.ribs import _live_web_fracs
from implicit_geometry_kernel.mesh.delaunay import mesh_rib_cells_many as native_many
from implicit_geometry_kernel.mesh.polygon import mesh_rib_cells_many as poly_many
from implicit_geometry_kernel.mesh.gmsh.plane import mesh_rib_cells as gmsh_one_rib
from implicit_geometry_kernel.mesh.gmsh.plane import mesh_rib_cells_many as gmsh_many


def _09_stations():
    s = 12.0 / 20.0
    root, foil = Airfoil.ellipse(1.0), Airfoil.naca4("0024")
    return [
        Station(0.0 * s, 1.80 * s, root, twist=12.0, pitch_axis=0.50),
        Station(1.0 * s, 1.88 * s, Airfoil.blend(root, foil, 0.15), twist=11.5, pitch_axis=0.47),
        Station(2.5 * s, 2.05 * s, Airfoil.blend(root, foil, 0.50), twist=10.0, pitch_axis=0.40),
        Station(4.0 * s, 2.28 * s, Airfoil.blend(root, foil, 0.85), twist=8.8, pitch_axis=0.33),
        Station(5.0 * s, 2.40 * s, foil, twist=8.0, pitch_axis=0.30),
        Station(6.5 * s, 2.20 * s, Airfoil.naca4("0022"), twist=6.5, pitch_axis=0.30),
        Station(8.5 * s, 1.85 * s, Airfoil.naca4("0020"), twist=4.8, pitch_axis=0.30),
        Station(11.0 * s, 1.50 * s, Airfoil.naca4("0018"), twist=2.8, pitch_axis=0.30),
        Station(14.0 * s, 1.10 * s, Airfoil.naca4("0015"), twist=1.0, pitch_axis=0.30),
        Station(17.0 * s, 0.90 * s, Airfoil.naca4("0013"), twist=-0.6, pitch_axis=0.30),
        Station(20.0 * s, 0.75 * s, Airfoil.naca4("0012"), twist=-2.0, pitch_axis=0.30),
    ]


def _mid_rib_job():
    sc = 12.0 / 20.0
    st = _09_stations()
    geo = BladeGeometry.from_concept(
        BladeConcept("5D-F-PA-10R"), st, 0.004 * sc,
        spar_cap_thickness=0.008 * sc, shear_web_thickness=0.004 * sc,
        pivot_fraction=([sta.z for sta in st], [sta.pitch_axis for sta in st]),
    )
    members = list(geo.transverse)
    assert members
    zm = float(members[len(members) // 2].z)
    cut = section_cut_near(geo, zm)
    assert cut is not None
    sg = cut.shell_graph()
    tmpl = ShellTopologyTemplate.from_region_graph(
        geo.region_graph, float(getattr(cut, "z", zm)), n_chord=24,
        live_web_fracs=_live_web_fracs(cut) or None, t_c=_cut_tc(cut),
        dropped_web_indices=tuple(cut.dropped_web_indices),
    )
    samples = sample_region_polylines(cut, sg, tmpl, geo.region_graph)
    rings = []
    for member in members:
        rings.extend(rib_cell_polygons_from_samples(tmpl, samples, member))
    return tmpl, samples, members, rings


def _ngon_area(ring: np.ndarray) -> float:
    xy = np.asarray(ring, float)[:, :2]
    return 0.5 * abs(float(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1))))


def _faces_area(faces) -> float:
    tot = 0.0
    for face in faces:
        xy = np.asarray(face, float)[:, :2]
        tot += 0.5 * abs(float(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1))))
    return tot


def _time(fn, repeats: int = 3):
    fn()
    times = []
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return min(times) * 1000.0, result


def test_rib_mesh_area_native():
    tmpl, samples, members, rings = _mid_rib_job()
    faces = native_many(tmpl, samples, members)
    parent = sum(_ngon_area(r) for _n, r in rings)
    got = _faces_area(f for _n, f in faces)
    assert faces
    assert abs(got - parent) / max(parent, 1e-12) < 0.08


@pytest.mark.skipif(importlib.util.find_spec("gmsh") is None, reason="gmsh not installed")
def test_rib_mesh_bench_print():
    tmpl, samples, members, rings = _mid_rib_job()
    parent = sum(_ngon_area(r) for _n, r in rings)
    n_rings = len(rings)

    def per_rib():
        out = []
        for member in members:
            out.extend(gmsh_one_rib(tmpl, samples, member))
        return out

    rows = []
    for label, fn in (
        ("gmsh one-session", lambda: gmsh_many(tmpl, samples, members)),
        ("native scipy", lambda: native_many(tmpl, samples, members, engine="scipy")),
        ("native numba", lambda: native_many(tmpl, samples, members, engine="numba")),
        ("native n-gon", lambda: poly_many(tmpl, samples, members)),
    ):
        ms, faces = _time(fn)
        ntri = len(faces)
        area = _faces_area(f for _n, f in faces)
        frac = abs(area - parent) / max(parent, 1e-12)
        rows.append((label, ms, ntri, frac))
        assert ntri > 0
        assert frac < 0.08

    print(f"\n5D-F-PA-10R mid-span ribs: {len(members)} members, {n_rings} cells, parent area {parent:.4f} m^2")
    print(f"{'backend':<20} {'ms':>10} {'tris':>8} {'area err':>10}")
    for label, ms, ntri, frac in rows:
        print(f"{label:<20} {ms:10.1f} {ntri:8d} {frac:10.3%}")
