"""Rib and tip-closure triangles cover each station plane once (G19).

A fan split at a reflex corner can put a triangle outside its cell, over the next
one, and rewinding it to positive area hides that. Here every rib and tip
triangle at a station must share one orientation, and their areas must add up
to the area of their union: no triangle is folded and none overlaps another.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.mesh import shell_mesh
from implicit_geometry_kernel.section.boundary import Polygon, unary_union


def _blade(layout: str, twist: float) -> BladeGeometry:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=twist),
        Station(10.0, 0.8, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, layout, 0.004, pivot_fraction=0.3,
                         spar_cap_thickness=0.008, shear_web_thickness=0.004)


@pytest.mark.parametrize("layout, twist", [("2D-F-LC-2R", 0.0), ("1B-CN-LC-6R", 8.0)])
def test_rib_and_tip_triangles_neither_fold_nor_overlap(layout, twist):
    mesh = shell_mesh(_blade(layout, twist), h=0.25)
    nodes = np.asarray(mesh.nodes)
    by_z: dict[float, list[np.ndarray]] = defaultdict(list)
    for surface, element in zip(mesh.surface, mesh.elements):
        if surface in ("rib", "tip_closure"):
            assert len(element) == 3
            by_z[round(float(nodes[element[0], 2]), 9)].append(nodes[list(element)][:, :2])
    assert len(by_z) >= 2
    for z, tris in by_z.items():
        t = np.asarray(tris)
        signed = 0.5 * ((t[:, 1, 0] - t[:, 0, 0]) * (t[:, 2, 1] - t[:, 0, 1])
                        - (t[:, 2, 0] - t[:, 0, 0]) * (t[:, 1, 1] - t[:, 0, 1]))
        assert np.all(signed > 0.0) or np.all(signed < 0.0), f"folded triangle at z = {z}"
        union = unary_union([Polygon(tri) for tri in tris]).area
        assert abs(np.abs(signed).sum() - union) <= 1e-9 * union, f"overlap at z = {z}"
