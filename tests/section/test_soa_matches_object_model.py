"""SoA shared quadtree is the production path: two builds match, properties are finite."""

from __future__ import annotations

import numpy as np
import pytest

from implicit_geometry_kernel.section.api.factory import (
    build_section_view,
    resolve_system_type,
)
from implicit_geometry_kernel.section.properties import area_mesh as am
from implicit_geometry_kernel.section.properties.backend import quadtree_spec_for_airfoil
from implicit_geometry_kernel.integration import quadtree as qt
from implicit_geometry_kernel.integration.quadtree.quadrature import section_properties_from_quadtree
from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF


def _sorted_bounds(arr):
    a = np.asarray(arr, dtype=float)
    order = np.lexsort((a[:, 3], a[:, 2], a[:, 1], a[:, 0]))
    return a[order]


@pytest.mark.parametrize("twist_deg", [0.0, 12.0])
def test_soa_leaf_set_is_deterministic(twist_deg):
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    section = build_section_view(
        airfoil, resolve_system_type("2D-CN"), twist_angle_rad=np.deg2rad(twist_deg)
    )
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=150000)

    ev = am._section_evaluator(section)
    labels = am._section_label_keys(section)
    samplers = ev.make_samplers(labels)
    supports = am._component_supports(samplers, spec, exprs=ev.exprs_for(labels))
    multi = ev.make_multi_sampler(labels)

    a = qt.build_shared_quadtree(spec, samplers, supports=supports, multi_sampler=multi)
    b = qt.build_shared_quadtree(spec, samplers, supports=supports, multi_sampler=multi)

    lb_a = qt.leaf_bounds(a)
    lb_b = qt.leaf_bounds(b)
    assert lb_a.shape[0] > 1
    assert np.array_equal(_sorted_bounds(lb_a), _sorted_bounds(lb_b))

    for lbl in labels:
        pa = section_properties_from_quadtree(a, a.phi_by_label, lbl, sampler=samplers[lbl])
        pb = section_properties_from_quadtree(b, b.phi_by_label, lbl, sampler=samplers[lbl])
        for key in ("area", "cx", "cy", "Ixx", "Iyy", "Ixy"):
            assert np.isfinite(pa[key])
            np.testing.assert_allclose(
                float(pa[key]), float(pb[key]), rtol=1e-12, atol=0.0, err_msg=f"{lbl}.{key}"
            )
