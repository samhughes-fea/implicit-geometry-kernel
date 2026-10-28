"""Struct-of-arrays quadtree vs the object-model reference.

``build_shared_quadtree`` builds the shared section quadtree with the SoA
refinement (`_refine_soa`). This must produce the *same leaf set* (identical
cell bounds, identical count) as the retained object-model path
(`build_square_root_forest` + `refine_quadtree_multi`). Per-label properties may
differ only by float summation round-off, because SoA emits leaves in
refinement order rather than DFS order (see PERFORMANCE_WORK_SUMMARY Phase 9).
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
# Real home, not the deprecated ``engine`` shim: that shim star-imports, and
# star-import never re-exports underscore names.
from blade_geometry_model.section_geometry.props_emit import area_mesh as am
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.area_2d import quadtree as qt
from blade_geometry_model.area_2d.quadrature import section_properties_from_quadtree
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF


def _object_model_mesh(spec, samplers, supports, multi_sampler):
    """Pre-SoA object-model build (square-cells branch)."""
    mesh = qt.build_square_root_forest(spec, samplers)
    mesh = qt.refine_quadtree_multi(
        mesh, samplers, supports=supports, multi_sampler=multi_sampler
    )
    qt._refresh_leaf_cache(mesh)
    mesh.phi_by_label = qt._phi_by_label_from_leaves(mesh.leaf_cells(), samplers)
    if supports is not None:
        mesh.component_supports = supports
    return mesh


def _sorted_bounds(arr):
    a = np.asarray(arr, dtype=float)
    order = np.lexsort((a[:, 3], a[:, 2], a[:, 1], a[:, 0]))
    return a[order]


@pytest.mark.parametrize("twist_deg", [0.0, 12.0])
def test_soa_leaf_set_and_properties_match_object_model(twist_deg):
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

    soa = qt.build_shared_quadtree(spec, samplers, supports=supports, multi_sampler=multi)
    obj = _object_model_mesh(spec, samplers, supports, multi)

    lb_soa = qt.leaf_bounds(soa)
    lb_obj = qt.leaf_bounds(obj)

    # Identical leaf set (count + exact bounds, order aside).
    assert lb_soa.shape[0] == lb_obj.shape[0]
    assert np.array_equal(_sorted_bounds(lb_soa), _sorted_bounds(lb_obj))

    # Properties match to summation round-off.
    for lbl in labels:
        ps = section_properties_from_quadtree(soa, soa.phi_by_label, lbl, sampler=samplers[lbl])
        po = section_properties_from_quadtree(obj, obj.phi_by_label, lbl, sampler=samplers[lbl])
        for key in ("area", "cx", "cy", "Ixx", "Iyy", "Ixy"):
            np.testing.assert_allclose(
                float(ps[key]), float(po[key]), rtol=1e-7, atol=1e-12,
                err_msg=f"{lbl}.{key}",
            )
