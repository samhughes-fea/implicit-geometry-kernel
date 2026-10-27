"""Subdivision phi-reuse parity.

``refine_quadtree_multi`` builds child-corner phi by reusing the decision-phase
corner cache + interior probes (``_subdivide_refined_reuse``) instead of
re-sampling. This must be bit-for-bit identical to the reference per-cell
re-evaluation (``_subdivide_cell_multi``), since the probe points coincide with
the children's new corners.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.area_2d.quadtree import (
    QuadCell,
    _cell_interior_points,
    _n_interior_probes,
    _subdivide_cell_multi,
    _subdivide_refined_reuse,
)
from blade_geometry_model.area_2d.spec import QuadtreeSpec


def _samplers():
    def make(cx, cy, r):
        def s(pts):
            a = np.asarray(pts, dtype=float).reshape(-1, 2)
            return np.hypot(a[:, 0] - cx, a[:, 1] - cy) - r

        return s

    return {"a": make(0.1, 0.0, 0.30), "b": make(-0.2, 0.1, 0.25)}


def _multi(samplers):
    def m(pts):
        return {k: np.asarray(s(pts), dtype=float).reshape(-1) for k, s in samplers.items()}

    return m


@pytest.mark.parametrize(
    "supports",
    [None, {"a": (-1.0, 1.0, -1.0, 1.0), "b": (0.0, 1.0, -1.0, 1.0)}],
)
def test_subdivide_reuse_matches_per_cell(supports):
    spec = QuadtreeSpec(-1.0, 1.0, -1.0, 1.0, max_depth=6, min_cell_size=0.001)
    samplers = _samplers()
    labels = ["a", "b"]
    msampler = _multi(samplers)

    # Two identical parents (one per code path), with cached corner phi as a
    # refined-mesh cell would carry.
    bounds = (-0.5, 0.5, -0.5, 0.5)
    p_ref = QuadCell(*bounds, level=2)
    p_new = QuadCell(*bounds, level=2)
    cphi = msampler(p_ref.corners())
    p_ref.phi_fields = {l: cphi[l].copy() for l in labels}
    p_new.phi_fields = {l: cphi[l].copy() for l in labels}

    # Reference: re-evaluate all child corners.
    ref_children = _subdivide_cell_multi(
        p_ref, samplers, supports=supports, spec=spec, multi_sampler=msampler
    )

    # Reuse: decision-phase corner cache + interior probes.
    n_probes = _n_interior_probes(spec)
    probe_phi = msampler(_cell_interior_points(p_new, spec))
    corners_by_label = {l: cphi[l].reshape(1, 4) for l in labels}
    interior_by_label = {l: probe_phi[l].reshape(1, n_probes) for l in labels}
    new_children = _subdivide_refined_reuse(
        [p_new], np.array([0]), labels, corners_by_label, interior_by_label, supports, spec
    )

    assert len(ref_children) == len(new_children) == 4
    for cr, cn in zip(ref_children, new_children):
        assert (cr.x_min, cr.x_max, cr.y_min, cr.y_max) == (
            cn.x_min,
            cn.x_max,
            cn.y_min,
            cn.y_max,
        )
        assert cr.level == cn.level
        for l in labels:
            assert np.array_equal(cr.phi_fields[l], cn.phi_fields[l]), l
