"""Centre-φ stash from refine must match a full re-eval attach."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.area_2d import QuadtreeSpec, build_shared_quadtree
from blade_geometry_model.area_2d.quadrature import (
    _attach_center_phi,
    section_properties_from_quadtree,
)


def _circle_sampler(radius: float = 0.25):
    def sampler(pts: np.ndarray) -> np.ndarray:
        arr = np.asarray(pts, dtype=float).reshape(-1, 2)
        return np.hypot(arr[:, 0], arr[:, 1]) - radius

    return sampler


def test_center_phi_stash_matches_full_reeval():
    spec = QuadtreeSpec(-0.6, 0.6, -0.6, 0.6, max_depth=6, min_cell_size=0.02, refine_band=0.03)
    samplers = {"a": _circle_sampler(0.15), "b": _circle_sampler(0.35)}

    def multi(pts: np.ndarray) -> dict[str, np.ndarray]:
        return {lbl: np.asarray(fn(pts), dtype=float).reshape(-1) for lbl, fn in samplers.items()}

    mesh = build_shared_quadtree(spec, samplers, multi_sampler=multi)
    stash = getattr(mesh, "_leaf_center_phi", None)
    assert stash is not None
    assert set(stash) == set(samplers)
    n_leaves = next(iter(stash.values())).shape[0]
    assert n_leaves > 0
    # Early-stop leaves should dominate; at least some centres must be stashed.
    known = sum(int(np.sum(~np.isnan(v))) for v in stash.values())
    assert known > 0

    corners = {lbl: np.asarray(mesh.phi_by_label[lbl], dtype=float)[:, :4] for lbl in samplers}
    with_reuse = _attach_center_phi(mesh, corners, samplers, multi_sampler=multi)

    # Full re-eval: drop stash and attach again.
    mesh._leaf_center_phi = None  # type: ignore[attr-defined]
    full = _attach_center_phi(mesh, corners, samplers, multi_sampler=multi)

    for lbl in samplers:
        np.testing.assert_array_equal(with_reuse[lbl][:, 4], full[lbl][:, 4])

    # Moments identical whether centres came from stash or re-eval.
    mesh._leaf_center_phi = stash  # type: ignore[attr-defined]
    props_reuse = {
        lbl: section_properties_from_quadtree(mesh, with_reuse, lbl) for lbl in samplers
    }
    props_full = {
        lbl: section_properties_from_quadtree(mesh, full, lbl) for lbl in samplers
    }
    for lbl in samplers:
        for key in ("area", "cx", "cy", "Ixx", "Iyy", "Ixy"):
            assert props_reuse[lbl][key] == props_full[lbl][key]


def test_center_phi_stash_present_on_shared_quadtree():
    spec = QuadtreeSpec(-0.5, 0.5, -0.5, 0.5, max_depth=5, min_cell_size=0.03, refine_band=0.05)
    mesh = build_shared_quadtree(spec, {"circle": _circle_sampler(0.2)})
    stash = getattr(mesh, "_leaf_center_phi", None)
    assert stash is not None
    assert "circle" in stash
    assert stash["circle"].shape[0] == mesh.phi_by_label["circle"].shape[0]
