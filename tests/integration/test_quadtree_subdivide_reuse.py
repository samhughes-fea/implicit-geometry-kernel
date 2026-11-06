"""SoA subdivision gathers child-corner φ from parent corners + interior probes."""

from __future__ import annotations

import numpy as np
import pytest

from implicit_geometry_kernel.integration.quadtree.soa import (
    _CHILD_CORNER_RECIPE,
    _expand_children_soa,
)
from implicit_geometry_kernel.integration.spec import QuadtreeSpec


@pytest.mark.parametrize(
    "supports",
    [None, {"a": (-1.0, 1.0, -1.0, 1.0)}],
)
def test_expand_children_soa_bounds_and_gathered_phi(supports):
    spec = QuadtreeSpec(-1.0, 1.0, -1.0, 1.0, max_depth=6, min_cell_size=0.001)
    bounds = np.array([[-0.5, 0.5, -0.5, 0.5]], dtype=float)
    levels = np.array([2], dtype=int)
    labels = ["a"]
    corners = {"a": np.array([[1.0, 2.0, 3.0, 4.0]], dtype=float)}
    probes = {"a": np.array([[10.0, 11.0, 12.0, 13.0, 14.0]], dtype=float)}

    child_bounds, child_levels, child_corners = _expand_children_soa(
        bounds, levels, corners, probes, labels, supports, spec
    )
    assert child_bounds.shape == (4, 4)
    assert np.array_equal(child_levels, np.array([3, 3, 3, 3]))
    np.testing.assert_allclose(child_bounds[0], [-0.5, 0.0, -0.5, 0.0])
    np.testing.assert_allclose(child_bounds[1], [0.0, 0.5, -0.5, 0.0])
    np.testing.assert_allclose(child_bounds[2], [0.0, 0.5, 0.0, 0.5])
    np.testing.assert_allclose(child_bounds[3], [-0.5, 0.0, 0.0, 0.5])

    c = corners["a"][0]
    p = probes["a"][0]
    for k, recipe in enumerate(_CHILD_CORNER_RECIPE):
        expected = np.array([c[idx] if kind == "c" else p[idx] for kind, idx in recipe])
        np.testing.assert_allclose(child_corners["a"][k], expected)
