"""CF-002 / SP-001: refuse silent constitutive invent."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from blade_structure_utilities.geometry_utilities.coordinate_frames.section_frames.stiffness import (
    transform_stiffness_B_to_Sprime,
)
from blade_structure_utilities.geometry_utilities.section_properties.stiffness.materials import (
    material_eg_rho,
)


def _rect_dims(*, n_points: int = 4, twist: float = 0.0):
    coords = np.array(
        [[-1.0, -0.5], [1.0, -0.5], [1.0, 0.5], [-1.0, 0.5], [-1.0, -0.5]],
        dtype=float,
    )
    padded = np.zeros((1, max(n_points, coords.shape[0]), 2))
    use = min(n_points, coords.shape[0])
    padded[0, :use] = coords[:use]
    return SimpleNamespace(
        outer_loop_geometry_S=padded,
        outer_loop_point_counts=np.array([n_points], dtype=int),
        zL_grid=np.array([0.0]),
        twist=np.array([twist]),
    )


def test_material_eg_rho_refuses_e_over_g_invent():
    with pytest.raises(ValueError, match="SP-001"):
        material_eg_rho(SimpleNamespace(E=70e9, density=1800.0))


def test_material_eg_rho_uses_supplied_g():
    E, G, rho, mode = material_eg_rho(
        SimpleNamespace(E=70e9, G=27e9, density=1800.0)
    )
    assert mode == "smeared"
    assert E == pytest.approx(70e9)
    assert G == pytest.approx(27e9)
    assert rho == pytest.approx(1800.0)


def test_stiffness_transform_refuses_degenerate_loop():
    dims = _rect_dims(n_points=2)
    with pytest.raises(ValueError, match="CF-002"):
        transform_stiffness_B_to_Sprime(
            70e9, np.array([1.0e6]), dims, EI_edge_B=None
        )


def test_stiffness_transform_with_supplied_edge_is_finite():
    dims = _rect_dims(n_points=5)
    ei_xx, ei_yy = transform_stiffness_B_to_Sprime(
        70e9, np.array([2.0e6]), dims, EI_edge_B=np.array([1.0e6])
    )
    assert np.isfinite(ei_xx).all() and np.isfinite(ei_yy).all()
    assert float(np.asarray(ei_yy).reshape(-1)[0]) > 0.0
