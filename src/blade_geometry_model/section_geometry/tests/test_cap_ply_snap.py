"""BL-019: post-clamp spar-cap heights floor onto the ply pitch."""

from __future__ import annotations

import pytest

from blade_geometry_model.section_geometry.sections.cap_feasibility import (
    cap_feasibility_clamp,
)
from blade_geometry_model.section_geometry.sections.cap_ply_snap import (
    snap_cap_height_at_most,
    snap_cap_heights,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)
from blade_structure_utilities.material_utilities.layup.ply_grid import PlyGrid


def test_snap_cap_height_floors_onto_pitch() -> None:
    pitch = 0.000125
    assert snap_cap_height_at_most(0.00770153, pitch) == pytest.approx(0.007625)
    assert snap_cap_height_at_most(0.00775, pitch) == pytest.approx(0.00775)
    # Sub-ply residual is geometric run-out (align with PlyGrid.largest_at_most).
    assert snap_cap_height_at_most(0.00005, pitch) == pytest.approx(0.0)
    assert PlyGrid(t_ply=pitch, unit_plies=1).largest_at_most(0.00005) == pytest.approx(0.0)


def test_clamp_then_snap_is_on_grid_and_at_most_clamped() -> None:
    af = AirfoilSDF.from_naca("0012", n_points=200, chord=0.24)
    pitch = 0.000125
    grid = PlyGrid(t_ply=pitch, unit_plies=1)
    # Request on-grid height that the cavity will reduce.
    requested = 0.020
    _lo, _hi, h_u, h_l, info = cap_feasibility_clamp(
        af,
        0.003,
        0.30 * af.chord,
        0.50 * af.chord,
        requested,
        requested,
        y_search=None,
        phi_inner_ref=None,
        aft_half_depth_rule=False,
    )
    assert info["scaled_heights"] or h_u < requested
    snapped_u, snapped_l = snap_cap_heights(h_u, h_l, pitch)
    assert snapped_u <= h_u + 1e-15
    assert snapped_l <= h_l + 1e-15
    assert grid.is_on_grid(snapped_u)
    assert grid.is_on_grid(snapped_l)
    # Laminate ceil of an on-grid built height is identity.
    assert grid.realize(snapped_u) == pytest.approx(snapped_u)
