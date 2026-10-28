"""Export shared-memo grid evaluation parity.

``export_section_json`` evaluates every component + skin-boundary field on the
contour grid through a single shared IR memo (``_shared_grid_phi``) so the shared
airfoil polygon is computed once instead of once per field. The resulting phi
must be bit-for-bit identical to a per-field ``grid.eval`` (only the evaluation
memo is shared), so the exported contours are unchanged.
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.section_geometry.props_emit.export import (
    _multicell_host,
    _shared_grid_phi,
)
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid


@pytest.mark.parametrize("twist_deg", [0.0, 15.0])
def test_shared_grid_phi_matches_grid_eval(twist_deg):
    airfoil = AirfoilSDF.from_naca("2412", n_points=80, chord=1.0)
    section = build_section_view(
        airfoil, resolve_system_type("2D-CN"), twist_angle_rad=np.deg2rad(twist_deg)
    )
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=5, max_cells=50000)
    grid = SDFGrid.from_bbox(spec.x_min, spec.x_max, spec.y_min, spec.y_max, 96, 48)

    shared = _shared_grid_phi(section, grid)

    fields = [section[label] for label in section.labels]
    host = _multicell_host(section)
    for attr in ("_skin_outer_boundary_sdf", "_skin_inner_boundary_sdf"):
        sdf = getattr(host, attr, None)
        if sdf is not None:
            fields.append(sdf)

    n_checked = 0
    for fn in fields:
        if id(fn) not in shared:
            continue
        n_checked += 1
        direct = np.asarray(grid.eval(fn), dtype=float)
        assert np.array_equal(shared[id(fn)], direct)

    # All components plus the two skin boundaries must be covered by the shared
    # eval (untwisted and twisted), so nothing falls back to an unshared eval.
    assert n_checked == len(list(section.labels)) + 2
