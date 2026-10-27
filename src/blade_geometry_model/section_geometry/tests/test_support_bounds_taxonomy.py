"""Generalization guard: analytic support bounds across the system-type taxonomy.

`docs/system_type_xyz_taxonomy.md` defines 40 registered system types (0A/0B
skin-only + 38 multicell across X=1..5 webs, Y=A/B/C/D families, Z=CN/F). The
analytic-IR-bounds path (`expr_bbox` / `_component_supports`) and the quadtree
property backend must work for all of them, not just the 2D-CN benchmark.

These tests assert, per representative system type:
- support boxes are conservative supersets of each component's true region;
- multicell families derive bounds analytically (zero SDF sampling), while the
  zero_web_catalogue 0A/0B layouts fall back to sampling without raising;
- the quadtree backend agrees with the dense-grid backend (correctness).
"""

from __future__ import annotations

import numpy as np
import pytest

from blade_geometry_model.section_geometry.api.factory import (
    resolve_system_type,
    build_section_view,
)
# Import the real home, not the deprecated ``engine`` shim: that shim is a
# ``from ..props_emit.area_mesh import *``, and star-import never re-exports
# underscore names, so ``_section_evaluator`` / ``_section_label_keys`` are
# invisible through it.
from blade_geometry_model.section_geometry.props_emit import area_mesh as am
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.props_emit.props_compute import compute_section_properties
from blade_geometry_model.section_geometry.sdf.grid import SDFGrid

# Representative subset: every family (A/B/C/D), both orientations (CN/F),
# a range of web counts, plus the two skin-only layouts.
REPRESENTATIVE = [
    "0A", "0B",
    "1A-CN", "1B-F", "1C-CN",
    "2B-CN", "2C-F", "2D-CN", "2D-F",
    "3A-F", "3C-CN",
    "5A-CN", "5D-F",
]
MULTICELL = [k for k in REPRESENTATIVE if k not in ("0A", "0B")]
PARITY = ["0A", "0B", "2D-CN", "2C-F", "3A-F"]


@pytest.fixture(scope="module")
def airfoil():
    return AirfoilSDF.from_naca("2412", n_points=100, chord=1.0)


def _twist_for(layout):
    return float(np.deg2rad(12.0)) if layout.web_orientation == "flapwise" else 0.0


@pytest.mark.parametrize("key", REPRESENTATIVE)
def test_supports_cover_true_region(key, airfoil):
    layout = resolve_system_type(key)
    section = build_section_view(airfoil, layout, twist_angle_rad=_twist_for(layout))
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=150000)
    ev = am._section_evaluator(section)
    labels = am._section_label_keys(section)
    samplers = ev.make_samplers(labels)

    supports = am._component_supports(samplers, spec, exprs=ev.exprs_for(labels))
    assert set(supports) == set(labels)

    xs = np.linspace(spec.x_min, spec.x_max, 220)
    ys = np.linspace(spec.y_min, spec.y_max, 220)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    for lbl in labels:
        phi = np.asarray(samplers[lbl](pts), dtype=float)
        inside = phi <= 0.0
        if not np.any(inside):
            continue
        x_in, y_in = pts[inside, 0], pts[inside, 1]
        x0, x1, y0, y1 = supports[lbl]
        tol = 1e-9
        assert x0 <= x_in.min() + tol and x1 >= x_in.max() - tol, f"{key}:{lbl} x"
        assert y0 <= y_in.min() + tol and y1 >= y_in.max() - tol, f"{key}:{lbl} y"


@pytest.mark.parametrize("key", MULTICELL)
def test_multicell_supports_zero_sampling(key, airfoil):
    """Multicell families have full IR exprs ⇒ analytic bounds, no SDF sampling."""
    layout = resolve_system_type(key)
    section = build_section_view(airfoil, layout, twist_angle_rad=_twist_for(layout))
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=150000)
    ev = am._section_evaluator(section)
    labels = am._section_label_keys(section)

    calls = {"n": 0}

    def counting(sampler):
        def wrapped(pts):
            calls["n"] += 1
            return sampler(pts)

        return wrapped

    samplers = {lbl: counting(s) for lbl, s in ev.make_samplers(labels).items()}
    am._component_supports(samplers, spec, exprs=ev.exprs_for(labels))
    assert calls["n"] == 0, f"{key}: expected analytic bounds, got {calls['n']} sampler calls"


@pytest.mark.parametrize("key", PARITY)
def test_quadtree_backend_matches_grid(key, airfoil):
    """Quadtree backend agrees with the dense-grid backend for every layout."""
    layout = resolve_system_type(key)
    twist = _twist_for(layout)
    section = build_section_view(airfoil, layout, twist_angle_rad=twist)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=6, max_cells=150000)
    grid = SDFGrid.from_bbox(spec.x_min, spec.x_max, spec.y_min, spec.y_max, 800, 480)

    pq, _ = compute_section_properties(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=airfoil
    )
    pg, _ = compute_section_properties(section, backend="grid", grid=grid)

    for lbl in pq:
        aq, ag = float(pq[lbl]["area"]), float(pg[lbl]["area"])
        # Area agreement to within a small fraction (quadtree vs grid discretization).
        assert abs(aq - ag) <= 5e-4 + 0.02 * max(ag, 1e-9), f"{key}:{lbl} area q={aq} g={ag}"
