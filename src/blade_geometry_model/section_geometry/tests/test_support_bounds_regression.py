"""Perf-regression guard for analytic IR support bounds.

Locks in the analytic-bounds wave on the canonical benchmark (2D-CN NACA 2412,
depth 7). Asserts:

1. ``_component_supports`` performs **zero** SDF evaluations (analytic path) and
   every support box is a conservative superset of the component's true region.
2. The deterministic profiling counters (leaf count, polygon/phi eval counts)
   and per-label section properties match the rebaselined snapshot in
   ``baseline_analytic_bounds.json`` within tolerance.

If a future change legitimately shifts these numbers, regenerate the baseline
(see the capture snippet in PERFORMANCE_WORK_SUMMARY.md) and review the diff.
"""

from __future__ import annotations

import json
from pathlib import Path

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
from blade_geometry_model.section_geometry.sdf.airfoil import AirfoilSDF
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.sdf.profiling import reset_stats

BASELINE = Path(__file__).with_name("baseline_analytic_bounds.json")


@pytest.fixture(scope="module")
def baseline():
    return json.loads(BASELINE.read_text())


@pytest.fixture(scope="module")
def canonical():
    airfoil = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    section = build_section_view(airfoil, resolve_system_type("2D-CN"), twist_angle_rad=0.0)
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=150000)
    return airfoil, section, spec


def test_supports_do_zero_sdf_eval_and_bound_region(canonical):
    _airfoil, section, spec = canonical
    ev = am._section_evaluator(section)
    labels = am._section_label_keys(section)

    calls = {"n": 0}

    def counting(sampler):
        def wrapped(pts):
            calls["n"] += 1
            return sampler(pts)

        return wrapped

    base = ev.make_samplers(labels)
    samplers = {lbl: counting(s) for lbl, s in base.items()}

    supports = am._component_supports(samplers, spec, exprs=ev.exprs_for(labels))

    # No fallback sampling: every component has a bounded analytic bbox.
    assert calls["n"] == 0
    assert set(supports) == set(labels)

    # Each analytic support must contain the component's true region (sampled).
    xs = np.linspace(spec.x_min, spec.x_max, 200)
    ys = np.linspace(spec.y_min, spec.y_max, 200)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    for lbl in labels:
        phi = np.asarray(base[lbl](pts), dtype=float)
        inside = phi <= 0.0
        if not np.any(inside):
            continue
        x_in, y_in = pts[inside, 0], pts[inside, 1]
        x0, x1, y0, y1 = supports[lbl]
        tol = 1e-9
        assert x0 <= x_in.min() + tol and x1 >= x_in.max() - tol, lbl
        assert y0 <= y_in.min() + tol and y1 >= y_in.max() - tol, lbl


def test_counts_and_properties_match_baseline(canonical, baseline, monkeypatch):
    _airfoil, section, spec = canonical
    monkeypatch.setenv("SECTION_GEOMETRY_PROFILE", "1")
    reset_stats()
    props, meta = compute_section_properties(
        section, backend="quadtree", quadtree_spec=spec, airfoil_for_spec=_airfoil
    )
    prof = meta["profiling"]

    # Deterministic counters: must match exactly.
    for key in ("sdf_polygon_count", "phi_eval_count", "quadtree_leaf_count"):
        assert prof[key] == baseline["counts"][key], (
            f"{key}: got {prof[key]}, baseline {baseline['counts'][key]}"
        )

    # Per-label section properties within tight tolerance (sub-discretization).
    for lbl, ref in baseline["properties"].items():
        for key, ref_v in ref.items():
            got = float(props[lbl][key])
            np.testing.assert_allclose(
                got, ref_v, rtol=1e-6, atol=1e-9,
                err_msg=f"{lbl}.{key}: got {got}, baseline {ref_v}",
            )
