"""Tests for SectionFieldEvaluator shared-batch evaluation."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry.sdf.ir import Circle, Intersect, Offset, Polygon, Shell
from blade_geometry_model.section_geometry.props_emit.section_evaluator import SectionFieldEvaluator


def test_eval_labels_shares_polygon_subtree():
    verts = tuple((float(np.cos(t)), float(np.sin(t))) for t in np.linspace(0, 2 * np.pi, 64, endpoint=False))
    poly = Polygon(vertices=verts)
    shared = Shell(poly, thickness=0.04)
    exprs = {
        "a": Intersect(shared, Offset(poly, amount=-0.01)),
        "b": Intersect(shared, Offset(poly, amount=-0.02)),
    }
    evaluator = SectionFieldEvaluator(exprs)
    pts = np.column_stack(
        np.meshgrid(np.linspace(-0.5, 0.5, 30), np.linspace(-0.3, 0.3, 20), indexing="xy")
    ).reshape(-1, 2)
    memo: dict = {}
    batch = evaluator.eval_labels(["a", "b"], pts, memo=memo)
    assert set(batch) == {"a", "b"}
    assert poly in memo
    assert len([k for k in memo if isinstance(k, Polygon)]) == 1
    single_a = evaluator.sample("a", pts)
    np.testing.assert_allclose(batch["a"], single_a, atol=1e-12, rtol=0.0)


def test_eval_union_uses_cached_expr():
    exprs = {
        "a": Circle(cx=-0.1, cy=0.0, r=0.3),
        "b": Circle(cx=0.1, cy=0.0, r=0.3),
    }
    evaluator = SectionFieldEvaluator(exprs)
    assert evaluator._union_expr is not None
    x = np.linspace(-0.5, 0.5, 20)
    y = np.linspace(-0.3, 0.3, 15)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    out = evaluator.eval_union(xx, yy)
    assert np.isfinite(out).all()
