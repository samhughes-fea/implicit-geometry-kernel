"""CSG shapes are SDF nodes: one class per shape, no compile step."""

import numpy as np
import pytest

from implicit_geometry_kernel.sdf import nodes
from implicit_geometry_kernel.sdf.compile import compile_expr_to_node
from implicit_geometry_kernel.sdf.field import SDFNode
from implicit_geometry_kernel.sdf.ir import (
    Box, Circle, Intersect, Offset, Polygon, Rotate, Shell, Slab, Subtract, Union, eval_expr_batch,
)
from implicit_geometry_kernel.sdf.nodes import ExprSDF, TagNode

SQUARE = Polygon(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)))
SHAPES = [
    Circle(0.2, 0.1, 0.5),
    Box(0.0, 0.0, 1.0, 0.5),
    SQUARE,
    Union(Circle(0, 0, 1), Circle(1.5, 0, 1)),
    Intersect(Circle(0, 0, 1), Box(0, 0, 0.5, 2)),
    Subtract(Box(0, 0, 1, 1), Circle(0, 0, 0.5)),
    Offset(SQUARE, -0.1),
    Rotate(Box(0, 0, 1, 0.2), 0.3, 0.5, 0.0),
    Shell(Circle(0, 0, 1), 0.1),
    Slab(SQUARE, 0.05),
]


@pytest.mark.parametrize("shape", SHAPES, ids=lambda s: type(s).__name__)
def test_shape_is_a_node_and_matches_memoised_evaluation(shape):
    assert isinstance(shape, SDFNode)
    assert compile_expr_to_node(shape) is shape
    rng = np.random.default_rng(3)
    xy = rng.uniform(-2.0, 2.0, size=(500, 2))
    direct = shape.evaluate(np.column_stack([xy, np.zeros(500)])).phi
    memo = eval_expr_batch(shape, xy[:, 0], xy[:, 1], {})
    plain = eval_expr_batch(shape, xy[:, 0], xy[:, 1])
    assert np.array_equal(direct, plain)
    assert np.allclose(direct, memo, rtol=0, atol=1e-15)
    assert hash(shape) == hash(shape) and shape == shape


def test_bbox_is_a_method():
    assert Circle(1.0, 2.0, 0.5).bbox() == pytest.approx((0.5, 1.5, 1.5, 2.5))


def test_feature_ids_follow_the_combinator_rules():
    a, b = TagNode(Circle(0, 0, 1), 1), TagNode(Circle(1.5, 0, 1), 2)
    pts = np.array([[-0.5, 0, 0], [2.0, 0, 0]])
    assert Union(a, b).evaluate(pts).meta.tolist() == [1, 2]
    assert Offset(a, 0.1).evaluate(pts).meta.tolist() == [1, 1]
    assert ExprSDF(Circle(0, 0, 1), 7).evaluate(pts).meta.tolist() == [7, 7]


def test_old_node_names_are_the_shapes():
    assert nodes.UnionNode is Union and nodes.OffsetNode is Offset and nodes.SlabNode is Slab
