"""Tests for analytic IR bounding boxes (``expr_bbox``).

``expr_bbox`` must return a conservative *superset* of the region ``expr <= 0``
so it is safe for support-based cell culling (no boundary/interior cell is ever
dropped). These tests check the closed-form bounds for each node type and a
property-based check that the box contains the empirical ``expr <= 0`` extent.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from blade_geometry_model.section_geometry.sdf.ir import (
    Box,
    Capsule,
    Circle,
    CallableField,
    HalfPlane,
    Intersect,
    Offset,
    OrientedBox,
    Polygon,
    Rotate,
    RoundedBox,
    Shell,
    Slab,
    Subtract,
    Union,
    eval_expr_batch,
    expr_bbox,
)


def _empirical_region_bbox(expr, n=400, pad=2.0):
    """Sample a dense grid and return the bbox of points with ``expr <= 0``."""
    xs = np.linspace(-pad, pad, n)
    ys = np.linspace(-pad, pad, n)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    phi = eval_expr_batch(expr, xx.ravel(), yy.ravel()).reshape(n, n)
    inside = phi <= 0.0
    if not np.any(inside):
        return None
    yi, xi = np.where(inside)
    return (float(xs[xi].min()), float(xs[xi].max()), float(ys[yi].min()), float(ys[yi].max()))


def _contains(outer, inner, tol=1e-9):
    """True if box ``inner`` is inside box ``outer`` (within ``tol``)."""
    return (
        outer[0] <= inner[0] + tol
        and outer[1] >= inner[1] - tol
        and outer[2] <= inner[2] + tol
        and outer[3] >= inner[3] - tol
    )


def test_primitive_boxes_exact():
    assert expr_bbox(Circle(1.0, -0.5, 0.25)) == (0.75, 1.25, -0.75, -0.25)
    assert expr_bbox(Box(0.0, 0.0, 2.0, 1.0)) == (-2.0, 2.0, -1.0, 1.0)
    assert expr_bbox(RoundedBox(0.0, 0.0, 1.0, 0.5, 0.1)) == (-1.1, 1.1, -0.6, 0.6)
    assert expr_bbox(Capsule(-1.0, 0.0, 1.0, 0.0, 0.3)) == (-1.3, 1.3, -0.3, 0.3)
    poly = Polygon(((0.0, 0.0), (2.0, 0.0), (1.0, 3.0)))
    assert expr_bbox(poly) == (0.0, 2.0, 0.0, 3.0)


def test_unbounded_and_opaque_return_none():
    assert expr_bbox(HalfPlane(0.0, 1.0, 0.0)) is None
    assert expr_bbox(CallableField(key="k", fn=lambda x, y: x)) is None


def test_offset_shell_slab_expand_or_keep():
    poly = Polygon(((0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)))
    # positive offset dilates
    assert expr_bbox(Offset(poly, amount=0.1)) == (-0.1, 2.1, -0.1, 1.1)
    # negative offset erodes -> base box still bounds it
    assert expr_bbox(Offset(poly, amount=-0.1)) == (0.0, 2.0, 0.0, 1.0)
    # shell band of half-thickness
    assert expr_bbox(Shell(poly, thickness=0.2)) == (-0.1, 2.1, -0.1, 1.1)
    # slab region is a sub-band of the base interior
    assert expr_bbox(Slab(poly, depth=0.3)) == (0.0, 2.0, 0.0, 1.0)


def test_intersect_subtract_use_finite_child():
    box = Box(0.0, 0.0, 1.0, 1.0)
    # intersect with an unbounded half-plane keeps the finite child's box
    assert expr_bbox(Intersect(box, HalfPlane(1.0, 0.0, 0.0))) == (-1.0, 1.0, -1.0, 1.0)
    # subtract uses the base box (cutter only removes interior)
    assert expr_bbox(Subtract(box, Circle(0.0, 0.0, 0.5))) == (-1.0, 1.0, -1.0, 1.0)


def test_rotate_box_corners():
    bb = expr_bbox(Rotate(Box(0.0, 0.0, 2.0, 1.0), angle=math.pi / 2.0, cx=0.0, cy=0.0))
    assert bb == pytest.approx((-1.0, 1.0, -2.0, 2.0), abs=1e-9)


@pytest.mark.parametrize(
    "expr",
    [
        Circle(0.3, -0.2, 0.4),
        Box(0.1, 0.0, 0.6, 0.3),
        OrientedBox(0.0, 0.0, 0.8, 0.3, angle=0.5),
        Capsule(-0.5, -0.2, 0.5, 0.3, 0.2),
        Polygon(((-0.6, -0.4), (0.7, -0.5), (0.5, 0.6), (-0.4, 0.5))),
        Shell(Polygon(((-0.6, -0.4), (0.7, -0.5), (0.5, 0.6), (-0.4, 0.5))), 0.1),
        Union(Circle(-0.4, 0.0, 0.3), Circle(0.4, 0.0, 0.3)),
        Intersect(Box(0.0, 0.0, 0.6, 0.6), HalfPlane(1.0, 0.0, 0.1)),
        Subtract(Box(0.0, 0.0, 0.6, 0.6), Circle(0.0, 0.0, 0.2)),
        Rotate(Box(0.2, 0.0, 0.5, 0.2), angle=0.7, cx=0.0, cy=0.0),
    ],
)
def test_bbox_is_conservative_superset(expr):
    box = expr_bbox(expr)
    emp = _empirical_region_bbox(expr)
    assert box is not None
    assert emp is not None
    assert _contains(box, emp, tol=1e-2), f"expr_bbox {box} does not contain region {emp}"
