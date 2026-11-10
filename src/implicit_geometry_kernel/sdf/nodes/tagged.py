"""Tagged SDF combinators: stamp and propagate integer feature ids."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from implicit_geometry_kernel.sdf.nodes.base import EvalCache, FieldSample, LeafNode, NO_FEATURE, SDFNode
from implicit_geometry_kernel.sdf.kernels.polygon_numba import sdf_polyline_unsigned_numba
from implicit_geometry_kernel.sdf.kernels.primitives import sdf_circle, sdf_polygon


class CircleSDF(LeafNode):
    """Exact circle in the section plane (z ignored)."""

    def __init__(self, cx: float = 0.0, cy: float = 0.0, r: float = 1.0, meta: int = NO_FEATURE) -> None:
        super().__init__(meta)
        self.cx, self.cy, self.r = float(cx), float(cy), float(r)

    def _phi(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(sdf_circle(points[:, 0], points[:, 1], self.cx, self.cy, self.r), dtype=np.float64)


class PolygonSDF(LeafNode):
    """Exact closed-polygon SDF in the section plane (z ignored)."""

    def __init__(self, vertices: np.ndarray, meta: int = NO_FEATURE) -> None:
        super().__init__(meta)
        self.vertices = np.asarray(vertices, dtype=np.float64)

    def _phi(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(sdf_polygon(points[:, 0], points[:, 1], self.vertices), dtype=np.float64)


class PolylineBandSDF(LeafNode):
    """Unsigned distance to an open polyline, minus ``radius`` (a wall strip)."""

    def __init__(self, vertices: np.ndarray, radius: float, meta: int = NO_FEATURE) -> None:
        super().__init__(meta)
        self.vertices = np.ascontiguousarray(vertices, dtype=np.float64)
        self.radius = float(radius)

    def _phi(self, points: np.ndarray) -> np.ndarray:
        d = sdf_polyline_unsigned_numba(points[:, 0], points[:, 1], self.vertices)
        return np.asarray(d, dtype=np.float64) - self.radius


def as_node(func: Callable, meta: int = NO_FEATURE) -> SDFNode:
    """Wrap ``f(x, y) -> phi`` as an ``SDFNode`` (z ignored)."""
    return CallableSDF(func, meta=meta)


class CallableSDF(LeafNode):
    """Wrap ``f(x, y) -> phi`` as an ``SDFNode`` (z ignored)."""

    def __init__(self, func: Callable, meta: int = NO_FEATURE) -> None:
        super().__init__(meta)
        self.func = func

    def _phi(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(self.func(points[:, 0], points[:, 1]), dtype=np.float64)


class TagNode(SDFNode):
    """Stamp a constant feature id; ``phi`` unchanged."""

    def __init__(self, node: SDFNode, meta: int) -> None:
        self.node = node
        self.meta_id = int(meta)

    def children(self) -> tuple[SDFNode, ...]:
        return (self.node,)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        s = self.node.evaluate(points, cache)
        return FieldSample(s.phi, np.full(s.phi.shape, self.meta_id, dtype=np.int32))


class ExprSDF(SDFNode):
    """A CSG shape stamped with a feature id (``meta``); ``phi`` unchanged.

    Every shape in :mod:`implicit_geometry_kernel.sdf.nodes` is already an ``SDFNode``;
    this wrapper only adds the id (without one it behaves as the shape itself).
    """

    def __init__(self, expr, meta: int = NO_FEATURE) -> None:
        self.expr = expr
        self.meta_id = int(meta)

    def children(self) -> tuple[SDFNode, ...]:
        return (self.expr,)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        s = self.expr.evaluate(points, cache)
        return FieldSample(s.phi, np.full(s.phi.shape, self.meta_id, dtype=np.int32))


class PriorityUnionNode(SDFNode):
    """Union whose metadata is the highest-priority child containing the point.

    ``phi`` is ``min``. Outside every child, meta is that of the nearest child.
    """

    def __init__(self, nodes: Sequence[SDFNode], priorities: Sequence[int]) -> None:
        if len(nodes) != len(priorities) or not nodes:
            raise ValueError("nodes and priorities must be non-empty and of equal length")
        self.nodes = tuple(nodes)
        self.priorities = np.asarray(priorities, dtype=np.int64)

    def children(self) -> tuple[SDFNode, ...]:
        return self.nodes

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        samples = [n.evaluate(points, cache) for n in self.nodes]
        phi = np.stack([s.phi for s in samples])
        meta = np.stack([s.meta for s in samples])
        cols = np.arange(phi.shape[1])
        inside = phi < 0.0
        score = np.where(inside, self.priorities[:, None], np.iinfo(np.int64).min)
        best = np.argmax(score, axis=0)
        nearest = np.argmin(phi, axis=0)
        k = np.where(inside.any(axis=0), best, nearest)
        return FieldSample(phi.min(axis=0), meta[k, cols].astype(np.int32))


class TrimNode(SDFNode):
    """Hard clip ``max(child, clip)``; metadata from the child."""

    def __init__(self, node: SDFNode, clip: SDFNode) -> None:
        self.node, self.clip = node, clip

    def children(self) -> tuple[SDFNode, ...]:
        return (self.node, self.clip)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        s = self.node.evaluate(points, cache)
        c = self.clip.evaluate(points, cache)
        return FieldSample(np.maximum(s.phi, c.phi), s.meta)


class FieldOffsetNode(SDFNode):
    """``phi = base + scale * f(points)``; meta from ``base``.

    Wall mid-surface: ``FieldOffsetNode(oml, thickness, 0.5)``. Cavity:
    ``FieldOffsetNode(oml, thickness, 1.0)``. ``thickness`` maps ``(N, 3) -> (N,)``.
    """

    def __init__(self, base: SDFNode, field, scale: float = 1.0) -> None:
        self.base = base
        self.field = field
        self.scale = float(scale)

    def children(self) -> tuple[SDFNode, ...]:
        return (self.base,)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        s = self.base.evaluate(points, cache)
        f = np.asarray(self.field(points), dtype=np.float64).reshape(-1)
        return FieldSample(s.phi + self.scale * f, s.meta)


class MidSurfaceNode(SDFNode):
    """Medial field ``(outer + cavity) / 2``; zero set is the wall mid-surface."""

    def __init__(self, outer: SDFNode, cavity: SDFNode) -> None:
        self.outer, self.cavity = outer, cavity

    def children(self) -> tuple[SDFNode, ...]:
        return (self.outer, self.cavity)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        o = self.outer.evaluate(points, cache)
        c = self.cavity.evaluate(points, cache)
        return FieldSample(0.5 * (o.phi + c.phi), o.meta)

    def thickness(self, points: np.ndarray) -> np.ndarray:
        cache = EvalCache(np.asarray(points, dtype=np.float64))
        return self.cavity.evaluate(cache.points, cache).phi - self.outer.evaluate(cache.points, cache).phi


__all__ = [
    "CallableSDF", "CircleSDF", "ExprSDF", "FieldOffsetNode", "MidSurfaceNode",
    "PolygonSDF", "PolylineBandSDF", "PriorityUnionNode",
    "TagNode", "TrimNode", "as_node",
]
