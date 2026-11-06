"""Combinators that structure assembly needs on top of the basic CSG nodes.

* :class:`ComplementNode` flips inside and outside.
* :class:`MaxNode` intersects any number of fields, taking ``meta`` from one child.
* :class:`LawOffsetNode` offsets a field by an amount that varies with ``z``.
* :class:`SpanIntervalNode` keeps a field only for ``z`` in a closed interval.

All of them keep the 2.5-D contract of :mod:`implicit_geometry_kernel.sdf.field`:
``phi`` is a distance in ``(x, y)`` at each point's own ``z``.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from implicit_geometry_kernel.sdf.field import EvalCache, FieldSample, SDFNode

SpanFunction = Callable[[np.ndarray], np.ndarray]
"""Any vectorised ``z -> value`` map, for example a ``structure`` span law."""


class ComplementNode(SDFNode):
    """``-phi``: the outside becomes the inside. ``meta`` is kept."""

    def __init__(self, node: SDFNode) -> None:
        self.node = node

    def children(self) -> tuple[SDFNode, ...]:
        return (self.node,)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        s = self.node.evaluate(points, cache)
        return FieldSample(-s.phi, s.meta)


class MaxNode(SDFNode):
    """Intersection of several fields, ``max`` over children.

    ``meta`` comes from child ``meta_from`` when given (a mask should not relabel the
    member it restricts), otherwise from the child that sets the maximum.
    """

    def __init__(self, *nodes: SDFNode, meta_from: int | None = None) -> None:
        if not nodes:
            raise ValueError("MaxNode needs at least one child")
        if meta_from is not None and not 0 <= meta_from < len(nodes):
            raise ValueError(f"meta_from={meta_from} is not a child index")
        self.nodes = tuple(nodes)
        self.meta_from = meta_from

    def children(self) -> tuple[SDFNode, ...]:
        return self.nodes

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        samples = [n.evaluate(points, cache) for n in self.nodes]
        phi = np.stack([s.phi for s in samples])
        if self.meta_from is not None:
            meta = samples[self.meta_from].meta
        else:
            k = np.argmax(phi, axis=0)
            meta = np.stack([s.meta for s in samples])[k, np.arange(phi.shape[1])]
        return FieldSample(phi.max(axis=0), meta)


class LawOffsetNode(SDFNode):
    """``phi - scale * law(z)``: a positive amount grows the inside."""

    def __init__(self, node: SDFNode, law: SpanFunction, scale: float = 1.0) -> None:
        self.node = node
        self.law = law
        self.scale = float(scale)

    def children(self) -> tuple[SDFNode, ...]:
        return (self.node,)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        s = self.node.evaluate(points, cache)
        amount = np.asarray(self.law(points[:, 2]), dtype=np.float64).reshape(-1)
        return FieldSample(s.phi - self.scale * amount, s.meta)


class SpanIntervalNode(SDFNode):
    """The child for ``z_start <= z <= z_end``; empty (``phi = +inf``) elsewhere.

    The member does not exist outside its interval, so there is no in-plane distance
    to it there; ``+inf`` keeps unions and priorities correct.
    """

    def __init__(self, node: SDFNode, z_start: float, z_end: float, tol: float = 1e-9) -> None:
        if not z_end >= z_start:
            raise ValueError(f"span interval needs z_end >= z_start, got ({z_start}, {z_end})")
        self.node = node
        self.z_start = float(z_start)
        self.z_end = float(z_end)
        self.tol = float(tol)

    def children(self) -> tuple[SDFNode, ...]:
        return (self.node,)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        s = self.node.evaluate(points, cache)
        z = points[:, 2]
        inside = (z >= self.z_start - self.tol) & (z <= self.z_end + self.tol)
        return FieldSample(np.where(inside, s.phi, np.inf), s.meta)


__all__ = ["ComplementNode", "LawOffsetNode", "MaxNode", "SpanFunction", "SpanIntervalNode"]
