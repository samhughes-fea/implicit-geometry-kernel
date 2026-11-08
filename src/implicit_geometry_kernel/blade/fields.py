"""Spanwise SDF DAG on a blade: OML offset, wall mid-surface, tagged wall_meta."""

from __future__ import annotations

from collections import OrderedDict

import numpy as np

from implicit_geometry_kernel.contract import InfeasibleSectionGeometry
from implicit_geometry_kernel.sdf.nodes.base import EvalCache, FeatureTable, FieldSample, SDFNode


def _polyline_distance(xy: np.ndarray, poly: np.ndarray) -> np.ndarray:
    xy = np.asarray(xy, float)
    poly = np.asarray(poly, float)
    if len(poly) < 2 or len(xy) == 0:
        return np.full(len(xy), np.inf)
    best = np.full(len(xy), np.inf)
    for a, b in zip(poly[:-1], poly[1:]):
        ab = b - a
        L2 = float(np.dot(ab, ab))
        if L2 < 1e-30:
            d = np.linalg.norm(xy - a, axis=1)
        else:
            t = np.clip(((xy - a) @ ab) / L2, 0.0, 1.0)
            d = np.linalg.norm(xy - (a + t[:, None] * ab), axis=1)
        best = np.minimum(best, d)
    return best


class StackThickness:
    """Skin thickness everywhere, plus cap thickness inside cap footprints.

    Nearest loop-run ``T`` (a step, not a ramp). Shared cap/skin vertices may
    evaluate to either side. Not a Newton joint spec.
    """

    def __init__(self, geo, cache_size: int = 256) -> None:
        self.geo = geo
        self._cache: OrderedDict[float, object] = OrderedDict()
        self._cache_size = int(cache_size)

    def _loop(self, z: float):
        key = round(float(z), 9)
        hit = self._cache.get(key)
        if hit is not None:
            self._cache.move_to_end(key)
            return hit
        try:
            graph = self.geo.section_at(key).shell_graph()
            payload = (graph.loop, float(self.geo._skin(key)))
        except InfeasibleSectionGeometry:
            payload = None
        self._cache[key] = payload
        while len(self._cache) > self._cache_size:
            self._cache.popitem(last=False)
        return payload

    def __call__(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, float)
        out = np.empty(points.shape[0], dtype=np.float64)
        if points.shape[0] == 0:
            return out
        z = np.round(points[:, 2], 9)
        order = np.argsort(z, kind="stable")
        groups = np.split(order, np.nonzero(np.diff(z[order]))[0] + 1)
        for group in groups:
            zkey = float(z[group[0]])
            payload = self._loop(zkey)
            if payload is None:
                out[group] = float(self.geo._skin(zkey))
                continue
            loop, skin_t = payload
            xy = points[group, :2]
            T = np.full(len(group), skin_t)
            best = np.full(len(group), np.inf)
            for run in loop:
                d = _polyline_distance(xy, run.points[:, :2])
                closer = d < best
                stack_t = float(sum(run.layer_thicknesses_m)) if run.layer_thicknesses_m else float(run.thickness)
                T[closer] = stack_t
                best[closer] = d[closer]
            out[group] = T
        return out


def span_thickness(fn):
    """Lift ``f(z) -> float`` to ``(N, 3) -> (N,)``."""

    def field(points: np.ndarray) -> np.ndarray:
        z = np.round(np.asarray(points, float)[:, 2], 12)
        uz, inv = np.unique(z, return_inverse=True)
        return np.array([float(fn(float(u))) for u in uz], dtype=np.float64)[inv]

    return field


class CutWallMetaSDF(SDFNode):
    """``wall_meta`` at each point's ``z`` from ``section_at``, cached per station.

    Infeasible cuts fall back to a skin band of the OML polygon at that ``z``.
    """

    def __init__(self, geo, table: FeatureTable | None = None, cache_size: int = 256) -> None:
        self.geo = geo
        self.table = table if table is not None else FeatureTable()
        self._cache: OrderedDict[float, SDFNode] = OrderedDict()
        self._cache_size = int(cache_size)

    def _at(self, z: float) -> SDFNode:
        from implicit_geometry_kernel.blade.network import wall_meta_from_cut

        key = round(float(z), 9)
        hit = self._cache.get(key)
        if hit is not None:
            self._cache.move_to_end(key)
            return hit
        try:
            node, _ = wall_meta_from_cut(self.geo.section_at(key), self.table)
        except InfeasibleSectionGeometry:
            from implicit_geometry_kernel.sdf.nodes.tagged import PolylineBandSDF, TagNode

            ring = self.geo.oml_polyline(key)[:, :2]
            if np.allclose(ring[0], ring[-1]):
                ring = ring[:-1]
            ring = np.vstack([ring, ring[0]])
            t = max(float(self.geo._skin(key)), 1e-4)
            spec = self.table.add("skin", "skin")
            node = TagNode(PolylineBandSDF(ring, 0.5 * t, spec.id), spec.id)
        self._cache[key] = node
        while len(self._cache) > self._cache_size:
            self._cache.popitem(last=False)
        return node

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        if points.shape[0] == 0:
            return FieldSample(np.zeros(0), np.zeros(0, dtype=np.int32))
        z = np.round(points[:, 2], 9)
        order = np.argsort(z, kind="stable")
        groups = np.split(order, np.nonzero(np.diff(z[order]))[0] + 1)
        phi = np.empty(points.shape[0])
        meta = np.empty(points.shape[0], dtype=np.int32)
        for group in groups:
            node = self._at(float(z[group[0]]))
            s = node.evaluate(points[group])
            phi[group] = s.phi
            meta[group] = s.meta
        return FieldSample(phi, meta)
