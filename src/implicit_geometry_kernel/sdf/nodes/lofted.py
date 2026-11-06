"""Exact polygon SDF at each point's own ``z`` (no distance interpolation)."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable

import numpy as np

from implicit_geometry_kernel.sdf.field import NO_FEATURE, LeafNode
from implicit_geometry_kernel.sdf.primitives import sdf_polygon

SectionPolygon = Callable[[float], np.ndarray]


class LoftedPolygonSDF(LeafNode):
    """OML as a function of ``z``: the exact 2-D polygon SDF of that station.

    Points are grouped by distinct ``z``. Distances are **not** interpolated between
    stations (interpolated distances are not distances; twist misplaces a thin TE).
    """

    def __init__(self, provider: SectionPolygon, meta: int = NO_FEATURE, cache_size: int = 4096) -> None:
        super().__init__(meta)
        self.provider = provider
        self._cache: OrderedDict[float, np.ndarray] = OrderedDict()
        self._cache_size = int(cache_size)

    def polygon_at(self, z: float) -> np.ndarray:
        key = round(float(z), 9)
        hit = self._cache.get(key)
        if hit is not None:
            self._cache.move_to_end(key)
            return hit
        poly = np.asarray(self.provider(key), dtype=np.float64)
        if poly.ndim != 2 or poly.shape[1] < 2:
            raise ValueError(f"section polygon at z={key} must be (N, 2)")
        self._cache[key] = poly
        while len(self._cache) > self._cache_size:
            self._cache.popitem(last=False)
        return poly

    def _phi(self, points: np.ndarray) -> np.ndarray:
        if points.shape[0] == 0:
            return np.zeros(0)
        z = np.round(points[:, 2], 9)
        order = np.argsort(z, kind="stable")
        zs = z[order]
        cuts = np.nonzero(np.diff(zs))[0] + 1
        groups = np.split(order, cuts)
        out = np.empty(points.shape[0])
        for group in groups:
            poly = self.polygon_at(float(z[group[0]]))
            xy = points[group]
            out[group] = sdf_polygon(xy[:, 0], xy[:, 1], poly)
        return out
