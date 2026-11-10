"""Exact polygon SDF at each point's own ``z`` (no distance interpolation)."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable, Mapping
from types import MappingProxyType

import numpy as np

from implicit_geometry_kernel.sdf.kernels.polygon_bvh import build_polygon_accel, sdf_polygon_bvh
from implicit_geometry_kernel.sdf.kernels.primitives import sdf_polygon
from implicit_geometry_kernel.sdf.nodes.base import NO_FEATURE, LeafNode

SectionPolygon = Callable[[float], np.ndarray]


class LoftedPolygonSDF(LeafNode):
    """OML as a function of ``z``: the exact 2-D polygon SDF of that station.

    Points are grouped by distinct ``z``. Distances are **not** interpolated between
    stations (interpolated distances are not distances; twist misplaces a thin TE).

    ``accelerate=True`` (default) caches a uniform-grid BVH / AABB of the station
    polygon for batched queries. ``phi`` matches the unaccelerated edge loop.
    Read the cache with :meth:`accel_at` or the :attr:`bvh` mapping.
    """

    def __init__(
        self,
        provider: SectionPolygon,
        meta: int = NO_FEATURE,
        cache_size: int = 4096,
        accelerate: bool = True,
    ) -> None:
        super().__init__(meta)
        self.provider = provider
        self.accelerate = bool(accelerate)
        self._cache: OrderedDict[float, np.ndarray] = OrderedDict()
        self._accel: OrderedDict[float, dict] = OrderedDict()
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
        if self.accelerate:
            self._accel[key] = build_polygon_accel(poly)
            self._accel.move_to_end(key)
        while len(self._cache) > self._cache_size:
            old, _ = self._cache.popitem(last=False)
            self._accel.pop(old, None)
        return poly

    def accel_at(self, z: float) -> dict | None:
        """Uniform-grid AABB for the station polygon, or ``None`` if not accelerated."""
        key = round(float(z), 9)
        self.polygon_at(key)
        return self._accel.get(key)

    @property
    def bvh(self) -> Mapping[float, dict]:
        """Read-only map ``round(z, 9) →`` accel dict (empty when ``accelerate`` is off)."""
        return MappingProxyType(self._accel)

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
            z0 = float(z[group[0]])
            poly = self.polygon_at(z0)
            xy = points[group]
            if self.accelerate:
                accel = self._accel[round(z0, 9)]
                out[group] = sdf_polygon_bvh(xy[:, 0], xy[:, 1], accel)
            else:
                out[group] = sdf_polygon(xy[:, 0], xy[:, 1], poly)
        return out


__all__ = ["LoftedPolygonSDF", "SectionPolygon"]
