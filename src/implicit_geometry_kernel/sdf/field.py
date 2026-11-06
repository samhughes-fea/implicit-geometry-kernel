"""Batched SDF evaluation with integer feature ids (``meta``).

Contract (2.5-D): ``evaluate`` takes ``(N, 3)`` points ``(x, y, z)``. ``phi`` is
1-Lipschitz in ``(x, y)`` at fixed ``z``. ``meta`` is a feature id or ``NO_FEATURE``.
This is not a Euclidean 3-D SDF: ``d phi / d z`` follows the planform.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

NO_FEATURE = -1


class FeatureKind(str, Enum):
    """Wall kind stamped into ``FeatureSpec`` (same vocabulary as section labels)."""

    CAP = "cap"
    WEB = "web"
    SKIN = "skin"


KIND_PRIORITY = {FeatureKind.CAP.value: 3, FeatureKind.WEB.value: 2, FeatureKind.SKIN.value: 1}


@dataclass(frozen=True)
class FieldSample:
    """``phi`` (signed distance) and ``meta`` (feature id), both shape ``(N,)``."""

    phi: np.ndarray
    meta: np.ndarray

    def __post_init__(self) -> None:
        phi = np.asarray(self.phi, dtype=np.float64).reshape(-1)
        meta = np.asarray(self.meta, dtype=np.int32).reshape(-1)
        if phi.shape != meta.shape:
            raise ValueError(f"phi {phi.shape} and meta {meta.shape} shapes differ")
        object.__setattr__(self, "phi", phi)
        object.__setattr__(self, "meta", meta)


@dataclass(frozen=True)
class FeatureSpec:
    """One named region; ``id`` is the integer stamped into ``FieldSample.meta``."""

    id: int
    name: str
    kind: str = "skin"
    priority: int = 0


@dataclass
class FeatureTable:
    """Stable integer ids for wall labels (geometry only; no layup)."""

    _specs: list[FeatureSpec] = field(default_factory=list)
    _by_name: dict[str, FeatureSpec] = field(default_factory=dict)

    def add(self, name: str, kind: str = "skin", priority: int | None = None) -> FeatureSpec:
        hit = self._by_name.get(name)
        if hit is not None:
            return hit
        pri = KIND_PRIORITY.get(kind, 0) if priority is None else int(priority)
        spec = FeatureSpec(len(self._specs), name, kind, pri)
        self._specs.append(spec)
        self._by_name[name] = spec
        return spec

    def __getitem__(self, fid: int) -> FeatureSpec:
        return self._specs[fid]

    def by_name(self, name: str) -> FeatureSpec:
        return self._by_name[name]

    def __len__(self) -> int:
        return len(self._specs)

    def __iter__(self):
        return iter(self._specs)


class EvalCache:
    """Memoises node results for one batch of points."""

    __slots__ = ("_store", "points")

    def __init__(self, points: np.ndarray) -> None:
        self.points = points
        self._store: dict[int, FieldSample] = {}

    def lookup(self, node: object) -> FieldSample | None:
        return self._store.get(id(node))

    def store(self, node: object, sample: FieldSample) -> None:
        self._store[id(node)] = sample


class SDFNode(ABC):
    """Implicit leaf or combinator: ``(N, 3) -> FieldSample``."""

    def evaluate(self, points: np.ndarray, cache: EvalCache | None = None) -> FieldSample:
        pts = np.asarray(points, dtype=np.float64)
        if pts.ndim != 2 or pts.shape[1] != 3:
            raise ValueError(f"points must have shape (N, 3), got {pts.shape}")
        if cache is None or cache.points is not pts:
            cache = EvalCache(pts)
        hit = cache.lookup(self)
        if hit is not None:
            return hit
        sample = self._evaluate(pts, cache)
        cache.store(self, sample)
        return sample

    def evaluate_plane(self, xy: np.ndarray, z: float) -> FieldSample:
        xy = np.asarray(xy, dtype=np.float64)
        pts = np.column_stack([xy[:, 0], xy[:, 1], np.full(xy.shape[0], float(z))])
        return self.evaluate(pts)

    @abstractmethod
    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        ...

    def children(self) -> Sequence[SDFNode]:
        return ()

    def as_xy(self, z: float = 0.0):
        """``f(x, y) -> phi`` adapter at fixed ``z`` (unifies ``SDFCallable`` and ``SDFNode``)."""
        node, zc = self, float(z)

        def _fn(x, y):
            xa = np.asarray(x, dtype=np.float64)
            ya = np.asarray(y, dtype=np.float64)
            xa, ya = np.broadcast_arrays(xa, ya)
            xy = np.column_stack([xa.ravel(), ya.ravel()])
            return node.evaluate_plane(xy, zc).phi.reshape(xa.shape)

        return _fn


class LeafNode(SDFNode):
    """A leaf with a constant metadata id."""

    def __init__(self, meta: int = NO_FEATURE) -> None:
        self.meta_id = int(meta)

    def _meta(self, n: int) -> np.ndarray:
        return np.full(n, self.meta_id, dtype=np.int32)

    def _evaluate(self, points: np.ndarray, cache: EvalCache) -> FieldSample:
        return FieldSample(self._phi(points), self._meta(points.shape[0]))

    @abstractmethod
    def _phi(self, points: np.ndarray) -> np.ndarray:
        ...
