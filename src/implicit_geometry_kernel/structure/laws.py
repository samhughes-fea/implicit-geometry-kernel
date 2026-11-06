"""Span laws: a scalar as a vectorised function of the spanwise coordinate ``z``.

Member dimensions (thicknesses, depths, positions) are laws so that one assembled
structure describes the whole span. Laws are small picklable objects, never lambdas,
so a structure can be sent to worker processes. All lengths are in metres.

:class:`TableLaw` interpolates knots with monotone cubic Hermite (PCHIP,
Fritsch-Carlson) or linearly. PCHIP keeps a continuous first derivative and never
overshoots the data, so a tapering thickness never dips below its smaller knot.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence

import numpy as np


class SpanLaw(ABC):
    """A scalar as a function of ``z``. Subclasses implement :meth:`_eval`."""

    def __call__(self, z) -> np.ndarray:
        zz = np.asarray(z, dtype=np.float64)
        if zz.size > 1:
            flat = zz.ravel()
            if flat.min() == flat.max():  # one section plane, the common case
                return np.full(zz.shape, self._eval(flat[:1])[0])
        return self._eval(zz)

    @abstractmethod
    def _eval(self, z: np.ndarray) -> np.ndarray:
        """Evaluate at arbitrary ``z`` (same shape out as in)."""

    def at(self, z: float) -> float:
        """Value at one station."""
        return float(self(np.asarray([z], dtype=np.float64))[0])

    def __mul__(self, factor: float) -> SpanLaw:
        return ScaledLaw(self, factor=float(factor))

    __rmul__ = __mul__


ThicknessLaw = SpanLaw
"""A span law whose value is a thickness in metres."""


class ConstantLaw(SpanLaw):
    """The same value everywhere."""

    def __init__(self, value: float) -> None:
        self.value = float(value)

    def _eval(self, z: np.ndarray) -> np.ndarray:
        return np.full(np.shape(z), self.value, dtype=np.float64)

    def __repr__(self) -> str:
        return f"ConstantLaw({self.value!r})"


def _pchip_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fritsch-Carlson knot slopes with the usual three-point end conditions."""
    h = np.diff(x)
    delta = np.diff(y) / h
    n = x.size
    d = np.zeros(n)
    if n == 2:
        d[:] = delta[0]
        return d
    for k in range(1, n - 1):
        if delta[k - 1] * delta[k] > 0.0:
            w1 = 2.0 * h[k] + h[k - 1]
            w2 = h[k] + 2.0 * h[k - 1]
            d[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])

    def end(h0, h1, m0, m1):
        s = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
        if np.sign(s) != np.sign(m0):
            return 0.0
        if np.sign(m0) != np.sign(m1) and abs(s) > abs(3.0 * m0):
            return 3.0 * m0
        return s

    d[0] = end(h[0], h[1], delta[0], delta[1])
    d[-1] = end(h[-1], h[-2], delta[-1], delta[-2])
    return d


class TableLaw(SpanLaw):
    """Interpolated knots, held constant beyond the first and last knot.

    ``kind`` is ``"pchip"`` (monotone cubic, the default) or ``"linear"``.
    """

    def __init__(self, z: Sequence[float], values: Sequence[float], kind: str = "pchip") -> None:
        zk = np.asarray(z, dtype=np.float64)
        vk = np.asarray(values, dtype=np.float64)
        if zk.ndim != 1 or zk.shape != vk.shape or zk.size == 0:
            raise ValueError("TableLaw needs matching 1-D knot arrays")
        if not (np.all(np.isfinite(zk)) and np.all(np.isfinite(vk))):
            raise ValueError("TableLaw knots must be finite")
        order = np.argsort(zk, kind="stable")
        zk, vk = zk[order], vk[order]
        if np.any(np.diff(zk) <= 0.0):
            raise ValueError("TableLaw knot positions must be distinct")
        if kind not in ("pchip", "linear"):
            raise ValueError(f"unknown interpolation kind {kind!r}")
        self.z, self.values, self.kind = zk, vk, kind
        self._slopes = _pchip_slopes(zk, vk) if kind == "pchip" and zk.size >= 2 else None

    def _eval(self, z: np.ndarray) -> np.ndarray:
        zc = np.clip(np.asarray(z, dtype=np.float64), self.z[0], self.z[-1])
        if self.z.size == 1:
            return np.full(zc.shape, self.values[0])
        if self._slopes is None:
            return np.interp(zc, self.z, self.values)
        k = np.clip(np.searchsorted(self.z, zc, side="right") - 1, 0, self.z.size - 2)
        h = self.z[k + 1] - self.z[k]
        t = (zc - self.z[k]) / h
        t2, t3 = t * t, t * t * t
        return ((2 * t3 - 3 * t2 + 1) * self.values[k] + (t3 - 2 * t2 + t) * h * self._slopes[k]
                + (-2 * t3 + 3 * t2) * self.values[k + 1] + (t3 - t2) * h * self._slopes[k + 1])

    def __repr__(self) -> str:
        return f"TableLaw(z={self.z.tolist()}, values={self.values.tolist()}, kind={self.kind!r})"


class ScaledLaw(SpanLaw):
    """``factor * law(z) + shift``."""

    def __init__(self, law: SpanLaw, factor: float = 1.0, shift: float = 0.0) -> None:
        self.law, self.factor, self.shift = law, float(factor), float(shift)

    def _eval(self, z: np.ndarray) -> np.ndarray:
        return self.factor * self.law(z) + self.shift


class SumLaw(SpanLaw):
    """``a(z) + b(z)``, for example skin plus cap for a partition's full thickness."""

    def __init__(self, a: SpanLaw, b: SpanLaw) -> None:
        self.a, self.b = a, b

    def _eval(self, z: np.ndarray) -> np.ndarray:
        return self.a(z) + self.b(z)


class ProductLaw(SpanLaw):
    """``a(z) * b(z)``, for example a chord fraction times the chord."""

    def __init__(self, a: SpanLaw, b: SpanLaw) -> None:
        self.a, self.b = a, b

    def _eval(self, z: np.ndarray) -> np.ndarray:
        return self.a(z) * self.b(z)


def as_law(value: float | SpanLaw) -> SpanLaw:
    """A law unchanged, or a number as a :class:`ConstantLaw`."""
    return value if isinstance(value, SpanLaw) else ConstantLaw(float(value))


__all__ = [
    "ConstantLaw", "ProductLaw", "ScaledLaw", "SpanLaw", "SumLaw", "TableLaw", "ThicknessLaw",
    "as_law",
]
