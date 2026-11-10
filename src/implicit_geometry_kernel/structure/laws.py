"""Span laws: a scalar as a vectorised function of the spanwise coordinate ``z``.

Member dimensions (thicknesses, depths, positions) are laws so that one assembled
structure describes the whole span. Laws are small picklable objects, never lambdas,
so a structure can be sent to worker processes. All lengths are in metres.

:class:`TableLaw` interpolates knots linearly and holds the end values outside the
table. That matches ply-drop gauges: constant, then a linear ramp, then constant.
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


class ThicknessLaw(ABC):
    """Thickness in metres at span ``z`` and chord fraction ``s``.

    A ply stack implements this; geometry only consumes metres. Span-only
    gauges wrap a :class:`SpanLaw` with :meth:`from_span` (``s`` ignored).
    """

    @abstractmethod
    def __call__(self, z, s) -> np.ndarray:
        """Broadcast ``z`` and ``s``; return thickness in metres."""

    @classmethod
    def from_span(cls, law: SpanLaw) -> ThicknessLaw:
        return _SpanThickness(law)


class _SpanThickness(ThicknessLaw):
    def __init__(self, law: SpanLaw) -> None:
        self.law = law

    def __call__(self, z, s) -> np.ndarray:
        zz, ss = np.broadcast_arrays(np.asarray(z, dtype=np.float64), np.asarray(s, dtype=np.float64))
        return np.reshape(self.law(zz.ravel()), zz.shape)


class ConstantLaw(SpanLaw):
    """The same value everywhere."""

    def __init__(self, value: float) -> None:
        self.value = float(value)

    def _eval(self, z: np.ndarray) -> np.ndarray:
        return np.full(np.shape(z), self.value, dtype=np.float64)

    def __repr__(self) -> str:
        return f"ConstantLaw({self.value!r})"


class TableLaw(SpanLaw):
    """Linear interpolation of knots, held constant beyond the first and last knot.

    Neighbouring equal values are plateaus. That is the ply-drop shape: constant,
    then a linear ramp, then constant.
    """

    def __init__(self, z: Sequence[float], values: Sequence[float]) -> None:
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
        self.z, self.values = zk, vk

    def _eval(self, z: np.ndarray) -> np.ndarray:
        zc = np.clip(np.asarray(z, dtype=np.float64), self.z[0], self.z[-1])
        if self.z.size == 1:
            return np.full(zc.shape, self.values[0])
        return np.interp(zc, self.z, self.values)

    def __repr__(self) -> str:
        return f"TableLaw(z={self.z.tolist()}, values={self.values.tolist()})"


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
