"""Span laws: a scalar as a vectorised function of the spanwise coordinate ``z``.

Member dimensions (thicknesses, depths, positions) are laws so that one assembled
structure describes the whole span. Laws are small picklable objects, never lambdas,
so a structure can be sent to worker processes. All lengths are in metres.

:class:`TableLaw` interpolates knots linearly and holds the end values outside the
table. That matches ply-drop gauges: constant, then a linear ramp, then constant.

:class:`PchipLaw` is the smooth alternative for authored laws (decision 0014): a
shape-preserving cubic (C1, no overshoot past neighbouring knots, flat plateaus kept
flat) with analytic first and second derivatives. A PCHIP loft adds curvature between
knots, so it needs denser knots than a linear one; :func:`check_pchip_sampling` tests
that with a leave-one-out prediction and :class:`PchipLaw` refuses knots that fail.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


class SpanLaw(ABC):
    """A scalar as a function of ``z``. Subclasses implement :meth:`_eval`."""

    def __call__(self, z: ArrayLike) -> np.ndarray:
        zz = np.asarray(z, dtype=np.float64)
        if zz.size > 1:
            flat = zz.ravel()
            if np.min(flat) == np.max(flat):  # one section plane, the common case
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
    def __call__(self, z: ArrayLike, s: ArrayLike) -> np.ndarray:
        """Broadcast ``z`` and ``s``; return thickness in metres."""

    @classmethod
    def from_span(cls, law: SpanLaw) -> ThicknessLaw:
        return _SpanThickness(law)


class _SpanThickness(ThicknessLaw):
    def __init__(self, law: SpanLaw) -> None:
        self.law = law

    def __call__(self, z: ArrayLike, s: ArrayLike) -> np.ndarray:
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


class UnderSampledSpanLaw(ValueError):
    """PCHIP knots too sparse to resolve the law (decision 0014).

    ``report`` is the :class:`SamplingReport` that failed. The fix is more knots
    (stations), not a looser tolerance.
    """

    def __init__(self, message: str, report: "SamplingReport") -> None:
        super().__init__(message)
        self.report = report


@dataclass(frozen=True)
class SamplingReport:
    """Leave-one-out check of a PCHIP law's knots.

    Attributes
    ----------
    name:
        The quantity checked (for messages), for example ``"chord"``.
    z:
        Knot positions, increasing.
    loo_error:
        Absolute error at each knot when that knot is dropped and predicted by
        PCHIP through the others. End knots and laws with fewer than three knots
        have zero error (nothing to predict).
    tolerance:
        The absolute tolerance used, in the quantity's own units.
    bay_gap:
        Largest ``|pchip - linear|`` inside each bay (``len(z) - 1`` values): the
        geometry PCHIP adds over a linear loft, in the quantity's units.
    """

    name: str
    z: np.ndarray
    loo_error: np.ndarray
    tolerance: float
    bay_gap: np.ndarray

    @property
    def ok(self) -> bool:
        """True when every leave-one-out error is within the tolerance."""
        return bool(np.all(self.loo_error <= self.tolerance))

    @property
    def worst(self) -> int:
        """Index of the knot with the largest leave-one-out error."""
        return int(np.argmax(self.loo_error)) if self.loo_error.size else 0

    def message(self) -> str:
        """One line naming the worst knot, its error and the tolerance."""
        i = self.worst
        return (
            f"{self.name}: PCHIP knots are too sparse at z = {self.z[i]:g} "
            f"(leave-one-out error {self.loo_error[i]:.3g} > tolerance {self.tolerance:.3g}); "
            "add stations near there or use linear interpolation"
        )


def _sorted_knots(z: Sequence[float], values: Sequence[float], name: str) -> tuple[np.ndarray, np.ndarray]:
    zk = np.asarray(z, dtype=np.float64)
    vk = np.asarray(values, dtype=np.float64)
    if zk.ndim != 1 or zk.shape != vk.shape or zk.size == 0:
        raise ValueError(f"{name} needs matching 1-D knot arrays")
    if not (np.all(np.isfinite(zk)) and np.all(np.isfinite(vk))):
        raise ValueError(f"{name} knots must be finite")
    order = np.argsort(zk, kind="stable")
    zk, vk = zk[order], vk[order]
    if np.any(np.diff(zk) <= 0.0):
        raise ValueError(f"{name} knot positions must be distinct")
    return zk, vk


def pchip_slopes(z: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Knot slopes of the Fritsch–Carlson shape-preserving cubic.

    Interior slopes are the weighted harmonic mean of the neighbouring secants,
    zero where the secants change sign or one is zero (a plateau or extremum
    stays flat). End slopes are the one-sided three-point estimate, limited so
    the end bay stays monotone. Two knots give the secant (a straight line).
    """
    z = np.asarray(z, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    n = z.size
    if n == 1:
        return np.zeros(1)
    h = np.diff(z)
    delta = np.diff(v) / h
    if n == 2:
        return np.array([delta[0], delta[0]])
    d = np.zeros(n)
    h0, h1 = h[:-1], h[1:]
    d0, d1 = delta[:-1], delta[1:]
    w1 = 2.0 * h1 + h0
    w2 = h1 + 2.0 * h0
    same = (d0 * d1) > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        harm = (w1 + w2) / (w1 / d0 + w2 / d1)
    d[1:-1] = np.where(same, harm, 0.0)
    d[0] = _pchip_end(h[0], h[1], delta[0], delta[1])
    d[-1] = _pchip_end(h[-1], h[-2], delta[-1], delta[-2])
    return d


def _pchip_end(h0: float, h1: float, m0: float, m1: float) -> float:
    d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    if np.sign(d) != np.sign(m0):
        return 0.0
    if np.sign(m0) != np.sign(m1) and abs(d) > abs(3.0 * m0):
        return 3.0 * m0
    return float(d)


def _hermite(z: np.ndarray, zk: np.ndarray, vk: np.ndarray, dk: np.ndarray, order: int) -> np.ndarray:
    """Cubic Hermite value (``order`` 0), slope (1) or curvature (2) at ``z`` inside the knots."""
    k = np.clip(np.searchsorted(zk, z, side="right") - 1, 0, zk.size - 2)
    h = zk[k + 1] - zk[k]
    t = (z - zk[k]) / h
    y0, y1 = vk[k], vk[k + 1]
    m0, m1 = dk[k] * h, dk[k + 1] * h
    if order == 0:
        t2, t3 = t * t, t * t * t
        return (
            (2 * t3 - 3 * t2 + 1) * y0 + (t3 - 2 * t2 + t) * m0
            + (-2 * t3 + 3 * t2) * y1 + (t3 - t2) * m1
        )
    if order == 1:
        t2 = t * t
        return (
            (6 * t2 - 6 * t) * y0 + (3 * t2 - 4 * t + 1) * m0
            + (-6 * t2 + 6 * t) * y1 + (3 * t2 - 2 * t) * m1
        ) / h
    return ((12 * t - 6) * y0 + (6 * t - 4) * m0 + (-12 * t + 6) * y1 + (6 * t - 2) * m1) / (h * h)


def _pchip_eval(z: np.ndarray, zk: np.ndarray, vk: np.ndarray, order: int = 0) -> np.ndarray:
    zz = np.asarray(z, dtype=np.float64)
    if zk.size == 1:
        return np.full(zz.shape, vk[0] if order == 0 else 0.0)
    dk = pchip_slopes(zk, vk)
    zc = np.clip(zz, zk[0], zk[-1])
    out = _hermite(zc, zk, vk, dk, order)
    if order > 0:
        out = np.where((zz < zk[0]) | (zz > zk[-1]), 0.0, out)
    return out


def check_pchip_sampling(
    z: Sequence[float],
    values: Sequence[float],
    *,
    name: str = "law",
    rtol: float = 5e-3,
    atol: float = 0.0,
    samples_per_bay: int = 16,
) -> SamplingReport:
    """Leave-one-out test of whether knots are dense enough for PCHIP.

    Each interior knot is dropped and predicted by PCHIP through the others. A
    knot the others cannot predict within the tolerance marks a shape the knots
    do not resolve, so PCHIP would invent geometry there. The test is
    conservative: dropping a knot doubles the local spacing.

    Parameters
    ----------
    z, values:
        Knots (any order; sorted here).
    name:
        Quantity name used in messages.
    rtol, atol:
        Tolerance ``atol + rtol * (max(values) - min(values))`` in the
        quantity's own units. The default is 0.5 % of the range over the span.
    samples_per_bay:
        Points per bay used to measure ``bay_gap``.

    Returns
    -------
    SamplingReport
        ``report.ok`` is False when any knot fails.
    """
    zk, vk = _sorted_knots(z, values, name)
    n = zk.size
    tol = float(atol) + float(rtol) * float(np.ptp(vk))
    loo = np.zeros(n)
    if n >= 3:
        for i in range(1, n - 1):
            keep = np.ones(n, dtype=bool)
            keep[i] = False
            loo[i] = abs(float(_pchip_eval(np.array([zk[i]]), zk[keep], vk[keep])[0]) - vk[i])
    if n >= 2:
        t = np.linspace(0.0, 1.0, max(int(samples_per_bay), 2) + 1)[1:-1]
        zs = (zk[:-1, None] + t[None, :] * np.diff(zk)[:, None])
        gap = np.abs(_pchip_eval(zs.ravel(), zk, vk) - np.interp(zs.ravel(), zk, vk)).reshape(zs.shape)
        bay_gap = gap.max(axis=1) if gap.size else np.zeros(n - 1)
    else:
        bay_gap = np.zeros(0)
    return SamplingReport(name=name, z=zk, loo_error=loo, tolerance=tol, bay_gap=bay_gap)


class PchipLaw(SpanLaw):
    """Shape-preserving cubic (PCHIP) through knots, held constant beyond the ends.

    C1 with analytic derivatives (:meth:`derivative`). No overshoot past
    neighbouring knots, and neighbouring equal values stay a flat plateau. On
    construction the knots are checked with :func:`check_pchip_sampling`; knots
    too sparse for PCHIP raise :class:`UnderSampledSpanLaw` unless
    ``sampling="flag"``, which keeps the report on :attr:`sampling` instead.

    Parameters
    ----------
    z, values:
        Knots.
    name:
        Quantity name for messages.
    rtol, atol:
        Sampling tolerance (see :func:`check_pchip_sampling`).
    sampling:
        ``"raise"`` (default) or ``"flag"``.
    """

    def __init__(
        self,
        z: Sequence[float],
        values: Sequence[float],
        *,
        name: str = "PchipLaw",
        rtol: float = 5e-3,
        atol: float = 0.0,
        sampling: str = "raise",
    ) -> None:
        if sampling not in ("raise", "flag"):
            raise ValueError("sampling must be 'raise' or 'flag'")
        self.z, self.values = _sorted_knots(z, values, name)
        self.slopes = pchip_slopes(self.z, self.values)
        self.name = str(name)
        self.sampling = check_pchip_sampling(
            self.z.tolist(), self.values.tolist(), name=name, rtol=rtol, atol=atol
        )
        if sampling == "raise" and not self.sampling.ok:
            raise UnderSampledSpanLaw(self.sampling.message(), self.sampling)

    def _eval(self, z: np.ndarray) -> np.ndarray:
        return _pchip_eval(z, self.z, self.values)

    def derivative(self, z: ArrayLike, order: int = 1) -> np.ndarray:
        """First (``order=1``) or second (``order=2``) derivative in z; zero outside the knots.

        The second derivative is piecewise linear and may jump at knots (C1, not C2).
        """
        if order not in (1, 2):
            raise ValueError("order must be 1 or 2")
        return _pchip_eval(np.asarray(z, dtype=np.float64), self.z, self.values, order)

    def __repr__(self) -> str:
        return f"PchipLaw(z={self.z.tolist()}, values={self.values.tolist()}, name={self.name!r})"


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
    "ConstantLaw", "PchipLaw", "ProductLaw", "SamplingReport", "ScaledLaw", "SpanLaw", "SumLaw",
    "TableLaw", "ThicknessLaw", "UnderSampledSpanLaw", "as_law", "check_pchip_sampling",
    "pchip_slopes",
]
