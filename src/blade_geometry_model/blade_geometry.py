"""Blade geometry: airfoils, spanwise stations, skin regions and shear webs.

These types describe the shape and the structural layout of a blade or wing, and
nothing about materials or meshes: an ``Airfoil`` (chord-normalised surfaces),
``Station`` (span position, chord, twist, pitch axis and offset), ``SkinRegion``
(a chordwise band of skin on both surfaces) and ``ShearWeb`` (a web at a chord
fraction). ``composite_laminate_model.BladeBuilder`` meshes them as shells, and the
section model cuts sections from them.

Axes: z runs along the span, x along the chord (leading edge to trailing edge) and
y through the thickness (suction side positive). The pitch axis sits at x = 0.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple, Union

import numpy as np


# ---------------------------------------------------------------------------
# Airfoils
# ---------------------------------------------------------------------------
class Airfoil:
    """Airfoil from coordinates normalised by the chord (x from 0 at the LE to 1 at the TE).

    ``points`` run in Selig order (TE, upper surface, LE, lower surface, TE) or in the
    reverse order; the surfaces are split at the leading edge (smallest x). The
    trailing edge can be sharp or blunt (a flatback).
    """

    def __init__(self, points: Sequence[Sequence[float]], name: str = "airfoil"):
        p = np.asarray(points, dtype=float)
        if p.ndim != 2 or p.shape[1] != 2 or len(p) < 5:
            raise ValueError("an airfoil needs at least 5 (x, y) points")
        self.name = name
        le = int(np.argmin(p[:, 0]))
        first, second = p[: le + 1][::-1], p[le:]
        if np.mean(first[:, 1]) < np.mean(second[:, 1]):
            first, second = second, first
        self._upper = self._monotone(first)
        self._lower = self._monotone(second)
        x0 = min(self._upper[0, 0], self._lower[0, 0])
        x1 = max(self._upper[-1, 0], self._lower[-1, 0])
        self._x0, self._span = x0, x1 - x0

    @staticmethod
    def _monotone(side: np.ndarray) -> np.ndarray:
        """Densely resampled surface, LE to TE, from a cubic spline through the points.

        The spline is parametrised by arc length, so blunt leading edges (where y(x)
        has a vertical tangent) are represented accurately.
        """
        if side[0, 0] > side[-1, 0]:
            side = side[::-1]
        keep = np.concatenate([[True], np.linalg.norm(np.diff(side, axis=0), axis=1) > 1e-12])
        side = side[keep]
        t = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(side, axis=0), axis=1))])
        try:
            from scipy.interpolate import CubicSpline

            u = t[-1] * (1.0 - np.cos(np.linspace(0.0, 0.5 * np.pi, 4001)))  # fine at the LE end
            dense = np.column_stack([CubicSpline(t, side[:, 0])(u), CubicSpline(t, side[:, 1])(u)])
        except ImportError:  # pragma: no cover - SciPy is an optional dependency
            dense = side
        order = np.argsort(dense[:, 0], kind="stable")
        dense = dense[order]
        keep = np.concatenate([[True], np.diff(dense[:, 0]) > 1e-14])
        return dense[keep]

    @classmethod
    def naca4(cls, code: str, n: int = 81, closed_te: bool = True) -> "Airfoil":
        """NACA 4-digit section (for example ``"0018"`` or ``"2412"``) with cosine spacing."""
        code = str(code)
        if len(code) != 4 or not code.isdigit():
            raise ValueError("a NACA 4-digit code has four digits, e.g. '0018'")
        m, p, t = int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0
        x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, n)))
        a4 = -0.1036 if closed_te else -0.1015
        yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 + a4 * x**4)
        if m > 0 and p > 0:
            yc = np.where(x < p, m / p**2 * (2 * p * x - x**2), m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x**2))
            dyc = np.where(x < p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
        else:
            yc, dyc = np.zeros_like(x), np.zeros_like(x)
        th = np.arctan(dyc)
        upper = np.column_stack([x - yt * np.sin(th), yc + yt * np.cos(th)])
        lower = np.column_stack([x + yt * np.sin(th), yc - yt * np.cos(th)])
        return cls(np.vstack([upper[::-1], lower[1:]]), f"NACA {code}")

    @classmethod
    def ellipse(cls, thickness_ratio: float = 1.0, n: int = 81) -> "Airfoil":
        """Elliptical section (a circle for ``thickness_ratio = 1``), as used at blade roots."""
        phi = np.linspace(0.0, 2.0 * np.pi, 2 * n - 1)
        x = 0.5 * (1.0 + np.cos(phi))
        y = 0.5 * thickness_ratio * np.sin(phi)
        return cls(np.column_stack([x, y]), f"ellipse t/c={thickness_ratio:g}")

    def surface(self, side: str, x_c: Union[float, np.ndarray]) -> np.ndarray:
        """y/c of the ``"upper"`` or ``"lower"`` surface at chord fractions ``x_c``."""
        pts = self._upper if side == "upper" else self._lower
        x = self._x0 + np.asarray(x_c, dtype=float) * self._span
        return np.interp(x, pts[:, 0], pts[:, 1])

    def thickness(self, x_c) -> np.ndarray:
        return self.surface("upper", x_c) - self.surface("lower", x_c)

    @property
    def max_thickness(self) -> float:
        x = np.linspace(0.0, 1.0, 401)
        return float(self.thickness(x).max())

    def __repr__(self) -> str:
        return f"Airfoil({self.name!r})"


# ---------------------------------------------------------------------------
# Blade definition
# ---------------------------------------------------------------------------
@dataclass
class Station:
    """A spanwise station: position ``z``, chord, twist (degrees, nose up about +z) and airfoil.

    ``pitch_axis`` is the chord fraction that lies on the blade axis (x = 0), and
    ``offset`` an (x, y) shift of that point (prebend, sweep).
    """

    z: float
    chord: float
    airfoil: Airfoil
    twist: float = 0.0
    pitch_axis: float = 0.25
    offset: Tuple[float, float] = (0.0, 0.0)

    def __post_init__(self):
        if self.chord <= 0.0:
            raise ValueError("chord must be positive")


@dataclass(frozen=True)
class SkinRegion:
    """Skin between chord fractions ``start`` and ``end`` on both surfaces, with ``n`` elements."""

    name: str
    start: float
    end: float
    n: int = 4

    def __post_init__(self):
        if not 0.0 <= self.start < self.end <= 1.0:
            raise ValueError(f"region '{self.name}' needs 0 <= start < end <= 1")
        if self.n < 1:
            raise ValueError("a region needs at least one element")


@dataclass(frozen=True)
class ShearWeb:
    """A flat web joining the surfaces at chord fraction ``position`` with ``n`` elements."""

    name: str
    position: float
    n: int = 4
