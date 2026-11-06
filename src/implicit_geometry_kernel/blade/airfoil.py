"""Airfoils: chord-normalised upper and lower surfaces.

An ``Airfoil`` is built from coordinates, as a NACA 4-digit section or as an
ellipse, and answers surface and thickness queries at chord fractions. Axes: x
along the chord (leading edge to trailing edge), y through the thickness (suction
side positive).
"""

from __future__ import annotations

from typing import Sequence, Union

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
        from implicit_geometry_kernel.sdf.shapes.naca import naca_four_digit_vertices

        # ``n`` is points per surface; the shared helper takes the total around the contour
        verts = naca_four_digit_vertices(m, p, t, 2 * int(n) - 1, 1.0, closed_te=closed_te)
        return cls(verts, f"NACA {code}")

    @classmethod
    def ellipse(cls, thickness_ratio: float = 1.0, n: int = 81) -> "Airfoil":
        """Elliptical section (a circle for ``thickness_ratio = 1``), as used at blade roots."""
        phi = np.linspace(0.0, 2.0 * np.pi, 2 * n - 1)
        x = 0.5 * (1.0 + np.cos(phi))
        y = 0.5 * thickness_ratio * np.sin(phi)
        return cls(np.column_stack([x, y]), f"ellipse t/c={thickness_ratio:g}")  # type: ignore

    def surface(self, side: str, x_c: Union[float, np.ndarray]) -> np.ndarray:
        """y/c of the ``"upper"`` or ``"lower"`` surface at chord fractions ``x_c``."""
        pts = self._upper if side == "upper" else self._lower
        x = self._x0 + np.asarray(x_c, dtype=float) * self._span
        return np.interp(x, pts[:, 0], pts[:, 1])

    def thickness(self, x_c) -> np.ndarray:
        return self.surface("upper", x_c) - self.surface("lower", x_c)

    def contour_points(self, n: int = 161) -> np.ndarray:
        """Closed-contour vertices in Selig order, chord-normalised, shape ``(N, 2)``."""
        x = np.linspace(0.0, 1.0, int(n))
        upper = np.column_stack([x[::-1], self.surface("upper", x[::-1])])
        lower = np.column_stack([x[1:], self.surface("lower", x[1:])])
        return np.vstack([upper, lower])

    def to_sdf(self, *, chord: float = 1.0, n: int = 161):
        """Implicit polygon adapter at ``chord`` metres."""
        from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF

        return AirfoilSDF.from_airfoil(self, chord=chord, n=n)

    @property
    def max_thickness(self) -> float:
        x = np.linspace(0.0, 1.0, 401)
        return float(self.thickness(x).max())

    def __repr__(self) -> str:
        return f"Airfoil({self.name!r})"
