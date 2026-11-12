"""Airfoils: chord-normalised upper and lower surfaces.

An ``Airfoil`` is built from coordinates, as a NACA 4-digit section or as an
ellipse, and answers surface and thickness queries at chord fractions. Axes: x
along the chord (leading edge to trailing edge), y through the thickness (suction
side positive).
"""

from __future__ import annotations

from pathlib import Path
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

    def __init__(self, points: Sequence[Sequence[float]] | np.ndarray, name: str = "airfoil"):
        p = np.asarray(points, dtype=float)
        if p.ndim != 2 or p.shape[1] != 2 or len(p) < 5:
            raise ValueError("an airfoil needs at least 5 (x, y) points")
        self.name = name
        self._analytic_t: float | None = None
        self._closed_te: bool = True
        te = 0.5 * (p[0] + p[-1])
        le = int(np.argmax(np.linalg.norm(p - te, axis=1)))
        first, second = p[: le + 1][::-1], p[le:]
        if np.mean(first[:, 1]) < np.mean(second[:, 1]):
            first, second = second, first
        closed = bool(np.linalg.norm(p[0] - p[-1]) <= 1e-12)
        self._upper = self._monotone(first, closed)
        self._lower = self._monotone(second, closed)
        x0 = min(self._upper[0, 0], self._lower[0, 0])
        x1 = max(self._upper[-1, 0], self._lower[-1, 0])
        self._x0, self._span = x0, x1 - x0

    @staticmethod
    def _monotone(side: np.ndarray, closed_te: bool = False) -> np.ndarray:
        """Densely resampled surface, LE to TE, from a cubic spline through the points.

        The spline is parametrised by arc length, so blunt leading edges (where y(x)
        has a vertical tangent) are represented accurately.

        A round trailing end (the surfaces meet at one point, ``closed_te``, and this
        one arrives there steeper than 45 degrees: a circle, an ellipse, a root
        transition section) gets a vertical end tangent and ends at its end point. With
        a free end the spline bulges 1e-5 of chord aft of that point, where y(x) is not
        single-valued and sorting the samples by x scrambles the surface.
        """
        if side[0, 0] > side[-1, 0]:
            side = side[::-1]
        keep = np.concatenate([[True], np.linalg.norm(np.diff(side, axis=0), axis=1) > 1e-12])
        side = side[keep]
        t = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(side, axis=0), axis=1))])
        try:
            from scipy.interpolate import CubicSpline

            u = t[-1] * (1.0 - np.cos(np.linspace(0.0, 0.5 * np.pi, 4001)))  # fine at the LE end
            # the given points are kept as samples, so the surface passes through them:
            # the cosine samples are 0.04 % of chord apart at the TE end and would cut a
            # corner sharper than that (the flare into a blunt trailing edge)
            u = np.union1d(u, t)
            last = side[-1] - side[-2]
            if closed_te and abs(last[1]) > abs(last[0]):
                ends = ("not-a-knot", (1, 0.0)), ("not-a-knot", (1, float(np.sign(last[1]))))
            else:
                ends = "not-a-knot", "not-a-knot"
            dense = np.column_stack([CubicSpline(t, side[:, 0], bc_type=ends[0])(u),
                                     CubicSpline(t, side[:, 1], bc_type=ends[1])(u)])
            if ends[0] != "not-a-knot":
                # y(x) cannot hold what little of the round end still lies aft of its
                # end point (under 1e-6 of chord): the surface is cut at that x
                inside = dense[:, 0] < side[-1, 0]
                inside[-1] = True
                dense = dense[inside]
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
        foil = cls(verts, f"NACA {code}")
        foil._analytic_t = t
        foil._closed_te = closed_te
        return foil

    @classmethod
    def naca5(cls, code: str, n: int = 81, closed_te: bool = True) -> "Airfoil":
        """NACA 5-digit section (for example ``"23012"``)."""
        code = str(code).strip()
        if len(code) != 5 or not code.isdigit():
            raise ValueError("a NACA 5-digit code has five digits, e.g. '23012'")
        from implicit_geometry_kernel.sdf.shapes.naca import naca_five_digit_vertices

        L, pq, tt = int(code[0]), int(code[1:3]), int(code[3:])
        verts = naca_five_digit_vertices(L, pq, tt, 2 * int(n) - 1, 1.0, closed_te=closed_te)
        return cls(verts, f"NACA {code}")

    @classmethod
    def naca6(cls, family: int, t_percent: float, n: int = 81, closed_te: bool = True) -> "Airfoil":
        """NACA 6-series family 63 or 64, thickness as a percent of chord."""
        from implicit_geometry_kernel.sdf.shapes.naca import naca_six_series_vertices

        verts = naca_six_series_vertices(int(family), float(t_percent), 2 * int(n) - 1, 1.0,
                                         closed_te=closed_te)
        return cls(verts, f"NACA {int(family)}-{int(t_percent):02d}")

    @classmethod
    def cst(cls, upper, lower, n: int = 81, te_thickness: float = 0.0, name: str = "CST") -> "Airfoil":
        """Kulfan CST section from upper and lower Bernstein weights (round LE, class ``√x(1-x)``)."""
        from math import comb

        au = np.asarray(upper, float).reshape(-1)
        al = np.asarray(lower, float).reshape(-1)
        if au.size < 1 or al.size < 1:
            raise ValueError("CST needs at least one weight per surface")
        x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, int(n))))
        klass = np.sqrt(np.maximum(x, 0.0)) * (1.0 - x)
        nu, nl = au.size - 1, al.size - 1
        su = sum(au[i] * comb(nu, i) * x**i * (1.0 - x) ** (nu - i) for i in range(au.size))
        sl = sum(al[i] * comb(nl, i) * x**i * (1.0 - x) ** (nl - i) for i in range(al.size))
        yu = klass * su + x * (0.5 * float(te_thickness))
        yl = -(klass * sl + x * (0.5 * float(te_thickness)))
        upper_xy = np.column_stack([x[::-1], yu[::-1]])
        lower_xy = np.column_stack([x[1:], yl[1:]])
        return cls(np.vstack([upper_xy, lower_xy]), name)

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
        t_frac = getattr(self, "_analytic_t", None)
        if t_frac is not None:
            from implicit_geometry_kernel.sdf.shapes.naca import naca_four_digit_thickness

            return naca_four_digit_thickness(x_c, t_frac, closed_te=getattr(self, "_closed_te", True))
        return self.surface("upper", x_c) - self.surface("lower", x_c)

    def with_thickness(self, t_c: float, n: int = 161) -> "Airfoil":
        """Same camber line, thickness scaled about it so max t/c is ``t_c``."""
        if t_c <= 0.0:
            raise ValueError("t/c must be positive")
        x = np.linspace(0.0, 1.0, int(n))
        yu, yl = self.surface("upper", x), self.surface("lower", x)
        camber = 0.5 * (yu + yl)
        thick = yu - yl
        peak = float(np.max(thick))
        if peak <= 0.0:
            raise ValueError("airfoil has no thickness to scale")
        scale = float(t_c) / peak
        yu2, yl2 = camber + 0.5 * thick * scale, camber - 0.5 * thick * scale
        return type(self)(self._selig_from_surfaces(x, yu2, yl2), f"{self.name} t/c={t_c:g}")

    def with_te_thickness(self, te: float, n: int = 161) -> "Airfoil":
        """Flatback: ramp trailing-edge thickness to ``te`` (chord units) from the LE."""
        if te < 0.0:
            raise ValueError("trailing-edge thickness must be >= 0")
        x = np.linspace(0.0, 1.0, int(n))
        yu, yl = self.surface("upper", x), self.surface("lower", x)
        add = 0.5 * (float(te) - float(yu[-1] - yl[-1]))
        return type(self)(self._selig_from_surfaces(x, yu + x * add, yl - x * add),
                          f"{self.name} TE={te:g}")

    @classmethod
    def blend(cls, a: "Airfoil", b: "Airfoil", s: float, n: int = 161, name: str | None = None) -> "Airfoil":
        """Linear surface blend; ``s = 0`` is ``a`` and ``s = 1`` is ``b``."""
        w = float(np.clip(s, 0.0, 1.0))
        x = np.linspace(0.0, 1.0, int(n))
        yu = (1.0 - w) * a.surface("upper", x) + w * b.surface("upper", x)
        yl = (1.0 - w) * a.surface("lower", x) + w * b.surface("lower", x)
        return cls(cls._selig_from_surfaces(x, yu, yl), name or f"{a.name}->{b.name} s={w:g}")

    @classmethod
    def cylinder_to_airfoil(cls, airfoil: "Airfoil", s: float, n: int = 161) -> "Airfoil":
        """Named root law: ``s = 0`` is a circle matching t/c, ``s = 1`` is ``airfoil``."""
        root = cls.ellipse(airfoil.max_thickness, n=n)
        return cls.blend(root, airfoil, s, n=n, name=f"cylinder→{airfoil.name} s={s:g}")

    @classmethod
    def interpolate_by_thickness(cls, a: "Airfoil", b: "Airfoil", t_c: float, n: int = 161) -> "Airfoil":
        """Blend two foils in t/c: surfaces lerp by thickness, then rescale to ``t_c``."""
        ta, tb = a.max_thickness, b.max_thickness
        w = 0.0 if abs(tb - ta) < 1e-15 else (float(t_c) - ta) / (tb - ta)
        x = np.linspace(0.0, 1.0, int(n))
        yu = (1.0 - w) * a.surface("upper", x) + w * b.surface("upper", x)
        yl = (1.0 - w) * a.surface("lower", x) + w * b.surface("lower", x)
        blended = cls(cls._selig_from_surfaces(x, yu, yl), f"{a.name}+{b.name}")
        return blended.with_thickness(t_c, n=n)

    @classmethod
    def from_selig(cls, source, name: str | None = None) -> "Airfoil":
        """Read Selig order (TE → upper → LE → lower → TE) from a ``.dat`` path or text."""
        text, label = _dat_text(source, name)
        pts = _parse_xy_pairs(text)
        return cls(pts, label or "selig")

    @classmethod
    def from_lednicer(cls, source, name: str | None = None) -> "Airfoil":
        """Read Lednicer order (upper LE→TE, then lower LE→TE) from a ``.dat`` path or text."""
        text, label = _dat_text(source, name)
        pts = _parse_lednicer(text)
        return cls(pts, label or "lednicer")

    @staticmethod
    def _selig_from_surfaces(x, yu, yl) -> np.ndarray:
        upper = np.column_stack([x[::-1], yu[::-1]])
        lower = np.column_stack([x[1:], yl[1:]])
        return np.vstack([upper, lower])

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


def _dat_text(source, name: str | None) -> tuple[str, str]:
    if isinstance(source, (str, Path)) and Path(source).is_file():
        path = Path(source)
        return path.read_text(encoding="utf-8"), name or path.stem
    return str(source), name or "airfoil"


def _parse_xy_pairs(text: str) -> np.ndarray:
    rows = []
    for line in text.splitlines():
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        rows.append((x, y))
    if len(rows) < 5:
        raise ValueError("airfoil file needs at least 5 coordinate pairs")
    return np.asarray(rows, float)


def _parse_lednicer(text: str) -> np.ndarray:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    counts = None
    start = 0
    for i, ln in enumerate(lines):
        parts = ln.replace(",", " ").split()
        if len(parts) == 2:
            try:
                nu, nl = int(float(parts[0])), int(float(parts[1]))
            except ValueError:
                continue
            if nu >= 3 and nl >= 3:
                counts = (nu, nl)
                start = i + 1
                break
    pts = _parse_xy_pairs("\n".join(lines[start:]))
    if counts is not None:
        nu, nl = counts
        if len(pts) < nu + nl:
            raise ValueError("Lednicer file has fewer points than the header counts")
        upper, lower = pts[:nu], pts[nu:nu + nl]
    else:
        le = int(np.argmin(pts[:, 0]))
        if le == 0:
            mid = len(pts) // 2
            upper, lower = pts[:mid], pts[mid:]
        else:
            return pts
    upper = upper[np.argsort(upper[:, 0])]
    lower = lower[np.argsort(lower[:, 0])]
    return np.vstack([upper[::-1], lower[1:]])
