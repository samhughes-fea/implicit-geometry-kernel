"""Span tables and airfoil morph helpers for BladeGeometry."""
from __future__ import annotations

from typing import Callable, Optional, Sequence, Union

import numpy as np

Thickness = Union[float, Callable[[float], float]]
SpanTable = Union[Callable, Sequence]
_FORK_TWIST_SIGN = -1.0


def _as_function(value: Thickness, name: str) -> Callable[[float], float]:
    if callable(value):
        return value
    v = float(value)
    if v < 0.0:
        raise ValueError(f"{name} must be >= 0")
    return lambda z: v


def _as_pivot_fraction(value: Union[float, SpanTable]):
    """Required piercing-point fraction: a constant in (0, 1) or a span table."""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        f = float(value)
        if not 0.0 < f < 1.0:
            raise ValueError("pivot_fraction must be a chord fraction in (0, 1)")
        return lambda z, _f=f: _f
    table = _as_span_table(value, "pivot_fraction")
    if table is None:
        raise TypeError("pivot_fraction is required (a chord fraction or span table)")
    return table


def _as_span_table(value: Optional[SpanTable], name: str, *, cols: int = 1):
    """Linear lookup ``f(z)``. ``None`` means use the bounding stations."""
    if value is None:
        return None
    if callable(value):
        return value
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError(f"{name} must be a callable or a (z, values) pair")
    z = np.asarray(value[0], float).reshape(-1)
    v = np.asarray(value[1], float)
    if z.size < 2 or np.any(np.diff(z) <= 0.0):
        raise ValueError(f"{name} needs strictly increasing z of length >= 2")
    if v.shape[0] != z.size:
        raise ValueError(f"{name} values must match z")
    if cols == 1:
        v = np.asarray(v, float).reshape(z.size)

        def lookup(zz: float) -> float:
            return float(np.interp(zz, z, v))

        return lookup
    if v.ndim != 2 or v.shape[1] != cols:
        raise ValueError(f"{name} values must have shape (n, {cols})")

    def lookup_vec(zz: float):
        return tuple(float(np.interp(zz, z, v[:, i])) for i in range(cols))

    return lookup_vec


def _smoothstep(w: float) -> float:
    w = float(np.clip(w, 0.0, 1.0))
    return w * w * (3.0 - 2.0 * w)


def _near_circle(foil) -> bool:
    name = str(getattr(foil, "name", "")).lower()
    return float(foil.max_thickness) >= 0.6 or "ellipse" in name


def _morph_surfaces(fa, fb, w: float, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """TE-preserving morph: smoothstep blend on circle roots, else thickness blend."""
    from implicit_geometry_kernel.blade.airfoil import Airfoil

    w = float(np.clip(w, 0.0, 1.0))
    n = max(int(len(x)), 5)
    if w <= 1e-15:
        return np.asarray(fa.surface("upper", x), float), np.asarray(fa.surface("lower", x), float)
    if w >= 1.0 - 1e-15:
        return np.asarray(fb.surface("upper", x), float), np.asarray(fb.surface("lower", x), float)
    if _near_circle(fa) or _near_circle(fb):
        s = _smoothstep(w) if _near_circle(fa) != _near_circle(fb) else w
        foil = Airfoil.blend(fa, fb, s, n=n)
    else:
        t_c = (1.0 - w) * float(fa.max_thickness) + w * float(fb.max_thickness)
        foil = Airfoil.interpolate_by_thickness(fa, fb, t_c, n=n)
    upper = np.asarray(foil.surface("upper", x), float)
    lower = np.asarray(foil.surface("lower", x), float)
    gap_a = abs(float(np.asarray(fa.surface("upper", 1.0)) - np.asarray(fa.surface("lower", 1.0))))
    gap_b = abs(float(np.asarray(fb.surface("upper", 1.0)) - np.asarray(fb.surface("lower", 1.0))))
    gap = (1.0 - w) * gap_a + w * gap_b
    mid = 0.5 * (float(upper[-1]) + float(lower[-1]))
    upper[-1] = mid + 0.5 * gap
    lower[-1] = mid - 0.5 * gap
    return upper, lower

def _match(polyline: np.ndarray, points: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Index of each of ``points`` among the vertices of ``polyline``."""
    P = np.asarray(polyline, float)
    out = []
    for p in np.asarray(points, float).reshape(-1, 2):
        d = np.linalg.norm(P - p, axis=1)
        k = int(np.argmin(d))
        if d[k] > tol * max(1.0, float(np.abs(P).max())):
            raise ValueError("a bond point is not a vertex of its strip")
        out.append(k)
    return np.array(out, dtype=int)


def _is_frac_sequence(value) -> bool:
    try:
        vals = [float(v) for v in value]
    except (TypeError, ValueError):
        return False
    return bool(vals) and all(0.0 < v < 1.0 for v in vals)


def _with_web_fractions(spec, fractions):
    import dataclasses

    return dataclasses.replace(spec, web_chord_fracs=tuple(fractions))
