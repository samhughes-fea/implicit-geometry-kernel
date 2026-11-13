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
    if isinstance(value, bool):
        raise TypeError("pivot_fraction is required (a chord fraction or span table)")
    if isinstance(value, (int, float)):
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


_MORPH_RULES = ("auto", "smoothstep", "thickness", "linear")


def _morph_surfaces(fa, fb, w: float, x: np.ndarray, rule: str = "auto") -> tuple[np.ndarray, np.ndarray]:
    """TE-preserving morph between two station airfoils.

    ``rule="auto"`` (station-table blades) uses a smoothstep blend when one
    parent is at least 60 % thick or named like an ellipse and the other is
    not, a linear blend when both are, and ``Airfoil.interpolate_by_thickness``
    otherwise. ``"smoothstep"``, ``"thickness"`` and ``"linear"`` force that
    one rule on every interval, so a consumer does not have to probe quarter
    stations to see which rule fired (G7).

    Surfaces are blended where they are sampled (``x``), from the parents' own
    splines. No intermediate airfoil is refitted, so a loft between equal
    airfoils reproduces them exactly at every ``z``.
    """
    w = float(np.clip(w, 0.0, 1.0))
    x = np.asarray(x, float)
    if rule not in _MORPH_RULES:
        raise ValueError(f"airfoil_morph must be one of {', '.join(_MORPH_RULES)}; got {rule!r}")
    if w <= 1e-15:
        return np.asarray(fa.surface("upper", x), float), np.asarray(fa.surface("lower", x), float)
    if w >= 1.0 - 1e-15:
        return np.asarray(fb.surface("upper", x), float), np.asarray(fb.surface("lower", x), float)

    def lerp(s: float, xs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        yu = (1.0 - s) * np.asarray(fa.surface("upper", xs), float) + s * np.asarray(fb.surface("upper", xs), float)
        yl = (1.0 - s) * np.asarray(fa.surface("lower", xs), float) + s * np.asarray(fb.surface("lower", xs), float)
        return yu, yl

    if rule == "linear":
        upper, lower = lerp(w, x)
    elif rule == "smoothstep":
        upper, lower = lerp(_smoothstep(w), x)
    elif rule == "thickness" or not (_near_circle(fa) or _near_circle(fb)):
        # Lerp the surfaces in t/c, then scale thickness about the camber line so the
        # peak t/c is the lerped one (Airfoil.interpolate_by_thickness, evaluated at x).
        # Peaks come from the sampled surfaces themselves, not Airfoil.max_thickness
        # (analytic for NACA sections), so equal parents give a scale of exactly 1.
        dense = np.linspace(0.0, 1.0, 401)  # the grid Airfoil.max_thickness uses
        ta = float(np.max(np.asarray(fa.surface("upper", dense)) - np.asarray(fa.surface("lower", dense))))
        tb = float(np.max(np.asarray(fb.surface("upper", dense)) - np.asarray(fb.surface("lower", dense))))
        t_c = (1.0 - w) * ta + w * tb
        s = 0.0 if abs(tb - ta) < 1e-15 else (t_c - ta) / (tb - ta)
        du, dl = lerp(s, dense)
        peak = float(np.max(du - dl))
        scale = t_c / peak if peak > 0.0 else 1.0
        yu, yl = lerp(s, x)
        camber, thick = 0.5 * (yu + yl), yu - yl
        upper, lower = camber + 0.5 * thick * scale, camber - 0.5 * thick * scale
    else:
        # auto, and the parents differ about the 60 % circle test
        s = _smoothstep(w) if _near_circle(fa) != _near_circle(fb) else w
        upper, lower = lerp(s, x)
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


def _metre_rows(raw, name: str):
    """One float per web, or ``(a, b)`` pairs (two feet, or a cap window)."""
    if raw is None:
        return None
    if isinstance(raw, (int, float)) and not isinstance(raw, bool):
        return (float(raw),)
    rows = []
    for item in raw:
        if isinstance(item, (int, float)) and not isinstance(item, bool):
            rows.append(float(item))
        else:
            pair = tuple(float(v) for v in item)
            if len(pair) != 2:
                raise ValueError(f"{name} rows are a metre or an (a, b) pair")
            rows.append(pair)
    if not rows:
        raise ValueError(f"{name} is empty")
    if any(isinstance(r, tuple) for r in rows) and any(not isinstance(r, tuple) for r in rows):
        raise ValueError(f"{name} mixes metres and pairs")
    return tuple(rows)


def _metres_law(value, name: str, *, pairs: bool = False):
    """Constant metres or a callable ``f(z)``.

    Webs: a float per web, or nested ``(suction_m, pressure_m)`` pairs.
    Caps (``pairs=True``): each window is ``(x0, x1)`` metres from the pitch
    axis; a single ``(x0, x1)`` is one window. Span variation is a callable,
    because a two-number sequence is already one window or two webs.
    """
    if value is None:
        return None
    if callable(value):
        return value
    if pairs and isinstance(value, (tuple, list)) and len(value) == 2 and all(
        isinstance(v, (int, float)) and not isinstance(v, bool) for v in value
    ):
        value = (value,)
    frozen = _metre_rows(value, name)
    return lambda z, _w=frozen: _w


def _is_frac_sequence(value) -> bool:
    try:
        vals = [float(v) for v in value]
    except (TypeError, ValueError):
        return False
    return bool(vals) and all(0.0 < v < 1.0 for v in vals)


def _with_web_fractions(spec, fractions):
    import dataclasses

    return dataclasses.replace(spec, web_chord_fracs=tuple(fractions))
