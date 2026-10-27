"""Typed clip algebra for section subcomponents (IR + Shapely).

Public clip types (:class:`ClipSpec`) are the SoT for both SDF IR and
constructive Shapely interpreters. Camber-side splits live here — not in
component modules.

Low-level callables still implemented in :mod:`clips` are re-exported here so
production code imports only this module (see ``test_clip_import_boundary``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Literal, Union

import numpy as np
from shapely.geometry import Polygon

from blade_geometry_model.section_geometry.sdf.ir import (
    HalfPlane,
    Intersect,
)

from ..constructive.helpers import (
    chordwise_box,
    ensure_polygon,
    flapwise_bay_clip,
    half_plane_clip,
)
from .clips import (
    _airfoil_offset_arc,
    _parallel_web_half_gt,
    _parallel_web_half_lt,
    _parallel_web_strip_clip,
    airfoil_expr,
    chordwise_clip_expr,
    flapwise_web_normal,
    lower_half_plane_expr,
    strip_clip_expr,
    upper_half_plane_expr,
    web_station_projection,
)


@dataclass(frozen=True)
class HalfPlaneClip:
    nx: float
    ny: float
    d: float


@dataclass(frozen=True)
class ChordBand:
    x_start: float
    x_end: float


@dataclass(frozen=True)
class FlapwiseStrip:
    nx: float
    ny: float
    c_lo: float | None = None
    c_hi: float | None = None
    clip_mode: Literal["strip", "half_lt", "half_gt"] = "strip"


@dataclass(frozen=True)
class CamberSide:
    surface: Literal["upper", "lower"]


ClipSpec = Union[HalfPlaneClip, ChordBand, FlapwiseStrip, CamberSide]


def camber_side_field(camber_xy, surface: str) -> Callable:
    """Signed split at the cavity midline: negative on the cap's / skin's own side.

    Splitting at the chord line ``y=0`` lets ring-band material from the
    opposite skin leak across on cambered sections.
    """
    xs_c = np.asarray(camber_xy, dtype=float)[:, 0]
    ys_c = np.asarray(camber_xy, dtype=float)[:, 1]
    sign = -1.0 if surface == "upper" else 1.0

    def _split(x, y):
        xa = np.asarray(x, dtype=float)
        ya = np.asarray(y, dtype=float)
        cam = np.interp(xa, xs_c, ys_c)
        return sign * (ya - cam)

    return _split


def to_ir(clip: ClipSpec):
    """Map a typed clip to SDF IR using existing ``clips`` helpers where possible."""
    if isinstance(clip, HalfPlaneClip):
        return HalfPlane(nx=float(clip.nx), ny=float(clip.ny), d=float(clip.d))
    if isinstance(clip, ChordBand):
        return chordwise_clip_expr(float(clip.x_start), float(clip.x_end))
    if isinstance(clip, CamberSide):
        if clip.surface == "upper":
            return upper_half_plane_expr()
        if clip.surface == "lower":
            return lower_half_plane_expr()
        raise ValueError(
            f"CamberSide.surface must be 'upper' or 'lower', got {clip.surface!r}"
        )
    if isinstance(clip, FlapwiseStrip):
        mode = str(clip.clip_mode or "strip")
        if mode == "strip":
            if clip.c_lo is None or clip.c_hi is None:
                raise ValueError("FlapwiseStrip strip mode requires c_lo and c_hi.")
            # Inside when c_lo <= n·p <= c_hi (phi <= 0).
            nx_n, ny_n = float(clip.nx), float(clip.ny)
            lo = float(min(clip.c_lo, clip.c_hi))
            hi = float(max(clip.c_lo, clip.c_hi))
            return Intersect(
                HalfPlane(nx=-nx_n, ny=-ny_n, d=lo),
                HalfPlane(nx=nx_n, ny=ny_n, d=-hi),
            )
        if mode == "half_lt":
            if clip.c_lo is None:
                raise ValueError("FlapwiseStrip half_lt requires c_lo.")
            return HalfPlane(nx=-float(clip.nx), ny=-float(clip.ny), d=float(clip.c_lo))
        if mode == "half_gt":
            c = clip.c_hi if clip.c_hi is not None else clip.c_lo
            if c is None:
                raise ValueError("FlapwiseStrip half_gt requires c_hi or c_lo.")
            return HalfPlane(nx=float(clip.nx), ny=float(clip.ny), d=-float(c))
        raise ValueError(f"Unknown FlapwiseStrip.clip_mode {mode!r}")
    raise TypeError(f"Unsupported clip type: {type(clip)!r}")


def apply_shapely(clip: ClipSpec, polygon, *, camber_xy: Any = None):
    """Apply ``clip`` to a Shapely polygon via constructive helpers."""
    if polygon is None or getattr(polygon, "is_empty", False):
        return Polygon()
    if isinstance(clip, HalfPlaneClip):
        nx_n, ny_n, plane_d = float(clip.nx), float(clip.ny), float(clip.d)
        n_norm = float(np.hypot(nx_n, ny_n))
        if n_norm < 1e-15:
            raise ValueError("HalfPlaneClip normal must be non-zero.")
        bay = flapwise_bay_clip(
            polygon,
            nx=nx_n / n_norm,
            ny=ny_n / n_norm,
            clip_mode="half_lt",
            c_lo=-plane_d / n_norm,
            c_hi=None,
        )
        return ensure_polygon(polygon.intersection(bay), allow_empty=True)
    if isinstance(clip, ChordBand):
        box = chordwise_box(polygon, float(clip.x_start), float(clip.x_end))
        return ensure_polygon(polygon.intersection(box), allow_empty=True)
    if isinstance(clip, FlapwiseStrip):
        bay = flapwise_bay_clip(
            polygon,
            nx=float(clip.nx),
            ny=float(clip.ny),
            clip_mode=str(clip.clip_mode or "strip"),
            c_lo=clip.c_lo,
            c_hi=clip.c_hi,
        )
        return ensure_polygon(polygon.intersection(bay), allow_empty=True)
    if isinstance(clip, CamberSide):
        cam = None if camber_xy is None else np.asarray(camber_xy, dtype=float)
        if cam is not None and cam.ndim == 2 and cam.shape[0] >= 2 and cam.shape[1] >= 2:
            from ..constructive.emitters import _camber_clip

            return ensure_polygon(
                _camber_clip(polygon, cam, surface=str(clip.surface)),
                allow_empty=True,
            )
        return half_plane_clip(polygon, surface=str(clip.surface))
    raise TypeError(f"Unsupported clip type: {type(clip)!r}")


__all__ = [
    "CamberSide",
    "ChordBand",
    "ClipSpec",
    "FlapwiseStrip",
    "HalfPlaneClip",
    "_airfoil_offset_arc",
    "_parallel_web_half_gt",
    "_parallel_web_half_lt",
    "_parallel_web_strip_clip",
    "airfoil_expr",
    "apply_shapely",
    "camber_side_field",
    "chordwise_clip_expr",
    "flapwise_web_normal",
    "lower_half_plane_expr",
    "strip_clip_expr",
    "to_ir",
    "upper_half_plane_expr",
    "web_station_projection",
]
