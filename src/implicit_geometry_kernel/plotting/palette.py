"""Shared colour tokens and field-drawing helpers for section_geometry figures.

Two jobs:

1. **One palette.** Region colours are keyed on the region label, so
   ``shear_web_0`` is the same green in the anatomy figure, the expression-tree
   figure and the fly-through. The accent tokens match the conference deck
   chrome (``docs/conferences/CDFAM-2026``), so a figure dropped on a slide does
   not read as a foreign object.

2. **Distance fields that look like distance fields.** Every SDF panel used to
   draw ``levels=[0.0]`` and nothing else, which throws away the one thing a
   signed-distance field shows for free: the concentric iso-distance banding
   away from the surface. :func:`draw_distance_field` paints the field with a
   norm clipped to the feature scale and overlays those bands.

Companion to :mod:`figure_io`, which owns save/DPI/rcParams.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# --- Deck-matched tokens -----------------------------------------------------
# Lifted from the CDFAM-2026 deck ``:root`` block so figures and slide chrome
# agree. Changing these changes every figure at once, which is the point.

NAVY = "#041E42"
RED = "#D50032"
TEAL = "#007288"
MID_BLUE = "#1a3a6b"

LIGHT_BLUE = "#E8F4F8"
LIGHT_RED = "#FDE8ED"
LIGHT_TEAL = "#e8f6f6"

INK = "#1a1a2e"
MUTED = "#4b5563"
SURFACE = "#f4f6f9"

def _truncated_diverging(name: str = "RdBu_r", *, keep: float = 0.82):
    """``RdBu_r`` with its darkest ends trimmed.

    The full map runs to near-black maroon and navy. On a field clipped to the
    feature scale, everything past the clip lands on those extremes, so a panel
    ends up carrying most of its ink as saturated far-field with the informative
    near-field squeezed into a thin pale band. Trimming the ends keeps the
    diverging structure and hands the contrast back to the geometry.
    """
    base = colormaps[name]
    margin = 0.5 * (1.0 - float(keep))
    samples = base(np.linspace(margin, 1.0 - margin, 256))
    return LinearSegmentedColormap.from_list(f"{name}_trimmed", samples)


#: Diverging map for signed distance. Blue inside, red outside, white at φ = 0.
SDF_COLORMAP = _truncated_diverging()

# --- Region colours ----------------------------------------------------------
# Keyed on the region label as it appears in ``MultiCellSection`` /
# ``properties_map``. Indexed families (``shear_web_0``, ``sandwich_core_2``)
# resolve through their stem, so a five-web section stays legible without
# inventing a colour per index.

REGION_PALETTE: dict[str, str] = {
    "skin_suction_surface": NAVY,
    "skin_pressure_surface": "#5B8DB8",
    "spar_cap_suction_surface": RED,
    "spar_cap_pressure_surface": "#F0883E",
    "shear_web": TEAL,
    "sandwich_core": "#9E86C8",
}

#: Fallback for a label the palette has never seen.
REGION_FALLBACK = MUTED


def region_color(label: str) -> str:
    """Resolve a region label to its palette colour.

    Exact match first, then the indexed-family stem, so ``shear_web_3`` lands on
    the ``shear_web`` entry. Unknown labels get :data:`REGION_FALLBACK` rather
    than raising — a figure with one odd-coloured region is a better failure
    than a demo that will not render.
    """
    text = str(label)
    if text in REGION_PALETTE:
        return REGION_PALETTE[text]
    stem = text.rstrip("0123456789").rstrip("_")
    return REGION_PALETTE.get(stem, REGION_FALLBACK)


# --- Distance fields ---------------------------------------------------------


def distance_field_scales(
    phi: Any,
    *,
    feature_scale: float | None = None,
) -> tuple[float, float]:
    """Return ``(inner_scale, outer_scale)`` for a signed-distance colour ramp.

    ``inner_scale`` is the deepest interior distance — the inradius of the
    region — so the whole interior spends the full blue half of the map rather
    than a sliver of it. ``outer_scale`` clips the exterior at twice that, which
    is what stops a panel drowning in saturated far-field colour: without the
    clip, a 4 mm web inside a 1 m domain renders as one white hairline on a
    uniformly dark ground.
    """
    values = np.asarray(phi, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 1.0, 1.0

    interior = finite[finite < 0.0]
    if interior.size:
        inner = float(-interior.min())
    else:
        # No interior on this panel (an empty intersection, say): fall back to a
        # robust fraction of the exterior spread.
        inner = float(np.percentile(np.abs(finite), 5.0))

    if inner <= 0.0 or not np.isfinite(inner):
        inner = float(np.abs(finite).max()) or 1.0

    outer = float(feature_scale) if feature_scale is not None else 2.0 * inner
    if outer <= 0.0 or not np.isfinite(outer):
        outer = inner

    return inner, outer


def distance_field_norm(
    phi: Any,
    *,
    feature_scale: float | None = None,
) -> TwoSlopeNorm:
    """Colour norm centred on the zero level set and clipped to the feature scale."""
    inner, outer = distance_field_scales(phi, feature_scale=feature_scale)
    return TwoSlopeNorm(vcenter=0.0, vmin=-inner, vmax=outer)


def draw_distance_field(
    axis: Any,
    X: Any,
    Y: Any,
    phi: Any,
    *,
    feature_scale: float | None = None,
    norm: Any = None,
    band_interval: float | None = None,
    band_count: int = 8,
    colormap: str = SDF_COLORMAP,
    band_color: str = INK,
    band_alpha: float = 0.28,
    band_line_width: float = 0.6,
    zero_color: str = "k",
    zero_line_width: float = 1.8,
    draw_zero: bool = True,
    draw_bands: bool = True,
    shading: str = "auto",
) -> Any:
    """Paint ``phi`` as a distance field with iso-distance banding.

    The bands are the whole point. A signed-distance field carries magnitude
    everywhere, not just a boundary, and drawing contours at regular distance
    intervals is what makes that legible — and what makes an implicit-geometry
    figure read as one.

    Pass ``norm`` to override the derived colour scale entirely; otherwise it
    comes from :func:`distance_field_norm`.

    Returns the ``QuadMesh`` so the caller can attach a colourbar.
    """
    values = np.asarray(phi, dtype=float)
    if norm is None:
        norm = distance_field_norm(values, feature_scale=feature_scale)

    mesh = axis.pcolormesh(X, Y, values, cmap=colormap, norm=norm, shading=shading)

    if draw_bands:
        interval = float(band_interval) if band_interval is not None else None
        if interval is None or interval <= 0.0:
            interval = float(norm.vmax) / max(int(band_count), 1)
        if interval > 0.0:
            lower = float(np.floor(norm.vmin / interval)) * interval
            upper = float(np.ceil(norm.vmax / interval)) * interval
            levels = np.arange(lower, upper + 0.5 * interval, interval)
            # Drop the zero level: it is drawn separately and heavier below.
            levels = levels[np.abs(levels) > 0.5 * interval]
            if levels.size:
                axis.contour(
                    X,
                    Y,
                    values,
                    levels=levels,
                    colors=band_color,
                    alpha=band_alpha,
                    linewidths=band_line_width,
                )

    if draw_zero:
        axis.contour(X, Y, values, levels=[0.0], colors=zero_color, linewidths=zero_line_width)

    return mesh


# --- Axes framing ------------------------------------------------------------


def tighten_to_section(
    axis: Any,
    X: Any = None,
    Y: Any = None,
    *,
    pad_fraction: float = 0.06,
    equal_aspect: bool = True,
) -> None:
    """Set equal aspect and crop to the drawn geometry.

    A blade station is a thin lens: left to itself, matplotlib gives it square
    axes and the section occupies a fifth of the panel. Cropping to the data
    bounds with a small pad routinely doubles the apparent size at no cost.

    With ``X``/``Y`` omitted, the current data limits are used.
    """
    if equal_aspect:
        axis.set_aspect("equal", adjustable="datalim")

    if X is not None and Y is not None:
        x_values = np.asarray(X, dtype=float)
        y_values = np.asarray(Y, dtype=float)
        x_lo, x_hi = float(np.nanmin(x_values)), float(np.nanmax(x_values))
        y_lo, y_hi = float(np.nanmin(y_values)), float(np.nanmax(y_values))
    else:
        axis.relim()
        axis.autoscale_view()
        x_lo, x_hi = axis.get_xlim()
        y_lo, y_hi = axis.get_ylim()

    x_pad = (x_hi - x_lo) * float(pad_fraction)
    y_pad = (y_hi - y_lo) * float(pad_fraction)
    if x_hi - x_lo > 0.0:
        axis.set_xlim(x_lo - x_pad, x_hi + x_pad)
    if y_hi - y_lo > 0.0:
        axis.set_ylim(y_lo - y_pad, y_hi + y_pad)


def strip_axes_chrome(axis: Any, *, keep_frame: bool = False) -> None:
    """Remove ticks, tick labels and (by default) spines.

    For pictorial panels only — anything a viewer reads a number off keeps its
    axes. On a slide, ``0.6 0.4 0.2 0.0`` under an aerofoil is ink that carries
    nothing and crowds out the drawing.
    """
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xlabel("")
    axis.set_ylabel("")
    if not keep_frame:
        for spine in axis.spines.values():
            spine.set_visible(False)


__all__ = [
    "INK",
    "LIGHT_BLUE",
    "LIGHT_RED",
    "LIGHT_TEAL",
    "MID_BLUE",
    "MUTED",
    "NAVY",
    "RED",
    "REGION_FALLBACK",
    "REGION_PALETTE",
    "SDF_COLORMAP",
    "SURFACE",
    "TEAL",
    "distance_field_norm",
    "distance_field_scales",
    "draw_distance_field",
    "region_color",
    "strip_axes_chrome",
    "tighten_to_section",
]
