"""Clamp spar-cap heights / chordwise extent so upper and lower caps never meet."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from .web_layout import inner_y_at_xs

# Minimum void kept between upper and lower spar caps so the shear core
# (ellipse void) is never consumed: absolute floor and fraction of local depth.
MIN_CAP_VOID_ABS_M = 0.001
MIN_CAP_VOID_FRAC = 0.15

# Solids camber-clip resolution (SG-010). Chord-approx error of the cavity
# midline scales ~ band / n² for similar sections. Budget is well below the
# airfoil-polygon LE floor (~62 µm at n_points=200 / 1 m chord, SG-008), so
# solid areas inherit the polygon contract rather than this knob.
CAMBER_CLIP_DEV_BUDGET_M = 1e-6  # 1 µm
_CAMBER_CLIP_N_MIN = 33
# Fitted to measured 18.56 µm deviation at n=33 on a 0.50 m band (x/c 0.15–0.65
# at 1 m chord, NACA 2412, 6 mm skin): err ≈ C * band / n².
_CAMBER_CLIP_C = 18.56e-6 * (33**2) / 0.50

ExtentAnchor = Literal["deepest", "band_center"]


def default_camber_clip_n_samples(x_lo: float, x_hi: float) -> int:
    """Sample count so cavity-midline chord error stays within the solid-clip budget."""
    band = max(abs(float(x_hi) - float(x_lo)), 1e-12)
    n = int(math.ceil(math.sqrt(_CAMBER_CLIP_C * band / CAMBER_CLIP_DEV_BUDGET_M)))
    return max(_CAMBER_CLIP_N_MIN, n)


def _scale_heights_to_central_half(
    *,
    allowed: np.ndarray,
    caps_sum: float,
    cap_h_upper: float,
    cap_h_lower: float,
    info: dict,
) -> tuple[float, float, float]:
    """Scale both cap heights to fit the central half of the band (or max allowed)."""
    i0, i1 = len(allowed) // 4, (3 * len(allowed)) // 4 + 1
    a_central = float(np.min(allowed[i0:i1]))
    target = max(a_central, 0.5 * float(np.max(allowed)))
    if target <= 0.0:
        target = float(np.max(allowed))
    if target <= 0.0:
        return float(cap_h_upper), float(cap_h_lower), float(caps_sum)
    scale = min(1.0, target / caps_sum) if caps_sum > 0.0 else 1.0
    cap_h_upper = float(cap_h_upper) * scale
    cap_h_lower = float(cap_h_lower) * scale
    caps_sum = cap_h_upper + cap_h_lower
    info["scaled_heights"] = True
    info["scale"] = scale
    return cap_h_upper, cap_h_lower, caps_sum


def _contiguous_extent_deepest(
    *,
    xs_q: np.ndarray,
    depth: np.ndarray,
    allowed: np.ndarray,
    caps_sum_q: float,
    x_lo: float,
    x_hi: float,
) -> tuple[float, float]:
    """Legacy extent: grow the feasible island from the deepest cavity sample."""
    feasible = caps_sum_q <= allowed + 1e-12
    if not np.any(feasible):
        return float(x_lo), float(x_hi)
    i_deep = int(np.nanargmax(np.where(np.isfinite(depth), depth, -np.inf)))
    if not feasible[i_deep]:
        i_deep = int(np.where(feasible)[0][0])
    lo_i = i_deep
    while lo_i > 0 and feasible[lo_i - 1]:
        lo_i -= 1
    hi_i = i_deep
    while hi_i < len(xs_q) - 1 and feasible[hi_i + 1]:
        hi_i += 1
    return float(xs_q[lo_i]), float(xs_q[hi_i])


def _symmetric_extent_about_center(
    *,
    xs_q: np.ndarray,
    allowed: np.ndarray,
    caps_sum_q: float,
    x_lo: float,
    x_hi: float,
    center_x: float,
) -> tuple[float, float]:
    """Largest ``[center-half_w, center+half_w]`` still feasible and inside the authored band.

    Equal trim on both sides keeps family-C pads centred on the web (SG-013).
    Never grows past ``[x_lo, x_hi]``.
    """
    feasible = caps_sum_q <= allowed + 1e-12
    if not np.any(feasible):
        return float(x_lo), float(x_hi)

    # Seed at the sample nearest the intended center that is still feasible.
    i_center = int(np.argmin(np.abs(xs_q - float(center_x))))
    if not feasible[i_center]:
        feasible_idx = np.where(feasible)[0]
        i_center = int(feasible_idx[np.argmin(np.abs(xs_q[feasible_idx] - float(center_x)))])

    half_authored = 0.5 * (float(x_hi) - float(x_lo))
    # Walk outward one sample at a time while both sides remain feasible and
    # the interval stays symmetric about center_x (clip to authored band).
    max_half = min(float(center_x) - float(x_lo), float(x_hi) - float(center_x), half_authored)
    if max_half <= 0.0:
        return float(center_x), float(center_x)

    # Discrete search: largest k such that every sample in [c-h, c+h] is feasible.
    # Binary search on half-width in sample index space via xs_q values.
    lo_h, hi_h = 0.0, float(max_half)
    best_half = 0.0
    for _ in range(48):
        mid_h = 0.5 * (lo_h + hi_h)
        x_a = float(center_x) - mid_h
        x_b = float(center_x) + mid_h
        # All samples strictly inside / on the candidate interval must be feasible.
        in_band = (xs_q >= x_a - 1e-15) & (xs_q <= x_b + 1e-15)
        if not np.any(in_band):
            # Interval thinner than sample spacing — accept if center seed is feasible.
            ok = bool(feasible[i_center])
        else:
            ok = bool(np.all(feasible[in_band]))
        if ok:
            best_half = mid_h
            lo_h = mid_h
        else:
            hi_h = mid_h
        if hi_h - lo_h < 1e-15 * max(half_authored, 1e-9):
            break

    x_lo_eff = float(center_x) - best_half
    x_hi_eff = float(center_x) + best_half
    # Clamp numerically into authored band.
    x_lo_eff = max(x_lo_eff, float(x_lo))
    x_hi_eff = min(x_hi_eff, float(x_hi))
    # Re-assert symmetry if floating clip broke it (prefer equal half about center).
    half = min(float(center_x) - x_lo_eff, x_hi_eff - float(center_x))
    return float(center_x) - half, float(center_x) + half


def cap_feasibility_clamp(
    airfoil_sdf,
    skin_thickness,
    x_lo: float,
    x_hi: float,
    cap_h_upper: float,
    cap_h_lower: float,
    *,
    y_search,
    phi_inner_ref,
    n_samples: int | None = None,
    aft_half_depth_rule: bool = False,
    extent_anchor: ExtentAnchor = "band_center",
    center_x: float | None = None,
) -> tuple[float, float, float, float, dict]:
    """Clamp cap heights / chordwise extent so the caps never meet.

    Samples the inner-cavity depth over ``[x_lo, x_hi]`` and:

    1. Restricts the cap band to a feasible chordwise interval where
       ``cap_h_upper + cap_h_lower + void_min(x) <= depth(x)``.
    2. If no such interval exists, scales both cap heights proportionally so
       the caps fit over the central half of the band, then retries extent.

    Extent policy (``extent_anchor``, ignored when ``aft_half_depth_rule``):

    * ``"band_center"`` (default, families B/C) — keep the pad **symmetric**
      about ``center_x`` (default midpoint of the authored band / web station).
      Shrink half-width equally on both sides; never slide the midpoint
      (SG-013).
    * ``"deepest"`` — legacy growth of the contiguous feasible island from
      ``argmax(depth)``. Kept for explicit comparisons only.

    With ``aft_half_depth_rule=True`` (continuous-spar / family D), the extent
    is never trimmed: each cap height is instead clamped to half the cavity
    depth at the aft end of the band (the aft web), so the caps span the full
    web-to-web band and at most butt at the cavity midline at the aft web.

    Returns ``(x_lo_eff, x_hi_eff, cap_h_upper_eff, cap_h_lower_eff, info)``.
    ``info["camber_xy"]`` holds sampled cavity-midline points for the cap split
    (solids / properties only — not shell midlines; see bug_log SG-010).

    Default ``n_samples`` is derived from band width so the clip polyline's
    chord-approximation stays within :data:`CAMBER_CLIP_DEV_BUDGET_M`.
    """
    caps_sum = float(cap_h_upper) + float(cap_h_lower)
    info: dict = {
        "clamped_extent": False,
        "scaled_heights": False,
        "scale": 1.0,
        "camber_xy": None,
        "extent_anchor": (
            "aft_half_depth" if aft_half_depth_rule else str(extent_anchor)
        ),
    }
    if caps_sum <= 0.0 or x_hi <= x_lo:
        return float(x_lo), float(x_hi), float(cap_h_upper), float(cap_h_lower), info

    n = (
        int(n_samples)
        if n_samples is not None
        else default_camber_clip_n_samples(x_lo, x_hi)
    )
    xs_q = np.linspace(float(x_lo), float(x_hi), n)
    depth = np.full(xs_q.shape, -np.inf)
    y_mid = np.full(xs_q.shape, np.nan)
    # One shared φ grid eval for all sample stations. Feasibility only needs
    # depth / camber to clamp accuracy; skip per-crossing bisection refine.
    probes = inner_y_at_xs(
        airfoil_sdf,
        skin_thickness,
        xs_q,
        y_search=y_search,
        phi_inner=phi_inner_ref,
        bisection_iters=0,
        raise_on_miss=False,
    )
    for i, pair in enumerate(probes):
        if pair is None:
            continue
        yt, yb = pair
        depth[i] = float(yt) - float(yb)
        y_mid[i] = 0.5 * (float(yt) + float(yb))
    if not np.any(np.isfinite(depth)):
        return float(x_lo), float(x_hi), float(cap_h_upper), float(cap_h_lower), info
    ok = np.isfinite(y_mid)
    if np.any(ok):
        info["camber_xy"] = np.column_stack([xs_q[ok], y_mid[ok]])

    if aft_half_depth_rule:
        # Continuous spar: cap height limited to half the cavity depth at the
        # aft web (the shallowest end of the band), full extent preserved.
        i_fin = np.where(np.isfinite(depth))[0]
        h_allow = 0.5 * float(depth[i_fin[-1]])
        if h_allow > 0.0:
            h_u = min(float(cap_h_upper), h_allow)
            h_l = min(float(cap_h_lower), h_allow)
            if h_u < float(cap_h_upper) or h_l < float(cap_h_lower):
                info["scaled_heights"] = True
                info["scale"] = (h_u + h_l) / caps_sum
            return float(x_lo), float(x_hi), h_u, h_l, info
        return float(x_lo), float(x_hi), float(cap_h_upper), float(cap_h_lower), info

    cx = float(center_x) if center_x is not None else 0.5 * (float(x_lo) + float(x_hi))
    # Keep center inside the authored band.
    cx = min(max(cx, float(x_lo)), float(x_hi))

    allowed = depth - np.maximum(MIN_CAP_VOID_ABS_M, MIN_CAP_VOID_FRAC * depth)
    if caps_sum > float(np.max(allowed)):
        # Caps do not fit anywhere at sized thickness: scale to fit the
        # central half of the band (falls back to the deepest point).
        cap_h_upper, cap_h_lower, caps_sum = _scale_heights_to_central_half(
            allowed=allowed,
            caps_sum=caps_sum,
            cap_h_upper=cap_h_upper,
            cap_h_lower=cap_h_lower,
            info=info,
        )
        if caps_sum <= 0.0:
            return float(x_lo), float(x_hi), float(cap_h_upper), float(cap_h_lower), info

    def _extent(caps_sum_q: float) -> tuple[float, float]:
        if extent_anchor == "deepest":
            return _contiguous_extent_deepest(
                xs_q=xs_q,
                depth=depth,
                allowed=allowed,
                caps_sum_q=caps_sum_q,
                x_lo=float(x_lo),
                x_hi=float(x_hi),
            )
        return _symmetric_extent_about_center(
            xs_q=xs_q,
            allowed=allowed,
            caps_sum_q=caps_sum_q,
            x_lo=float(x_lo),
            x_hi=float(x_hi),
            center_x=cx,
        )

    authored_width = float(x_hi) - float(x_lo)

    def _scale_to_fit_full_band() -> None:
        nonlocal cap_h_upper, cap_h_lower, caps_sum
        finite = allowed[np.isfinite(allowed)]
        if finite.size == 0 or caps_sum <= 0.0:
            return
        target = float(np.min(finite))
        if target <= 0.0 or caps_sum <= target + 1e-15:
            return
        scale = target / caps_sum
        cap_h_upper = float(cap_h_upper) * scale
        cap_h_lower = float(cap_h_lower) * scale
        caps_sum = cap_h_upper + cap_h_lower
        info["scaled_heights"] = True
        info["scale"] = float(info["scale"]) * scale

    def _center_feasible(caps_sum_q: float) -> bool:
        feasible = caps_sum_q <= allowed + 1e-12
        if not np.any(feasible):
            return False
        i_c = int(np.argmin(np.abs(xs_q - cx)))
        return bool(feasible[i_c])

    x_lo_eff, x_hi_eff = _extent(caps_sum)
    if (x_hi_eff - x_lo_eff) < 0.5 * authored_width:
        # Barely-feasible / collapsed band: thin the caps rather than emit a
        # needle-thin pad (stiffness must not decrease as height grows).
        # band_center (SG-013): scale to the worst sample in the full authored
        # band so the symmetric +/- half-width about the web can recover.
        # deepest: keep the historic central-half rule.
        if extent_anchor == "band_center":
            _scale_to_fit_full_band()
        elif not info["scaled_heights"]:
            cap_h_upper, cap_h_lower, caps_sum = _scale_heights_to_central_half(
                allowed=allowed,
                caps_sum=caps_sum,
                cap_h_upper=cap_h_upper,
                cap_h_lower=cap_h_lower,
                info=info,
            )
        x_lo_eff, x_hi_eff = _extent(caps_sum)

    # Pad run-out: the web station itself cannot carry a centred pad (no inner
    # intercept / zero half-width). Prefer height 0 over a zero-width strip
    # that later fails midline emit.
    if extent_anchor == "band_center" and (
        (x_hi_eff - x_lo_eff) < 1e-12 or not _center_feasible(caps_sum)
    ):
        info["clamped_extent"] = True
        info["pad_run_out"] = True
        return float(cx), float(cx), 0.0, 0.0, info

    if x_lo_eff > float(x_lo) + 1e-12 or x_hi_eff < float(x_hi) - 1e-12:
        info["clamped_extent"] = True
    return x_lo_eff, x_hi_eff, float(cap_h_upper), float(cap_h_lower), info


# Historic private alias used by factory / tests / examples.
_cap_feasibility_clamp = cap_feasibility_clamp

__all__ = [
    "MIN_CAP_VOID_ABS_M",
    "MIN_CAP_VOID_FRAC",
    "CAMBER_CLIP_DEV_BUDGET_M",
    "ExtentAnchor",
    "default_camber_clip_n_samples",
    "cap_feasibility_clamp",
    "_cap_feasibility_clamp",
]
