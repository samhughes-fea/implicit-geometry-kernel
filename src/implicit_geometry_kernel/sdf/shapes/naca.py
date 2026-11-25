"""
Parametric NACA section generators (4-, 5-, and 6-series).

``Airfoil.naca4`` and ``AirfoilSDF.from_naca`` both call these helpers; keep the
coefficients and vertex construction here so the two types stay in lock-step.
``Airfoil`` is the chord-normalised shape; ``AirfoilSDF`` is its implicit field
(``Airfoil.to_sdf()``). This module must not import ``blade``.

Five-digit camber/thickness conventions follow Abbott & von Doenhoff / NASA summaries.
Six-series profiles for families 63 and 64 use embedded UIUC Airfoil Data Site
coordinates (unit chord), with linear thickness scaling in y when t_percent differs
from the reference section thickness.
"""

from __future__ import annotations

import math
import re
from functools import lru_cache

import numpy as np
from numpy.typing import NDArray

# Published NACA 5-digit normal (Q = 0) camber parameters for C_l,i = 0.3
# (Abbott & von Doenhoff, Theory of Wing Sections, section 6.6): columns [p, r, k1],
# p the position of maximum camber (P / 20), r the join of the cubic and the line.
_NACA5_NORMAL: NDArray[np.float64] = np.array(
    [
        [0.05, 0.0580, 361.400],
        [0.10, 0.1260, 51.640],
        [0.15, 0.2025, 15.957],
        [0.20, 0.2900, 6.643],
        [0.25, 0.3910, 3.230],
    ],
    dtype=np.float64,
)

# Published NACA 5-digit reflex (S=1) camber parameters for Cli=0.3
# (Abbott & von Doenhoff / NASA TM 4741): columns [p, r, k1, k2/k1].
_NACA5_REFLEX: NDArray[np.float64] = np.array(
    [
        [0.10, 0.1300, 51.990, 0.000764],
        [0.15, 0.2170, 15.793, 0.006770],
        [0.20, 0.3180, 6.520, 0.030300],
        [0.25, 0.4410, 3.191, 0.135500],
    ],
    dtype=np.float64,
)


def naca_four_digit_thickness(xc, t_frac, *, closed_te=True):
    """Analytic NACA 4-digit thickness ``2 yt(x/c)`` (chord-normalised)."""
    return 2.0 * _thickness_yt_4digit(np.asarray(xc, float), float(t_frac), closed_te=closed_te)


def _thickness_yt_4digit(xc, t_frac, *, closed_te):
    # Classic Abbott coefficients: closed TE → yt(1)=0 uses a4=-0.1036;
    # open TE (finite trailing-edge thickness) uses a4=-0.1015.
    a4 = -0.1036 if closed_te else -0.1015
    a = np.array([0.2969, -0.1260, -0.3516, 0.2843, a4], dtype=np.float64)
    return (t_frac / 0.2) * (
        a[0] * np.sqrt(np.maximum(xc, 0.0))
        + a[1] * xc + a[2] * xc**2 + a[3] * xc**3 + a[4] * xc**4
    )


def _wrap_upper_lower(xc, yc, dyc, yt, chord):
    theta = np.arctan(dyc)
    xu = (xc - yt * np.sin(theta)) * chord
    yu = (yc + yt * np.cos(theta)) * chord
    xl = (xc + yt * np.sin(theta)) * chord
    yl = (yc - yt * np.cos(theta)) * chord
    upper = np.column_stack([xu[::-1], yu[::-1]])
    lower = np.column_stack([xl[1:], yl[1:]])
    return np.vstack([upper, lower])


def naca_four_digit_vertices(m, p, t_frac, n_points, chord, *, closed_te=True):
    """Return (N,2) vertices for classic NACA 4-digit section."""
    beta = np.linspace(0.0, math.pi, max(4, n_points // 2 + 1), dtype=np.float64)
    xc = 0.5 * (1.0 - np.cos(beta))
    yt = _thickness_yt_4digit(xc, float(t_frac), closed_te=closed_te)
    m = float(m)
    p = float(p)
    if m <= 0.0 or p <= 0.0:
        yc = np.zeros_like(xc)
        dyc = np.zeros_like(xc)
    else:
        yc = np.where(
            xc <= p,
            (m / p**2) * (2 * p * xc - xc**2),
            (m / (1 - p) ** 2) * ((1 - 2 * p) + 2 * p * xc - xc**2),
        )
        dyc = np.where(
            xc <= p,
            (2 * m / p**2) * (p - xc),
            (2 * m / (1 - p) ** 2) * (p - xc),
        )
    verts = _wrap_upper_lower(xc, yc, dyc, yt, float(chord))
    if closed_te:
        # Collapse residual floating-point / camber TE split to one tip.
        tip = 0.5 * (verts[0] + verts[-1])
        verts[0] = tip
        verts[-1] = tip
    return verts


def _naca5_join(p):
    """Join ``r`` of the normal 5-digit mean line: the root in (p, 1) of 3p^2 - 6pr + r^2(3 - r) = 0 (R2.12)."""
    roots = np.roots([-1.0, 3.0, -6.0 * float(p), 3.0 * float(p) ** 2])
    real = [float(z.real) for z in roots if abs(z.imag) < 1e-12 and float(p) < z.real < 1.0]
    if not real:
        raise ValueError(f"no 5-digit mean-line join for p = {p}")
    return min(real)


def _interp_naca5_table(p):
    """Return (r, k1) for maximum-camber position ``p`` at C_l,i = 0.3.

    The published values are used at their five tabulated positions; between them
    ``r`` comes from (R2.12) and ``k1`` is interpolated in log space.
    """
    p = float(np.clip(p, _NACA5_NORMAL[0, 0], _NACA5_NORMAL[-1, 0]))
    pp = _NACA5_NORMAL[:, 0]
    hit = np.isclose(pp, p, rtol=0.0, atol=1e-12)
    if np.any(hit):
        row = _NACA5_NORMAL[int(np.argmax(hit))]
        return float(row[1]), float(row[2])
    k1 = float(np.exp(np.interp(p, pp, np.log(_NACA5_NORMAL[:, 2]))))
    return _naca5_join(p), k1


def _interp_naca5_reflex(p):
    """Return (r, k1, k2/k1) for design max-camber location ``p`` (Cli=0.3 table)."""
    p = float(np.clip(p, _NACA5_REFLEX[0, 0], _NACA5_REFLEX[-1, 0]))
    pp = _NACA5_REFLEX[:, 0]
    r = float(np.interp(p, pp, _NACA5_REFLEX[:, 1]))
    k1 = float(np.interp(p, pp, _NACA5_REFLEX[:, 2]))
    k2_over_k1 = float(np.interp(p, pp, _NACA5_REFLEX[:, 3]))
    return r, k1, k2_over_k1


def _reflex_k2_over_k1(r, p):
    """Closed-form k2/k1 from join ``r`` and max-camber station ``p``."""
    r = float(r)
    p = float(p)
    denom = (1.0 - r) ** 3
    if abs(denom) < 1e-18:
        raise ValueError("Reflex camber join r too close to trailing edge.")
    return (3.0 * (r - p) ** 2 - r**3) / denom


def naca_five_digit_vertices(L, PQ, TT, n_points, chord, *, closed_te=True, reflex=False):
    """NACA five-digit section ``L P Q TT`` (Abbott & von Doenhoff; R2, Section 3.3).

    ``PQ`` is the two middle digits: ``P`` sets the position of maximum camber
    ``p = P / 20`` and ``Q = 1`` selects the reflexed mean line (as does
    ``reflex=True``). ``L`` sets the design lift ``C_l,i = 3L / 20``.
    """
    if L < 1 or L > 9:
        raise ValueError(f"Five-digit L must be 1..9, got {L}.")
    if PQ < 0 or PQ > 99:
        raise ValueError(f"Five-digit PQ must be 0..99, got {PQ}.")
    if TT < 1 or TT > 99:
        raise ValueError(f"Five-digit TT must be 1..99, got {TT}.")
    cli = (3.0 / 20.0) * float(L)
    lift_scale = cli / 0.3
    beta = np.linspace(0.0, math.pi, max(4, n_points // 2 + 1), dtype=np.float64)
    xc = 0.5 * (1.0 - np.cos(beta))
    yt = _thickness_yt_4digit(xc, float(TT) / 100.0, closed_te=closed_te)

    P, Q = int(PQ) // 10, int(PQ) % 10
    if Q not in (0, 1):
        raise ValueError(f"Five-digit Q (third digit) must be 0 or 1, got {Q}.")
    if not 1 <= P <= 5:
        raise ValueError(f"Five-digit P (second digit) must be 1..5, got {P}.")
    p_design = P / 20.0
    if reflex or Q == 1:
        r, k1_tab, k2_ratio = _interp_naca5_reflex(p_design)
        # Prefer tabulated k2/k1; fall back to closed form if needed.
        if not np.isfinite(k2_ratio):
            k2_ratio = _reflex_k2_over_k1(r, p_design)
        k1 = k1_tab * lift_scale
        yc = np.zeros_like(xc)
        dyc = np.zeros_like(xc)
        mask0 = xc <= r
        # Standard reflex camber (NASA TM 4741 / Abbott & von Doenhoff).
        yc = np.where(
            mask0,
            (k1 / 6.0)
            * (
                (xc - r) ** 3
                - k2_ratio * (1.0 - r) ** 3 * xc
                - r**3 * xc
                + r**3
            ),
            yc,
        )
        yc = np.where(
            ~mask0,
            (k1 / 6.0)
            * (
                k2_ratio * (xc - r) ** 3
                - k2_ratio * (1.0 - r) ** 3 * xc
                - r**3 * xc
                + r**3
            ),
            yc,
        )
        dyc = np.where(
            mask0,
            (k1 / 6.0)
            * (3.0 * (xc - r) ** 2 - k2_ratio * (1.0 - r) ** 3 - r**3),
            dyc,
        )
        dyc = np.where(
            ~mask0,
            (k1 / 6.0)
            * (3.0 * k2_ratio * (xc - r) ** 2 - k2_ratio * (1.0 - r) ** 3 - r**3),
            dyc,
        )
        return _wrap_upper_lower(xc, yc, dyc, yt, float(chord))

    r, k1 = _interp_naca5_table(p_design)
    k1 = k1 * lift_scale
    yc = np.zeros_like(xc)
    dyc = np.zeros_like(xc)
    mask0 = xc <= r
    yc = np.where(mask0, k1 / 6.0 * (xc**3 - 3.0 * r * xc**2 + r**2 * (3.0 - r) * xc), yc)
    yc = np.where(~mask0, k1 * r**3 / 6.0 * (1.0 - xc), yc)
    dyc = np.where(mask0, k1 / 6.0 * (3.0 * xc**2 - 6.0 * r * xc + r**2 * (3.0 - r)), dyc)
    dyc = np.where(~mask0, -k1 * r**3 / 6.0, dyc)
    return _wrap_upper_lower(xc, yc, dyc, yt, float(chord))


def _parse_uiuc_ads_flat_dat(text):
    lines = [ln.strip() for ln in text.strip().splitlines() if ln.strip()]
    nums = []
    for ln in lines:
        if re.match(r"^NACA", ln, re.I):
            m = re.search(r"0\.0+\s+0\.0+", ln)
            if m:
                ln = ln[m.start():]
            else:
                continue
        parts = ln.replace(",", " ").split()
        for p in parts:
            try:
                nums.append(float(p))
            except ValueError:
                pass
    if len(nums) < 8 or len(nums) % 2 != 0:
        raise ValueError("UIUC flat dat: expected an even count of numeric tokens.")
    return np.asarray(nums, dtype=np.float64).reshape(-1, 2)


def _embed_n63415():
    return """NACA 63-415 AIRFOIL 26. 26. 0.000000 0.000000 0.003000 0.012870 0.005250 0.015850 0.009910 0.020740 0.021980 0.029640 0.046600 0.042640 0.071470 0.052610 0.096470 0.060770 0.146690 0.073480 0.197050 0.082790 0.247500 0.089410 0.298000 0.093620 0.348520 0.095590 0.399050 0.095270 0.449550 0.092890 0.500000 0.088710 0.550390 0.082980 0.600700 0.075950 0.650930 0.067800 0.701060 0.058770 0.751090 0.049070 0.801020 0.039000 0.850850 0.028850 0.900590 0.018840 0.950280 0.009310 1.000000 0.000000 0.000000 0.000000 0.007000 -0.010870 0.009750 -0.013050 0.015090 -0.016460 0.028020 -0.022200 0.053400 -0.030000 0.078530 -0.035650 0.103530 -0.040090 0.153310 -0.046560 0.202950 -0.050950 0.252500 -0.053610 0.302000 -0.054740 0.351480 -0.054390 0.400950 -0.052430 0.450450 -0.049090 0.500000 -0.044590 0.549610 -0.039180 0.599300 -0.033110 0.649070 -0.026600 0.698940 -0.019890 0.748910 -0.013270 0.799890 -0.007160 0.849150 -0.001930 0.899410 0.001840 0.949720 0.003330 1.000000 0.000000"""


def _embed_n64212():
    return """NACA 64(1)-212 26.0 26.0 0.0000000 0.0000000 0.0041800 0.0102500 0.0065900 0.0124500 0.0114700 0.0159300 0.0238200 0.0221800 0.0486800 0.0312300 0.0736400 0.0381500 0.0986500 0.0438600 0.1487200 0.0529100 0.1988600 0.0596800 0.2490300 0.0647000 0.2992100 0.0681500 0.3494100 0.0700800 0.3996100 0.0705200 0.4498200 0.0689300 0.5000000 0.0658300 0.5501600 0.0615100 0.6002900 0.0561900 0.6503900 0.0500400 0.7004500 0.0432200 0.7504700 0.0359000 0.8004500 0.0282500 0.8503800 0.0205400 0.9002700 0.0130300 0.9501300 0.0060400 1.0000000 0.0000000 0.0000000 0.0000000 0.0058200 -0.0092500 0.0084100 -0.0110500 0.0135300 -0.0137900 0.0261800 -0.0184600 0.0513200 -0.0249100 0.0763600 -0.0296700 0.1013500 -0.0335200 0.1512800 -0.0394500 0.2011400 -0.0437600 0.2509700 -0.0468000 0.3007900 -0.0487100 0.3505900 -0.0494800 0.4003900 -0.0491000 0.4501800 -0.0470300 0.5000000 -0.0437700 0.5498400 -0.0396100 0.5997100 -0.0347700 0.6496100 -0.0294400 0.6995500 -0.0237800 0.7495300 -0.0180000 0.7995500 -0.0123300 0.8496200 -0.0070800 0.8997300 -0.0026900 0.9498700 0.0002800 1.0000000 0.0000000"""


def _six_series_reference_xy(family):
    if family == 63:
        raw = _parse_uiuc_ads_flat_dat(_embed_n63415())
    elif family == 64:
        raw = _parse_uiuc_ads_flat_dat(_embed_n64212())
    else:
        raise ValueError(f"Unsupported 6-series family: {family}.")
    n = raw.shape[0] // 2
    upper = raw[:n]
    lower = raw[n:]
    upper_o = upper[np.argsort(upper[:, 0])][::-1]
    lower_o = lower[np.argsort(lower[:, 0])]
    return np.vstack([upper_o, lower_o[1:]])


#: design lift coefficient of each embedded reference section (63-415, 64(1)-212)
_SIX_SERIES_REF_CL = {63: 0.4, 64: 0.2}


def six_series_camber_thickness(family, x):
    """Vertical mid-line and vertical thickness of the reference section at chord fractions ``x``."""
    ref = _six_series_reference_xy(int(family))
    n = len(ref) // 2
    upper = ref[: n + 1][::-1]
    lower = ref[n:]
    yu = np.interp(x, upper[:, 0], upper[:, 1])
    yl = np.interp(x, lower[:, 0], lower[:, 1])
    return 0.5 * (yu + yl), yu - yl


def naca_six_series_vertices(family, t_percent, n_points, chord, *, closed_te=True, design_cl_tenth=None):
    """Six-series section from embedded UIUC coordinates (families 63 and 64).

    The reference section (NACA 63-415 or 64(1)-212) is split into its vertical
    mid-line and thickness (R2.19); only the thickness is scaled to ``t_percent``,
    so the camber stays the reference's. ``design_cl_tenth`` (the design lift in
    tenths, as in the code 63-4xx) scales the mid-line linearly from the
    reference's design lift; ``None`` keeps it (R2, Section 3.4).
    """
    _ = closed_te
    fam = int(family)
    if fam not in _SIX_SERIES_REF_CL:
        raise ValueError(f"Unsupported 6-series family: {family}.")
    beta = np.linspace(0.0, math.pi, max(4, int(n_points) // 2 + 1), dtype=np.float64)
    x = 0.5 * (1.0 - np.cos(beta))
    camber, thick = six_series_camber_thickness(fam, x)
    t_ref = float(np.max(thick))
    if t_ref <= 1e-9:
        raise ValueError("Reference six-series airfoil has near-zero thickness.")
    if design_cl_tenth is not None:
        camber = camber * (float(design_cl_tenth) / 10.0) / _SIX_SERIES_REF_CL[fam]
    half = 0.5 * thick * (float(t_percent) / 100.0) / t_ref
    upper = np.column_stack([x[::-1], (camber + half)[::-1]])
    lower = np.column_stack([x[1:], (camber - half)[1:]])
    return np.vstack([upper, lower]) * float(chord)


def airfoil_vertices_from_spanwise(series, m, p, xx, n_points, chord, *, closed_te=True):
    """Dispatch NACA geometry from precompute naca_series + (naca_m, naca_p, naca_xx)."""
    s = int(series)
    key = (
        s,
        float(round(float(m), 10)),
        float(round(float(p), 10)),
        float(round(float(xx), 10)),
        int(n_points),
        float(round(float(chord), 10)),
        bool(closed_te),
    )
    return _airfoil_vertices_from_spanwise_cached(*key).copy()


@lru_cache(maxsize=512)
def _airfoil_vertices_from_spanwise_cached(series, m, p, xx, n_points, chord, closed_te):
    s = int(series)
    if s == 4:
        mi = int(np.clip(int(round(float(m))), 0, 9))
        pi = int(np.clip(int(round(float(p))), 0, 9))
        xxi = int(np.clip(int(round(float(xx))), 0, 99))
        return naca_four_digit_vertices(mi / 100.0, pi / 10.0, xxi / 100.0, n_points, chord, closed_te=closed_te)
    if s == 5:
        L = int(np.clip(int(round(float(m))), 1, 9))
        PQ = int(np.clip(int(round(float(p))), 0, 99))
        TT = int(np.clip(int(round(float(xx))), 1, 99))
        return naca_five_digit_vertices(L, PQ, TT, n_points, chord, closed_te=closed_te)
    if s == 6:
        fam = int(round(float(m)))
        if fam not in (63, 64):
            raise ValueError(f"Six-series family {fam} is not supported.")
        return naca_six_series_vertices(fam, float(xx), n_points, chord, closed_te=closed_te, design_cl_tenth=float(p))
    raise ValueError(f"Unsupported NACA series {series!r} (expected 4, 5, or 6).")


def spanwise_airfoil_label(series, m, p, xx):
    """Short label for logs / plots."""
    s = int(series)
    if s == 4:
        mi = int(np.clip(int(round(float(m))), 0, 9))
        pi = int(np.clip(int(round(float(p))), 0, 9))
        xxi = int(np.clip(int(round(float(xx))), 0, 99))
        return f"{mi:d}{pi:d}{xxi:02d}"
    if s == 5:
        L = int(np.clip(int(round(float(m))), 1, 9))
        PQ = int(np.clip(int(round(float(p))), 0, 99))
        TT = int(np.clip(int(round(float(xx))), 1, 99))
        return f"{L:d}{PQ:02d}{TT:02d}"
    if s == 6:
        fam = int(round(float(m)))
        clt = int(np.clip(int(round(float(p))), 0, 9))
        t = int(np.clip(int(round(float(xx))), 1, 99))
        return f"{fam}-{clt}{t:02d}"
    return f"series{s}"
