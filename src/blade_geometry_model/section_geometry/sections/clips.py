"""Private clip IR helpers — import via :mod:`clip_algebra` in production code."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_contract.errors import (
    InfeasibleSectionGeometry,
)
from blade_geometry_model.section_geometry.sdf.ir import CallableField, HalfPlane, Intersect, Polygon


def _ensure_f64(a):
    """Return float64 ndarray with a fast-path for ndarray inputs."""
    if isinstance(a, np.ndarray) and a.dtype == np.float64:
        return a
    return np.asarray(a, dtype=float)


def airfoil_expr(airfoil_sdf):
    """Return polygon IR for AirfoilSDF-like inputs, else None."""
    verts = getattr(airfoil_sdf, "vertices", None)
    if verts is None:
        return None
    arr = np.asarray(verts, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2 or arr.shape[0] < 3:
        return None
    return Polygon(tuple((float(x), float(y)) for x, y in arr))


def chordwise_clip_expr(x_start: float, x_end: float):
    return Intersect(
        HalfPlane(nx=-1.0, ny=0.0, d=float(x_start)),
        HalfPlane(nx=1.0, ny=0.0, d=-float(x_end)),
    )


def upper_half_plane_expr():
    return HalfPlane(nx=0.0, ny=-1.0, d=0.0)


def lower_half_plane_expr():
    return HalfPlane(nx=0.0, ny=1.0, d=0.0)


def strip_clip_expr(strip_clip, key: str):
    """Wrap a callable strip clip as IR (non-shared geometry only)."""
    return CallableField(key=key, fn=strip_clip)


def _chordwise_clip(x_start, x_end):
    """Return an SDF callable that is negative between x_start and x_end."""

    def _clip(x, y):
        x = _ensure_f64(x)
        left = x - x_start
        right = x_end - x
        return np.maximum(-left, -right)

    return _clip


def _parallel_web_strip_clip(nx, ny, c_a, c_b):
    """Strip between two parallel lines n·p=c_a and n·p=c_b."""
    lo = float(min(c_a, c_b))
    hi = float(max(c_a, c_b))
    nx = float(nx)
    ny = float(ny)

    def _clip(x, y):
        xa = _ensure_f64(x)
        ya = _ensure_f64(y)
        nd = nx * xa + ny * ya
        return np.maximum(lo - nd, nd - hi)

    return _clip


def _parallel_web_half_lt(nx, ny, c_bound):
    """Interior: n·p < c_bound (e.g. fore of first flapwise web)."""
    nx = float(nx)
    ny = float(ny)
    c = float(c_bound)

    def _clip(x, y):
        xa = _ensure_f64(x)
        ya = _ensure_f64(y)
        nd = nx * xa + ny * ya
        return nd - c

    return _clip


def _parallel_web_half_gt(nx, ny, c_bound):
    """Interior: n·p > c_bound (e.g. aft of last flapwise web)."""
    nx = float(nx)
    ny = float(ny)
    c = float(c_bound)

    def _clip(x, y):
        xa = _ensure_f64(x)
        ya = _ensure_f64(y)
        nd = nx * xa + ny * ya
        return c - nd

    return _clip


def flapwise_web_normal(twist_angle_rad):
    """Unit normal **n** (in chord frame) for lines parallel to flapwise web run."""
    w = float(twist_angle_rad)
    return -np.cos(w), np.sin(w)


def web_station_projection(nx, ny, x_web, y_mid):
    """Scalar n·p at a web chord station (mid-thickness point on the web axis)."""
    return float(nx) * float(x_web) + float(ny) * float(y_mid)


def _airfoil_offset_branch(airfoil_sdf, *, depth, surface):
    """Open, strictly x-ascending branch of the inward normal-offset locus.

    Splits the closed offset ring at its chordwise extremes — the direct
    analogue of ``AirfoilSDF.upper_surface()`` / ``lower_surface()``, which split
    the vertex ring at the LE with the TE as the array seam.
    """
    from ..constructive.offset_locus import airfoil_offset_ring

    ring = airfoil_offset_ring(airfoil_sdf, depth)
    # Seam at the TE apex so the LE split yields two clean branches.
    ring = np.roll(ring, -int(np.argmax(ring[:, 0])), axis=0)
    i_le = int(np.argmin(ring[:, 0]))
    branches = (ring[: i_le + 1], np.vstack([ring[i_le:], ring[:1]]))
    pick = max if surface == "upper" else min
    br = np.asarray(pick(branches, key=lambda b: float(np.mean(b[:, 1]))), dtype=float)
    if br.shape[0] >= 2 and br[0, 0] > br[-1, 0]:
        br = br[::-1]
    if br.shape[0] < 2:
        raise InfeasibleSectionGeometry(
            f"Offset branch for surface={surface!r} at depth={float(depth):g} "
            "has fewer than 2 vertices."
        )
    dx = np.diff(br[:, 0])
    if np.any(dx <= 0.0):
        bad = int(np.argmax(dx <= 0.0))
        raise InfeasibleSectionGeometry(
            f"Offset branch for surface={surface!r} at depth={float(depth):g} is "
            f"not strictly x-ascending (first reversal at vertex {bad}: "
            f"x={br[bad, 0]:g} -> {br[bad + 1, 0]:g}). A non-graph locus cannot "
            "satisfy the ShellMidlineStrip 'ascending in chord-frame x' contract."
        )
    return br


def _longest_inside_run(mask):
    """``(i_lo, i_hi)`` of the longest run of True in ``mask``, or ``None``.

    Preferred over first/last-inside-index, which would silently bridge a gap
    between two disjoint inside runs and fabricate material across it.
    """
    best = None
    start = None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if best is None or (i - 1 - start) > (best[1] - best[0]):
                best = (start, i - 1)
            start = None
    if start is not None and (
        best is None or (len(mask) - 1 - start) > (best[1] - best[0])
    ):
        best = (start, len(mask) - 1)
    return best


def _airfoil_offset_arc(
    airfoil_sdf,
    *,
    skin_thickness,
    cap_height,
    x_start,
    x_end,
    surface,
    n=80,
    strip_clip=None,
):
    """Sample the cap mid-surface in the chord frame.

    The locus is the **normal** offset at depth ``skin_thickness + cap_height/2``
    — the mid-surface of the ``Slab(Offset(poly, -skin_thickness), cap_height)``
    solid that ``SparCap`` builds. It was previously a *vertical* displacement of
    the outer surface, which lands at normal depth ``d*cos(theta)`` and so ran
    shallow wherever the surface is steep, in the worst measured case sitting
    inside the skin band rather than in the cap.
    """
    if surface not in ("upper", "lower"):
        raise ValueError(f"surface must be 'upper' or 'lower', got {surface!r}")
    depth = float(skin_thickness) + float(cap_height) / 2.0
    br = _airfoil_offset_branch(airfoil_sdf, depth=depth, surface=surface)
    xb, yb = br[:, 0], br[:, 1]

    if strip_clip is None:
        x0, x1 = float(x_start), float(x_end)
    else:
        phi = np.asarray(strip_clip(xb, yb), dtype=float)
        run = _longest_inside_run(phi <= 0.0)
        if run is None:
            return br[:0]
        i_lo, i_hi = run

        def _cross(i, j):
            a, b = float(phi[i]), float(phi[j])
            if a == b:
                return float(xb[i])
            return float(xb[i] + (a / (a - b)) * (xb[j] - xb[i]))

        x0 = _cross(i_lo, i_lo - 1) if i_lo > 0 else float(xb[i_lo])
        x1 = _cross(i_hi, i_hi + 1) if i_hi < xb.size - 1 else float(xb[i_hi])

    # Clamp to where the offset locus actually exists.
    x0 = max(float(x0), float(xb[0]))
    x1 = min(float(x1), float(xb[-1]))
    if not (x1 > x0):
        return br[:0]

    # Merge the locus' own vertices into the sample grid. Without them the
    # polyline chords across its corners and sags by several percent of
    # ``depth`` between samples, even though the samples themselves are exact.
    xs = np.linspace(x0, x1, int(max(2, n)))
    knots = xb[(xb > x0) & (xb < x1)]
    xs = np.unique(np.concatenate([xs, knots]))
    tol = 1e-9 * (x1 - x0)
    keep = [float(xs[0])]
    for v in xs[1:]:
        if float(v) - keep[-1] > tol:
            keep.append(float(v))
    keep[-1] = x1
    xs = np.asarray(keep, dtype=float)
    return np.column_stack([xs, np.interp(xs, xb, yb)])
