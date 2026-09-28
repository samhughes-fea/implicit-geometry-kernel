"""Web anchor geometry helpers for MultiCellSection."""

from __future__ import annotations

from blade_structure_utilities.geometry_utilities.section_contract.errors import (
    InfeasibleSectionGeometry,
)

import numpy as np

from blade_structure_utilities.geometry_utilities.section_geometry.sdf.csg import offset


def _default_y_search(airfoil_sdf, y_search=None):
    """Resolve a y-window wide enough for tip chords (see ``inner_y_at_x``)."""
    if y_search is not None:
        return tuple(y_search)
    verts = getattr(airfoil_sdf, "vertices", None)
    if verts is not None:
        v = np.asarray(verts, dtype=float)
        if v.ndim == 2 and v.shape[0] > 0 and v.shape[1] >= 2:
            ymin = float(v[:, 1].min())
            ymax = float(v[:, 1].max())
            # Tip chords are O(1e-2) thick; a fixed (-0.5, 0.5) window with
            # ~240 samples undersamples and falsely reports zero crossings.
            pad = max(0.05 * max(ymax - ymin, 1e-9), 0.01)
            return (ymin - pad, ymax + pad)
    return (-0.5, 0.5)


def _interp_crossing_1d(y_vals, phi_vals, idx, x_query, phi_inner, bisection_iters):
    y0, y1 = float(y_vals[idx]), float(y_vals[idx + 1])
    p0, p1 = float(phi_vals[idx]), float(phi_vals[idx + 1])
    if p0 == 0.0:
        return y0
    if p1 == 0.0:
        return y1
    if p0 * p1 < 0.0:
        for _ in range(int(bisection_iters)):
            ym = 0.5 * (y0 + y1)
            pm = float(phi_inner(np.array([x_query]), np.array([ym]))[0])
            if pm == 0.0:
                return ym
            if p0 * pm < 0.0:
                y1, p1 = ym, pm
            else:
                y0, p0 = ym, pm
    return y0 - p0 * (y1 - y0) / (p1 - p0 + 1e-30)


def inner_y_at_x(
    airfoil_sdf,
    skin_thickness,
    x_query,
    y_search=None,
    n_samples=240,
    *,
    phi_inner=None,
    bisection_iters=10,
):
    """Estimate the upper and lower inner-skin y-coordinates at a given x."""
    y_search = _default_y_search(airfoil_sdf, y_search)
    y_vals = np.linspace(y_search[0], y_search[1], int(n_samples))
    x_vals = np.full_like(y_vals, x_query)

    if phi_inner is None:
        phi_inner = offset(airfoil_sdf, -skin_thickness)
    phi_vals = phi_inner(x_vals, y_vals)

    crossings = np.where(np.diff(np.sign(phi_vals)))[0]

    if len(crossings) < 2:
        raise InfeasibleSectionGeometry(
            "Could not determine inner skin intercepts at "
            f"x={x_query:.6g} within y_search={tuple(y_search)}. "
            "Provide explicit web_y_coords or widen y_search."
        )

    y_crossings = sorted(
        [
            _interp_crossing_1d(y_vals, phi_vals, i, x_query, phi_inner, bisection_iters)
            for i in crossings
        ]
    )
    y_bot = y_crossings[0]
    y_top = y_crossings[-1]
    return y_top, y_bot


def inner_y_at_xs(
    airfoil_sdf,
    skin_thickness,
    x_queries,
    y_search=None,
    n_samples=240,
    *,
    phi_inner=None,
    bisection_iters=10,
    raise_on_miss: bool = False,
):
    """Batch :func:`inner_y_at_x` over many ``x`` with one shared φ grid eval.

    Returns a list of ``(y_top, y_bot)`` or ``None`` (when ``raise_on_miss`` is
    false and that column has fewer than two sign crossings).
    """
    xs = [float(x) for x in x_queries]
    if not xs:
        return []

    y_search = _default_y_search(airfoil_sdf, y_search)
    y_vals = np.linspace(y_search[0], y_search[1], int(n_samples))
    if phi_inner is None:
        phi_inner = offset(airfoil_sdf, -skin_thickness)

    # One (n_x, n_y) φ evaluation instead of n_x independent columns.
    x_arr = np.asarray(xs, dtype=float)
    X = np.repeat(x_arr[:, None], y_vals.size, axis=1)
    Y = np.broadcast_to(y_vals[None, :], X.shape)
    phi = np.asarray(phi_inner(X.ravel(), Y.ravel()), dtype=float).reshape(X.shape)

    out: list[tuple[float, float] | None] = []
    for j, x_query in enumerate(xs):
        phi_vals = phi[j]
        crossings = np.where(np.diff(np.sign(phi_vals)))[0]
        if len(crossings) < 2:
            if raise_on_miss:
                raise InfeasibleSectionGeometry(
                    "Could not determine inner skin intercepts at "
                    f"x={x_query:.6g} within y_search={tuple(y_search)}. "
                    "Provide explicit web_y_coords or widen y_search."
                )
            out.append(None)
            continue
        y_crossings = sorted(
            [
                _interp_crossing_1d(
                    y_vals, phi_vals, i, x_query, phi_inner, bisection_iters
                )
                for i in crossings
            ]
        )
        out.append((float(y_crossings[-1]), float(y_crossings[0])))
    return out


def resolve_web_anchors(
    airfoil_sdf,
    skin_thickness,
    web_x_positions,
    y_search=None,
    *,
    phi_inner=None,
    min_web_gap_frac=0.02,
    min_cavity_depth_m: float = 0.0,
    allow_le_ward_migration: bool = False,
    return_dropped: bool = False,
):
    """Resolve ``(x, y_top, y_bot)`` for each web; drop those that do not fit.

    Near the tip / TE, local thickness can be smaller than the inner-skin offset
    (``skin_thickness/2``), so sampling yields zero SDF crossings and the web
    cannot be anchored there.

    ``allow_le_ward_migration`` selects what happens then:

    * ``False`` (default) — the web is **dropped**. Web chordwise positions are
      anchored spanwise, so a web that does not fit at its anchored station is
      one the section can no longer carry; sliding it forward would move the
      web centroid station to station, which is exactly what the anchored model
      exists to prevent.
    * ``True`` — legacy behaviour: walk the web LE-ward in up to 12 steps until
      intercepts exist, dropping it only if none are found. Retained because it
      is reachable for comparison, not because it is recommended.

    Webs closer together than ``min_web_gap_frac × chord`` collapse onto the
    previous one under either setting, rather than stacking duplicates.

    ``min_cavity_depth_m`` drops a web whose local cavity is too shallow to be
    worth anchoring, even where an intercept still exists. Having an intercept
    is a topological test, not a structural one: as the chord shrinks, an
    anchored web migrates aft in ``x/c`` into ever-shallower cavity, and on the
    reference rotor one sits at 0.85c in a 1.0 mm cavity while still probing
    two crossings. Family D then keys its cap height on that web —
    ``h_allow = 0.5 x depth`` — and reports a 0.5 mm allowance against a 2 mm
    minimum gauge, so the section carries a web it cannot build a cap against.
    Defaults to ``0.0``, which is exactly the previous behaviour; callers that
    know their minimum gauge opt in.

    With ``return_dropped`` the result is ``(resolved, dropped)`` where
    ``dropped`` is a list of :class:`DroppedWeb` — callers that step the
    layout down need to know *which* webs went, not just how many.
    """
    from blade_structure_utilities.geometry_utilities.section_geometry.api.worker import (
        DroppedWeb,
    )

    xs_in = [float(x) for x in web_x_positions]
    if not xs_in:
        return ([], []) if return_dropped else []

    verts = getattr(airfoil_sdf, "vertices", None)
    xmin = xmax = None
    if verts is not None:
        v = np.asarray(verts, dtype=float)
        if v.ndim == 2 and v.shape[0] > 0 and v.shape[1] >= 2:
            xmin = float(v[:, 0].min())
            xmax = float(v[:, 0].max())
            if y_search is None:
                ymin = float(v[:, 1].min())
                ymax = float(v[:, 1].max())
                pad = max(0.05 * max(xmax - xmin, 1e-9), 0.01)
                y_search = (ymin - pad, ymax + pad)

    chord = float(getattr(airfoil_sdf, "chord", 0.0) or 0.0)
    if chord <= 0.0 and xmin is not None and xmax is not None:
        chord = max(xmax - xmin, 1e-12)
    if chord <= 0.0:
        chord = 1.0
    min_gap = float(min_web_gap_frac) * chord
    min_depth = max(0.0, float(min_cavity_depth_m))
    x_le = xmin if xmin is not None else min(xs_in)
    x_floor = x_le + 0.05 * chord

    resolved = []
    dropped: list[DroppedWeb] = []

    # Common path: no LE-ward walk — batch all intercept probes in one φ eval.
    if not allow_le_ward_migration:
        probes = inner_y_at_xs(
            airfoil_sdf,
            skin_thickness,
            sorted(xs_in),
            y_search=y_search,
            phi_inner=phi_inner,
            raise_on_miss=False,
        )
        for x0, anchors_yb in zip(sorted(xs_in), probes):
            if anchors_yb is None:
                dropped.append(DroppedWeb(x_m=float(x0), reason="no_inner_skin_intercept"))
                continue
            y_top, y_bot = anchors_yb
            if float(y_top) - float(y_bot) < min_depth:
                dropped.append(DroppedWeb(x_m=float(x0), reason="cavity_too_shallow"))
                continue
            anchors = (float(x0), float(y_top), float(y_bot))
            if resolved and abs(anchors[0] - resolved[-1][0]) < min_gap:
                dropped.append(DroppedWeb(x_m=float(x0), reason="merged_into_previous_web"))
                continue
            resolved.append(anchors)
        if return_dropped:
            return resolved, dropped
        return resolved

    for x0 in sorted(xs_in):
        x_try = float(x0)
        anchors = None
        n_steps = 12
        for step in range(n_steps + 1):
            try:
                y_top, y_bot = inner_y_at_x(
                    airfoil_sdf,
                    skin_thickness,
                    x_try,
                    y_search=y_search,
                    phi_inner=phi_inner,
                )
                anchors = (x_try, float(y_top), float(y_bot))
                break
            except ValueError:
                # Pull LE-ward; last step forces near leading-edge floor.
                frac = (step + 1) / float(n_steps)
                x_try = float(x0) - frac * (float(x0) - x_floor)
                if x_try <= x_floor + 1e-12:
                    x_try = x_floor
        if anchors is None:
            dropped.append(DroppedWeb(x_m=float(x0), reason="no_inner_skin_intercept"))
            continue
        if anchors[1] - anchors[2] < min_depth:
            # Applied after the LE-ward walk, not during it: the walk exists to
            # find *an* intercept, and stopping it early at a shallow station
            # would drop webs this path is meant to rescue.
            dropped.append(DroppedWeb(x_m=float(x0), reason="cavity_too_shallow"))
            continue
        if resolved and abs(anchors[0] - resolved[-1][0]) < min_gap:
            # Collapse onto previous web rather than stacking duplicates.
            dropped.append(DroppedWeb(x_m=float(x0), reason="merged_into_previous_web"))
            continue
        resolved.append(anchors)

    if return_dropped:
        return resolved, dropped
    return resolved


def flapwise_web_sframe_inputs(
    x: float,
    y_top: float,
    y_bot: float,
    twist_angle: float,
    phi_inner_b,
    y_search_b,
):
    """Compute S-frame capsule endpoints for a flapwise web at twist."""
    cos_t = float(np.cos(float(twist_angle)))
    sin_t = float(np.sin(float(twist_angle)))
    mid_y_init = 0.5 * (y_top + y_bot)
    x_b_web = cos_t * x - sin_t * mid_y_init
    try:
        y_b_top, y_b_bot = inner_y_at_x(
            None,
            None,
            x_b_web,
            y_search=y_search_b,
            phi_inner=phi_inner_b,
        )
    except (ValueError, InfeasibleSectionGeometry) as exc:
        # Do not rotate stale S-frame IML anchors into B: under twist that lands
        # the web feet off the skin mid-surface by ~t_skin/2 (SG-002 tip mono).
        # Callers that froze web_y at a thinner skin must re-resolve or step down.
        raise InfeasibleSectionGeometry(
            "flapwise web cannot re-anchor on the twisted inner skin at "
            f"x_B={x_b_web:.6g} (S-frame x={float(x):.6g}); frozen web_y "
            f"({float(y_top):.6g},{float(y_bot):.6g}) is not usable at this "
            f"skin gauge under twist. Re-resolve the station layout. ({exc})"
        ) from exc
    x_top_S = cos_t * x_b_web + sin_t * y_b_top
    y_top_S = -sin_t * x_b_web + cos_t * y_b_top
    x_bot_S = cos_t * x_b_web + sin_t * y_b_bot
    y_bot_S = -sin_t * x_b_web + cos_t * y_b_bot
    return x_top_S, y_top_S, x_bot_S, y_bot_S
