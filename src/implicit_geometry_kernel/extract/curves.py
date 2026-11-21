"""Spanwise feature-curve locators and trailing-edge classification."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from implicit_geometry_kernel.sdf.solve import GeometrySolver


def locate_shell_vertex(geo: Any, vid: Any, zs: ArrayLike, n_chord: int = 8) -> np.ndarray:
    """XYZ of a :data:`~implicit_geometry_kernel.blade.topology.ShellVertexId` at each ``z``.

    Uses the authored keys on the live section (OML TE snap, web feet), not a
    straight chord between stations. Infeasible stations are omitted.
    """
    from implicit_geometry_kernel.blade.topology import (
        place_station_keys,
        section_cut_near,
    )
    from implicit_geometry_kernel.blade.topology.map import _key_name
    from implicit_geometry_kernel.blade.topology.ribs import _template_from_cut

    out = []
    for z in np.asarray(zs, float).ravel():
        cut = section_cut_near(geo, float(z))
        if cut is None:
            continue
        sg = cut.shell_graph()
        tmpl = _template_from_cut(geo.region_graph, cut, float(z), n_chord)
        wrap_on = any(b.wrap for b in tmpl.loop)
        keys = place_station_keys(cut, sg, geo.region_graph, float(z), te_blunt=wrap_on)
        name = _key_name(vid)
        side = "suction"
        if len(vid) > 1 and vid[0] == "TE" and vid[1] in ("suction", "pressure"):
            side = vid[1]
        elif len(vid) > 2 and vid[2] in ("suction", "pressure"):
            side = vid[2]
        pt = keys.get((name, side))
        if pt is None:
            pt = keys.get((name, "suction")) or keys.get((name, "pressure"))
        if pt is not None:
            out.append(np.asarray(pt, float)[:3])
    return np.vstack(out) if out else np.zeros((0, 3), dtype=float)


def feature_curve(a: Any, b: Any, zs: ArrayLike, guess_xy: ArrayLike | None = None,
                  solver: GeometrySolver | None = None) -> np.ndarray:
    """Trace a spanwise feature.

    Newton form: ``feature_curve(phi_a, phi_b, zs, guess_xy)`` intersects
    ``phi_a = phi_b = 0`` in each ``z`` plane (no silent lerp).

    Named form: ``feature_curve(geo, name_a, name_b, zs)`` delegates to
    :func:`~implicit_geometry_kernel.blade.shell_graph.partition_edge_curve`
    (graph vertices; stations without the pair are omitted).

    Vertex form: ``feature_curve(geo, vid, zs)`` locates a shell vertex id.
    """
    if isinstance(b, tuple) and not isinstance(b, np.ndarray):
        return locate_shell_vertex(a, b, zs)
    if isinstance(b, str):
        from implicit_geometry_kernel.blade.shell_graph import partition_edge_curve

        if not isinstance(zs, str):
            raise TypeError("named feature_curve needs (geo, name_a, name_b, zs)")
        return partition_edge_curve(a, b, zs, guess_xy)
    solver = solver or GeometrySolver(tolerance=1e-8, max_step=0.05)
    g = np.asarray(guess_xy, float).reshape(2)
    out = []
    for z in np.asarray(zs, float):
        start = np.array([[g[0], g[1], float(z)]])
        hit = solver.intersect(a, b, start)
        p = hit.points[0]
        out.append(p)
        g = p[:2]
    return np.asarray(out)


_TE_CORNER_DEG = 150.0
"""A convex-hull vertex sharper than this is a trailing-edge corner; a smooth nose has none."""


def _hull(pts: np.ndarray) -> np.ndarray:
    """Counter-clockwise convex hull of 2D points (monotone chain), no repeated end."""
    order = np.lexsort((pts[:, 1], pts[:, 0]))
    p = pts[order]

    def half(seq: np.ndarray) -> list[np.ndarray]:
        out: list[np.ndarray] = []
        for q in seq:
            while len(out) >= 2:
                a, b = out[-2], out[-1]
                if (b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0]) > 0.0:
                    break
                out.pop()
            out.append(q)
        return out

    lower, upper = half(p), half(p[::-1])
    return np.asarray(lower[:-1] + upper[:-1], float)


def classify_trailing_edge(ring_xy: Any, blunt_tol: float = 1e-4) -> str:
    """``sharp`` if the TE is a point; ``blunt`` if a finite edge remains.

    Accepts an OML ring ``(N, 2|3)`` or a :class:`~implicit_geometry_kernel.blade.shell_graph.ShellGraph`
    (stacked skin loop via ``loop_polygon``). The test needs no frame: the
    trailing edge is the sharpest convex-hull corner, and it is blunt when the
    next hull corner is also sharper than ``_TE_CORNER_DEG`` and lies more than
    ``blunt_tol`` from it (the two ends of a flatback). A twisted section
    classifies the same as an untwisted one.
    """
    if hasattr(ring_xy, "loop_polygon"):
        pts = np.asarray(ring_xy.loop_polygon(), float)
    else:
        pts = np.asarray(ring_xy, float)
    pts = pts[:, :2]
    if np.linalg.norm(pts[0] - pts[-1]) < 1e-12:
        pts = pts[:-1]
    hull = _hull(pts)
    n = len(hull)
    if n < 3:
        return "sharp"
    into = hull - np.roll(hull, 1, axis=0)
    out = np.roll(hull, -1, axis=0) - hull
    cos = -(into * out).sum(1) / np.maximum(np.linalg.norm(into, axis=1) * np.linalg.norm(out, axis=1), 1e-300)
    angle = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    i = int(np.argmin(angle))
    for j in ((i - 1) % n, (i + 1) % n):
        if angle[j] < _TE_CORNER_DEG and np.linalg.norm(hull[j] - hull[i]) > blunt_tol:
            return "blunt"
    return "sharp"
