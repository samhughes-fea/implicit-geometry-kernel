"""Rib faces from mapped loft samples."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.stations import TransverseMember
from implicit_geometry_kernel.blade.topology.graph import WallRegionGraph, _cut_tc
from implicit_geometry_kernel.blade.topology.map import (
    _key_name,
    sample_region_polylines,
)
from implicit_geometry_kernel.blade.topology.map_loft import (
    _append_ring,
    _web_pts_from_su,
)
from implicit_geometry_kernel.blade.topology.map_walk import _chord_fracs
from implicit_geometry_kernel.blade.topology.template import ShellTopologyTemplate

def rib_cell_polygons_from_samples(
    template: ShellTopologyTemplate,
    samples: Mapping[str, np.ndarray],
    member: TransverseMember | str,
) -> list[tuple[str, np.ndarray]]:
    """One n-gon per cell from loft skin/web nodes (no Steiner vertices)."""
    label = member if isinstance(member, str) else str(member.name)
    su = [b for b in template.loop if b.side == "suction" and not b.wrap]
    pr = [b for b in template.loop if b.side == "pressure" and not b.wrap]
    wrap = next((b for b in template.loop if b.wrap), None)
    webs = list(template.webs)
    by_key = {_key_name(w.suction): w for w in webs}
    n_cells = min(len(su), len(pr))
    out: list[tuple[str, np.ndarray]] = []
    for k, (sb, pb) in enumerate(zip(su, pr)):
        if sb.name not in samples or pb.name not in samples:
            continue
        su_pts = np.asarray(samples[sb.name], float)
        pr_pts = np.asarray(samples[pb.name], float)
        if len(su_pts) < 2 or len(pr_pts) < 2:
            continue
        aft = by_key.get(_key_name(sb.end))
        fwd = by_key.get(_key_name(sb.start))
        aft_pts = _web_pts_from_su(samples, aft, su_pts[-1])
        fwd_pts = _web_pts_from_su(samples, fwd, su_pts[0])
        wrap_pts = None
        if wrap is not None and wrap.name in samples and _key_name(sb.end) == "te":
            wrap_pts = np.asarray(samples[wrap.name], float)
            if len(wrap_pts) and np.linalg.norm(wrap_pts[0] - su_pts[-1]) > np.linalg.norm(
                wrap_pts[-1] - su_pts[-1]
            ):
                wrap_pts = wrap_pts[::-1].copy()
        parts: list[np.ndarray] = []
        _append_ring(parts, su_pts)
        _append_ring(parts, aft_pts)
        _append_ring(parts, wrap_pts)
        _append_ring(parts, pr_pts[::-1])
        if fwd_pts is not None:
            _append_ring(parts, fwd_pts[::-1])
        if not parts:
            continue
        ring = np.vstack(parts)
        kept = [ring[0]]
        for p in ring[1:]:
            if not np.allclose(kept[-1], p, atol=1e-12, rtol=0.0):
                kept.append(p)
        ring = np.asarray(kept, float)
        if len(ring) > 2 and np.allclose(ring[0], ring[-1], atol=1e-12, rtol=0.0):
            ring = ring[:-1]
        if len(ring) < 3:
            continue
        name = label if n_cells == 1 else f"{label}:cell_{k}"
        out.append((name, ring))
    return out


def rib_face_polygons(geo, *, n_chord: int = 12) -> list[tuple[str, np.ndarray]]:
    """Rib/tip n-gons at each ``TransverseMember.z``, vertices from the mapped loft."""
    out: list[tuple[str, np.ndarray]] = []
    graph = getattr(geo, "region_graph", None)
    if graph is None:
        return out
    for member in getattr(geo, "transverse", ()) or ():
        cut = section_cut_near(geo, float(member.z))
        if cut is None:
            continue
        sg = cut.shell_graph()
        tmpl = _template_from_cut(graph, cut, float(member.z), int(n_chord))
        samples = sample_region_polylines(cut, sg, tmpl, graph)
        out.extend(rib_cell_polygons_from_samples(tmpl, samples, member))
    return out


def rib_face_triangles(geo, *, n_chord: int = 12) -> list[tuple[str, np.ndarray]]:
    """STL-only tessellation: fan n-gons from the first loft node (no centroid)."""
    faces = []
    for name, ring in rib_face_polygons(geo, n_chord=n_chord):
        body = ring[:-1] if len(ring) > 3 and np.allclose(ring[0], ring[-1]) else ring
        if len(body) < 3:
            continue
        for i in range(1, len(body) - 1):
            faces.append((name, np.stack([body[0], body[i], body[i + 1]])))
    return faces


def section_cut_near(geo, z: float):
    """``section_at(z)``, or a neighbour if that station is infeasible."""
    from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

    z = float(z)
    try:
        return geo.section_at(z)
    except InfeasibleSectionGeometry:
        pass
    for dz in (1e-6, -1e-6, 1e-4, -1e-4, 1e-3, -1e-3):
        try:
            return geo.section_at(z + dz)
        except (InfeasibleSectionGeometry, ValueError):
            continue
    return None


def _live_web_fracs(cut) -> list[float]:
    """Twist-aware live web chord fractions from the shell graph."""
    live = []
    for run in cut.shell_graph().webs:
        pts = np.asarray(run.points, float)[:, :2]
        live.append(float(np.clip(_chord_fracs(cut, pts[:1])[0], 0.0, 1.0)))
    return live


def _template_from_cut(graph: WallRegionGraph, cut, z: float, n_chord: int) -> ShellTopologyTemplate:
    dropped = getattr(cut, "dropped_web_indices", None)
    return ShellTopologyTemplate.from_region_graph(
        graph, float(z), n_chord=int(n_chord),
        live_web_fracs=_live_web_fracs(cut) or None, t_c=_cut_tc(cut),
        dropped_web_indices=None if dropped is None else tuple(int(i) for i in dropped),
    )


def templates_along(geo, zs: Sequence[float], n_chord: int = 12) -> list[ShellTopologyTemplate]:
    graph = geo.region_graph
    out = []
    for z in zs:
        cut = section_cut_near(geo, float(z))
        if cut is None:
            continue
        out.append(_template_from_cut(graph, cut, float(z), n_chord))
    return out
