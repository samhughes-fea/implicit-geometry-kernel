"""Sample authored regions at a station and map loft strips between stations."""
from __future__ import annotations

from typing import Sequence

import numpy as np

from implicit_geometry_kernel.blade.locked_loft import arc_resample_open
from implicit_geometry_kernel.blade.shell_graph import ShellGraph, ShellRun
from implicit_geometry_kernel.blade.topology.graph import (
    ShellVertexId,
    WallRegionGraph,
)
from implicit_geometry_kernel.blade.topology.template import (
    ShellTopologyTemplate,
)
from implicit_geometry_kernel.blade.topology.map_loft import (
    _append_ring as _append_ring,
    _quad_bowtie as _quad_bowtie,
    _web_pts_from_su as _web_pts_from_su,
    loft_mapped_strips as loft_mapped_strips,
)
from implicit_geometry_kernel.blade.topology.map_walk import (
    _as_xyz,
    _chord_fracs as _chord_fracs,
    _oml_te_xyz as _oml_te_xyz,
    _side_arcs,
    _skin_mid_loop,
    _te_corner_indices,
    _te_frac_drop,
    _te_wrap_arc as _te_wrap_arc,
    _xyz_at_frac,
)
from implicit_geometry_kernel.blade.topology.map_loft import (
    _resample_keys,
)


def place_station_keys(cut, graph: ShellGraph, region_graph: WallRegionGraph, z: float
                       ) -> dict[tuple[str, str], np.ndarray]:
    """Project authored anchors onto the TE-started skin-mid loop.

    Keys are ``(anchor_name, side)`` with ``side`` in ``suction`` / ``pressure``.
    TE is the start of each TE→LE side arc (not a nearest-OML snap).
    """
    pts = _skin_mid_loop(graph, z)
    keys: dict[tuple[str, str], np.ndarray] = {}
    if len(pts) < 3:
        return keys
    su, pr = _side_arcs(pts, cut)
    su_f, pr_f = _chord_fracs(cut, su), _chord_fracs(cut, pr)
    anchors = dict(region_graph.evaluate_anchors(z))
    drop = _te_frac_drop(cut)
    for side, arc, fr in (("suction", su, su_f), ("pressure", pr, pr_f)):
        f0 = float(fr[0]) if len(fr) else 1.0
        f_c = f0 - drop
        if len(arc) and f_c > float(fr[-1]) + 1e-12:
            keys[("te", side)] = _xyz_at_frac(arc, fr, f_c)
        else:
            keys[("te", side)] = arc[0].copy()
        keys[("le", side)] = arc[-1].copy()
        for name, frac in anchors.items():
            if name in ("le", "te"):
                continue
            keys[(name, side)] = _xyz_at_frac(arc, fr, frac)
    return keys


def _key_name(vid: ShellVertexId) -> str:
    if vid == ("TE",) or (len(vid) > 0 and vid[0] == "TE") or (len(vid) > 1 and vid[1] == "te"):
        return "te"
    if vid[0] == "anchor":
        return str(vid[1])
    return "te"


def sample_region_polylines(cut, graph: ShellGraph, template: ShellTopologyTemplate,
                            region_graph: WallRegionGraph) -> dict[str, np.ndarray]:
    """Resample authored regions as walks between shared loop keys (held ``n_chord``)."""
    n = max(int(template.n_chord), 2)
    z = float(template.z)
    keys = place_station_keys(cut, graph, region_graph, z)
    out: dict[str, np.ndarray] = {}
    pts = _skin_mid_loop(graph, z)
    if len(pts) < 3 or not keys:
        return out
    su, pr = _side_arcs(pts, cut)
    arcs = {"suction": su, "pressure": pr}

    chord = max(float(np.linalg.norm(su[0] - su[-1])), 1e-9)
    for branch in template.loop:
        if branch.wrap:
            a = keys.get(("te", "suction"), su[0])
            b = keys.get(("te", "pressure"), pr[0])
            i_le, i_te, _ = _te_corner_indices(pts, cut)
            if float(np.linalg.norm(a - b)) < 1e-6 * chord:
                ridge = pts[int(i_te)]
                if float(np.linalg.norm(a[:2] - ridge[:2])) < 1e-6 * chord:
                    continue
            wrap_arc = _te_wrap_arc(pts, a, b, i_le, i_te)
            if wrap_arc is None:
                continue
            if float(np.linalg.norm(wrap_arc[0] - a)) > float(np.linalg.norm(wrap_arc[0] - b)):
                wrap_arc = wrap_arc[::-1]
            sampled = arc_resample_open(wrap_arc, max(n, 4))
            sampled[0], sampled[-1] = a, b
            out[branch.name] = sampled
            continue
        side = branch.side if branch.side in arcs else "suction"
        n0, n1 = _key_name(branch.start), _key_name(branch.end)
        p0 = keys.get((n0, side))
        p1 = keys.get((n1, side))
        if p0 is None or p1 is None:
            continue
        if float(np.linalg.norm(p0 - p1)) < 1e-6 * chord:
            continue
        n_use = n + n // 2 if branch.name.startswith("te_") else n
        out[branch.name] = _resample_keys(
            p0, p1, arcs[side], n_use, bunch_end=branch.name.startswith("te_"),
        )

    by_frac = {w.index: w.at(z) for w in region_graph.webs_at(z)}
    for web in template.webs:
        su_name = _key_name(web.suction)
        pr_name = _key_name(web.pressure)
        su_ft = keys.get((su_name, "suction"))
        pr_ft = keys.get((pr_name, "pressure"))
        if su_ft is None or pr_ft is None:
            continue
        match = _web_run_at_frac(graph.webs, cut, by_frac.get(web.index))
        ends = np.vstack([su_ft, pr_ft])
        if match is None:
            sampled = arc_resample_open(ends, n)
        else:
            mid = _as_xyz(np.asarray(match.points, float), z)
            d0 = min(np.linalg.norm(mid[0] - su_ft), np.linalg.norm(mid[-1] - su_ft))
            d1 = min(np.linalg.norm(mid[0] - pr_ft), np.linalg.norm(mid[-1] - pr_ft))
            span = max(float(np.linalg.norm(su_ft - pr_ft)), 1e-12)
            if d0 > 0.05 * span or d1 > 0.05 * span:
                sampled = arc_resample_open(ends, n)
            else:
                sampled = arc_resample_open(mid, n)
        if np.linalg.norm(sampled[0] - su_ft) <= np.linalg.norm(sampled[0] - pr_ft):
            sampled[0], sampled[-1] = su_ft, pr_ft
        else:
            sampled[0], sampled[-1] = pr_ft, su_ft
        out[web.name] = sampled
        _hold_key_on_side(out, template, "suction", su_ft)
        _hold_key_on_side(out, template, "pressure", pr_ft)
    for side, tag in (("suction", "su"), ("pressure", "pr")):
        te_xyz = keys.get(("te", side))
        name = f"te_{tag}"
        if te_xyz is not None and name in out and len(out[name]):
            out[name][-1] = np.asarray(te_xyz, float)
    skin = [out[b.name] for b in template.loop if b.name in out]
    if skin:
        wall = np.vstack(skin)
        chord = max(float(np.linalg.norm(su[0] - su[-1])), 1e-9)
        for web in template.webs:
            if web.name not in out:
                continue
            pts = np.asarray(out[web.name], float)
            g0 = float(np.min(np.linalg.norm(wall - pts[0], axis=1)))
            g1 = float(np.min(np.linalg.norm(wall - pts[-1], axis=1)))
            if g0 > 1e-3 * chord or g1 > 1e-3 * chord:
                del out[web.name]
    return out


def _hold_key_on_side(out: dict[str, np.ndarray], template: ShellTopologyTemplate,
                      side: str, xyz: np.ndarray) -> None:
    """Replace the nearest *interior* skin sample on ``side`` with the web key."""
    best_name, best_i, best_d = None, 0, float("inf")
    for branch in template.loop:
        if branch.wrap or branch.side != side or branch.name not in out:
            continue
        pts = np.asarray(out[branch.name], float)
        if len(pts) < 3:
            continue
        d = np.linalg.norm(pts - xyz, axis=1)
        d[0] = d[-1] = np.inf
        i = int(np.argmin(d))
        if np.isfinite(d[i]) and float(d[i]) < best_d:
            best_name, best_i, best_d = branch.name, i, float(d[i])
    if best_name is not None:
        out[best_name][best_i] = np.asarray(xyz, float)


def _web_run_at_frac(runs: Sequence[ShellRun], cut, target: float | None) -> ShellRun | None:
    if target is None or not runs:
        return None
    best, best_d = None, float("inf")
    for run in runs:
        pts = np.asarray(run.points, float)
        if pts.ndim != 2 or len(pts) < 2:
            continue
        f = float(np.mean(_chord_fracs(cut, pts[:, :2])))
        d = abs(f - float(target))
        if d < best_d:
            best, best_d = run, d
    return best


