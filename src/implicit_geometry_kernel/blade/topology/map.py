"""Sample authored regions at a station and map loft strips between stations."""
from __future__ import annotations

from typing import Sequence

import numpy as np

from implicit_geometry_kernel.blade.locked_loft import arc_resample_open
from implicit_geometry_kernel.blade.topology.sizing import resample_open_h, station_live_thickness
from implicit_geometry_kernel.blade.shell_graph import ShellGraph, ShellRun
from implicit_geometry_kernel.blade.topology.graph import (
    ShellVertexId,
    WallRegionGraph,
    _match_web_indices,
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
    _te_wrap_arc as _te_wrap_arc,
    _xyz_at_frac,
    te_is_blunt,
)
from implicit_geometry_kernel.blade.topology.map_loft import (
    _resample_keys,
)


def place_station_keys(cut, graph: ShellGraph, region_graph: WallRegionGraph, z: float,
                       te_blunt: bool | None = None,
                       ) -> dict[tuple[str, str], np.ndarray]:
    """Project authored anchors onto the TE-started skin-mid loop.

    Keys are ``(anchor_name, side)`` with ``side`` in ``suction`` / ``pressure``.
    Blunt TE keys are the real edge corners; sharp TE is one shared ridge.
    """
    pts = _skin_mid_loop(graph, z)
    keys: dict[tuple[str, str], np.ndarray] = {}
    if len(pts) < 3:
        return keys
    i_le, i_su, i_pr = _te_corner_indices(pts, cut)
    blunt = te_is_blunt(cut, pts) if te_blunt is None else bool(te_blunt)
    if not blunt:
        i_ridge = int(np.argmax(_chord_fracs(cut, pts)))
        i_su = i_pr = i_ridge
    su, pr = _side_arcs(pts, cut, i_le, i_su, i_pr)
    su_f, pr_f = _chord_fracs(cut, su), _chord_fracs(cut, pr)
    anchors = dict(region_graph.evaluate_anchors(z))
    keys[("te", "suction")] = pts[int(i_su)].copy()
    keys[("te", "pressure")] = pts[int(i_pr)].copy()
    for side, arc, fr in (("suction", su, su_f), ("pressure", pr, pr_f)):
        keys[("le", side)] = arc[-1].copy()
        for name, frac in anchors.items():
            if name in ("le", "te"):
                continue
            keys[(name, side)] = _xyz_at_frac(arc, fr, frac)
    for name, run in _web_name_runs(graph.webs, cut, anchors):
        mid = _as_xyz(np.asarray(run.points, float), float(z))
        if len(mid) < 2:
            continue
        su_ft, pr_ft = _web_feet_on_sides(mid, su, pr)
        keys[(name, "suction")] = np.asarray(su_ft, float)
        keys[(name, "pressure")] = np.asarray(pr_ft, float)
    return keys


def _key_name(vid: ShellVertexId) -> str:
    if vid == ("TE",) or (len(vid) > 0 and vid[0] == "TE") or (len(vid) > 1 and vid[1] == "te"):
        return "te"
    if vid[0] == "anchor":
        return str(vid[1])
    return "te"


def sample_region_polylines(cut, graph: ShellGraph, template: ShellTopologyTemplate,
                            region_graph: WallRegionGraph, *,
                            size_field: bool = False) -> dict[str, np.ndarray]:
    """Resample authored regions as walks between shared loop keys.

    Default holds ``n_chord``. ``size_field=True`` spaces nodes by in-plane ``h``
    (CSG viz); template ``n_chord`` is unused then.
    """
    n = max(int(template.n_chord), 2)
    z = float(template.z)
    wrap_on = any(b.wrap for b in template.loop)
    keys = place_station_keys(cut, graph, region_graph, z, te_blunt=wrap_on)
    out: dict[str, np.ndarray] = {}
    pts = _skin_mid_loop(graph, z)
    if len(pts) < 3 or not keys:
        return out
    i_le, i_su, i_pr = _te_corner_indices(pts, cut)
    if not wrap_on:
        i_ridge = int(np.argmax(_chord_fracs(cut, pts)))
        i_su = i_pr = i_ridge
    su, pr = _side_arcs(pts, cut, i_le, i_su, i_pr)
    arcs = {"suction": su, "pressure": pr}

    chord = max(float(np.linalg.norm(su[0] - su[-1])), 1e-9)
    c_h = max(float(getattr(cut, "chord", 0.0) or 0.0), chord)
    t_live = station_live_thickness(cut) if size_field else 0.0
    for branch in template.loop:
        if branch.wrap:
            a = keys.get(("te", "suction"), su[0])
            b = keys.get(("te", "pressure"), pr[0])
            if float(np.linalg.norm(a - b)) < 1e-6 * chord:
                continue
            i_le, i_te, _ = _te_corner_indices(pts, cut)
            wrap_arc = _te_wrap_arc(pts, a, b, i_le, i_te)
            if wrap_arc is None:
                continue
            if float(np.linalg.norm(wrap_arc[0] - a)) > float(np.linalg.norm(wrap_arc[0] - b)):
                wrap_arc = wrap_arc[::-1]
            if size_field:
                sampled = resample_open_h(wrap_arc, c_h, t_live)
            else:
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
            p0, p1, arcs[side], n_use,
            bunch_end=False if size_field else branch.name.startswith("te_"),
            size_field=size_field, chord=c_h, t=t_live,
        )

    name_to_run = dict(_web_name_runs(graph.webs, cut, dict(region_graph.evaluate_anchors(z))))
    for web in template.webs:
        su_name = _key_name(web.suction)
        pr_name = _key_name(web.pressure)
        su_ft = keys.get((su_name, "suction"))
        pr_ft = keys.get((pr_name, "pressure"))
        if su_ft is None or pr_ft is None:
            continue
        match = name_to_run.get(f"web_{web.index}") or name_to_run.get(su_name)
        ends = np.vstack([su_ft, pr_ft])
        if match is None:
            sampled = resample_open_h(ends, c_h, t_live) if size_field else arc_resample_open(ends, n)
        else:
            mid = _as_xyz(np.asarray(match.points, float), z)
            sampled = resample_open_h(mid, c_h, t_live) if size_field else arc_resample_open(mid, n)
        if np.linalg.norm(sampled[0] - su_ft) <= np.linalg.norm(sampled[0] - pr_ft):
            sampled[0], sampled[-1] = su_ft, pr_ft
        else:
            sampled[0], sampled[-1] = pr_ft, su_ft
        out[web.name] = sampled
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


def _web_feet_on_sides(mid: np.ndarray, su: np.ndarray, pr: np.ndarray
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Assign emitted web endpoints to suction / pressure by nearest skin arc."""
    e0, e1 = mid[0], mid[-1]

    def _d(arc: np.ndarray, p: np.ndarray) -> float:
        return float(np.min(np.linalg.norm(arc[:, :2] - np.asarray(p, float)[:2], axis=1)))

    if _d(su, e0) + _d(pr, e1) <= _d(su, e1) + _d(pr, e0):
        return e0, e1
    return e1, e0


def _web_run_frac(cut, run: ShellRun) -> float | None:
    pts = np.asarray(run.points, float)
    if pts.ndim != 2 or len(pts) < 2:
        return None
    return float(np.mean(_chord_fracs(cut, pts[:, :2])))


def _web_name_runs(
    runs: Sequence[ShellRun],
    cut,
    anchors: dict[str, float],
) -> list[tuple[str, ShellRun]]:
    """Pair ``web_*`` anchors to emitted runs, unique and order-preserving in chord fraction."""
    names = [n for n in anchors if str(n).startswith("web_")]
    live: list[tuple[float, ShellRun]] = []
    for run in runs:
        frac = _web_run_frac(cut, run)
        if frac is not None:
            live.append((frac, run))
    if not names or not live:
        return []
    names_sorted = sorted(names, key=lambda n: float(anchors[n]))
    live.sort(key=lambda t: t[0])
    authored = [float(anchors[n]) for n in names_sorted]
    keep = set(_match_web_indices([f for f, _ in live], authored, list(range(len(names_sorted)))))
    kept = [names_sorted[i] for i in sorted(keep)]
    n = min(len(kept), len(live))
    return [(kept[i], live[i][1]) for i in range(n)]


