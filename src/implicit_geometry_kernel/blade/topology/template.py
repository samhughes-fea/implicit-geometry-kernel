"""Station shell-topology templates and compatible span zones."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from implicit_geometry_kernel.blade.topology.graph import (
    ShellVertexId,
    WallRegionGraph,
    _cluster_web_members,
    _match_web_indices,
    live_cap_windows,
)
from implicit_geometry_kernel.blade.topology.layout_compile import _regions_from_ordered
from implicit_geometry_kernel.blade.topology.signatures import signatures_compatible

def _regions_per_side(graph: WallRegionGraph, fracs: Mapping[str, float], ordered: Sequence[str],
                      side_fracs: Mapping[str, Mapping[str, float]], *, z: float,
                      active_webs: Sequence) -> tuple:
    """TE wrap, then each side's regions in that side's own anchor order."""
    out: list = []
    for side in ("suction", "pressure"):
        fr = dict(fracs)
        fr.update({k: float(v) for k, v in side_fracs.get(side, {}).items() if k in fr})
        side_order = sorted(ordered, key=lambda name: float(fr[name]))
        windows = live_cap_windows(graph, fr, side_order, z=z, active_webs=active_webs)
        regions = _regions_from_ordered(side_order, fracs=fr, windows=windows)
        if not out:
            out.extend(r for r in regions if r.wrap)
        out.extend(r for r in regions if not r.wrap and r.side == side)
    return tuple(out)


@dataclass(frozen=True)
class LoopBranch:
    name: str
    start: ShellVertexId
    end: ShellVertexId
    n: int
    side: str
    wrap: bool = False


def _loop_vid(name: str, side: str = "") -> ShellVertexId:
    """Skin-loop key id. Web ends reuse these (SHARED); they are not extra vertices."""
    if name in ("te", "TE"):
        if side in ("suction", "pressure"):
            return ("TE", side)
        return ("TE",)
    if name == "le":
        return ("anchor", "le")
    if side in ("suction", "pressure"):
        return ("anchor", name, side)
    return ("anchor", name)


@dataclass(frozen=True)
class WebBranch:
    """Web midline chord. ``suction`` / ``pressure`` are the skin-loop keys (SHARED)."""

    name: str
    index: int
    suction: ShellVertexId
    pressure: ShellVertexId
    n: int


@dataclass(frozen=True)
class ShellTopologyTemplate:
    """Station template: loop branches = regions; webs = members active at z."""

    z: float
    signature: tuple
    loop: tuple[LoopBranch, ...]
    webs: tuple[WebBranch, ...]
    vertices: tuple[ShellVertexId, ...]
    n_chord: int = 2

    @classmethod
    def from_region_graph(cls, graph: WallRegionGraph, z: float, n_chord: int = 2,
                          live_web_fracs: Sequence[float] | None = None,
                          t_c: float | None = None,
                          dropped_web_indices: Sequence[int] | None = None,
                          te_blunt: bool | None = None,
                          side_fracs: Mapping[str, Mapping[str, float]] | None = None,
                          ) -> "ShellTopologyTemplate":
        """Station template from the region graph at ``z``.

        ``side_fracs`` (optional) gives each side's own chord fraction of an
        anchor (``{"suction": {name: frac}, "pressure": {...}}``), as the station
        keys sit on the cut. A cap edge and a web foot can then take different
        orders on the two sides (upright webs on a twisted section), and each
        side's regions follow its own order.
        """
        n = max(int(n_chord), 2)
        te_arc = float(graph.te_arc)
        authored = tuple(graph.webs_at(z))
        if dropped_web_indices is not None:
            dropped = {int(i) for i in dropped_web_indices}
            active_webs = tuple(w for w in authored if w.index not in dropped)
        else:
            authored = tuple(w for w in authored if te_arc < w.at(z) < 1.0 - te_arc)
            authored = _cluster_web_members(
                authored, z, min_sep=max(0.08, 2.0 * float(graph.min_gap)), t_c=t_c)
            live = None
            if live_web_fracs is not None:
                live = [float(f) for f in live_web_fracs if te_arc < float(f) < 1.0 - te_arc]
            if live:
                keep = set(_match_web_indices(live, [w.at(z) for w in authored],
                                              [w.index for w in authored]))
                active_webs = tuple(w for w in authored if w.index in keep)
            else:
                active_webs = authored
        fracs = dict(graph.evaluate_anchors(z))
        if dropped_web_indices is not None:
            keep_web = {f"web_{w.index}" for w in active_webs}
            ordered = [n for n in fracs if not str(n).startswith("web_") or n in keep_web]
        else:
            ordered = list(fracs)
        if side_fracs is None:
            windows = live_cap_windows(
                graph, fracs, ordered, z=float(z), active_webs=active_webs,
            )
            regions = _regions_from_ordered(ordered, fracs=fracs, windows=windows)
        else:
            regions = _regions_per_side(graph, fracs, ordered, side_fracs,
                                        z=float(z), active_webs=active_webs)
        vertices: list[ShellVertexId] = [
            ("TE",), ("TE", "suction"), ("TE", "pressure"), ("anchor", "le"),
        ]
        loop = []
        wrap_on = True if te_blunt is None else bool(te_blunt)
        for region in regions:
            start: ShellVertexId
            end: ShellVertexId
            if region.wrap:
                if not wrap_on:
                    continue
                start, end = ("TE", "suction"), ("TE", "pressure")
            else:
                start = _loop_vid(region.start, region.side)
                end = _loop_vid(region.end, region.side)
                if not wrap_on:
                    if region.start in ("te", "TE"):
                        start = ("TE",)
                    if region.end in ("te", "TE"):
                        end = ("TE",)
            loop.append(LoopBranch(region.name, start, end, n, region.side, region.wrap))
            for vid in (start, end):
                if vid not in vertices:
                    vertices.append(vid)
        webs = []
        for w in active_webs:
            su = _loop_vid(f"web_{w.index}", "suction")
            pr = _loop_vid(f"web_{w.index}", "pressure")
            for vid in (su, pr):
                if vid not in vertices:
                    vertices.append(vid)
            webs.append(WebBranch(w.name or f"web_{w.index}", w.index, su, pr, n))
        sig = (tuple(b.name for b in loop), tuple(w.index for w in webs))
        return cls(float(z), sig, tuple(loop), tuple(webs), tuple(dict.fromkeys(vertices)), n)

    def map_vertex(self, vid: ShellVertexId, other: "ShellTopologyTemplate") -> ShellVertexId:
        if vid in other.vertices:
            return vid
        if vid == ("TE",) or (len(vid) > 0 and vid[0] == "TE") or (
                vid[0] == "anchor" and len(vid) > 1 and vid[1] == "te"):
            side = vid[1] if len(vid) > 1 and vid[0] == "TE" and vid[1] in ("suction", "pressure") else ""
            te_cands: tuple[ShellVertexId | None, ...] = (
                ("TE", side) if side else None, ("TE", "suction"), ("TE", "pressure"), ("TE",),
            )
            for cand in te_cands:
                if cand is not None and cand in other.vertices:
                    return cand
            return ("TE",)
        if vid[0] == "anchor":
            name = vid[1]
            side = vid[2] if len(vid) > 2 else ""
            if str(name).startswith("web_"):
                try:
                    idx = int(str(name).split("_")[-1])
                except ValueError:
                    idx = 0
                others = [v for v in other.vertices
                          if len(v) > 2 and v[0] == "anchor" and str(v[1]).startswith("web_")
                          and v[2] == side]
                if others:
                    def _idx(v):
                        try:
                            return int(str(v[1]).split("_")[-1])
                        except ValueError:
                            return 0
                    rest = [_idx(v) for v in others]
                    if idx < min(rest) and ("anchor", "le") in other.vertices:
                        return ("anchor", "le")
                    if idx > max(rest):
                        te = ("TE", side) if side else ("TE",)
                        if te in other.vertices:
                            return te
                    return min(others, key=lambda v: abs(_idx(v) - idx))
            if name.startswith("cap_"):
                for wname in ("web_0", "web_1", "le"):
                    cap_cand = _loop_vid(wname, side)
                    if cap_cand in other.vertices:
                        return cap_cand
            if side and ("TE", side) in other.vertices:
                return ("TE", side)
            if ("TE",) in other.vertices:
                return ("TE",)
        if not signatures_compatible(self.signature, other.signature):
            raise KeyError(vid)
        return ("TE",)


@dataclass(frozen=True)
class ZonedShellTopology:
    """Zones of compatible templates; 1-station blips are transition bands."""

    templates: tuple[ShellTopologyTemplate, ...]
    zones: tuple[slice, ...]


def zone_region_graphs(templates: Sequence[ShellTopologyTemplate]) -> ZonedShellTopology:
    tpls = tuple(templates)
    if not tpls:
        return ZonedShellTopology((), ())
    sigs = [t.signature for t in tpls]
    n = len(sigs)
    flags = [True] * n
    for i in range(1, n - 1):
        if sigs[i] != sigs[i - 1] and sigs[i - 1] == sigs[i + 1]:
            flags[i] = False
    zones: list[slice] = []
    start = 0
    cur = _zone_sig(sigs, flags, 0)
    for i in range(1, n):
        nxt = _zone_sig(sigs, flags, i)
        if nxt != cur and flags[i]:
            if not signatures_compatible(cur, nxt):
                zones.append(slice(start, i))
                start = i
                cur = nxt
            else:
                cur = nxt
        elif flags[i]:
            cur = nxt
    zones.append(slice(start, n))
    return ZonedShellTopology(tpls, tuple(zones))


def _zone_sig(sigs, flags, i):
    if flags[i]:
        return sigs[i]
    if i > 0:
        return sigs[i - 1]
    return sigs[min(i + 1, len(sigs) - 1)]

