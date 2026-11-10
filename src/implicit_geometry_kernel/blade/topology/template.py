"""Station shell-topology templates and compatible span zones."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from implicit_geometry_kernel.blade.topology.graph import (
    ShellVertexId,
    WallRegionGraph,
    _cluster_web_members,
    _match_web_indices,
    _regions_from_ordered,
    live_cap_windows,
    signatures_compatible,
)

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
                          ) -> "ShellTopologyTemplate":
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
        windows = live_cap_windows(
            graph, fracs, ordered, z=float(z), active_webs=active_webs,
        )
        regions = _regions_from_ordered(ordered, fracs=fracs, windows=windows)
        vertices: list[ShellVertexId] = [
            ("TE",), ("TE", "suction"), ("TE", "pressure"), ("anchor", "le"),
        ]
        loop = []
        wrap_on = True if te_blunt is None else bool(te_blunt)
        for region in regions:
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
            for cand in (("TE", side) if side else None, ("TE", "suction"), ("TE", "pressure"), ("TE",)):
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
                    cand = _loop_vid(wname, side)
                    if cand in other.vertices:
                        return cand
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

