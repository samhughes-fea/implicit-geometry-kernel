"""Authored wall-region graph and mapped loft (tidalshell identity, no mesher).

Loft identity is a :class:`WallRegionGraph` (anchors, named loop regions, webs),
not inferred ``ShellRun`` names. ``shared_loft_names`` stays the Wave 2 event
helper; figure loft and ``write_mid_surfaces`` join with ``map_vertex``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.locked_loft import arc_resample_open
from implicit_geometry_kernel.blade.shell_graph import ShellGraph, ShellRun
from implicit_geometry_kernel.blade.stations import TransverseMember

TE_ARC = 0.03
MIN_GAP = 0.002

FracLaw = Callable[[float], float]
ShellVertexId = tuple


def _interp_table(z: float, zs: Sequence[float], vs: Sequence[float]) -> float:
    zz = np.asarray(zs, float).reshape(-1)
    vv = np.asarray(vs, float).reshape(-1)
    if zz.size == 0:
        return 0.0
    if zz.size == 1:
        return float(vv[0])
    return float(np.interp(float(z), zz, vv))


def _as_law(value) -> FracLaw:
    if callable(value):
        return value
    if isinstance(value, (tuple, list)) and len(value) == 2 and np.ndim(value[0]) == 1:
        zs, vs = np.asarray(value[0], float), np.asarray(value[1], float)
        return lambda z, _z=zs, _v=vs: _interp_table(z, _z, _v)
    c = float(value)
    return lambda z, _c=c: _c


@dataclass(frozen=True)
class Anchor:
    """A feature line on the skin-mid loop, as chord fraction from the LE."""

    name: str
    frac: object
    z0: float | None = None
    z1: float | None = None

    def at(self, z: float) -> float:
        return float(np.clip(_as_law(self.frac)(z), 0.0, 1.0))

    def active(self, z: float) -> bool:
        if self.z0 is not None and float(z) < float(self.z0) - 1e-12:
            return False
        if self.z1 is not None and float(z) > float(self.z1) + 1e-12:
            return False
        return True


@dataclass(frozen=True)
class Region:
    """Coverage between two anchors on one skin side (or the TE wrap)."""

    name: str
    start: str
    end: str
    is_cap: bool = False
    side: str = "suction"
    wrap: bool = False


@dataclass(frozen=True)
class WebMember:
    """A shear web: published index, chord-fraction law, optional span window."""

    index: int
    frac: object
    z0: float | None = None
    z1: float | None = None
    name: str = ""

    def at(self, z: float) -> float:
        return float(np.clip(_as_law(self.frac)(z), 0.0, 1.0))

    def active(self, z: float) -> bool:
        if self.z0 is not None and float(z) < float(self.z0) - 1e-12:
            return False
        if self.z1 is not None and float(z) > float(self.z1) + 1e-12:
            return False
        return True


@dataclass(frozen=True)
class WallRegionGraph:
    """Authored feature lines for one blade / concept instance."""

    anchors: tuple[Anchor, ...]
    regions: tuple[Region, ...]
    webs: tuple[WebMember, ...]
    te_arc: float = TE_ARC
    min_gap: float = MIN_GAP

    def evaluate_anchors(self, z: float) -> tuple[tuple[str, float], ...]:
        """Active anchors at ``z``, TE-merged and ``min_gap`` regularised."""
        raw = [(a.name, a.at(z)) for a in self.anchors if a.active(z)]
        return _regularise_anchors(raw, te_arc=self.te_arc, min_gap=self.min_gap)

    def regions_at(self, z: float) -> tuple[Region, ...]:
        """Authored regions (constant names). Span filters apply only to webs."""
        del z
        return self.regions

    def webs_at(self, z: float) -> tuple[WebMember, ...]:
        return tuple(w for w in self.webs if w.active(z) and 0.0 < w.at(z) < 1.0)

    def signature_at(self, z: float, web_fracs: Sequence[float] | None = None) -> tuple:
        regions = tuple(r.name for r in self.regions_at(z))
        if web_fracs is None:
            webs = tuple(w.index for w in self.webs_at(z))
        else:
            webs = tuple(_match_web_indices(web_fracs, [w.at(z) for w in self.webs if w.active(z)],
                                            [w.index for w in self.webs if w.active(z)]))
        return (regions, webs)

    @classmethod
    def from_section_layout(cls, layout, *, te_arc: float = TE_ARC,
                            min_gap: float = MIN_GAP) -> "WallRegionGraph":
        if isinstance(layout, str):
            from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
            layout = SectionLayout.from_catalogue_key(layout)
        elif not hasattr(layout, "cap_placement"):
            from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
            from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
            cap = CapPlacementRule.from_structural_family(
                getattr(layout, "structural_family", "D"),
                spar_cap_positioning=getattr(layout, "spar_cap_positioning", None),
            )
            layout = SectionLayout(
                web_chord_fracs=tuple(float(f) for f in getattr(layout, "web_chord_fracs", ()) or ()),
                orientation=getattr(layout, "orientation", getattr(layout, "web_orientation", "flapwise")),
                cap_placement=cap,
                catalogue_key=getattr(layout, "key", None),
                structural_family=str(getattr(layout, "structural_family", "") or ""),
            )
        webs_frac = tuple(float(f) for f in getattr(layout, "web_chord_fracs", ()) or ())
        cap = getattr(layout, "cap_placement", None)
        anchors, regions, webs = _layout_parts(webs_frac, cap)
        return cls(anchors, regions, webs, te_arc=float(te_arc), min_gap=float(min_gap))

    @classmethod
    def from_windio_structure(cls, structure: Mapping, z_of_s=None, *,
                              te_arc: float = TE_ARC,
                              min_gap: float = MIN_GAP) -> "WallRegionGraph":
        from implicit_geometry_kernel.blade.windio import (
            cap_window_table_from_windio_structure,
            web_position_table_from_windio_structure,
        )

        from implicit_geometry_kernel.blade.windio import web_span_windows_from_windio_structure

        z_table = None if z_of_s is None else np.asarray(z_of_s, float)
        anchors: list[Anchor] = [Anchor("le", 0.0), Anchor("te", 1.0)]
        web_table = web_position_table_from_windio_structure(structure, z_table)
        windows_span = web_span_windows_from_windio_structure(structure, z_table)
        webs: list[WebMember] = []
        if web_table is not None:
            zg, cols = web_table
            names = [str(w.get("name", f"web_{j}"))
                     for j, w in enumerate(structure.get("webs") or ())
                     if isinstance(w, Mapping)]
            for j in range(cols.shape[1]):
                z0 = z1 = None
                if j < len(windows_span) and windows_span[j] is not None:
                    z0, z1 = windows_span[j]
                law = (zg, cols[:, j])
                name = names[j] if j < len(names) else f"web_{j}"
                webs.append(WebMember(int(j), law, z0, z1, name=name))
                anchors.append(Anchor(f"web_{j}", law, z0, z1))
        windows = cap_window_table_from_windio_structure(structure, z_table)
        if windows is not None:
            zg, vals = windows
            anchors.append(Anchor("cap_lo", (zg, vals[:, 0])))
            anchors.append(Anchor("cap_hi", (zg, vals[:, 1])))
            web_fracs = tuple(float(np.median(web_table[1][:, j])) for j in range(web_table[1].shape[1])) if web_table is not None else ()
            dummy_cap = _CapWindow(tuple(float(np.median(windows[1][:, k])) for k in (0, 1)))
            _, regions, _ = _layout_parts(web_fracs, dummy_cap)
        else:
            _, regions, _ = _layout_parts(tuple(w.at(0.0) for w in webs) if webs else (), None)
        return cls(tuple(anchors), regions, tuple(webs), te_arc=float(te_arc),
                   min_gap=float(min_gap))


class _CapWindow:
    kind = "fixed_band"
    chord_window = None
    half_width = None

    def __init__(self, window: tuple[float, float]):
        self.chord_window = (min(window), max(window))
        self.half_width = 0.5 * (self.chord_window[1] - self.chord_window[0])


def _layout_parts(web_fracs: Sequence[float], cap) -> tuple[tuple[Anchor, ...], tuple[Region, ...], tuple[WebMember, ...]]:
    webs_frac = tuple(sorted(float(f) for f in web_fracs if 0.0 < float(f) < 1.0))
    points: list[tuple[str, float]] = [("le", 0.0), ("te", 1.0)]
    webs = []
    for i, f in enumerate(webs_frac):
        points.append((f"web_{i}", f))
        webs.append(WebMember(i, f, name=f"web_{i}"))
    cap_lo = cap_hi = None
    kind = getattr(cap, "kind", "none") if cap is not None else "none"
    window = getattr(cap, "chord_window", None) if cap is not None else None
    if kind == "none" or cap is None:
        pass
    elif window is not None:
        cap_lo, cap_hi = float(min(window)), float(max(window))
    elif kind == "continuous_box" and len(webs_frac) >= 2:
        cap_lo, cap_hi = float(webs_frac[0]), float(webs_frac[-1])
    elif kind == "per_web_pads" and webs_frac:
        pad = 0.04
        cap_lo = max(0.02, float(webs_frac[0]) - pad)
        cap_hi = min(0.98, float(webs_frac[-1]) + pad)
    elif kind == "fixed_band":
        hw = float(getattr(cap, "half_width", None) or 0.12)
        cap_lo, cap_hi = max(0.02, 0.30 - hw), min(0.98, 0.30 + hw)
    # Distinct cap edges only when they are not already web feet (Family D box).
    if cap_lo is not None and cap_hi is not None and cap_hi - cap_lo > 1e-6:
        existing = [p[1] for p in points]
        if min(abs(cap_lo - e) for e in existing) > MIN_GAP:
            points.append(("cap_lo", cap_lo))
        if min(abs(cap_hi - e) for e in existing) > MIN_GAP:
            points.append(("cap_hi", cap_hi))
    points.sort(key=lambda kv: (kv[1], kv[0]))
    merged: list[tuple[str, float]] = []
    for name, frac in points:
        if merged and abs(frac - merged[-1][1]) < MIN_GAP:
            if name in ("le", "te"):
                merged[-1] = (name, frac)
            continue
        merged.append((name, frac))
    anchors = [Anchor(name, frac) for name, frac in merged]
    ordered = [name for name, _ in merged]
    regions = _regions_from_ordered(ordered, cap_lo=cap_lo, cap_hi=cap_hi)
    return tuple(anchors), regions, tuple(webs)


def _regions_from_ordered(ordered: Sequence[str], *, cap_lo: float | None,
                          cap_hi: float | None) -> tuple[Region, ...]:
    out: list[Region] = [Region("te", "te", "te", False, "te", wrap=True)]
    has_cap = cap_lo is not None and cap_hi is not None
    for side in ("suction", "pressure"):
        bay = 0
        cap_i = 0
        tag = side[:2]
        for a, b in zip(ordered[:-1], ordered[1:]):
            is_cap = False
            if a == "le":
                name = f"le_{tag}"
            elif b == "te":
                name = f"te_{tag}"
            elif has_cap and (
                (a.startswith("web_") and b.startswith("web_"))
                or {a, b} == {"cap_lo", "cap_hi"}
                or (a in ("cap_lo", "web_0") and b in ("cap_hi", "web_1"))
            ):
                name = f"cap_{tag}" if cap_i == 0 else f"cap_{tag}_{cap_i}"
                cap_i += 1
                is_cap = True
            else:
                name = f"bay_{tag}_{bay}"
                bay += 1
            out.append(Region(name, a, b, is_cap, side, wrap=False))
    return tuple(out)


def _regularise_anchors(raw: Sequence[tuple[str, float]], *, te_arc: float,
                        min_gap: float) -> tuple[tuple[str, float], ...]:
    """Merge TE-near and closer-than-min_gap anchors. LE/TE always kept."""
    by: dict[str, float] = {}
    for name, frac in raw:
        by[name] = float(frac)
    by.setdefault("le", 0.0)
    by.setdefault("te", 1.0)
    drop = []
    for name, frac in list(by.items()):
        if name in ("le", "te") or str(name).startswith("web"):
            continue
        if min(frac, 1.0 - frac) < te_arc and frac >= 0.5:
            drop.append(name)
        elif frac < te_arc:
            drop.append(name)
    for name in drop:
        by.pop(name, None)
    items = sorted(by.items(), key=lambda kv: (kv[1], kv[0]))
    kept: list[tuple[str, float]] = []
    for name, frac in items:
        if kept and abs(frac - kept[-1][1]) < min_gap:
            prev = kept[-1][0]
            web_new, web_prev = str(name).startswith("web"), str(prev).startswith("web")
            cap_new, cap_prev = str(name).startswith("cap"), str(prev).startswith("cap")
            if name in ("le", "te"):
                kept[-1] = (name, frac if name == "te" else 0.0)
            elif web_new and cap_prev:
                kept[-1] = (name, frac)
            elif cap_new and web_prev:
                continue
            elif web_new and web_prev:
                kept.append((name, frac))
            elif prev in ("le", "te"):
                continue
            else:
                continue
        else:
            kept.append((name, frac))
    names = {n for n, _ in kept}
    if "le" not in names:
        kept.insert(0, ("le", 0.0))
    if "te" not in names:
        kept.append(("te", 1.0))
    return tuple(sorted(kept, key=lambda kv: (kv[1], kv[0])))


def _cluster_web_members(webs: Sequence[WebMember], z: float, *, min_sep: float,
                         t_c: float | None) -> tuple[WebMember, ...]:
    """One web per coincident chord-fraction cluster; omit TE plugs on a cylinder.

    Cylinder stations still drop TE-named plugs. Cluster radius stays ``min_sep``
    so a PA lattice at 0.2 chord (``0.6 - 0.4`` is slightly under 0.20 in float)
    is not merged into the neighbouring spar.
    """
    items = list(webs)
    if not items:
        return ()
    circular = t_c is not None and float(t_c) >= 0.6
    sep = float(min_sep)
    order = sorted(range(len(items)), key=lambda i: items[i].at(z))
    used = [False] * len(items)
    clustered: list[WebMember] = []
    for i in order:
        if used[i]:
            continue
        group = [items[i]]
        used[i] = True
        fi = items[i].at(z)
        for j in order:
            if used[j]:
                continue
            if abs(items[j].at(z) - fi) < sep:
                group.append(items[j])
                used[j] = True
        spars = [w for w in group if "te" not in (w.name or "").lower()]
        clustered.append(spars[0] if spars else group[0])
    spars = [w for w in clustered if "te" not in (w.name or "").lower()]
    max_spar = max((w.at(z) for w in spars), default=0.0)
    out = []
    for w in clustered:
        if "te" in (w.name or "").lower():
            if circular or w.at(z) < max_spar + min_sep:
                continue
        out.append(w)
    return tuple(out)


def _cut_tc(cut) -> float | None:
    skins = getattr(cut, "walls_of", lambda *_: ())("skin")
    if not skins:
        return None
    pts = np.vstack([np.asarray(w.points, float)[:, :2] for w in skins])
    return float(np.ptp(pts[:, 1])) / max(float(cut.chord), 1e-9)


def _match_web_indices(live: Sequence[float], authored: Sequence[float],
                       indices: Sequence[int]) -> tuple[int, ...]:
    """Order-preserving min-cost assignment of live chord-fractions to authored webs.

    Greedy nearest-frac matching can map a surviving mid-chord web onto a dropped
    LE-side PA/LE web. Live and authored are both sorted by fraction; a subsequence
    of authored of length ``len(live)`` minimises the L1 gap.
    """
    if not live or not authored:
        return ()
    order = sorted(range(len(authored)), key=lambda i: float(authored[i]))
    live_s = sorted(float(f) for f in live)
    m, n = len(order), len(live_s)
    if n > m:
        live_s = live_s[:m]
        n = m
    inf = float("inf")
    dp = [[inf] * (n + 1) for _ in range(m + 1)]
    take = [[False] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = 0.0
    for i in range(1, m + 1):
        af = float(authored[order[i - 1]])
        for j in range(1, n + 1):
            dp[i][j] = dp[i - 1][j]
            take[i][j] = False
            cost = dp[i - 1][j - 1] + abs(af - live_s[j - 1])
            if cost < dp[i][j]:
                dp[i][j] = cost
                take[i][j] = True
    chosen: list[int] = []
    i, j = m, n
    while j > 0 and i > 0:
        if take[i][j]:
            chosen.append(int(indices[order[i - 1]]))
            i, j = i - 1, j - 1
        else:
            i -= 1
    return tuple(sorted(chosen))


def role_atom(name: str) -> tuple:
    """Paint-stable role for compatibility (cap flicker shares a cap atom)."""
    if name == "te" or name.startswith("te_"):
        return ("TE",)
    if name.startswith("le_"):
        return ("LE", name[-2:])
    if name.startswith("cap_"):
        return ("cap", name[4:6])
    if name.startswith("bay_"):
        parts = name.split("_")
        side = parts[1] if len(parts) > 1 else ""
        return ("bay", side)
    if name.startswith("web"):
        return ("web", name)
    return ("panel", name)


def signatures_compatible(a, b) -> bool:
    """True if role atoms match after TE wrap / cap-paint flicker.

    Web add/drop is **not** compatible (new zone) even if loop atoms match.
    """
    if a == b:
        return True
    if not isinstance(a, tuple) or not isinstance(b, tuple):
        return False
    if len(a) == 2 and len(b) == 2 and isinstance(a[0], tuple):
        ra, wa = a
        rb, wb = b
        if wa != wb:
            return False
        return _loop_atoms(ra) == _loop_atoms(rb)
    return _loop_atoms(a) == _loop_atoms(b)


def _loop_atoms(names: Sequence[str]) -> frozenset:
    atoms = []
    for n in names:
        atom = role_atom(n)
        if atom[0] == "bay":
            atoms.append(("panel", atom[1]))
        elif atom[0] == "cap":
            atoms.append(("panel", atom[1]))
        else:
            atoms.append(atom)
    return frozenset(atoms)


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
        cap_lo, cap_hi = fracs.get("cap_lo"), fracs.get("cap_hi")
        regions = _regions_from_ordered(ordered, cap_lo=cap_lo, cap_hi=cap_hi)
        vertices: list[ShellVertexId] = [
            ("TE",), ("TE", "suction"), ("TE", "pressure"), ("anchor", "le"),
        ]
        loop = []
        for region in regions:
            if region.wrap:
                start, end = ("TE", "suction"), ("TE", "pressure")
            else:
                start = _loop_vid(region.start, region.side)
                end = _loop_vid(region.end, region.side)
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


def _chord_fracs(cut, xy: np.ndarray) -> np.ndarray:
    t = np.radians(float(cut.twist))
    c, s = np.cos(t), np.sin(t)
    x = xy[:, 0] - float(cut.offset[0])
    y = xy[:, 1] - float(cut.offset[1])
    local_x = c * x - s * y
    return (local_x + float(cut.pitch_axis) * float(cut.chord)) / max(float(cut.chord), 1e-12)


def _as_xyz(pts: np.ndarray, z: float) -> np.ndarray:
    pts = np.asarray(pts, float)
    if pts.ndim != 2 or len(pts) == 0:
        return np.zeros((0, 3))
    if pts.shape[1] == 2:
        return np.column_stack([pts, np.full(len(pts), float(z))])
    return pts[:, :3]


def _skin_mid_loop(graph: ShellGraph, z: float) -> np.ndarray:
    """Closed skin-mid polyline chained by shared endpoints (closing vertex omitted)."""
    chain = graph.closed_loop_chain(rel_tol=1e-9) if hasattr(graph, "closed_loop_chain") else None
    if chain:
        by_name = {r.name: r for r in graph.loop}
        chunks: list[np.ndarray] = []
        for name, _xy, rev in chain:
            run = by_name.get(name)
            if run is None:
                continue
            pts = _as_xyz(np.asarray(run.points, float), z)
            if rev:
                pts = pts[::-1]
            chunks.append(pts)
        if chunks:
            pts = chunks[0]
            for nxt in chunks[1:]:
                if np.allclose(pts[-1], nxt[0]):
                    pts = np.vstack([pts, nxt[1:]])
                else:
                    pts = np.vstack([pts, nxt])
            if len(pts) > 2 and np.allclose(pts[0], pts[-1]):
                pts = pts[:-1]
            return pts
    try:
        xy = np.asarray(graph.loop_polygon(), float)
    except (ValueError, AttributeError):
        xy = np.zeros((0, 2))
    if len(xy) == 0:
        return np.zeros((0, 3))
    if len(xy) > 2 and np.allclose(xy[0], xy[-1]):
        xy = xy[:-1]
    return np.column_stack([xy[:, :2], np.full(len(xy), float(z))])


def _te_frac_drop(cut) -> float:
    """Chord-fraction length of the TE wrap: full ``TE_ARC`` on a cylinder, a
    fraction of that on a thin foil, so the cap width follows t/c instead of
    switching on a y-span threshold (that switch wrinkled the loft)."""
    tc = _cut_tc(cut)
    if tc is None:
        tc = 0.2
    return float(TE_ARC) * float(np.clip(float(tc) / 0.45, 0.25, 1.0))


def _te_corner_indices(pts: np.ndarray, cut) -> tuple[int, int, int]:
    """LE and the OML-TE ridge index (both corners). Side arcs start at the ridge;
    wrap corners are interpolated along those arcs at ``_te_frac_drop``.
    """
    frac = _chord_fracs(cut, pts)
    i_le = int(np.argmin(frac))
    oml = _oml_te_xyz(cut)
    if oml is not None:
        i_te = int(np.argmin(np.linalg.norm(pts[:, :2] - np.asarray(oml, float)[:2], axis=1)))
        if i_te != i_le:
            return i_le, i_te, i_te
    i_te = int(np.argmax(frac))
    return i_le, i_te, i_te


def _te_wrap_arc(pts: np.ndarray, a: np.ndarray, b: np.ndarray, i_le: int,
                 i_te: int | None = None) -> np.ndarray | None:
    """Skin-mid walk from TE suction to TE pressure that does not pass the LE.

    The walk must include the OML TE ridge when ``i_te`` is given so the loft
    keys sit on the trailing edge rather than a 0.65-chord shortcut.
    """
    if len(pts) < 2:
        return None
    ia = int(np.argmin(np.linalg.norm(pts - np.asarray(a, float), axis=1)))
    ib = int(np.argmin(np.linalg.norm(pts - np.asarray(b, float), axis=1)))
    le = pts[int(i_le)]
    ridge_i = None if i_te is None else int(i_te)

    def has_le(arc: np.ndarray) -> bool:
        return any(np.allclose(p, le) for p in arc)

    def has_ridge(arc: np.ndarray) -> bool:
        if ridge_i is None or ridge_i < 0 or ridge_i >= len(pts):
            return True
        ridge = pts[ridge_i]
        return any(np.allclose(p, ridge) for p in arc)

    def length(arc: np.ndarray) -> float:
        return float(np.sum(np.linalg.norm(np.diff(arc, axis=0), axis=1)))

    if ia == ib:
        if ridge_i is None or ridge_i == ia:
            return None
        via = np.vstack((pts[ia], pts[ridge_i], pts[ib]))
        return via if not has_le(via) else None

    cands = [_walk_arc_dir(pts, ia, ib, 1), _walk_arc_dir(pts, ia, ib, -1)]
    with_ridge = [arc for arc in cands if len(arc) > 1 and not has_le(arc) and has_ridge(arc)]
    if with_ridge:
        return min(with_ridge, key=length)
    if ridge_i is not None and ridge_i != ia and ridge_i != ib:
        via = np.vstack((_walk_arc(pts, ia, ridge_i), _walk_arc(pts, ridge_i, ib)[1:]))
        if len(via) > 1 and not has_le(via):
            return via
    ok = [arc for arc in cands if len(arc) > 1 and not has_le(arc)]
    if not ok:
        return None
    return min(ok, key=length)


def _orient_loop_te(pts: np.ndarray, cut) -> np.ndarray:
    if len(pts) == 0:
        return pts
    _i_le, i_su, _i_pr = _te_corner_indices(pts, cut)
    mid = 0.5 * (pts[i_su] + pts[_i_pr])
    i = int(np.argmin(np.linalg.norm(pts - mid, axis=1)))
    return np.roll(pts, -i, axis=0)


def _walk_arc_dir(pts: np.ndarray, i0: int, i1: int, step: int) -> np.ndarray:
    n = len(pts)
    if n == 0:
        return pts
    out = [pts[i0]]
    if i0 == i1:
        return np.asarray(out, float)
    i = i0
    for _ in range(n):
        i = (i + int(step)) % n
        out.append(pts[i])
        if i == i1:
            break
    return np.asarray(out, float)


def _walk_arc(pts: np.ndarray, i0: int, i1: int) -> np.ndarray:
    n = len(pts)
    if n == 0:
        return pts
    if i0 == i1:
        return pts[i0 : i0 + 1]
    d1, d2 = (i1 - i0) % n, (i0 - i1) % n
    step = 1 if d1 <= d2 else -1
    return _walk_arc_dir(pts, i0, i1, step)


def _xyz_at_frac(arc: np.ndarray, fracs: np.ndarray, f: float) -> np.ndarray:
    """Point on an open arc whose chord fraction is ``f`` (linear between vertices)."""
    if len(arc) == 0:
        return np.zeros(3)
    f = float(np.clip(f, 0.0, 1.0))
    fracs = np.asarray(fracs, float)
    for i in range(len(arc) - 1):
        a, b = float(fracs[i]), float(fracs[i + 1])
        if (a - f) * (b - f) <= 0.0 and abs(b - a) > 1e-15:
            u = (f - a) / (b - a)
            return arc[i] + u * (arc[i + 1] - arc[i])
    return arc[int(np.argmin(np.abs(fracs - f)))].copy()


def _side_arcs(pts: np.ndarray, cut) -> tuple[np.ndarray, np.ndarray]:
    """Suction and pressure arcs, each TE-corner → LE (open; TE corners may differ)."""
    i_le, i_su, i_pr = _te_corner_indices(pts, cut)
    n = len(pts)
    if i_su == i_pr:
        cw = _walk_arc_dir(pts, i_su, i_le, 1)
        ccw = _walk_arc_dir(pts, i_su, i_le, -1)
        y0 = float(cut.offset[1])
        interior = cw[1:-1] if len(cw) > 2 else cw
        if float(np.mean(interior[:, 1])) >= y0:
            return cw, ccw
        return ccw, cw

    def _away(i_from: int, i_other: int) -> int:
        if (i_from + 1) % n == i_other:
            return -1
        if (i_from - 1) % n == i_other:
            return 1
        d_le = (i_le - i_from) % n
        return 1 if d_le <= n - d_le else -1

    return (
        _walk_arc_dir(pts, i_su, i_le, _away(i_su, i_pr)),
        _walk_arc_dir(pts, i_pr, i_le, _away(i_pr, i_su)),
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


def _oml_te_xyz(cut) -> np.ndarray | None:
    """Trailing-edge point of the skin walls (nearest stand-in for OML TE)."""
    key_pts = getattr(cut, "key_points", None)
    if callable(key_pts):
        te = key_pts().get("TE")
        if te is not None:
            return np.asarray(te, float)[:3]
    skins = getattr(cut, "walls_of", lambda *_: ())("skin")
    if not skins:
        return None
    oml = np.vstack([np.asarray(w.points, float) for w in skins])
    fr = _chord_fracs(cut, oml)
    return oml[int(np.argmax(fr))].copy()


def _resample_keys(p0: np.ndarray, p1: np.ndarray, arc: np.ndarray, n: int,
                   *, bunch_end: bool = False) -> np.ndarray:
    """Arc-length resample between two keys that already lie on ``arc`` (TE→LE)."""
    n = max(int(n), 2)
    if len(arc) < 2:
        out = np.repeat(np.asarray(p0, float).reshape(1, -1), n, axis=0)
        out[-1] = p1
        return out
    i0 = int(np.argmin(np.linalg.norm(arc - p0, axis=1)))
    i1 = int(np.argmin(np.linalg.norm(arc - p1, axis=1)))
    if i0 == i1:
        sub = np.vstack([p0, p1])
    elif i0 < i1:
        sub = np.vstack([p0, arc[i0 + 1 : i1], p1])
    else:
        sub = np.vstack([p0, arc[i0 - 1 : i1 : -1], p1])
    if bunch_end and n > 2:
        sampled = _arc_resample_bunched(sub, n)
    else:
        sampled = arc_resample_open(sub, n)
    sampled[0] = np.asarray(p0, float)
    sampled[-1] = np.asarray(p1, float)
    return sampled


def _arc_resample_bunched(pts: np.ndarray, n: int) -> np.ndarray:
    """Arc-length resample clustered toward the last vertex (TE)."""
    pts = np.asarray(pts, float)
    n = max(int(n), 2)
    if len(pts) < 2:
        return arc_resample_open(pts, n)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 1e-15:
        return arc_resample_open(pts, n)
    u = np.sin(np.linspace(0.0, 0.5 * np.pi, n)) ** 2
    t = u * s[-1]
    out = np.empty((n, pts.shape[1]))
    for j, tj in enumerate(t):
        k = int(np.searchsorted(s, tj, side="right") - 1)
        k = min(max(k, 0), len(seg) - 1)
        span = seg[k] if seg[k] > 1e-15 else 1.0
        w = (tj - s[k]) / span
        out[j] = pts[k] + w * (pts[k + 1] - pts[k])
    return out


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
        chord = np.vstack([su_ft, pr_ft])
        if match is None:
            sampled = arc_resample_open(chord, n)
        else:
            mid = _as_xyz(np.asarray(match.points, float), z)
            d0 = min(np.linalg.norm(mid[0] - su_ft), np.linalg.norm(mid[-1] - su_ft))
            d1 = min(np.linalg.norm(mid[0] - pr_ft), np.linalg.norm(mid[-1] - pr_ft))
            span = max(float(np.linalg.norm(su_ft - pr_ft)), 1e-12)
            if d0 > 0.05 * span or d1 > 0.05 * span:
                sampled = arc_resample_open(chord, n)
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


def loft_mapped_strips(s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
                       t0: ShellTopologyTemplate, t1: ShellTopologyTemplate
                       ) -> list[tuple[str, np.ndarray, bool]]:
    """Quad (or collapsing) strips between two templates. Omit one-sided webs.

    Each t1 skin branch is ruled from the t0 skin walk between its mapped
    endpoints (web onset/drop splits a panel along the mould, not a chord).

    Returns ``(name, face, collapsing)`` with ``face`` shape ``(4, 3)`` or ``(3, 3)``.
    """
    faces: list[tuple[str, np.ndarray, bool]] = []
    used: set[str] = set()
    for branch in t1.loop:
        if branch.wrap or branch.name not in s1:
            continue
        pb = np.asarray(s1[branch.name], float)
        pa = _skin_parent(branch, s0, t0, t1, pb)
        if pa is None:
            continue
        birth = bool(np.allclose(pa[0], pa[-1]))
        faces.extend(_strip_faces(branch.name, pa, pb, t0, t1, branch, birth=birth))
        used.add(branch.name)
    for branch in t0.loop:
        if not branch.wrap or branch.name not in s0:
            continue
        dest = _strip_dest(branch, s1, t0, t1)
        if dest is None:
            continue
        pa, pb = np.asarray(s0[branch.name], float), np.asarray(s1[dest], float)
        faces.extend(_strip_faces(dest, pa, pb, t0, t1, branch))
        used.add(dest)
    names1 = {w.name for w in t1.webs}
    for web in t0.webs:
        if web.name not in s0 or web.name not in s1 or web.name not in names1:
            continue
        pa, pb = np.asarray(s0[web.name], float), np.asarray(s1[web.name], float)
        n = min(len(pa), len(pb))
        if n < 2:
            continue
        for i in range(n - 1):
            quad = np.stack([pa[i], pa[i + 1], pb[i + 1], pb[i]])
            if _quad_bowtie(quad):
                a, b, c, d = quad
                faces.append((web.name, np.stack([a, b, c]), True))
                faces.append((web.name, np.stack([a, c, d]), True))
                continue
            faces.append((web.name, quad, False))
    return faces


def _strip_faces(name: str, pa: np.ndarray, pb: np.ndarray, t0, t1, branch,
                 *, birth: bool = False) -> list[tuple[str, np.ndarray, bool]]:
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    if len(pa) != len(pb) or len(pa) < 2:
        n = min(len(pa), len(pb))
        if n < 2:
            return []
        pa, pb = pa[:n], pb[:n]
    collapsing = birth or np.allclose(pa[0], pa[-1]) or np.allclose(pb[0], pb[-1])
    if not birth:
        try:
            start_m = t0.map_vertex(branch.start, t1) if branch.start in t0.vertices else branch.start
            end_m = t0.map_vertex(branch.end, t1) if branch.end in t0.vertices else branch.end
            collapsing = collapsing or (start_m == end_m and not branch.wrap)
        except KeyError:
            pass
    if branch.wrap and np.allclose(pa[0], pa[-1]) and np.allclose(pb[0], pb[-1]):
        return []
    out = []
    if collapsing:
        for i in range(len(pa) - 1):
            out.append((name, np.stack([pa[i], pa[i + 1], pb[min(i + 1, len(pb) - 1)]]), True))
        return out
    for i in range(len(pa) - 1):
        quad = np.stack([pa[i], pa[i + 1], pb[i + 1], pb[i]])
        if _quad_bowtie(quad):
            a, b, c, d = quad
            out.append((name, np.stack([a, b, c]), True))
            out.append((name, np.stack([a, c, d]), True))
        else:
            out.append((name, quad, False))
    return out


def _xyz_of_vid(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate, vid) -> np.ndarray | None:
    for b in tmpl.loop:
        if b.name not in samples:
            continue
        pts = np.asarray(samples[b.name], float)
        if len(pts) == 0:
            continue
        if b.start == vid:
            return pts[0].copy()
        if b.end == vid:
            return pts[-1].copy()
    return None


def _concat_side_samples(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate,
                         side: str) -> np.ndarray | None:
    """LE→TE skin polyline on one side from sampled loop branches."""
    parts: list[np.ndarray] = []
    for b in tmpl.loop:
        if b.wrap or b.side != side or b.name not in samples:
            continue
        _append_ring(parts, np.asarray(samples[b.name], float))
    if not parts:
        return None
    return np.vstack(parts)


def _on_chain(chain: np.ndarray, p: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(np.asarray(chain, float)[:, :2] - np.asarray(p, float)[:2], axis=1)
    return np.asarray(chain, float)[int(np.argmin(d))].copy()


def _skin_parent(branch: LoopBranch, s0: Mapping[str, np.ndarray],
                 t0: ShellTopologyTemplate, t1: ShellTopologyTemplate, pb: np.ndarray
                 ) -> np.ndarray | None:
    """t0 skin arc between the mapped endpoints of a t1 branch (not a chord).

    A web that only exists on t1 is placed on the t0 mould by nearest XY, so the
    TE panel splits at the new foot instead of collapsing to the TE point.
    """
    pb = np.asarray(pb, float)
    n = max(len(pb), 2)
    chain = _concat_side_samples(s0, t0, branch.side)
    p0 = _xyz_of_vid(s0, t0, branch.start) if branch.start in t0.vertices else None
    p1 = _xyz_of_vid(s0, t0, branch.end) if branch.end in t0.vertices else None
    if chain is None:
        return _collapsed_parent(branch, s0, t0, t1, n)
    if p0 is None:
        p0 = _on_chain(chain, pb[0])
    if p1 is None:
        p1 = _on_chain(chain, pb[-1])
    return _resample_keys(p0, p1, chain, n)


def _collapsed_parent(branch: LoopBranch, s0: Mapping[str, np.ndarray],
                      t0: ShellTopologyTemplate, t1: ShellTopologyTemplate, n: int) -> np.ndarray | None:
    try:
        start0 = t1.map_vertex(branch.start, t0)
        end0 = t1.map_vertex(branch.end, t0)
    except KeyError:
        return None
    p0 = _xyz_of_vid(s0, t0, start0)
    p1 = _xyz_of_vid(s0, t0, end0)
    if p0 is None:
        p0 = p1
    if p1 is None:
        p1 = p0
    if p0 is None:
        return None
    n = max(int(n), 2)
    if np.allclose(p0, p1):
        return np.repeat(np.asarray(p0, float).reshape(1, -1), n, axis=0)
    return np.linspace(p0, p1, n)


def _same_web_set(t0: ShellTopologyTemplate, t1: ShellTopologyTemplate) -> bool:
    s0, s1 = t0.signature, t1.signature
    if isinstance(s0, tuple) and isinstance(s1, tuple) and len(s0) == 2 and len(s1) == 2:
        return s0[1] == s1[1]
    return s0 == s1


def _strip_dest(branch: LoopBranch, s1: Mapping[str, np.ndarray],
                t0: ShellTopologyTemplate, t1: ShellTopologyTemplate) -> str | None:
    if _same_web_set(t0, t1) and branch.name in s1:
        return branch.name
    try:
        start_m = t0.map_vertex(branch.start, t1)
        end_m = t0.map_vertex(branch.end, t1)
    except KeyError:
        return branch.name if branch.name in s1 else None
    if start_m == end_m and not branch.wrap:
        return None
    for other in t1.loop:
        if other.name not in s1 or other.side != branch.side:
            continue
        if other.start == start_m and other.end == end_m:
            return other.name
    return None


def _quad_bowtie(quad: np.ndarray) -> bool:
    """True if the two triangles of a loft quad have opposing normals (a fold)."""
    a, b, c, d = np.asarray(quad, float)
    n1 = np.cross(b - a, c - a)
    n2 = np.cross(c - a, d - a)
    if np.linalg.norm(n1) < 1e-18 or np.linalg.norm(n2) < 1e-18:
        return False
    return float(np.dot(n1, n2)) < -0.2 * float(np.linalg.norm(n1) * np.linalg.norm(n2))


def _append_ring(parts: list[np.ndarray], chain: np.ndarray | None) -> None:
    if chain is None:
        return
    c = np.asarray(chain, float)
    if len(c) == 0:
        return
    if parts and np.allclose(parts[-1][-1], c[0], atol=1e-12, rtol=0.0):
        c = c[1:]
    if len(c):
        parts.append(c)


def _web_pts_from_su(samples: Mapping[str, np.ndarray], web: WebBranch | None,
                     su_end: np.ndarray) -> np.ndarray | None:
    if web is None or web.name not in samples:
        return None
    pts = np.asarray(samples[web.name], float)
    if len(pts) < 2:
        return None
    if np.linalg.norm(pts[-1] - su_end) < np.linalg.norm(pts[0] - su_end):
        pts = pts[::-1].copy()
    return pts


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


__all__ = [
    "Anchor",
    "LoopBranch",
    "MIN_GAP",
    "Region",
    "ShellTopologyTemplate",
    "ShellVertexId",
    "TE_ARC",
    "WallRegionGraph",
    "WebBranch",
    "WebMember",
    "ZonedShellTopology",
    "loft_mapped_strips",
    "place_station_keys",
    "rib_cell_polygons_from_samples",
    "rib_face_polygons",
    "rib_face_triangles",
    "role_atom",
    "sample_region_polylines",
    "section_cut_near",
    "signatures_compatible",
    "templates_along",
    "zone_region_graphs",
]
