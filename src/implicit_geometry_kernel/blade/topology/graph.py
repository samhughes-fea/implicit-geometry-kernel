"""Authored wall-region graph (anchors, named loop regions, webs)."""

from __future__ import annotations

import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementKind, CapPlacementRule

if TYPE_CHECKING:
    from implicit_geometry_kernel.catalogue.system_layout import SystemLayoutSpec
    from implicit_geometry_kernel.section.layout.section_layout import SectionLayout

TE_ARC = 0.03
MIN_GAP = 0.002
PER_WEB_PAD = 0.04

_CAP_MODE_LEGACY: dict[CapPlacementKind, str] = {
    "none": "none",
    "fixed_band": "band",
    "per_web_pads": "pads",
    "continuous_box": "box",
}

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

        def _table(z: float, *, _z=zs, _v=vs) -> float:
            return _interp_table(z, _z, _v)

        return _table
    c = float(value)

    def _const(z: float, *, _c: float = c) -> float:
        return _c

    return _const


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
    """A shear web: published index, chord-fraction law, optional span window.

    Default ``name`` is ``web_{index}``, matching catalogue / station labels.
    """

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
    cap_windows: tuple[tuple[float, float], ...] = ()
    cap_placement: CapPlacementRule = field(default_factory=CapPlacementRule.none)
    cap_pad: float | None = None

    @property
    def cap_mode(self) -> str:
        """Deprecated alias of :attr:`cap_placement` kind (``band``/``pads``/``box``)."""
        warnings.warn(
            "WallRegionGraph.cap_mode is deprecated; use cap_placement.kind "
            "(CapPlacementKind). Removed in the next minor.",
            DeprecationWarning,
            stacklevel=2,
        )
        return _CAP_MODE_LEGACY[self.cap_placement.kind]

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
    def from_section_layout(
        cls,
        layout: SectionLayout | str | SystemLayoutSpec,
        *,
        te_arc: float = TE_ARC,
        min_gap: float = MIN_GAP,
    ) -> WallRegionGraph:
        from implicit_geometry_kernel.blade.topology.layout_compile import (
            graph_from_section_layout,
        )

        return graph_from_section_layout(cls, layout, te_arc=te_arc, min_gap=min_gap)

    @classmethod
    def from_windio_structure(
        cls, structure: Mapping, z_of_s=None, *,
        te_arc: float = TE_ARC,
        min_gap: float = MIN_GAP,
    ) -> WallRegionGraph:
        from implicit_geometry_kernel.blade.topology.layout_compile import (
            graph_from_windio_structure,
        )

        return graph_from_windio_structure(
            cls, structure, z_of_s, te_arc=te_arc, min_gap=min_gap
        )


def _indexed_cap_name(name: str) -> bool:
    parts = str(name).split("_")
    return len(parts) >= 3 and parts[0] == "cap" and parts[1].isdigit() and parts[2] in ("lo", "hi")


def live_cap_windows(
    graph: WallRegionGraph,
    fracs: Mapping[str, float],
    ordered: Sequence[str],
    *,
    z: float,
    active_webs: Sequence[WebMember],
) -> tuple[tuple[float, float], ...]:
    """Cap chord windows at ``z``: per-web pads, live box, or authored band."""
    from implicit_geometry_kernel.blade.topology.layout_compile import _merge_chord_windows

    kind = graph.cap_placement.kind
    if kind == "per_web_pads":
        pad = float(graph.cap_pad if graph.cap_pad is not None else PER_WEB_PAD)
        return _merge_chord_windows(
            [(max(0.02, float(w.at(z)) - pad), min(0.98, float(w.at(z)) + pad))
             for w in active_webs]
        )
    if kind == "continuous_box":
        webs = [float(fracs[n]) for n in ordered if str(n).startswith("web_") and n in fracs]
        if len(webs) >= 2:
            return ((min(webs), max(webs)),)
        return ()
    live = cap_windows_from_fracs(fracs)
    return live if live else tuple(graph.cap_windows or ())


def cap_windows_from_fracs(fracs: Mapping[str, float]) -> tuple[tuple[float, float], ...]:
    """Live ``cap_lo``/``cap_hi`` or indexed ``cap_{i}_lo``/``cap_{i}_hi`` pairs."""
    if "cap_lo" in fracs and "cap_hi" in fracs:
        return ((float(fracs["cap_lo"]), float(fracs["cap_hi"])),)
    out: list[tuple[float, float]] = []
    i = 0
    while f"cap_{i}_lo" in fracs or f"cap_{i}_hi" in fracs:
        lo, hi = fracs.get(f"cap_{i}_lo"), fracs.get(f"cap_{i}_hi")
        if lo is not None and hi is not None:
            out.append((float(lo), float(hi)))
        i += 1
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
