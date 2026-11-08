"""Authored wall-region graph (anchors, named loop regions, webs)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np

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
        cap_rule = getattr(layout, "cap_placement", None)
        anchors, regions, webs = _layout_parts(webs_frac, cap_rule)
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
                win = windows_span[j] if j < len(windows_span) else None
                if win is not None:
                    z0, z1 = win
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
            lo = float(np.median(windows[1][:, 0]))
            hi = float(np.median(windows[1][:, 1]))
            dummy_cap = _CapWindow((lo, hi))
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

