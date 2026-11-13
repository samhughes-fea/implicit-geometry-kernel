"""Compile a :class:`SectionLayout` into authored loft anchors and regions."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import TYPE_CHECKING

import numpy as np

from implicit_geometry_kernel.catalogue.system_layout import SystemLayoutSpec
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout

from implicit_geometry_kernel.blade.topology.graph import (
    MIN_GAP,
    PER_WEB_PAD,
    TE_ARC,
    Anchor,
    Region,
    WebMember,
)

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.topology.graph import WallRegionGraph


def coerce_section_layout(layout: SectionLayout | str | SystemLayoutSpec) -> SectionLayout:
    """Normalise catalogue keys and specs to :class:`SectionLayout`."""
    if isinstance(layout, str):
        return SectionLayout.from_catalogue_key(layout)
    if isinstance(layout, SectionLayout):
        return layout
    if isinstance(layout, SystemLayoutSpec):
        return SectionLayout.from_system_layout_spec(layout)
    raise TypeError(
        "from_section_layout expects SectionLayout, catalogue key, or "
        f"SystemLayoutSpec; got {type(layout)!r}"
    )


def _merge_chord_windows(
    windows: Sequence[tuple[float, float]], *, gap: float = MIN_GAP,
) -> tuple[tuple[float, float], ...]:
    """Union overlapping or min-gap-adjacent chord windows, LE to TE."""
    items = sorted(
        (float(min(a, b)), float(max(a, b)))
        for a, b in windows
        if abs(float(b) - float(a)) > 1e-6
    )
    out: list[tuple[float, float]] = []
    for lo, hi in items:
        if out and lo <= out[-1][1] + gap:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return tuple(out)


def _append_cap_edge(points: list[tuple[str, float]], name: str, frac: float) -> None:
    existing = [p[1] for p in points]
    if min(abs(float(frac) - e) for e in existing) > MIN_GAP:
        points.append((name, float(frac)))


def _segment_in_windows(
    a: str, b: str, fracs: Mapping[str, float],
    windows: Sequence[tuple[float, float]],
) -> bool:
    if not windows or a not in fracs or b not in fracs:
        return False
    fa, fb = float(fracs[a]), float(fracs[b])
    lo, hi = (fa, fb) if fa <= fb else (fb, fa)
    for wlo, whi in windows:
        if lo >= float(wlo) - 1e-9 and hi <= float(whi) + 1e-9:
            return True
    return False


def _regions_from_ordered(
    ordered: Sequence[str],
    *,
    cap_lo: float | None = None,
    cap_hi: float | None = None,
    fracs: Mapping[str, float] | None = None,
    windows: Sequence[tuple[float, float]] | None = None,
) -> tuple[Region, ...]:
    out: list[Region] = [Region("te", "te", "te", False, "te", wrap=True)]
    if windows is None:
        windows = () if cap_lo is None or cap_hi is None else ((float(cap_lo), float(cap_hi)),)
    frac_map: dict[str, float] = dict(fracs or {})
    if cap_lo is not None:
        frac_map.setdefault("cap_lo", float(cap_lo))
    if cap_hi is not None:
        frac_map.setdefault("cap_hi", float(cap_hi))
    has_cap = bool(windows)
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
                _segment_in_windows(a, b, frac_map, windows)
                or (
                    not frac_map
                    and (
                        (a.startswith("web_") and b.startswith("web_"))
                        or {a, b} == {"cap_lo", "cap_hi"}
                        or (a in ("cap_lo", "web_0") and b in ("cap_hi", "web_1"))
                    )
                )
            ):
                name = f"cap_{tag}" if cap_i == 0 else f"cap_{tag}_{cap_i}"
                cap_i += 1
                is_cap = True
            else:
                name = f"bay_{tag}_{bay}"
                bay += 1
            out.append(Region(name, a, b, is_cap, side, wrap=False))
    return tuple(out)


def layout_parts(
    web_fracs: Sequence[float], cap: CapPlacementRule | None
) -> tuple[
    tuple[Anchor, ...], tuple[Region, ...], tuple[WebMember, ...],
    tuple[tuple[float, float], ...], CapPlacementRule, float | None,
]:
    """Anchors, regions, webs, windows, resolved cap rule, and pad."""
    rule = CapPlacementRule.none() if cap is None else cap
    webs_frac = tuple(sorted(float(f) for f in web_fracs if 0.0 < float(f) < 1.0))
    points: list[tuple[str, float]] = [("le", 0.0), ("te", 1.0)]
    webs = []
    for i, f in enumerate(webs_frac):
        points.append((f"web_{i}", f))
        webs.append(WebMember(i, f, name=f"web_{i}"))
    kind = rule.kind
    window = rule.chord_window
    pad: float | None = None
    windows: tuple[tuple[float, float], ...] = ()
    if kind == "none":
        pass
    elif kind == "per_web_pads" and webs_frac:
        pad = float(rule.half_width if rule.half_width is not None else PER_WEB_PAD)
        windows = _merge_chord_windows(
            [(max(0.02, float(f) - pad), min(0.98, float(f) + pad)) for f in webs_frac]
        )
        if rule.half_width is None:
            rule = replace(rule, half_width=pad)
    elif window is not None:
        windows = ((float(min(window)), float(max(window))),)
    elif kind == "continuous_box" and len(webs_frac) >= 2:
        windows = ((float(webs_frac[0]), float(webs_frac[-1])),)
    elif kind == "fixed_band":
        hw = float(rule.half_width) if rule.half_width is not None else 0.12
        windows = ((max(0.02, 0.30 - hw), min(0.98, 0.30 + hw)),)
    for i, (lo, hi) in enumerate(windows):
        if hi - lo <= 1e-6:
            continue
        if len(windows) == 1:
            lo_name, hi_name = "cap_lo", "cap_hi"
        else:
            lo_name, hi_name = f"cap_{i}_lo", f"cap_{i}_hi"
        _append_cap_edge(points, lo_name, lo)
        _append_cap_edge(points, hi_name, hi)
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
    fracs = {name: frac for name, frac in merged}
    regions = _regions_from_ordered(ordered, fracs=fracs, windows=windows)
    return tuple(anchors), regions, tuple(webs), windows, rule, pad


def graph_from_section_layout(
    cls: type[WallRegionGraph],
    layout: SectionLayout | str | SystemLayoutSpec,
    *,
    te_arc: float = TE_ARC,
    min_gap: float = MIN_GAP,
) -> WallRegionGraph:
    sl = coerce_section_layout(layout)
    webs_frac = tuple(float(f) for f in sl.web_chord_fracs)
    anchors, regions, webs, windows, rule, pad = layout_parts(webs_frac, sl.cap_placement)
    return cls(
        anchors, regions, webs, te_arc=float(te_arc), min_gap=float(min_gap),
        cap_windows=windows, cap_placement=rule, cap_pad=pad,
    )


def graph_from_windio_structure(
    cls: type[WallRegionGraph],
    structure: Mapping,
    z_of_s=None,
    *,
    te_arc: float = TE_ARC,
    min_gap: float = MIN_GAP,
    geo=None,
) -> WallRegionGraph:
    if geo is not None:
        return _graph_from_windio_on_sections(cls, structure, z_of_s, geo, te_arc=te_arc,
                                              min_gap=min_gap)
    from implicit_geometry_kernel.blade.windio import (
        cap_window_table_from_windio_structure,
        web_position_table_from_windio_structure,
        web_span_windows_from_windio_structure,
    )

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
        dummy_cap = CapPlacementRule.from_chord_windows(((lo, hi),))
        _, regions, _, cap_windows, rule, pad = layout_parts(web_fracs, dummy_cap)
    else:
        _, regions, _, cap_windows, rule, pad = layout_parts(
            tuple(w.at(0.0) for w in webs) if webs else (), None)
    return cls(
        tuple(anchors), regions, tuple(webs), te_arc=float(te_arc),
        min_gap=float(min_gap), cap_windows=cap_windows, cap_placement=rule, cap_pad=pad,
    )


def _graph_from_windio_on_sections(cls: type[WallRegionGraph], structure: Mapping, z_of_s, geo,
                                   *, te_arc: float, min_gap: float) -> WallRegionGraph:
    """As :func:`graph_from_windio_structure`, with the file's arcs placed on ``geo``'s
    sections (:mod:`implicit_geometry_kernel.blade.windio_arc`): each web anchor is the
    mid-point of its two feet, the cap anchors the union of the two sides' windows."""
    from implicit_geometry_kernel.blade.windio_arc import (
        cap_side_windows_fn_from_windio_structure,
        cap_union_table_from_windio_structure,
        web_feet_series_from_windio_structure,
        web_midpoint_laws_from_windio_structure,
    )

    anchors: list[Anchor] = [Anchor("le", 0.0), Anchor("te", 1.0)]
    webs: list[WebMember] = []
    # Laws follow the arc placed on the section at z. A chord-fraction table of the
    # file's grid freezes at the last sample, and the cap edge then passes that
    # frozen web near the tip (IEA 15, about 113.3 m; BG-037).
    laws = web_midpoint_laws_from_windio_structure(structure, geo, z_of_s) or []
    series = web_feet_series_from_windio_structure(structure, geo, z_of_s) or []
    by_name = {name: (zs, feet) for name, zs, feet, _window in series}
    order = sorted(range(len(laws)), key=lambda j: float(np.median(by_name[laws[j][0]][1].sum(axis=1)))
                   if laws[j][0] in by_name and len(by_name[laws[j][0]][1]) else j)
    for j, k in enumerate(order):
        name, law, window = laws[k]
        z0, z1 = (None, None) if window is None else window
        webs.append(WebMember(int(j), law, z0, z1, name=name))
        anchors.append(Anchor(f"web_{j}", law, z0, z1))
    web_fracs = tuple(
        float(np.median(0.5 * by_name[laws[k][0]][1].sum(axis=1)))
        for k in order if laws[k][0] in by_name
    )
    windows = cap_union_table_from_windio_structure(structure, geo, z_of_s)
    side_windows = cap_side_windows_fn_from_windio_structure(structure, geo, z_of_s)
    cap = None
    if windows is not None and side_windows is not None:
        zg, vals = windows

        def _union(z: float, end: int, _fn=side_windows, _zg=zg, _vals=vals) -> float:
            live = _fn(float(z))
            if live and set(live) == {"suction", "pressure"}:
                lo = min(live["suction"][0], live["pressure"][0])
                hi = max(live["suction"][1], live["pressure"][1])
                return float(lo if end == 0 else hi)
            return float(np.interp(float(z), _zg, _vals[:, end]))

        z0, z1 = float(zg[0]), float(zg[-1])
        anchors.append(Anchor("cap_lo", lambda z, _u=_union: _u(z, 0), z0, z1))
        anchors.append(Anchor("cap_hi", lambda z, _u=_union: _u(z, 1), z0, z1))
        cap = CapPlacementRule.from_chord_windows(((float(np.median(vals[:, 0])),
                                                    float(np.median(vals[:, 1]))),))
    _, regions, _, cap_windows, rule, pad = layout_parts(web_fracs, cap)
    return cls(
        tuple(anchors), regions, tuple(webs), te_arc=float(te_arc),
        min_gap=float(min_gap), cap_windows=cap_windows, cap_placement=rule, cap_pad=pad,
    )


__all__ = [
    "coerce_section_layout",
    "graph_from_section_layout",
    "graph_from_windio_structure",
    "layout_parts",
]
