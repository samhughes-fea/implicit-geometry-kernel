"""The blade's conforming shell mesh: quads on the skin and webs, triangles on ribs and the tip.

:func:`shell_mesh` is the kernel's mesh product. Skin and web mid-surfaces are
all-quad mapped strips between span stations. That quad sheet is the node
source: each interior rib and the tip closure is a triangle mesh whose boundary
nodes are indices into the sheet. The fill may add interior nodes only. Every element carries the partition it
lies on and that partition's shell data from the cut
(:class:`~implicit_geometry_kernel.blade.shell_graph.ShellRun`): mesh
thickness, layer names and thicknesses, and the cap eccentricity. Nodes stay on
the skin mid-surface; a consumer applies the eccentricity in its constitutive
law.

A call with no size uses :data:`TIP_H_OVER_C` of the local chord. An explicit
``h`` or ``n_span`` leaves that default off.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field, replace
from typing import Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology import (
    place_station_keys,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
    section_cut_near,
)
from implicit_geometry_kernel.blade.topology.crossings import key_crossing_stations
from implicit_geometry_kernel.blade.topology.map_walk import (
    _skin_mid_loop, te_face_open, te_wrap_for_mesh,
)
from implicit_geometry_kernel.blade.topology.region_counts import region_counts
from implicit_geometry_kernel.blade.topology.ribs import _template_from_cut
from implicit_geometry_kernel.blade.topology.sizing import EPS_OVER_C
from implicit_geometry_kernel.mesh.fill import triangulate_ring
from implicit_geometry_kernel.mesh.quality import (
    QUALITY_DTYPE, MeshRunQuality, _measure, quality_rollup, volume_ratio_from_edges,
)
from implicit_geometry_kernel.mesh.tables import (
    CORNER_DTYPE, EDGE_DTYPE, StationDecision, StripDecision, assemble,
)
from implicit_geometry_kernel.profiling import report_stage, stage

SURFACES = ("skin", "web", "rib", "tip_closure")

TIP_H_OVER_C = 0.05
"""Default shell edge as a fraction of the local chord. ``h`` or ``n_span`` leaves it off."""


@dataclass(frozen=True)
class MeshRecipe:
    """The controls that built a shell mesh, and the edge they resolved.

    ``tip_h_over_c`` is the fraction in force (the default or a caller's
    ``h_over_c``). It is ``None`` when the mesh was sized by ``h`` alone or by
    a station count. ``edge`` is that target length at each station, and empty
    in count mode. ``h_over_c`` records only a fraction the caller passed.
    """

    tip_h_over_c: float | None
    h: float | None
    h_over_c: float | None
    n_span: int | None
    region_h: tuple[tuple[str, float], ...]
    refine: tuple[tuple[float, float, float], ...]
    edge: tuple[float, ...]


@dataclass(frozen=True)
class ShellMesh:
    """Nodes and elements of the blade shell, with per-element shell data.

    ``elements[i]`` lists node indices (3 for a rib or tip triangle, 4 for a quad)
    in loop order. The per-element tuples share that index. Normals (right-hand
    rule on that order) point out of the section on the skin, towards the
    trailing edge on webs, and along +z on ribs and the tip closure; a cap
    stack's mid-plane lies ``eccentricity`` inside the skin, along -normal.
    """

    nodes: np.ndarray
    elements: tuple[tuple[int, ...], ...]
    name: tuple[str, ...]
    """Mesh region: a skin region (``cap_su``, ``bay_pr_0``, …), a web, a rib cell or the tip closure."""
    surface: tuple[str, ...]
    """``skin``, ``web``, ``rib`` or ``tip_closure``."""
    partition: tuple[str, ...]
    """The cut's partition the element lies on (``spar_cap_suction_surface``, a rib name, …)."""
    thickness: np.ndarray
    eccentricity: np.ndarray
    """Cap eccentricity of the partition (``ShellRun.eccentricity_m``, m): half the cap stack; 0 off the caps."""
    layers: tuple[tuple[str, ...], ...]
    layer_thicknesses: tuple[tuple[float, ...], ...]
    stations: tuple[float, ...] = field(default=())
    quality: tuple[MeshRunQuality, ...] = field(default=())
    """Per-partition rollup of :attr:`element_quality`."""
    recipe: MeshRecipe | None = None
    """The size controls and the resolved edge at each station."""
    element_quality: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=QUALITY_DTYPE))
    """One quality row per element: aspect, skew, angles, Jacobian, area, shell volume, volume ratio, warp."""
    n_source: int = 0
    """Skin and web nodes occupy ``nodes[:n_source]``. Interior rib and tip nodes follow."""
    bay: tuple[int, ...] = ()
    """Span bay index of each element. Ribs and the tip use the station index."""
    z0: tuple[float, ...] = ()
    z1: tuple[float, ...] = ()
    constructor: tuple[str, ...] = ()
    """How the element was built: ``ruled``, ``midspan_split``, ``step_patch``, ``transition``, ``fan``, ``interior``."""
    n0: tuple[int, ...] = ()
    """Node count of the region on the root side of the bay. Zero on a rib or the tip."""
    n1: tuple[int, ...] = ()
    flipped: tuple[int, ...] = ()
    """1 when orientation reversed the element's loop."""
    node_origin: tuple[str, ...] = ()
    """Per node: ``station``, ``span_mid``, ``centroid``, ``station_as_midpoint``, or ``interior``."""
    section_distance: np.ndarray = field(default_factory=lambda: np.empty(0))
    """Distance of each node from the section at its span. Station samples are zero."""
    valence: tuple[int, ...] = ()
    edges: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=EDGE_DTYPE))
    """One :data:`~implicit_geometry_kernel.mesh.tables.EDGE_DTYPE` row per undirected pair."""
    loops: tuple = ()
    counts: tuple = ()
    station_log: tuple = ()
    strips: tuple = ()
    near_miss: tuple = ()
    junctions: tuple = ()
    valence_defects: tuple[int, ...] = ()
    edge_sagitta: dict = field(default_factory=dict)
    corners: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=CORNER_DTYPE))
    """One :data:`~implicit_geometry_kernel.mesh.tables.CORNER_DTYPE` row per element corner. Not recomputed on export."""
    diagonal_min: tuple[float, ...] = ()
    """Shortest distance between non-adjacent nodes of each element, metres. A collapsed fan is a tiny diagonal."""

    @property
    def quads(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) == 4)

    @property
    def triangles(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) == 3)

    @property
    def polygons(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) > 4)

    def edge_use(self) -> dict[tuple[int, int], int]:
        """How many elements use each edge (2 on a smooth sheet, more at a junction)."""
        out: dict[tuple[int, int], int] = {}
        for e in self.elements:
            for a, b in zip(e, e[1:] + e[:1]):
                if a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                out[key] = out.get(key, 0) + 1
        return out


def _unique_z(zs: Sequence[float], *, abs_tol: float = 1e-3) -> tuple[float, ...]:
    """Sorted stations merged when closer than ``abs_tol`` (metres).

    A key-crossing can land within a millimetre of an anchor; those must
    collapse to one station so neighbouring bays weld the same nodes.
    """
    ordered = sorted(float(z) for z in zs)
    if not ordered:
        return ()
    out = [ordered[0]]
    for z in ordered[1:]:
        if abs(z - out[-1]) > abs_tol:
            out.append(z)
    return tuple(out)


def _chord_at(geo, z: float) -> float:
    return float(geo._span_scalar(getattr(geo, "_chord", None),
                                  [s.chord for s in geo.stations], float(z)))


def _target_edge(geo, z: float, h: float | None, h_over_c: float | None) -> float:
    """Local quad edge: ``h`` in metres, ``h_over_c`` times the local chord, or the smaller of the two."""
    s = None
    if h_over_c is not None:
        if float(h_over_c) <= 0.0:
            raise ValueError("h_over_c must be positive")
        s = float(h_over_c) * _chord_at(geo, z)
    if h is not None:
        if float(h) <= 0.0:
            raise ValueError("element size h must be positive")
        s = float(h) if s is None else min(float(h), s)
    if s is None or s <= 0.0:
        raise ValueError("element size must be positive")
    return s


def _n_bays(distance: float, s: float) -> int:
    """Equal bays whose step is closest to ``s`` and still inside ``[s/2, s]``.

    A gap of at most ``s`` is one bay, including a short gap between two anchors.
    """
    if s <= 0.0 or distance <= s * (1.0 + 1e-9):
        return 1
    lo = max(int(np.ceil(distance / s - 1e-9)), 1)
    hi = int(np.floor(2.0 * distance / s + 1e-9))
    if hi < lo:
        return lo
    target = distance / s
    return min(range(lo, hi + 1), key=lambda n: (abs(n - target), n))


def _fill_gaps(stations: Sequence[float], s_at) -> tuple[float, ...]:
    """Fit regular stations into the gaps. ``s_at(z)`` is the local target edge."""
    zs = sorted(float(z) for z in stations)
    if len(zs) < 2:
        return tuple(zs)
    out = [zs[0]]
    for z in zs[1:]:
        start = out[-1]
        s = min(s_at(start), s_at(z))
        n = _n_bays(z - start, s)
        if n > 1:
            step = (z - start) / n
            out.extend(start + step * k for k in range(1, n))
        out.append(z)
    return _unique_z(out)


def _drop_stations(geo) -> tuple[float, ...]:
    """Span where a web with no stored end leaves the section.

    Leading-edge and pitch-axis webs are dropped, not slid, when they miss the
    inner skin. The station sits on that change.
    """
    webs = [w for w in getattr(geo.region_graph, "webs", ()) or () if w.z1 is None]
    if not webs:
        return ()
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    tip = section_cut_near(geo, z1)
    if tip is None:
        return ()
    dropped_tip = {int(i) for i in (getattr(tip, "dropped_web_indices", ()) or ())}
    need = [w for w in webs if int(w.index) in dropped_tip]
    if not need:
        return ()
    root = section_cut_near(geo, z0)
    dropped_root = {int(i) for i in (getattr(root, "dropped_web_indices", ()) or ())} if root is not None else set()
    found: list[float] = []
    span = max(z1 - z0, 1.0)
    for web in need:
        if int(web.index) in dropped_root:
            continue
        lo, hi = z0, z1
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            cut = section_cut_near(geo, mid)
            dropped = {int(i) for i in (getattr(cut, "dropped_web_indices", ()) or ())} if cut is not None else set()
            if int(web.index) in dropped:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-6 * span:
                break
        if z0 + 1e-9 < hi < z1 - 1e-9:
            found.append(hi)
    return tuple(found)


def _anchor_stations(geo) -> tuple[float, ...]:
    """Spanwise stations the section sample must contain.

    The root, the tip and every defining station. The birth and death of each
    web, including a web that leaves the section with no authored end. Every
    rib and other transverse member, the ends of stringers and ply drops, and
    a key crossing. The section sample is built from this list, and the mesh
    uses the same stations as its anchors.
    """
    cached = geo.__dict__.get("_mesh_anchor_z")
    if cached is not None:
        return cached
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = [float(station.z) for station in geo.stations]
    zs.append(float(geo.tip_closure_z))
    zs.extend(float(z) for z in geo.forced_stations())
    zs.extend(float(member.z) for member in geo.transverse)
    for anchor in getattr(geo.region_graph, "anchors", ()) or ():
        for z in (getattr(anchor, "z0", None), getattr(anchor, "z1", None)):
            if z is not None and z0 < float(z) < z1:
                zs.append(float(z))
    zs.extend(_drop_stations(geo))
    found = _unique_z(z for z in zs if z0 - 1e-9 <= z <= z1 + 1e-9)
    found = _unique_z(list(found) + list(key_crossing_stations(geo, found)))
    geo.__dict__["_mesh_anchor_z"] = found
    return found


def _refine_stations(zs: Sequence[float], z0: float, z1: float,
                     refine: Sequence[tuple[float, float, float]]) -> tuple[float, ...]:
    out = set(zs)
    for a, b, dz in refine:
        a, b = max(float(min(a, b)), z0), min(float(max(a, b)), z1)
        if float(dz) <= 0.0:
            raise ValueError("refinement spacing dz must be positive")
        if b > a:
            count = int(np.ceil((b - a) / float(dz) - 1e-9)) + 1
            out.update(float(z) for z in np.linspace(a, b, count))
    return _unique_z(out)


def _median_skin_edge(tmpl, smp) -> float | None:
    """Median chordwise edge of the skin regions. Webs and the blunt wrap are left out."""
    lengths = []
    for branch in tmpl.loop:
        if branch.wrap or branch.name not in smp:
            continue
        pts = np.asarray(smp[branch.name], float)
        edges = len(pts) - 1
        if edges <= 0:
            continue
        length = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        if length > 0.0:
            lengths.append(length / edges)
    if not lengths:
        return None
    return float(np.median(lengths))


def _grade_stations(geo, zs: Sequence[float], tmpls, smps, h: float | None,
                    h_over_c: float | None):
    """Split a bay longer than about 1.25 times the outboard skin edge. One pass.

    The new step stays at least half the local target, so a station is not
    placed on top of an anchor. A skin edge already shorter than that half is
    the geometry (a narrow cap, a short web), and the span step does not chase it.
    A bay left long is returned as a refused grade, with the reason.
    """
    anchors = {round(z, 9) for z in _anchor_stations(geo)}
    out = [float(zs[0])]
    changed = False
    refused: list[StationDecision] = []
    for i in range(len(zs) - 1):
        za, zb = float(zs[i]), float(zs[i + 1])
        med = _median_skin_edge(tmpls[i + 1], smps[i + 1])
        bay = zb - za
        mid = 0.5 * (za + zb)
        if med is not None and bay > 1.25 * med * (1.0 + 1e-12):
            step_min = 0.5 * _target_edge(geo, zb, h, h_over_c)
            if med < step_min:
                refused.append(StationDecision(
                    mid, "grade", True, "skin edge already under half the target"))
            else:
                n = max(int(np.ceil(bay / med - 1e-9)), 1)
                step = bay / n
                if n > 1 and step + 1e-9 >= step_min:
                    clear = True
                    for k in range(1, n):
                        p = za + step * k
                        for end in (za, zb):
                            if round(end, 9) in anchors and abs(p - end) < step_min - 1e-8:
                                clear = False
                    if clear:
                        out.extend(za + step * k for k in range(1, n))
                        changed = True
                    else:
                        refused.append(StationDecision(
                            mid, "grade", True, "step would land on an anchor"))
                else:
                    refused.append(StationDecision(
                        mid, "grade", True, "split step under half the target"))
        out.append(zb)
    return _unique_z(out), changed, tuple(refused)


def mesh_stations(geo, n_span: int | None = None, *, h: float | None = None,
                  h_over_c: float | None = None,
                  refine: Sequence[tuple[float, float, float]] = ()) -> tuple[float, ...]:
    """Span stations of the mesh.

    With ``h`` (metres) or ``h_over_c`` (a fraction of the local chord), the
    anchors are the defining stations, the root, the tip and the topology
    events, and regular stations fill the gaps so each bay is in ``[s/2, s]``
    of the local target edge.
    With ``n_span`` and no size, the stations are that count, plus the same
    anchors. ``refine`` holds ``(z0, z1, dz)`` windows. A key crossing found
    between stations is kept, and that bay is filled again.
    """
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    if h is not None and float(h) <= 0.0:
        raise ValueError("element size h must be positive")
    if h_over_c is not None and float(h_over_c) <= 0.0:
        raise ValueError("h_over_c must be positive")
    anchors = _anchor_stations(geo)
    if h is None and h_over_c is None:
        if n_span is None:
            raise ValueError("give n_span, h, or h_over_c")
        zs = _unique_z(list(np.linspace(z0, z1, int(n_span))) + list(anchors))
        zs = _refine_stations(zs, z0, z1, refine)
        base = [z for z in zs if z0 - 1e-9 <= z <= z1 + 1e-9]
        return _unique_z(list(base) + list(key_crossing_stations(geo, base)))
    s_at = lambda z, _h=h, _f=h_over_c: _target_edge(geo, z, _h, _f)
    filled = _fill_gaps(anchors, s_at)
    filled = _refine_stations(filled, z0, z1, refine)
    crossings = key_crossing_stations(geo, filled)
    if crossings:
        filled = _fill_gaps(tuple(filled) + tuple(crossings), s_at)
    return filled


def _run_distances(run, pts: np.ndarray) -> np.ndarray:
    """Distance of each ``(q, 2)`` point to one shell run."""
    xy = np.asarray(run.points, float)[:, :2]
    query = np.asarray(pts, float).reshape(-1, 2)
    if len(xy) == 0:
        return np.full(len(query), np.inf)
    if len(xy) == 1:
        return np.linalg.norm(query - xy[0], axis=1)
    a = xy[:-1]
    ab = np.diff(xy, axis=0)
    offset = query[:, None, :] - a[None, :, :]
    denom = np.maximum((ab * ab).sum(axis=1), 1e-30)
    t = np.clip(np.sum(offset * ab[None, :, :], axis=2) / denom, 0.0, 1.0)
    proj = a[None, :, :] + t[:, :, None] * ab[None, :, :]
    return np.linalg.norm(proj - query[:, None, :], axis=2).min(axis=1)


def _nearest_runs(runs, pts: np.ndarray):
    """The run nearest each point. Earlier runs win a tie, matching the one-point scan."""
    query = np.asarray(pts, float).reshape(-1, 2)
    best = [None] * len(query)
    best_d = np.full(len(query), np.inf)
    for run in runs:
        distance = _run_distances(run, query)
        take = distance < best_d
        best_d = np.where(take, distance, best_d)
        for i in np.flatnonzero(take):
            best[int(i)] = run
    return best


def _nearest_run(runs, p: np.ndarray):
    return _nearest_runs(runs, np.asarray(p, float).reshape(1, 2))[0]


def _area_vector(p: np.ndarray) -> np.ndarray:
    """Newell normal of a polygon: its area times its unit normal (right-hand rule on the node order)."""
    return 0.5 * np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)


def _plain_ruled(pairs) -> bool:
    """True when every strip is a run of four-corner quads with no fold.

    A bowtie or a collapsed corner stays on the face loft, which already splits
    those strips.
    """
    from implicit_geometry_kernel.blade.topology.map_loft import _ruled_quad_array, _ruled_quad_masks

    for _name, pa, pb in pairs:
        pa = np.asarray(pa, float)
        pb = np.asarray(pb, float)
        if len(pa) < 2 or len(pa) != len(pb):
            return False
        keep, bow = _ruled_quad_masks(_ruled_quad_array(pa, pb))
        if len(keep) == 0 or not np.all(keep) or np.any(bow):
            return False
    return True


def _ruled_station_queries(quads: np.ndarray, z_root: float, z_tip: float) -> tuple[np.ndarray, np.ndarray]:
    """Query point and tip-station flag for each ruled quad.

    The query is the mean of the corners that lie on the nearer station, which
    is the same point the per-face scan uses. A tie stays on the root station.
    """
    zc = quads[:, :, 2].mean(axis=1)
    tip = np.abs(zc - z_root) > np.abs(z_tip - zc)
    z_sel = np.where(tip, z_tip, z_root)
    on = np.abs(quads[:, :, 2] - z_sel[:, None]) <= 1e-12
    use = on.copy()
    empty = ~on.any(axis=1)
    use[empty] = True
    weight = use.astype(float)
    xy = (quads[:, :, :2] * weight[:, :, None]).sum(axis=1) / np.maximum(weight.sum(axis=1)[:, None], 1.0)
    return xy, tip


def _areas_means(nodes: np.ndarray, elements: list[tuple[int, ...]],
                 indices: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """Area vectors and vertex means of ``elements[indices]``, one row per index."""
    areas = np.zeros((len(indices), 3), float)
    means = np.zeros((len(indices), 3), float)
    grouped: dict[int, list[int]] = {}
    for local, element_id in enumerate(indices):
        grouped.setdefault(len(elements[element_id]), []).append(local)
    for n, locals_ in grouped.items():
        if n < 1:
            continue
        loc = np.asarray(locals_, np.int64)
        chosen = [indices[i] for i in locals_]
        block = nodes[np.asarray([list(elements[i]) for i in chosen], np.int64)]
        means[loc] = block.mean(axis=1)
        if n >= 3:
            areas[loc] = 0.5 * np.sum(np.cross(block, np.roll(block, -1, axis=1)), axis=1)
    return areas, means


def _oriented(nodes: np.ndarray, elements: list[tuple[int, ...]], surface: Sequence[str],
              zs: Sequence[float], te_at: Mapping[float, np.ndarray]
              ) -> tuple[list[tuple[int, ...]], list[int]]:
    """Elements reordered so that their normals follow one rule per surface.

    Skin: out of the section. The skin is made consistent across every edge two
    skin elements share, then each connected sheet takes the sign for which
    ``sum(A_e . (x_e - c))`` is positive (``A_e`` area vector, ``c`` any point
    of the station): for a closed section that sum is twice its area with
    outward normals, whatever the section's concavity. Webs: towards the
    trailing edge of their station. Ribs and the tip closure: along +z.
    """
    out = list(elements)
    flipped = [0] * len(out)
    z_st = np.asarray(sorted(zs), float)

    def flip(i: int) -> None:
        e = out[i]
        out[i] = (e[0],) + tuple(reversed(e[1:]))
        flipped[i] ^= 1

    def station(x: np.ndarray) -> float:
        return float(z_st[int(np.argmin(np.abs(z_st - x[2])))])

    skin = [i for i, s in enumerate(surface) if s == "skin"]
    by_edge: dict[tuple[int, int], list[int]] = {}
    for i in skin:
        e = out[i]
        for a, b in zip(e, e[1:] + e[:1]):
            by_edge.setdefault((min(a, b), max(a, b)), []).append(i)

    def directed(i: int, a: int, b: int) -> bool:
        e = out[i]
        return any(u == a and v == b for u, v in zip(e, e[1:] + e[:1]))

    seen: set[int] = set()
    for start in skin:
        if start in seen:
            continue
        seen.add(start)
        sheet, stack = [start], [start]
        while stack:
            i = stack.pop()
            e = out[i]
            for a, b in zip(e, e[1:] + e[:1]):
                users = by_edge[(min(a, b), max(a, b))]
                if len(users) != 2:
                    continue
                j = users[0] if users[1] == i else users[1]
                if j in seen:
                    continue
                if directed(j, a, b):  # a consistent neighbour runs the shared edge the other way
                    flip(j)
                seen.add(j)
                sheet.append(j)
                stack.append(j)
        centre: dict[float, np.ndarray] = {}
        for i in sheet:
            for k in out[i]:
                centre.setdefault(station(nodes[k]), nodes[k])
        areas, means = _areas_means(nodes, out, sheet)
        station_z = z_st[np.argmin(np.abs(means[:, 2:3] - z_st), axis=1)]
        # A face between two stations can sit nearer a third, so none of its
        # corners were filed under that station. Any section point fixes the
        # sign: on a closed sheet the constant cancels.
        keys = list(centre)
        key_z = np.asarray(keys, float)
        picked = []
        for z in station_z:
            key = float(z)
            point = centre.get(key)
            if point is None:
                point = centre[keys[int(np.argmin(np.abs(key_z - key)))]]
            picked.append(point)
        centres = np.asarray(picked, float)
        flux = float(np.sum(areas[:, :2] * (means[:, :2] - centres[:, :2])))
        if flux < 0.0:
            for i in sheet:
                flip(i)
    webs = [i for i, s in enumerate(surface) if s == "web"]
    areas, means = _areas_means(nodes, out, webs)
    if len(webs):
        station_z = z_st[np.argmin(np.abs(means[:, 2:3] - z_st), axis=1)]
        for i, area, mean, z in zip(webs, areas, means, station_z):
            te = te_at.get(float(z))
            if te is not None and float(area[:2] @ (te - mean)[:2]) < 0.0:
                flip(i)
    ribs = [i for i, s in enumerate(surface) if s in ("rib", "tip_closure")]
    areas, _means = _areas_means(nodes, out, ribs)
    for i, area in zip(ribs, areas):
        if area[2] < 0.0:
            flip(i)
    return out, flipped


def orient_shell_elements(
    nodes: np.ndarray,
    elements: list[tuple[int, ...]],
    surface: Sequence[str],
    zs: Sequence[float],
    te_at: Mapping[float, np.ndarray],
) -> list[tuple[int, ...]]:
    """Public wrapper: reorder element corners so normals follow the shell rule.

    Skin normals point out of the section, webs toward the trailing edge of
    their station, ribs and tip closures along +z. Used by comparison remeshers
    that build their own connectivity.
    """
    out, _flipped = _oriented(nodes, elements, surface, zs, te_at)
    return out


def _size_fraction(n_span: int | None, h: float | None, h_over_c: float | None,
                   explicit: bool) -> float | None:
    """The chord fraction in force, or ``None`` when the caller sized by metres or by count.

    A bare call uses :data:`TIP_H_OVER_C`. ``h`` alone, ``n_span``, and an
    explicit station list leave that default off. A passed ``h_over_c`` is the
    fraction and overrides the default.
    """
    if h_over_c is not None:
        if float(h_over_c) <= 0.0:
            raise ValueError("h_over_c must be positive")
        return float(h_over_c)
    if h is not None or n_span is not None or explicit:
        return None
    return TIP_H_OVER_C


def _point_key(p: np.ndarray) -> tuple[float, float, float]:
    return (float(p[0]), float(p[1]), float(p[2]))


def _split_crossed_triangles(elements: list, rows: list, meta: list,
                             nodes: list[np.ndarray], tol: float = 1e-8) -> None:
    """Split a rib or tip triangle whose edge passes through another node.

    A straight trailing edge puts a station node on the chord between two
    corners. The unstructured fill can draw an interior edge along that chord,
    and the station node is then a T-junction. Splitting the triangle at the
    node makes the node a corner. Rib and tip edges use one batched segment
    query; the area reject applies to those triangles.
    """
    from implicit_geometry_kernel.mesh.tables import _cell_groups, _segment_hits

    pts = np.asarray(nodes, float)
    if len(pts) == 0:
        return
    cell = 0.05
    for _pass in range(8):
        tri = [
            i for i, element in enumerate(elements)
            if len(element) == 3 and rows[i][1] in ("rib", "tip_closure")
        ]
        if not tri:
            break
        edges_a: list[int] = []
        edges_b: list[int] = []
        edges_c: list[int] = []
        edges_ei: list[int] = []
        for ei in tri:
            a, b, c = (int(v) for v in elements[ei])
            for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
                edges_a.append(u)
                edges_b.append(v)
                edges_c.append(w)
                edges_ei.append(ei)
        a = np.asarray(edges_a, np.int64)
        b = np.asarray(edges_b, np.int64)
        pa = pts[a]
        pb = pts[b]
        ab = pb - pa
        span2 = np.maximum(np.sum(ab * ab, axis=1), 1e-30)
        length = np.sqrt(span2)
        live = np.flatnonzero(length > tol)
        if len(live) == 0:
            break
        a, b = a[live], b[live]
        pa, ab, span2 = pa[live], ab[live], span2[live]
        ei_of = np.asarray(edges_ei, np.int64)[live]
        c_of = np.asarray(edges_c, np.int64)[live]
        lo = np.minimum(pa, pts[b]) - tol
        hi = np.maximum(pa, pts[b]) + tol
        groups = _cell_groups(np.floor(pts / cell).astype(np.int64))
        klo = np.floor(lo / cell).astype(np.int64)
        khi = np.floor(hi / cell).astype(np.int64)
        by_cell: dict[tuple[int, int, int], list[int]] = {}
        for edge in range(len(a)):
            for ix in range(int(klo[edge, 0]), int(khi[edge, 0]) + 1):
                for iy in range(int(klo[edge, 1]), int(khi[edge, 1]) + 1):
                    for iz in range(int(klo[edge, 2]), int(khi[edge, 2]) + 1):
                        by_cell.setdefault((ix, iy, iz), []).append(edge)
        edge_of = {(int(a[i]), int(b[i])): int(i) for i in range(len(a))}
        found = None
        for key, edge_ids in by_cell.items():
            nodes_here = groups.get(key)
            if nodes_here is None:
                continue
            hits = _segment_hits(
                pts, np.asarray(edge_ids, np.int64), nodes_here, a, b, pa, ab, span2, lo, hi, tol,
            )
            for ia, ib, node, _gap in hits:
                edge = edge_of.get((ia, ib))
                if edge is None:
                    continue
                owned = {int(a[edge]), int(b[edge]), int(c_of[edge])}
                if node in owned:
                    continue
                found = (int(ei_of[edge]), int(a[edge]), int(b[edge]), int(c_of[edge]), int(node))
                break
            if found is not None:
                break
        if found is None:
            break
        ei, ua, ub, uc, node = found
        elements[ei] = (ua, node, uc)
        elements.append((node, ub, uc))
        rows.append(rows[ei])
        meta.append(meta[ei])
    kept_e: list = []
    kept_r: list = []
    kept_m: list = []
    for element, row, rec in zip(elements, rows, meta):
        if len(element) == 3 and row[1] in ("rib", "tip_closure"):
            poly = pts[list(element)]
            area = float(np.linalg.norm(np.cross(poly[1] - poly[0], poly[2] - poly[0])))
            if area <= 1e-12:
                continue
        kept_e.append(element)
        kept_r.append(row)
        kept_m.append(rec)
    elements[:] = kept_e
    rows[:] = kept_r
    meta[:] = kept_m


def _weld_triangles(faces: np.ndarray, index: dict[tuple[float, float, float], int],
                    nodes: list[np.ndarray], origins: list[str]) -> np.ndarray:
    """Map every corner of ``faces`` onto a node id.

    Coordinates already in ``index`` keep that id. New points are appended in
    the order of one exact ``(x, y, z)`` unique, the same key as :func:`_point_key`.
    A triangle that is not three distinct nodes is marked with ``-1``.
    """
    flat = np.ascontiguousarray(np.asarray(faces, float).reshape(-1, 3))
    if len(flat) == 0:
        return np.zeros((0, 3), np.int64)
    void = flat.view(np.dtype((np.void, flat.dtype.itemsize * flat.shape[1]))).ravel()
    _unique, first, inverse = np.unique(void, return_index=True, return_inverse=True)

    def _coarse(p) -> tuple[float, float, float]:
        return tuple(round(float(p[k]), 8) for k in range(3))

    near = {_coarse(p): i for i, p in enumerate(nodes)}
    ids = np.empty(len(first), np.int64)
    for i, p in enumerate(flat[first]):
        key = _point_key(p)
        k = index.get(key)
        if k is None:
            k = near.get(_coarse(p))
        if k is None:
            k = len(nodes)
            index[key] = k
            nodes.append(np.asarray(key, float))
            origins.append("interior")
            near[_coarse(p)] = k
        ids[i] = k
    corners = ids[inverse].reshape(-1, 3)
    collapsed = (
        (corners[:, 0] == corners[:, 1])
        | (corners[:, 1] == corners[:, 2])
        | (corners[:, 0] == corners[:, 2])
    )
    corners[collapsed] = -1
    return corners


def shell_mesh(geo, n_span: int | None = None, *, h: float | None = None,
               h_over_c: float | None = None, n_chord: int = 8,
               zs: Sequence[float] | None = None, polygon_size: float | None = None,
               region_h: Mapping[str, float] | None = None,
               refine: Sequence[tuple[float, float, float]] = ()) -> ShellMesh:
    """The blade's conforming shell mesh.

    Stations are ``zs`` or :func:`mesh_stations`. With no size and no station
    count the edge is :data:`TIP_H_OVER_C` of the local chord. ``h`` (metres)
    alone leaves that default off. ``h_over_c`` is the same fraction and
    overrides the default; passing it with ``h`` uses the smaller edge at each
    station. ``n_span`` with no size stays a station count. Each region's node
    count then follows its length
    (:func:`~implicit_geometry_kernel.blade.topology.region_counts.region_counts`);
    a count with no size gives every region ``n_chord`` nodes (odd).
    ``region_h`` sets the element size of named regions
    (``{"le_su": h / 3, "le_pr": h / 3}`` refines the leading edge) and
    ``refine`` adds span stations in ``(z0, z1, dz)`` windows (see
    :func:`mesh_stations`). ``polygon_size`` is the rib and tip triangle spacing
    (default: the median edge of the station ring).

    The skin and web quads are built first. Every interior rib and the tip
    closure indexes those nodes. The triangle fill adds interior nodes only.
    """
    with stage("implicit_geometry_kernel", "shell_mesh") as rec:
        mesh = _shell_mesh_body(
            geo, n_span, h=h, h_over_c=h_over_c, n_chord=n_chord, zs=zs,
            polygon_size=polygon_size, region_h=region_h, refine=refine)
        rec.set("elements", len(mesh.elements))
        if h is not None:
            rec.set("edge_length", float(h))
        if mesh.recipe is not None and mesh.recipe.tip_h_over_c is not None:
            rec.set("h_over_c", float(mesh.recipe.tip_h_over_c))
        return mesh


def _closest_on_loop_many(loop: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Closest point of a closed polyline to each row of ``pts``."""
    query = np.asarray(pts, float).reshape(-1, 3)
    ring = np.asarray(loop, float)
    if len(ring) < 2 or len(query) == 0:
        return query.copy()
    closed = np.vstack([ring, ring[:1]])
    a = closed[:-1]
    ab = np.diff(closed, axis=0)
    denom = np.maximum((ab * ab).sum(axis=1), 1e-30)
    offset = query[:, None, :] - a[None, :, :]
    t = np.clip(np.sum(offset * ab[None, :, :], axis=2) / denom, 0.0, 1.0)
    proj = a[None, :, :] + t[:, :, None] * ab[None, :, :]
    pick = np.argmin(np.linalg.norm(proj - query[:, None, :], axis=2), axis=1)
    return proj[np.arange(len(query)), pick]


def _closest_on_loop(loop: np.ndarray, xyz: np.ndarray) -> np.ndarray:
    """Closest point of a closed polyline to ``xyz``."""
    return _closest_on_loop_many(loop, np.asarray(xyz, float).reshape(1, 3))[0]


def _le_te_points(geo, zs: Sequence[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Placed leading and trailing edges at each ``z``, and the local chord."""
    zs = np.asarray([float(z) for z in zs], float)
    empty = np.zeros((0, 3), float)
    if len(zs) == 0:
        return empty, empty, np.zeros(0, float)
    chord, twist, pitch, offset, x, upper, lower = geo.interpolate_many(zs)
    signed = np.array([geo._placement_twist(float(t)) for t in twist], float)
    i_le = int(np.argmin(x))
    i_te = int(np.argmax(x))

    def at(index: int) -> np.ndarray:
        y = 0.5 * (upper[:, index] + lower[:, index])
        local_x = float(x[index]) * chord - pitch * chord
        local_y = y * chord
        angle = np.radians(signed)
        c, s = np.cos(angle), np.sin(angle)
        return np.column_stack([
            c * local_x + s * local_y + offset[:, 0],
            -s * local_x + c * local_y + offset[:, 1],
            zs,
        ])

    return at(i_le), at(i_te), np.asarray(chord, float)


def _fill_le_te_cache(geo, zs: Sequence[float], cache: dict) -> None:
    """Store leading edge, trailing edge, and chord for each missing ``z``."""
    need = [float(z) for z in zs if round(float(z), 9) not in cache]
    if not need:
        return
    le, te, chord = _le_te_points(geo, need)
    for i, z in enumerate(need):
        cache[round(float(z), 9)] = (le[i], te[i], float(chord[i]))


def _span_edges(geo, zs: Sequence[float], *, cache: dict | None = None):
    """Placed leading and trailing edges, and each bay's sagitta and budget.

    The edge point is the unit-chord extreme of the interpolated airfoil, placed
    with the station's chord, twist and offset. That is the point ``_le_te_point``
    used to build one station at a time. Sagitta is the bow of the straight
    ruling, so it does not build a section. ``cache`` keeps evaluated edges
    across freeze rounds.
    """
    zs = np.asarray([float(z) for z in zs], float)
    empty = np.zeros((0, 3), float)
    if len(zs) == 0:
        return empty, empty, np.zeros(0), np.zeros(0), np.zeros(0)
    mids = 0.5 * (zs[:-1] + zs[1:]) if len(zs) > 1 else np.zeros(0, float)
    store = cache if cache is not None else {}
    query = np.concatenate([zs, mids]) if len(mids) else zs
    _fill_le_te_cache(geo, query, store)
    le = np.stack([store[round(float(z), 9)][0] for z in zs])
    te = np.stack([store[round(float(z), 9)][1] for z in zs])
    n = len(zs)
    if n < 2:
        return le, te, np.zeros(0), np.zeros(0), np.zeros(0)
    mid_le = np.stack([store[round(float(z), 9)][0] for z in mids])
    mid_te = np.stack([store[round(float(z), 9)][1] for z in mids])
    budget = EPS_OVER_C * np.asarray([store[round(float(z), 9)][2] for z in mids], float)

    def sag(a: np.ndarray, mid: np.ndarray, b: np.ndarray) -> np.ndarray:
        ab = b - a
        length = np.linalg.norm(ab, axis=1)
        out = np.zeros(n - 1, float)
        ok = length >= 1e-12
        if np.any(ok):
            step = np.sum((mid[ok] - a[ok]) * ab[ok], axis=1) / (length[ok] * length[ok])
            out[ok] = np.linalg.norm(mid[ok] - (a[ok] + step[:, None] * ab[ok]), axis=1)
        return out

    return le, te, sag(le[:-1], mid_le, le[1:]), sag(te[:-1], mid_te, te[1:]), budget


def _le_te_point(geo, z: float, which: str) -> np.ndarray:
    """Leading or trailing edge of the shell placement at ``z``, from the loft interpolation."""
    le, te, *_rest = _span_edges(geo, (float(z),))
    return le[0] if which == "le" else te[0]


def _span_sagitta(geo, za: float, zb: float) -> tuple[float, float, float]:
    """``(le, te, budget)`` sagitta of the straight ruling against the section at mid-span."""
    _le, _te, le, te, budget = _span_edges(geo, (float(za), float(zb)))
    return float(le[0]), float(te[0]), float(budget[0])


def _sagitta_pending(geo, zs: Sequence[float], h: float | None, fraction: float | None,
                     present: set[float], *, cache: dict | None = None
                     ) -> tuple[list[tuple[float, str]], list[StationDecision]]:
    """Mid-bay stations whose leading or trailing edge bows past the budget.

    One evaluation covers every bay. A bay shorter than half the local target,
    or a mid-point that is already a station, is refused and not inserted.
    """
    _le, _te, le, te, budget = _span_edges(geo, zs, cache=cache)
    pending: list[tuple[float, str]] = []
    refused: list[StationDecision] = []
    for i in range(len(zs) - 1):
        if max(float(le[i]), float(te[i])) <= float(budget[i]):
            continue
        mid = 0.5 * (float(zs[i]) + float(zs[i + 1]))
        bay = float(zs[i + 1]) - float(zs[i])
        target = _target_edge(geo, float(zs[i + 1]), h, fraction)
        if bay < 0.5 * target:
            refused.append(StationDecision(mid, "sagitta", True, "short bay before an apex"))
            continue
        key = round(mid, 9)
        if key in present or any(abs(mid - float(s)) <= 1e-6 for s in zs):
            refused.append(StationDecision(mid, "sagitta", True, "step would land on a station"))
        else:
            pending.append((mid, "sagitta"))
            present.add(key)
    return pending, refused


def _freeze_sagitta(geo, zs: Sequence[float], h: float | None, fraction: float | None) -> tuple[float, ...]:
    """Insert mid-bay stations until the leading and trailing edges hold the budget.

    This runs before any section is built. Eight rounds is the same depth the
    mesh used when it rebuilt the whole station list each pass. Edge points
    already evaluated stay in the cache; only new midpoints are interpolated.
    """
    chosen = [float(z) for z in zs]
    inserted: dict[float, str] = {}
    cache: dict = {}
    for _round in range(8):
        present = {round(z, 9) for z in chosen}
        pending, _refused = _sagitta_pending(geo, chosen, h, fraction, present, cache=cache)
        if not pending:
            break
        for z, origin in pending:
            inserted[round(float(z), 9)] = origin
        chosen = list(_unique_z(chosen + [z for z, _origin in pending]))
    return tuple(chosen), inserted


def _web_fraction_pairs(geo, zs: Sequence[float]) -> np.ndarray:
    """``(n, w, 2)`` suction and pressure chord fractions. Empty when the layout has no webs."""
    rows: list[tuple[tuple[float, float], ...]] = []
    width = 0
    for z in zs:
        feet = geo._web_feet_at(float(z))
        if feet:
            pairs = tuple((float(a), float(b)) for a, b in feet)
        else:
            layout = geo._layout_at(float(z))
            fracs = tuple(float(f) for f in getattr(layout, "web_chord_fracs", ()) or ())
            pairs = tuple((f, f) for f in fracs)
        width = max(width, len(pairs))
        rows.append(pairs)
    out = np.zeros((len(rows), width, 2), float)
    for i, pairs in enumerate(rows):
        if not pairs:
            continue
        out[i, :len(pairs)] = pairs
        if len(pairs) < width:
            out[i, len(pairs):] = pairs[-1]
    return out


def _stacked_shell(geo, zs: Sequence[float]):
    """Skin midline and web lines of every station, from one placed airfoil stack.

    The skin is the outer mould offset along its vertex normals by half the
    skin thickness. A web is the chord fraction on the upper and lower surfaces.
    Cap and web fractions then split that skin by arc length. The named region
    polylines stay on the section until this skin matches the section midline.
    """
    from implicit_geometry_kernel.blade.loft.span_place import offset_rings, place_oml, split_at_arc, web_lines

    zs = tuple(float(z) for z in zs)
    chord, twist, pitch, offset, x, upper, lower = geo.interpolate_many(zs)
    signed = np.array([geo._placement_twist(float(t)) for t in twist], float)
    rings = place_oml(chord, signed, pitch, offset, x, upper, lower, np.asarray(zs, float))
    half = np.array([0.5 * float(geo._skin(float(z))) for z in zs], float)
    skin = offset_rings(rings, half)
    pairs = _web_fraction_pairs(geo, zs)
    webs = web_lines(x, upper, lower, chord, pitch, signed, offset, np.asarray(zs, float), pairs)
    fractions = pairs[:, :, 0] if pairs.shape[1] else np.zeros((len(zs), 0), float)
    splits = split_at_arc(skin, fractions) if fractions.size else np.zeros((len(zs), 0, 3), float)
    le, te, le_sag, te_sag, budget = _span_edges(geo, zs)
    return skin, webs, splits, le, te, le_sag, te_sag, budget


def _shared_count_gap(tmpls) -> list[int]:
    """Bays where a region present on both stations changes its node count by more than two."""
    gaps = []
    for i in range(len(tmpls) - 1):
        a = {b.name: int(b.n) for b in (*tmpls[i].loop, *tmpls[i].webs)}
        b = {b.name: int(b.n) for b in (*tmpls[i + 1].loop, *tmpls[i + 1].webs)}
        if any(abs(a[name] - b[name]) > 2 for name in set(a) & set(b)):
            gaps.append(i)
    return gaps


def _station_kind(geo, z: float, inserted: Mapping[float, str], explicit: bool) -> str:
    """Anchor kind of a kept station. Inserted stations keep the reason they were added."""
    key = round(float(z), 9)
    if key in inserted:
        return inserted[key]
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    if abs(z - z0) <= 1e-8:
        return "root"
    if abs(z - z1) <= 1e-8 or abs(z - float(geo.tip_closure_z)) <= 1e-8:
        return "tip"
    for member in geo.transverse:
        if abs(float(member.z) - z) <= 1e-8:
            return "rib"
    for web in getattr(geo.region_graph, "webs", ()) or ():
        for end in (getattr(web, "z0", None), getattr(web, "z1", None)):
            if end is not None and abs(float(end) - z) <= 1e-8:
                return "web_end"
    for stringer in getattr(geo, "stringers", ()) or ():
        extent = getattr(stringer, "extent", None)
        if extent is not None and any(abs(float(end) - z) <= 1e-8 for end in extent):
            return "stringer"
    for drop in getattr(geo, "ply_drops", ()) or ():
        if any(abs(float(end) - z) <= 1e-8 for end in (drop.z0, drop.z1)):
            return "ply_drop"
    return "caller" if explicit else "grade"


def _nearest_key(keys: Mapping[tuple[str, str], np.ndarray], xyz: np.ndarray) -> np.ndarray | None:
    """The station key nearest ``xyz`` in the section, with ``xyz``'s span coordinate."""
    query = np.asarray(xyz, float)
    best = None
    best_d = 1e300
    for point in keys.values():
        p = np.asarray(point, float)
        d = float((p[0] - query[0]) ** 2 + (p[1] - query[1]) ** 2)
        if d < best_d:
            best, best_d = p, d
    if best is None:
        return None
    out = np.asarray(best, float).copy()
    out[2] = float(query[2])
    return out


def _loop_wind(unit: np.ndarray) -> float:
    """Signed area of an open chord-frame loop. The sign is the winding."""
    if len(unit) < 3:
        return 0.0
    x, y = unit[:, 0], unit[:, 1]
    return float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _open_unit_loop(loop: np.ndarray, keys: Mapping, cut) -> tuple[np.ndarray, dict]:
    """Chord-frame loop and keys, opened so the leading-edge key is the first point."""
    from implicit_geometry_kernel.blade.loft.span_place import to_chord_normal

    pts = np.asarray(loop, float)
    if len(pts) > 1 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    if len(pts) == 0:
        return np.zeros((0, 2), float), {}
    chord, twist, pitch, offset = cut.chord, cut.twist, cut.pitch_axis, cut.offset
    unit = to_chord_normal(pts[:, :2], chord, twist, pitch, offset)
    le = None
    for side in ("suction", "pressure"):
        hit = keys.get(("le", side))
        if hit is not None:
            le = to_chord_normal(
                np.asarray(hit, float)[:2].reshape(1, 2), chord, twist, pitch, offset)[0]
            break
    start = int(np.argmin(unit[:, 0]) if le is None else np.argmin(np.linalg.norm(unit - le, axis=1)))
    unit = np.roll(unit, -start, axis=0)
    key_unit = {}
    for name, xyz in keys.items():
        key_unit[name] = to_chord_normal(
            np.asarray(xyz, float)[:2].reshape(1, 2), chord, twist, pitch, offset)[0]
    return unit, key_unit


def _frame_at(geo, z: float) -> tuple[float, float, float, tuple[float, float]]:
    """Chord, signed twist, pitch and offset that place a point at ``z``."""
    chord, twist, pitch, offset, *_rest = geo.interpolate(float(z))
    offset = (float(offset[0]), float(offset[1]))
    return float(chord), float(geo._placement_twist(twist)), float(pitch), offset


def _blend_frames(left_unit, right_unit, left_keys, right_keys, weight: float,
                  chord: float, twist: float, pitch: float, offset, z: float):
    """Place the chord-frame blend of two station midlines at ``z``."""
    left_unit = np.asarray(left_unit, float)
    right_unit = np.asarray(right_unit, float)
    if len(left_unit) < 2 and len(right_unit) < 2:
        return None, {}
    if len(left_unit) < 2:
        unit = right_unit
    elif len(right_unit) < 2:
        unit = left_unit
    else:
        other = right_unit
        wind_l, wind_r = _loop_wind(left_unit), _loop_wind(other)
        if wind_l * wind_r < 0.0 and abs(wind_l) > 1e-12 and abs(wind_r) > 1e-12:
            other = np.roll(other[::-1], 1, axis=0)
        unit = _blend_open(left_unit, other, weight, orient=False)
    loop = _placed_unit(unit, chord, twist, pitch, offset, z)
    names = set(left_keys) & set(right_keys)
    keys = {}
    for name in names:
        blended = (1.0 - weight) * np.asarray(left_keys[name], float) + weight * np.asarray(right_keys[name], float)
        placed = _placed_unit(blended.reshape(1, 2), chord, twist, pitch, offset, z)
        keys[name] = np.asarray(placed, float).reshape(-1, 3)[0]
    return loop, keys


def _section_projector(geo, web_names: set[str]):
    """``on_section(name, xyz, at_end)``: skin nodes and strip ends on the shell midline.

    A web node that is not an end of its strip stays the straight average, so
    the web is not pulled onto the outer mould. An end is the nearest station
    key at its span, so a short trailing-edge face cannot snap the seam onto
    one corner and then the other. Any other skin node is the closest point of
    the midline, and the same coordinates always project to the same point.

    A mesh station uses the cut already stored for it. A span between stations
    uses the chord-frame blend of those two station midlines, so a seam does
    not build a section.
    """
    slots: dict[float, tuple[np.ndarray | None, dict]] = {}
    prepared: dict[float, tuple] = {}
    station_rows: list[tuple[float, object]] | None = None

    def stations() -> list[tuple[float, object]]:
        nonlocal station_rows
        if station_rows is None:
            built = getattr(geo, "_mesh_z_cuts", {})
            station_rows = sorted(
                (float(cut.z), cut) for cut in built.values()
            )
        return station_rows

    def load(cut) -> tuple:
        key = round(float(cut.z), 9)
        found = prepared.get(key)
        if found is not None:
            return found
        graph = cut.shell_graph()
        loop = _skin_mid_loop(graph, float(cut.z))
        region_graph = getattr(geo, "region_graph", None)
        keys: dict = {}
        if region_graph is not None and len(loop) >= 3:
            keys = place_station_keys(cut, graph, region_graph, float(cut.z))
        unit, key_unit = _open_unit_loop(loop, keys, cut)
        found = (loop, keys, unit, key_unit)
        prepared[key] = found
        return found

    def slot(z: float) -> tuple[np.ndarray | None, dict]:
        key = round(float(z), 9)
        if key in slots:
            return slots[key]
        rows = stations()
        if not rows:
            slots[key] = (None, {})
            return slots[key]
        zs = [row[0] for row in rows]
        i = int(np.searchsorted(zs, float(z)))
        exact = None
        if i < len(zs) and abs(zs[i] - float(z)) <= 1e-6:
            exact = i
        elif i > 0 and abs(zs[i - 1] - float(z)) <= 1e-6:
            exact = i - 1
        if exact is not None:
            loop, keys, _unit, _key_unit = load(rows[exact][1])
            slots[key] = (loop, keys)
            return slots[key]
        if i <= 0 or i >= len(rows):
            _loop, _keys, unit, key_unit = load(rows[0][1] if i <= 0 else rows[-1][1])
            chord, twist, pitch, offset = _frame_at(geo, z)
            loop, keys = _blend_frames(unit, unit, key_unit, key_unit, 0.0, chord, twist, pitch, offset, float(z))
            slots[key] = (loop, keys)
            return slots[key]
        left_z, left = rows[i - 1]
        right_z, right = rows[i]
        _l_loop, _l_keys, left_unit, left_keys = load(left)
        _r_loop, _r_keys, right_unit, right_keys = load(right)
        chord, twist, pitch, offset = _frame_at(geo, z)
        loop, keys = _blend_frames(
            left_unit, right_unit, left_keys, right_keys,
            _bay_weight(left_z, right_z, float(z)),
            chord, twist, pitch, offset, float(z),
        )
        slots[key] = (loop, keys)
        return slots[key]

    def loop_at(z: float) -> np.ndarray | None:
        return slot(z)[0]

    def project(name: str, xyz: np.ndarray, at_end: bool) -> np.ndarray:
        xyz = np.asarray(xyz, float)
        if xyz.ndim == 2:
            return _project_many(name, xyz, bool(at_end))
        if name in web_names and not at_end:
            return xyz
        loop, keys = slot(float(xyz[2]))
        if at_end and keys:
            snapped = _nearest_key(keys, xyz)
            if snapped is not None:
                return snapped
        if loop is None or len(loop) < 2:
            return xyz
        return _closest_on_loop(loop, xyz)

    def _project_many(name: str, xyz: np.ndarray, at_end: bool) -> np.ndarray:
        if len(xyz) == 0 or (name in web_names and not at_end):
            return np.array(xyz, dtype=float, copy=True)
        out = np.array(xyz, dtype=float, copy=True)
        zkey = np.round(xyz[:, 2], 9)
        for key in dict.fromkeys(zkey.tolist()):
            chosen = np.flatnonzero(zkey == key)
            loop, keys = slot(float(xyz[int(chosen[0]), 2]))
            pending = list(chosen)
            if at_end and keys:
                pending = []
                for i in chosen:
                    snapped = _nearest_key(keys, xyz[int(i)])
                    if snapped is None:
                        pending.append(int(i))
                    else:
                        out[int(i)] = snapped
            if not pending or loop is None or len(loop) < 2:
                continue
            idx = np.asarray(pending, np.int64)
            out[idx] = _closest_on_loop_many(loop, xyz[idx])
        return out

    return project, loop_at


def _edge_sagitta(geo, nodes: np.ndarray, origins: Sequence[str], zs: Sequence[float],
                  elements, surface: Sequence[str], span=None) -> dict:
    """Chord sagitta of each station edge, and span sagitta of the LE and TE rulings.

    A station edge's sagitta is the distance from its chord to the shell midline.
    A spanwise ruling's sagitta is the bow of the leading or trailing edge.
    """
    sag: dict[tuple[int, int], float] = {}
    chord: dict[float, list[tuple[tuple[int, int], np.ndarray]]] = {}
    cuts = getattr(geo, "_mesh_z_cuts", {}) or {}
    seen: set[tuple[int, int]] = set()
    for element, surf in zip(elements, surface):
        if surf not in ("skin", "web"):
            continue
        for a, b in zip(element, element[1:] + element[:1]):
            if origins[a] != "station" or origins[b] != "station":
                continue
            if abs(float(nodes[a, 2]) - float(nodes[b, 2])) > 1e-9:
                continue
            key = (a, b) if a < b else (b, a)
            if key in seen:
                continue
            seen.add(key)
            if surf != "skin":
                sag[key] = 0.0
                continue
            mid = 0.5 * (nodes[a] + nodes[b])
            chord.setdefault(round(float(mid[2]), 9), []).append((key, mid))
    for zkey, group in chord.items():
        mids = np.stack([mid for _key, mid in group])
        cut = cuts.get(zkey)
        loop = None if cut is None else _skin_mid_loop(cut.shell_graph(), float(mids[0, 2]))
        on = mids if loop is None or len(loop) < 2 else _closest_on_loop_many(loop, mids)
        for (key, mid), hit in zip(group, on):
            sag[key] = float(np.linalg.norm(mid - hit))
    if span is None:
        le_p, te_p, le_s, te_s, _budget = _span_edges(geo, zs)
    else:
        le_p, te_p, le_s, te_s = span
    if len(zs) > 1:
        za = np.asarray(zs[:-1], float)
        zb = np.asarray(zs[1:], float)
        query_z = np.concatenate([za, zb, za, zb])
        query_xy = np.concatenate([le_p[:-1], le_p[1:], te_p[:-1], te_p[1:]])
        ids = _nearest_stations(nodes, origins, query_z, query_xy)
        n = len(za)
        pairs = ((ids[:n], ids[n:2 * n], le_s), (ids[2 * n:3 * n], ids[3 * n:], te_s))
        for left, right, value in pairs:
            for i in range(n):
                a, b = int(left[i]), int(right[i])
                if a < 0 or b < 0 or a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                sag[key] = float(value[i])
    return sag


def _nearest_stations(nodes: np.ndarray, origins: Sequence[str], zs, targets) -> np.ndarray:
    """Nearest station-node index for each query. ``-1`` when the station is empty.

    Queries that share a span station are one distance against that station's nodes.
    """
    nodes = np.asarray(nodes, float)
    zs = np.asarray(zs, float).reshape(-1)
    xy = np.asarray(targets, float).reshape(len(zs), -1)[:, :2]
    out = np.full(len(zs), -1, np.int64)
    if len(zs) == 0 or len(nodes) == 0:
        return out
    station = np.fromiter((origin == "station" for origin in origins), dtype=bool, count=len(origins))
    if not np.any(station):
        return out
    zkey = np.round(zs, 9)
    node_z = nodes[:, 2]
    for key in dict.fromkeys(zkey.tolist()):
        chosen = np.flatnonzero(zkey == key)
        z = float(zs[int(chosen[0])])
        idx = np.flatnonzero(station & (np.abs(node_z - z) <= 1e-7))
        if len(idx) == 0:
            continue
        delta = nodes[idx, :2][None, :, :] - xy[chosen, None, :]
        out[chosen] = idx[np.argmin((delta * delta).sum(axis=2), axis=1)]
    return out


def _nearest_at(nodes: np.ndarray, origins: Sequence[str], z: float, target: np.ndarray) -> int | None:
    found = int(_nearest_stations(nodes, origins, [z], np.asarray(target, float).reshape(1, -1))[0])
    return None if found < 0 else found


def _gauge_key(value, chord: float) -> tuple[float, ...]:
    from implicit_geometry_kernel.blade.loft.loft_cut import _gauge_values

    chord = float(chord) if float(chord) else 1.0
    return tuple(round(float(v) / chord, 6) for v in _gauge_values(value))


def _window_key(geo, z: float) -> tuple[float, ...]:
    fn = getattr(geo, "_cap_windows_fn", None)
    if fn is None:
        return ()
    raw = fn(float(z))
    if raw is None:
        return ()
    if isinstance(raw, Mapping):
        flat = []
        for side in sorted(raw):
            window = raw[side]
            if window is None:
                continue
            flat.extend(float(v) for v in window)
        return tuple(round(v, 6) for v in flat)
    flat = np.atleast_1d(np.asarray(raw, float)).ravel()
    return tuple(round(float(v), 6) for v in flat)


def _shape_keys(geo, zs: Sequence[float]) -> tuple:
    """One hashable key per station. Stations of one shape share a key.

    The key is the unit airfoil, the thickness-to-chord ratios, the web feet,
    the cap windows and the manufacturing scale. Twist, pitch and offset are
    placement, so they do not split a shape.
    """
    from implicit_geometry_kernel.blade.loft.loft_cut import _min_thickness_scale

    zs = tuple(float(z) for z in zs)
    if not zs:
        return ()
    _chord, _twist, _pitch, _offset, _x, upper, lower = geo.interpolate_many(zs)
    chord = np.asarray(_chord, float)
    surfaces = np.round(np.concatenate([np.asarray(upper, float), np.asarray(lower, float)], axis=1), 6)
    pairs = _web_fraction_pairs(geo, zs)
    keys = []
    for i, z in enumerate(zs):
        skin = geo._skin(z)
        cap = geo._cap(z) if getattr(geo, "_cap", None) is not None else 0.0
        web = geo._web(z) if getattr(geo, "_web", None) is not None else 0.0
        k = _min_thickness_scale(skin, cap, web)
        feet = tuple(round(float(v), 6) for v in np.asarray(pairs[i], float).ravel())
        keys.append((
            surfaces[i].tobytes(),
            _gauge_key(skin, chord[i]),
            _gauge_key(cap, chord[i]),
            _gauge_key(web, chord[i]),
            round(float(k), 6),
            feet,
            _window_key(geo, z),
        ))
    return tuple(keys)


def _sibling_cut(geo, rep, units: Sequence[np.ndarray], z: float):
    """The representative section placed at ``z`` by chord, twist and offset."""
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal

    chord, twist, pitch, offset, *_rest = geo.interpolate(float(z))
    signed = float(geo._placement_twist(twist))
    pitch = float(pitch)
    chord = float(chord)
    scale = chord / float(rep.chord) if float(rep.chord) else 1.0
    walls = tuple(
        replace(wall, points=place_chord_normal(unit, chord, signed, pitch, offset, z)[0],
                thickness=float(wall.thickness) * scale)
        for wall, unit in zip(rep.walls, units)
    )
    cells = []
    for cell in rep.cells:
        from implicit_geometry_kernel.blade.loft.span_place import to_chord_normal

        unit = to_chord_normal(np.asarray(cell.points, float)[:, :2], rep.chord, rep.twist,
                               rep.pitch_axis, rep.offset)
        cells.append(replace(
            cell,
            points=place_chord_normal(unit, chord, signed, pitch, offset, z)[0],
            area=float(cell.area) * scale * scale,
        ))
    cut = replace(
        rep,
        z=float(z),
        chord=chord,
        twist=signed,
        pitch_axis=pitch,
        offset=(float(offset[0]), float(offset[1])),
        walls=walls,
        cells=tuple(cells),
        transverse=geo._transverse_at(z),
        stringers=geo._stringers_at(z),
        tip_closure=geo._tip_closure_at(z),
        host=None,
    )
    reference = getattr(rep, "_wall_reference", None)
    if reference is not None:
        object.__setattr__(cut, "_wall_reference", reference)
    graph = getattr(rep, "_partition_graph", None)
    if graph is not None:
        object.__setattr__(cut, "_partition_graph", _place_graph(
            graph, rep, chord, signed, pitch, offset, z, scale))
    return cut


def _place_runs(runs, rep, chord, signed, pitch, offset, z, scale):
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal, to_chord_normal

    placed = []
    for run in runs:
        pts = np.asarray(run.points, float)
        if pts.ndim != 2 or len(pts) == 0:
            placed.append(run)
            continue
        unit = to_chord_normal(pts[:, :2], rep.chord, rep.twist, rep.pitch_axis, rep.offset)
        points = place_chord_normal(unit, chord, signed, pitch, offset, z)[0]
        placed.append(replace(
            run,
            points=points,
            thickness=float(run.thickness) * scale,
            eccentricity_m=float(run.eccentricity_m) * scale,
            layer_thicknesses_m=tuple(float(t) * scale for t in run.layer_thicknesses_m),
        ))
    return tuple(placed)


def _place_graph(graph, rep, chord, signed, pitch, offset, z, scale):
    return replace(
        graph,
        z=float(z),
        loop=_place_runs(graph.loop, rep, chord, signed, pitch, offset, z, scale),
        webs=_place_runs(graph.webs, rep, chord, signed, pitch, offset, z, scale),
        faces=(),
        branches=_place_runs(graph.branches, rep, chord, signed, pitch, offset, z, scale),
        solids=(),
    )


def _named_window(windows, prefix: str) -> dict[str, float]:
    """One chord fraction per window entry, named so the two bay ends can share it."""
    values = np.asarray(windows, float).ravel()
    out = {}
    for i, value in enumerate(values):
        name = f"{prefix}_{i}" if prefix else f"cap_{i}"
        out[name] = float(value)
    return out


def _live_fractions(geo, z: float) -> dict[str, float]:
    """Cap windows and web feet that are live at ``z``.

    A web outside its span window is absent. Feet come from the region-graph
    law, replaced by ``_web_positions_fn`` for a web that function still covers.
    """
    z = float(z)
    fracs: dict[str, float] = {}
    graph = getattr(geo, "region_graph", None)
    webs = tuple(getattr(graph, "webs", ()) or ()) if graph is not None else ()
    positions = getattr(geo, "_web_positions_fn", None)
    raw = None
    if positions is not None:
        raw = np.atleast_1d(np.asarray(positions(z), float)).ravel()
    width = len(raw) // len(webs) if webs and raw is not None and len(raw) % len(webs) == 0 else 0
    if webs and width:
        table = raw.reshape(len(webs), width)
        for i, web in enumerate(webs):
            if not web.active(z):
                continue
            name = str(web.name or f"web_{web.index}")
            for k, value in enumerate(table[i]):
                fracs[name if width == 1 else f"{name}_{k}"] = float(value)
    elif webs:
        for web in webs:
            if web.active(z):
                fracs[str(web.name or f"web_{web.index}")] = float(web.at(z))
    elif raw is not None:
        for i, frac in enumerate(raw):
            fracs[f"web_{i}"] = float(frac)
    caps = getattr(geo, "_cap_windows_fn", None)
    if caps is not None:
        windows = caps(z)
        if isinstance(windows, dict):
            for side, side_windows in windows.items():
                if side_windows is not None:
                    fracs.update(_named_window(side_windows, str(side)))
        else:
            fracs.update(_named_window(windows, ""))
    elif graph is not None:
        for name, value in graph.evaluate_anchors(z):
            if str(name).startswith("cap_"):
                fracs[str(name)] = float(value)
    return fracs


def _loft_gap(geo, za: float, zb: float) -> float:
    """How far, in metres, a straight loft misses the section inside ``(za, zb)``.

    The samples are the airfoil surface and the cap and web feet live at both
    ends of the bay. A positive result is the worst excess over the sagitta budget.
    """
    weights = np.array((0.25, 0.5, 0.75))
    query = float(za) + weights * (float(zb) - float(za))
    chord, _twist, _pitch, _offset, _x, upper, lower = geo.interpolate_many(
        (float(za), *query.tolist(), float(zb)))
    chord = np.asarray(chord, float)
    blend_u = (1.0 - weights)[:, None] * upper[0] + weights[:, None] * upper[-1]
    blend_l = (1.0 - weights)[:, None] * lower[0] + weights[:, None] * lower[-1]
    gaps = np.maximum(
        np.max(np.abs(blend_u - upper[1:-1]), axis=1),
        np.max(np.abs(blend_l - lower[1:-1]), axis=1),
    ) * chord[1:-1]
    left, right = _live_fractions(geo, float(za)), _live_fractions(geo, float(zb))
    shared = [name for name in left if name in right]
    if shared:
        for i, z in enumerate(query):
            live = _live_fractions(geo, float(z))
            errs = [abs(live[name] - ((1.0 - weights[i]) * left[name] + weights[i] * right[name]))
                    for name in shared if name in live]
            if errs:
                gaps[i] = max(float(gaps[i]), max(errs) * float(chord[i + 1]))
    return float(np.max(gaps - EPS_OVER_C * chord[1:-1]))


def _refine_shape_sites(geo, sites: Sequence[float]) -> tuple[float, ...]:
    """The section sample, built from ``sites``.

    Every site is kept. A bay is split while a straight loft of the skin and of
    the feet live at both ends bows past ``EPS_OVER_C`` of the local chord.
    """
    seeds = tuple(round(float(z), 9) for z in _unique_z(sites))
    cached = geo.__dict__.get("_mesh_shape_sites")
    if cached is not None and cached[0] == seeds:
        return cached[1]
    chosen = _unique_z(sites)
    for _round in range(8):
        extra = []
        for za, zb in zip(chosen, chosen[1:]):
            if float(zb) - float(za) <= 1e-6:
                continue
            if _loft_gap(geo, float(za), float(zb)) > 0.0:
                extra.append(0.5 * (float(za) + float(zb)))
        if not extra:
            break
        chosen = _unique_z(list(chosen) + extra)
    found = tuple(chosen)
    geo.__dict__["_mesh_shape_sites"] = (seeds, found)
    return found


def _section_sites(geo) -> tuple[float, ...]:
    """Section sample built from the anchors.

    Every anchor stays. A bay between two anchors is split only while the skin
    loft bows past the sagitta budget. Mesh stations between the finished
    samples are a placement of that loft.
    """
    return _refine_shape_sites(geo, _anchor_stations(geo))


def _cache_graph(cut):
    graph = getattr(cut, "_partition_graph", None)
    if graph is None:
        graph = cut.shell_graph()
        object.__setattr__(cut, "_partition_graph", graph)
    return graph


def _ensure_site_cuts(geo) -> list[tuple[float, object]]:
    """Execute one section at each site and keep its shell graph."""
    cached = geo.__dict__.get("_mesh_site_list")
    if cached is not None:
        return cached
    store = {}
    sites = []
    for z in _section_sites(geo):
        cut = section_cut_near(geo, float(z))
        if cut is None:
            continue
        _cache_graph(cut)
        key = round(float(cut.z), 9)
        if key in store:
            continue
        store[key] = cut
        sites.append((float(cut.z), cut))
    sites.sort(key=lambda item: item[0])
    geo.__dict__["_mesh_site_cuts"] = store
    geo.__dict__["_mesh_site_list"] = sites
    return sites


def _wall_units(cut) -> tuple:
    from implicit_geometry_kernel.blade.loft.span_place import to_chord_normal

    return tuple(
        to_chord_normal(np.asarray(wall.points, float)[:, :2], cut.chord, cut.twist,
                        cut.pitch_axis, cut.offset)
        for wall in cut.walls
    )


def _blend_open(left: np.ndarray, right: np.ndarray, weight: float, *, orient: bool = True) -> np.ndarray:
    """Arc-length blend of two open polylines.

    ``orient`` turns the right-hand polyline around when its first point is
    nearer the left-hand end than the left-hand start. A named region already
    runs from the branch start to the branch end, so that loft leaves the
    stored order alone.
    """
    return _blend_open_many(left, right, (float(weight),), orient=orient)[0]


def _blend_open_many(left: np.ndarray, right: np.ndarray, weights, *, orient: bool = True) -> np.ndarray:
    """``_blend_open`` for many weights. Returns ``(n, p, d)``."""
    from implicit_geometry_kernel.blade.loft.span_place import resample_open_many

    a = np.asarray(left, float)
    b = np.asarray(right, float)
    w = np.asarray(weights, float).reshape(-1)
    if orient and len(a) >= 2 and len(b) >= 2 and np.linalg.norm(a[0] - b[0]) > np.linalg.norm(a[0] - b[-1]):
        b = b[::-1]
    count = max(len(a), len(b), 2)
    if len(a) != count:
        a = resample_open_many(a[None, :, :], count)[0]
    if len(b) != count:
        b = resample_open_many(b[None, :, :], count)[0]
    if len(w) == 0:
        return np.zeros((0, count, a.shape[-1]), float)
    return (1.0 - w)[:, None, None] * a + w[:, None, None] * b


def _placed_unit(unit: np.ndarray, chord: float, twist: float, pitch: float, offset, z: float) -> np.ndarray:
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal

    return place_chord_normal(unit, chord, twist, pitch, offset, z)[0]


def _bay_weight(left_z: float, right_z: float, z: float) -> float:
    span = float(right_z) - float(left_z)
    if abs(span) <= 1e-15:
        return 0.0
    return float(np.clip((float(z) - float(left_z)) / span, 0.0, 1.0))


def _run_units(runs, cut) -> list[tuple[object, np.ndarray | None]]:
    """Chord-normal polyline of each run, or ``None`` when the run has no points."""
    from implicit_geometry_kernel.blade.loft.span_place import to_chord_normal

    out: list[tuple[object, np.ndarray | None]] = []
    for run in runs:
        pts = np.asarray(run.points, float)
        if pts.ndim != 2 or len(pts) == 0:
            out.append((run, None))
            continue
        out.append((run, to_chord_normal(
            pts[:, :2], cut.chord, cut.twist, cut.pitch_axis, cut.offset)))
    return out


def _lerp_runs(left_runs, right_runs, left, right, chord, twist, pitch, offset, z, weight):
    many = _lerp_runs_many(
        left_runs, right_runs, left, right,
        np.asarray([chord], float), np.asarray([twist], float), np.asarray([pitch], float),
        np.asarray([offset], float).reshape(1, 2), np.asarray([z], float), np.asarray([weight], float),
    )
    return None if many is None else many[0]


def _lerp_runs_many(left_runs, right_runs, left, right, chord, twist, pitch, offset, zs, weights):
    """Blend the site runs onto every station of a bay. ``None`` on a name mismatch."""
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal

    other = {run.name: run for run in right_runs}
    left_units = _run_units(left_runs, left)
    right_by_name = {run.name: unit for run, unit in _run_units(right_runs, right)}
    n = len(weights)
    if n == 0:
        return []
    columns: list[list] = [[] for _ in range(n)]
    for run, left_unit in left_units:
        match = other.get(run.name)
        if match is None:
            return None
        right_unit = right_by_name.get(run.name)
        if left_unit is None or right_unit is None:
            for column in columns:
                column.append(run)
            continue
        units = _blend_open_many(left_unit, right_unit, weights)
        points = place_chord_normal(units, chord, twist, pitch, offset, zs)
        left_layers = tuple(float(t) for t in run.layer_thicknesses_m)
        right_layers = tuple(float(t) for t in match.layer_thicknesses_m)
        t0, t1 = float(run.thickness), float(match.thickness)
        e0, e1 = float(run.eccentricity_m), float(match.eccentricity_m)
        for i, w in enumerate(weights):
            if len(left_layers) == len(right_layers):
                layers = tuple((1.0 - float(w)) * a + float(w) * b for a, b in zip(left_layers, right_layers))
            else:
                layers = left_layers
            columns[i].append(replace(
                run,
                points=points[i],
                thickness=(1.0 - float(w)) * t0 + float(w) * t1,
                eccentricity_m=(1.0 - float(w)) * e0 + float(w) * e1,
                layer_thicknesses_m=layers,
            ))
    if any(len(column) != len(right_runs) for column in columns):
        return None
    return [tuple(column) for column in columns]


def _blend_wall_points(left, right, wall, other, weight, chord, signed, pitch, offset, z):
    from implicit_geometry_kernel.blade.loft.span_place import to_chord_normal

    unit = _blend_open(
        to_chord_normal(np.asarray(wall.points, float)[:, :2], left.chord, left.twist,
                        left.pitch_axis, left.offset),
        to_chord_normal(np.asarray(other.points, float)[:, :2], right.chord, right.twist,
                        right.pitch_axis, right.offset),
        weight,
    )
    return _placed_unit(unit, chord, signed, pitch, offset, z)


def _lerp_cut(geo, left, right, z: float, weight: float, *, place_walls: bool = False):
    """Loft the two site sections onto ``z``. A topology change keeps the nearer section.

    Mesh stations read the partition runs. Wall polylines stay empty unless
    ``place_walls`` is set, and the thickness is still the blend of the two ends.
    """
    return _lerp_cuts_bay(geo, left, right, (float(z),), (float(weight),), place_walls=place_walls)[0]


def _sibling_each(geo, left, right, zs, weights):
    return [
        _sibling_cut(
            geo,
            left if float(w) < 0.5 else right,
            _wall_units(left if float(w) < 0.5 else right),
            float(z),
        )
        for z, w in zip(zs, weights)
    ]


def _lerp_cuts_bay(geo, left, right, zs, weights, *, place_walls: bool = False):
    """Loft every interior station of one site bay. Site runs convert once."""
    zs = np.asarray([float(z) for z in zs], float)
    weights = np.asarray([float(w) for w in weights], float)
    if len(zs) == 0:
        return []
    right_walls = {wall.label: wall for wall in right.walls}
    if len(left.walls) != len(right.walls) or any(wall.label not in right_walls for wall in left.walls):
        return _sibling_each(geo, left, right, zs, weights)
    chord, twist, pitch, offset, *_rest = geo.interpolate_many(zs)
    signed = np.array([geo._placement_twist(float(t)) for t in twist], float)
    pitch = np.asarray(pitch, float)
    chord = np.asarray(chord, float)
    offset = np.asarray(offset, float).reshape(-1, 2)
    left_graph = getattr(left, "_partition_graph", None)
    right_graph = getattr(right, "_partition_graph", None)
    if left_graph is None or right_graph is None:
        return _sibling_each(geo, left, right, zs, weights)
    loops = _lerp_runs_many(
        left_graph.loop, right_graph.loop, left, right, chord, signed, pitch, offset, zs, weights)
    webs = _lerp_runs_many(
        left_graph.webs, right_graph.webs, left, right, chord, signed, pitch, offset, zs, weights)
    branches = _lerp_runs_many(
        left_graph.branches, right_graph.branches, left, right, chord, signed, pitch, offset, zs, weights)
    if loops is None or webs is None or branches is None:
        return _sibling_each(geo, left, right, zs, weights)
    out = []
    for i, (z, weight) in enumerate(zip(zs, weights)):
        walls = []
        for wall in left.walls:
            other = right_walls[wall.label]
            thickness = (1.0 - float(weight)) * float(wall.thickness) + float(weight) * float(other.thickness)
            points = (
                _blend_wall_points(
                    left, right, wall, other, float(weight),
                    float(chord[i]), float(signed[i]), float(pitch[i]),
                    (float(offset[i, 0]), float(offset[i, 1])), float(z),
                )
                if place_walls else np.zeros((0, 3))
            )
            walls.append(replace(wall, points=points, thickness=thickness))
        cut = replace(
            left,
            z=float(z),
            chord=float(chord[i]),
            twist=float(signed[i]),
            pitch_axis=float(pitch[i]),
            offset=(float(offset[i, 0]), float(offset[i, 1])),
            walls=tuple(walls),
            cells=(),
            transverse=geo._transverse_at(float(z)),
            stringers=geo._stringers_at(float(z)),
            tip_closure=geo._tip_closure_at(float(z)),
            host=None,
            scale=(1.0 - float(weight)) * float(left.scale) + float(weight) * float(right.scale),
        )
        object.__setattr__(cut, "_partition_graph", replace(
            left_graph, z=float(z), loop=loops[i], webs=webs[i], faces=(),
            branches=branches[i], solids=()))
        if not place_walls:
            object.__setattr__(cut, "_wall_sources", (left, right, float(weight)))
        out.append(cut)
    return out


def _bay_of(sites: list[tuple[float, object]], z: float) -> tuple[float, float]:
    if not sites:
        key = round(float(z), 9)
        return (key, key)
    zs = [site[0] for site in sites]
    i = int(np.searchsorted(zs, z))
    if i < len(sites) and abs(zs[i] - z) <= 1e-6:
        key = round(zs[i], 9)
        return (key, key)
    if i > 0 and abs(zs[i - 1] - z) <= 1e-6:
        key = round(zs[i - 1], 9)
        return (key, key)
    if i <= 0:
        key = round(zs[0], 9)
        return (key, key)
    if i >= len(sites):
        key = round(zs[-1], 9)
        return (key, key)
    return (round(zs[i - 1], 9), round(zs[i], 9))


def _cut_at_z(geo, sites: list[tuple[float, object]], z: float):
    if not sites:
        cut = section_cut_near(geo, z)
        if cut is None:
            raise ValueError(f"no section at z = {z}")
        _cache_graph(cut)
        key = round(float(cut.z), 9)
        return cut, (key, key)
    zs = [site[0] for site in sites]
    i = int(np.searchsorted(zs, z))
    if i < len(sites) and abs(zs[i] - z) <= 1e-6:
        key = round(zs[i], 9)
        return sites[i][1], (key, key)
    if i > 0 and abs(zs[i - 1] - z) <= 1e-6:
        key = round(zs[i - 1], 9)
        return sites[i - 1][1], (key, key)
    if i <= 0:
        cut = _sibling_cut(geo, sites[0][1], _wall_units(sites[0][1]), z)
        key = round(zs[0], 9)
        return cut, (key, key)
    if i >= len(sites):
        cut = _sibling_cut(geo, sites[-1][1], _wall_units(sites[-1][1]), z)
        key = round(zs[-1], 9)
        return cut, (key, key)
    left_z, left = sites[i - 1]
    right_z, right = sites[i]
    weight = _bay_weight(left_z, right_z, z)
    return _lerp_cut(geo, left, right, z, weight), (round(left_z, 9), round(right_z, 9))


def _placed_cuts(geo, zs: Sequence[float]):
    """One executed section per defining station. Mesh stations loft between them.

    A station already built for a ``z`` is reused, including across a later
    count-gap pass. The key of a station is the pair of sites it lies between.
    Interior stations of one site bay blend the site runs once.
    """
    sites = _ensure_site_cuts(geo)
    built = geo.__dict__.setdefault("_mesh_z_cuts", {})
    site_zs = [site[0] for site in sites]
    cuts: list = [None] * len(zs)
    keys: list = [None] * len(zs)
    by_bay: dict[int, list[int]] = {}
    for index, z in enumerate(zs):
        z = float(z)
        cached = built.get(round(z, 9))
        if cached is not None:
            cuts[index] = cached
            keys[index] = _bay_of(sites, z)
            continue
        if not sites:
            cut, key = _cut_at_z(geo, sites, z)
            built[round(z, 9)] = cut
            cuts[index] = cut
            keys[index] = key
            continue
        i = int(np.searchsorted(site_zs, z))
        if i < len(sites) and abs(site_zs[i] - z) <= 1e-6:
            cuts[index] = sites[i][1]
            keys[index] = (round(site_zs[i], 9), round(site_zs[i], 9))
            built[round(z, 9)] = cuts[index]
            continue
        if i > 0 and abs(site_zs[i - 1] - z) <= 1e-6:
            cuts[index] = sites[i - 1][1]
            keys[index] = (round(site_zs[i - 1], 9), round(site_zs[i - 1], 9))
            built[round(z, 9)] = cuts[index]
            continue
        if i <= 0:
            cut = _sibling_cut(geo, sites[0][1], _wall_units(sites[0][1]), z)
            built[round(z, 9)] = cut
            cuts[index] = cut
            keys[index] = (round(site_zs[0], 9), round(site_zs[0], 9))
            continue
        if i >= len(sites):
            cut = _sibling_cut(geo, sites[-1][1], _wall_units(sites[-1][1]), z)
            built[round(z, 9)] = cut
            cuts[index] = cut
            keys[index] = (round(site_zs[-1], 9), round(site_zs[-1], 9))
            continue
        by_bay.setdefault(i - 1, []).append(index)
    for bay, indices in by_bay.items():
        left_z, left = sites[bay]
        right_z, right = sites[bay + 1]
        batch_z = [float(zs[i]) for i in indices]
        batch_w = [_bay_weight(left_z, right_z, z) for z in batch_z]
        for index, cut in zip(indices, _lerp_cuts_bay(geo, left, right, batch_z, batch_w)):
            built[round(float(zs[index]), 9)] = cut
            cuts[index] = cut
            keys[index] = (round(left_z, 9), round(right_z, 9))
    return cuts, tuple(keys)


def _branch_names(tmpl) -> tuple[str, ...]:
    return tuple(branch.name for branch in tmpl.loop) + tuple(web.name for web in tmpl.webs)


def _pin_station_skin_ends(smp: dict, tmpl) -> dict:
    """Share adjacent skin endpoints after place/lerp (BG-049)."""
    from implicit_geometry_kernel.blade.topology.map import _share_adjacent_skin_ends

    _share_adjacent_skin_ends(smp, tmpl)
    return smp


def _place_polylines(samples: Mapping[str, np.ndarray], rep, cut) -> dict[str, np.ndarray]:
    """Carry the representative's region polylines onto ``cut``'s placement."""
    from implicit_geometry_kernel.blade.loft.span_place import place_chord_normal, to_chord_normal

    out = {}
    for name, pts in samples.items():
        pts = np.asarray(pts, float)
        unit = to_chord_normal(pts[:, :2], rep.chord, rep.twist, rep.pitch_axis, rep.offset)
        out[name] = place_chord_normal(
            unit, cut.chord, cut.twist, cut.pitch_axis, cut.offset, cut.z)[0]
    return out


def _lerp_polylines(left_samples, right_samples, left, right, cut, weight: float) -> dict[str, np.ndarray]:
    """Blend two stations' region polylines in chord-normal coordinates and place them."""
    from implicit_geometry_kernel.blade.loft.span_place import to_chord_normal

    out = {}
    for name, left_pts in left_samples.items():
        right_pts = right_samples.get(name)
        left_pts = np.asarray(left_pts, float)
        if right_pts is None or len(left_pts) == 0:
            if len(left_pts):
                out[name] = _place_polylines({name: left_pts}, left, cut)[name]
            continue
        right_pts = np.asarray(right_pts, float)
        if len(right_pts) == 0:
            continue
        unit = _blend_open(
            to_chord_normal(left_pts[:, :2], left.chord, left.twist, left.pitch_axis, left.offset),
            to_chord_normal(right_pts[:, :2], right.chord, right.twist, right.pitch_axis, right.offset),
            weight,
            orient=False,
        )
        out[name] = _placed_unit(unit, cut.chord, cut.twist, cut.pitch_axis, cut.offset, cut.z)
    return out


def _station_wrap(cut, z: float, h: float | None, h_over_c: float | None, geo) -> bool:
    skin = _skin_mid_loop(cut.shell_graph(), float(z))
    if h is None and h_over_c is None:
        return te_face_open(cut, skin)
    return te_wrap_for_mesh(cut, skin, _target_edge(geo, float(z), h, h_over_c))


def _template_at(geo, cut, z: float, n_chord: int, h: float | None, h_over_c: float | None, donor=None):
    wrap = _station_wrap(cut, z, h, h_over_c, geo)
    donor_wrap = any(branch.wrap for branch in donor.loop) if donor is not None else None
    if donor is not None and donor_wrap == bool(wrap):
        return replace(donor, z=float(z))
    sources = getattr(cut, "_wall_sources", None)
    if sources is not None and (not cut.walls or len(np.asarray(cut.walls[0].points)) == 0):
        cut = _lerp_cut(geo, sources[0], sources[1], float(z), float(sources[2]), place_walls=True)
    return _template_from_cut(geo.region_graph, cut, z, n_chord, te_blunt=wrap)


def _resample_counts(smps: list, tmpls) -> None:
    """Stations that share a node count are resampled by one interpolation.

    A sized branch is always re-arced even when the vertex count already
    matches: the first sample may have bunched TE spacing, and the length/``h``
    count assumes equal edges so max chordwise length stays within the ceiling.
    """
    from implicit_geometry_kernel.blade.loft.span_place import resample_open_many

    groups: dict[tuple, list] = {}
    for i, (smp, tmpl) in enumerate(zip(smps, tmpls)):
        for branch in (*tmpl.loop, *tmpl.webs):
            pts = smp.get(branch.name)
            if pts is None:
                continue
            count = int(branch.n)
            sized = bool(getattr(branch, "sized", False))
            if len(pts) == count and not sized:
                continue
            bunch = (not sized) and str(branch.name).startswith("te_")
            groups.setdefault((count, bunch, len(pts), pts.shape[-1]), []).append((i, branch.name, pts))
    for (count, bunch, _length, _dim), items in groups.items():
        stack = np.stack([np.asarray(pts, float) for _i, _name, pts in items])
        resampled = resample_open_many(stack, count, bunch=bunch)
        for row, (i, name, _pts) in enumerate(items):
            smps[i][name] = resampled[row]


def _mark(name: str, started: float, **measures: int | float) -> None:
    """Record one step of the mesh. Wall time is since ``started``."""
    report_stage(
        "implicit_geometry_kernel", name, time.perf_counter() - started, None, started,
        measures=measures or None,
    )


def _sample_stations(geo, zs: Sequence[float], n_chord: int, h: float | None,
                     h_over_c: float | None, region_h: Mapping[str, float] | None,
                     log: list | None = None):
    with stage("implicit_geometry_kernel", "shell_stations") as rec:
        started = time.perf_counter()
        cuts, keys = _placed_cuts(geo, zs)
        _mark("shell_place", started, stations=len(zs))
        started = time.perf_counter()
        sites = {round(z, 9): cut for z, cut in _ensure_site_cuts(geo)}
        pack = geo.__dict__.setdefault("_mesh_site_pack", {})
        for key, cut in sites.items():
            if key in pack:
                continue
            tmpl = _template_at(geo, cut, float(cut.z), n_chord, h, h_over_c)
            pack[key] = (tmpl, sample_region_polylines(
                cut, cut.shell_graph(), tmpl, geo.region_graph))
        tmpls = []
        for z, cut, key in zip(zs, cuts, keys):
            donor = pack.get(key[0], (None, None))[0]
            tmpls.append(_template_at(geo, cut, float(z), n_chord, h, h_over_c, donor))
        _mark("shell_templates", started, stations=len(zs))
        started = time.perf_counter()
        smps: list = []
        for cut, key, tmpl in zip(cuts, keys, tmpls):
            left = pack.get(key[0])
            right = pack.get(key[1])
            if left is None:
                smps.append(sample_region_polylines(
                    cut, cut.shell_graph(), tmpl, geo.region_graph))
                continue
            left_tmpl, left_smp = left
            if key[0] == key[1] or right is None:
                smp = _place_polylines(left_smp, sites[key[0]], cut)
                smps.append(_pin_station_skin_ends(smp, tmpl))
                continue
            right_tmpl, right_smp = right
            names = _branch_names(tmpl)
            if _branch_names(left_tmpl) != names or _branch_names(right_tmpl) != names:
                smp = sample_region_polylines(
                    cut, cut.shell_graph(), tmpl, geo.region_graph)
                smps.append(_pin_station_skin_ends(smp, tmpl))
                continue
            weight = _bay_weight(float(sites[key[0]].z), float(sites[key[1]].z), float(cut.z))
            smp = _lerp_polylines(left_smp, right_smp, sites[key[0]], sites[key[1]], cut, weight)
            smps.append(_pin_station_skin_ends(smp, tmpl))
        _mark("shell_polylines", started, stations=len(zs))
        started = time.perf_counter()
        if h is not None or h_over_c is not None:
            h_at = [_target_edge(geo, z, h, h_over_c) for z in zs]
            chords = [_chord_at(geo, z) for z in zs]
            tmpls = region_counts(
                tmpls, smps, h_at, region_h=region_h, zs=zs, chords=chords, log=log)
            _resample_counts(smps, tmpls)
            _mark("shell_counts", started, stations=len(zs))
            started = time.perf_counter()
        for smp, tmpl in zip(smps, tmpls):
            _pin_station_skin_ends(smp, tmpl)
        from implicit_geometry_kernel.blade.topology.map_loft import align_event_polylines
        bands = align_event_polylines(smps, tmpls)
        geo.__dict__["_mesh_bands"] = bands
        _mark("shell_align", started, stations=len(zs))
        rec.set("stations", len(zs))
        return cuts, tmpls, smps


def _shell_mesh_body(geo, n_span: int | None = None, *, h: float | None = None,
                     h_over_c: float | None = None, n_chord: int = 8,
                     zs: Sequence[float] | None = None, polygon_size: float | None = None,
                     region_h: Mapping[str, float] | None = None,
                     refine: Sequence[tuple[float, float, float]] = ()) -> ShellMesh:
    explicit = zs is not None
    fraction = _size_fraction(n_span, h, h_over_c, explicit)
    if region_h and h is None and fraction is None:
        raise ValueError("region_h needs h or h_over_c")
    started = time.perf_counter()
    if explicit:
        zs = tuple(sorted(float(z) for z in zs))
    else:
        zs = mesh_stations(geo, n_span, h=h, h_over_c=fraction, refine=refine)
    inserted: dict[float, str] = {}
    if (h is not None or fraction is not None) and not explicit:
        zs, inserted = _freeze_sagitta(geo, zs, h, fraction)
    _mark("shell_span", started, stations=len(zs))
    count_log: list = []
    _cuts, tmpls, smps = _sample_stations(geo, zs, n_chord, h, fraction, region_h, count_log)
    refused: list[StationDecision] = []
    started = time.perf_counter()
    passes = 0
    if (h is not None or fraction is not None) and not explicit:
        for _pass in range(8):
            passes += 1
            pending: list[tuple[float, str]] = []
            for i in _shared_count_gap(tmpls):
                pending.append((0.5 * (float(zs[i]) + float(zs[i + 1])), "count_step"))
            present = {round(float(z), 9) for z in zs}
            sag_pending, sag_refused = _sagitta_pending(geo, zs, h, fraction, present)
            pending.extend(sag_pending)
            refused.extend(sag_refused)
            graded, changed, grade_refused = _grade_stations(geo, zs, tmpls, smps, h, fraction)
            refused.extend(grade_refused)
            fresh: list[float] = []
            present = {round(float(z), 9) for z in zs}
            for z, origin in pending:
                key = round(float(z), 9)
                if key in present or any(abs(float(z) - float(s)) <= 1e-6 for s in zs):
                    refused.append(StationDecision(
                        float(z), origin, True, "step would land on a station"))
                else:
                    fresh.append(float(z))
                    inserted[key] = origin
                    present.add(key)
            if changed:
                for z in graded:
                    key = round(float(z), 9)
                    if key not in present:
                        fresh.append(float(z))
                        inserted[key] = "grade"
                        present.add(key)
            if not fresh:
                break
            zs = _unique_z(list(zs) + fresh)
            count_log.clear()
            _cuts, tmpls, smps = _sample_stations(geo, zs, n_chord, h, fraction, region_h, count_log)
    _mark("shell_insert", started, passes=passes, stations=len(zs))
    cuts = _cuts
    station_log = tuple(
        StationDecision(float(z), _station_kind(geo, float(z), inserted, explicit), False, "")
        for z in zs
    ) + tuple(refused)
    # A tip station that only carries LE branches is too collapsed to loft onto:
    # fanning every dying panel onto those two arcs leaves open fan edges. Drop
    # it from the loft and close the tip on the last full station instead.
    tip_z = float(geo.tip_closure_z)
    tip_closure_at = tip_z
    if (
        len(zs) >= 2
        and abs(float(zs[-1]) - tip_z) <= 1e-8
        and {b.name for b in tmpls[-1].loop} <= {"le_su", "le_pr"}
    ):
        zs = tuple(zs[:-1])
        cuts = cuts[:-1]
        tmpls = tmpls[:-1]
        smps = smps[:-1]
        tip_closure_at = float(zs[-1])
        if hasattr(geo, "__dict__") and "_mesh_bands" in geo.__dict__:
            geo.__dict__["_mesh_bands"] = list(geo.__dict__["_mesh_bands"][:-1])
    started = time.perf_counter()
    graphs = [c.shell_graph() for c in cuts]
    _mark("shell_graphs", started, stations=len(zs))

    index: dict[tuple[float, float, float], int] = {}
    nodes: list[np.ndarray] = []
    origins: list[str] = []

    def add_node(p: np.ndarray, origin: str = "station") -> int:
        key = _point_key(p)
        k = index.get(key)
        if k is None:
            k = index[key] = len(nodes)
            nodes.append(np.array(key))
            origins.append(origin)
        elif origin == "station":
            origins[k] = "station"
        return k

    station_ids: list[dict[str, np.ndarray]] = []
    for tmpl, smp in zip(tmpls, smps):
        ids_at: dict[str, np.ndarray] = {}
        for branch in (*tmpl.loop, *tmpl.webs):
            pts = smp.get(branch.name)
            if pts is None:
                continue
            pts = np.asarray(pts, float)
            row = np.empty(len(pts), np.int64)
            for k, point in enumerate(pts):
                row[k] = add_node(point, "station")
            ids_at[branch.name] = row
        station_ids.append(ids_at)

    elements: list[tuple[int, ...]] = []
    rows: list[tuple[str, str, str, float, float, tuple[str, ...], tuple[float, ...]]] = []
    meta: list[tuple[int, float, float, str, int, int]] = []

    def add(ids: list[int], row, record) -> None:
        if len(ids) > 1 and ids[0] == ids[-1]:
            ids = ids[:-1]
        # Consecutive weld of loft corners (collapsed edge) can leave three
        # distinct indices. Keep an all-quad sheet by repeating the last
        # corner; ``edge_use`` ignores the zero-length self-edge.
        if row[1] not in ("rib", "tip_closure"):
            compact: list[int] = []
            for k in ids:
                if not compact or compact[-1] != k:
                    compact.append(k)
            if len(compact) == 3:
                ids = [compact[0], compact[1], compact[2], compact[2]]
            elif len(compact) < 3:
                return
            else:
                ids = compact
        least = 3 if row[1] in ("rib", "tip_closure") else 4
        if len(ids) < least or (row[1] in ("rib", "tip_closure") and len(ids) != 3):
            raise ValueError(
                f"element of {row[0]!r} has {len(ids)} nodes; "
                "skin and webs are quads, ribs and the tip are triangles"
            )
        elements.append(tuple(ids))
        rows.append(row)
        meta.append(record)

    def ids_of(face: np.ndarray, face_origins: Sequence[str]) -> list[int]:
        out: list[int] = []
        for p, origin in zip(np.asarray(face, float), face_origins):
            k = add_node(p, origin)
            if not out or out[-1] != k:
                out.append(k)
        return out

    web_names = {w.name for t in tmpls for w in t.webs}
    project, loop_at = _section_projector(geo, web_names)
    strips: list[StripDecision] = []
    from implicit_geometry_kernel.blade.topology.map_loft import (
        _ruled_quad_array,
        _ruled_quad_masks,
        loft_mapped_strips,
    )

    def weld(pts: np.ndarray) -> np.ndarray:
        out = np.empty(len(pts), np.int64)
        for k, point in enumerate(np.asarray(pts, float)):
            out[k] = add_node(point, "station")
        return out

    def ruled_ids(station: int, name: str, pts: np.ndarray) -> np.ndarray:
        """Node ids for a ruled strip end.

        The station cache is keyed by region name, but an event bay's parent
        walk can reuse another region's polyline under a new name (or the
        reverse). Only reuse the cache when every point already maps to the
        cached id; otherwise weld from ``pts`` so neighbouring bays share
        the same station nodes.
        """
        pts = np.asarray(pts, float)
        known = station_ids[station].get(name)
        if known is not None and len(known) == len(pts):
            if all(index.get(_point_key(p)) == int(known[k]) for k, p in enumerate(pts)):
                return known
        return weld(pts)

    def record_strips(i: int, bay_log: list[dict]) -> None:
        c0 = _chord_at(geo, float(zs[i]))
        c1 = _chord_at(geo, float(zs[i + 1]))
        den = max(c0, c1)
        chord_ratio = abs(c1 - c0) / den if den > 0.0 else 0.0
        z0, z1 = float(zs[i]), float(zs[i + 1])
        for row in bay_log:
            strips.append(StripDecision(
                bay=i, z0=z0, z1=z1,
                region=row["region"], n0=int(row["n0"]), n1=int(row["n1"]),
                constructor=row["constructor"], trigger=row["trigger"],
                adv_a=int(row["adv_a"]), adv_b=int(row["adv_b"]),
                ties=int(row["ties"]), note=row["note"],
                chord_ratio=chord_ratio,
                site0=int(row.get("site0", -1)), site1=int(row.get("site1", -1)),
            ))

    def emit_faces(i: int, faces, by_region: dict) -> None:
        pending_faces = []
        z0, z1 = float(zs[i]), float(zs[i + 1])
        for face in faces:
            info = by_region.get(face.name, {})
            xyz = np.asarray(face.face, float)
            zc = float(np.mean(xyz[:, 2]))
            j = i if abs(zc - z0) <= abs(z1 - zc) else i + 1
            on = xyz[np.abs(xyz[:, 2] - zs[j]) <= 1e-12]
            p = (on if len(on) else xyz)[:, :2].mean(axis=0)
            surface = "web" if face.name in web_names else "skin"
            pending_faces.append((j, surface, p, face, info, xyz))
        buckets: dict[tuple[int, str], list] = {}
        for item in pending_faces:
            buckets.setdefault((item[0], item[1]), []).append(item)
        for (j, surface), items in buckets.items():
            runs = graphs[j].webs if surface == "web" else graphs[j].loop
            chosen = _nearest_runs(runs, np.stack([item[2] for item in items]))
            for run, (_j, _surface, _p, face, info, xyz) in zip(chosen, items):
                if run is None:
                    raise ValueError(f"no {surface} partition at z = {zs[j]} for {face.name!r}")
                add(ids_of(xyz, face.origins),
                    (face.name, surface, run.name, float(run.thickness), float(run.eccentricity_m),
                     tuple(run.layers), tuple(float(t) for t in run.layer_thicknesses_m)),
                    (i, z0, z1, face.constructor, int(info.get("n0", 0)), int(info.get("n1", 0))))

    def emit_ruled(i: int, pairs) -> None:
        z0, z1 = float(zs[i]), float(zs[i + 1])
        bay_log: list[dict] = []
        blocks = []
        for name, pa, pb in pairs:
            pa = np.asarray(pa, float)
            pb = np.asarray(pb, float)
            note = ""
            if len(pa) >= 2 and len(pa) == len(pb):
                quads = _ruled_quad_array(pa, pb)
                keep, bow = _ruled_quad_masks(quads)
                if np.any(bow):
                    note = "bowtie"
                if np.any(keep):
                    blocks.append((name, pa, pb, quads, keep, bow))
            bay_log.append({
                "region": name, "n0": int(len(pa)), "n1": int(len(pb)),
                "constructor": "ruled", "trigger": "",
                "adv_a": 0, "adv_b": 0, "ties": 0, "note": note,
                "site0": -1, "site1": -1,
            })
        names1 = {w.name for w in tmpls[i + 1].webs}
        names0 = {w.name for w in tmpls[i].webs}
        for web in tmpls[i].webs:
            if web.name not in smps[i + 1] or web.name not in names1:
                bay_log.append({
                    "region": web.name, "n0": int(len(smps[i].get(web.name, ()))), "n1": 0,
                    "constructor": "omitted", "trigger": "",
                    "adv_a": 0, "adv_b": 0, "ties": 0, "note": "web omitted",
                    "site0": -1, "site1": -1,
                })
        for web in tmpls[i + 1].webs:
            if web.name not in smps[i] or web.name not in names0:
                bay_log.append({
                    "region": web.name, "n0": 0, "n1": int(len(smps[i + 1].get(web.name, ()))),
                    "constructor": "omitted", "trigger": "",
                    "adv_a": 0, "adv_b": 0, "ties": 0, "note": "web omitted",
                    "site0": -1, "site1": -1,
                })
        record_strips(i, bay_log)
        if not blocks:
            return
        queries = []
        for block_i, (name, _pa, _pb, quads, keep, _bow) in enumerate(blocks):
            xy, tip = _ruled_station_queries(quads, z0, z1)
            surface = "web" if name in web_names else "skin"
            for k in np.flatnonzero(keep):
                queries.append((int(tip[k]), surface, xy[k], name, int(k), block_i))
        grouped: dict[tuple[int, str], list] = {}
        for item in queries:
            grouped.setdefault((item[0], item[1]), []).append(item)
        run_of: dict[tuple[int, int], object] = {}
        for (tip, surface), items in grouped.items():
            station = i + 1 if tip else i
            runs = graphs[station].webs if surface == "web" else graphs[station].loop
            chosen = _nearest_runs(runs, np.stack([item[2] for item in items]))
            for run, item in zip(chosen, items):
                if run is None:
                    raise ValueError(f"no {surface} partition at z = {zs[station]} for {item[3]!r}")
                run_of[(item[5], item[4])] = run
        for block_i, (name, pa, pb, _quads, keep, bow) in enumerate(blocks):
            left, right = ruled_ids(i, name, pa), ruled_ids(i + 1, name, pb)
            quad_ids = np.stack([left[:-1], left[1:], right[1:], right[:-1]], axis=1)
            n0, n1 = int(len(pa)), int(len(pb))
            for k in np.flatnonzero(keep):
                run = run_of[(block_i, int(k))]
                surface = "web" if name in web_names else "skin"
                row = (name, surface, run.name, float(run.thickness), float(run.eccentricity_m),
                       tuple(run.layers), tuple(float(t) for t in run.layer_thicknesses_m))
                ids = [int(v) for v in quad_ids[k]]
                if bow[k]:
                    add(ids[:3], row, (i, z0, z1, "bowtie", n0, n1))
                    add([ids[0], ids[2], ids[3]], row, (i, z0, z1, "bowtie", n0, n1))
                else:
                    add(ids, row, (i, z0, z1, "ruled", n0, n1))

    with stage("implicit_geometry_kernel", "shell_strips"):
        bands = geo.__dict__.get("_mesh_bands") or ()
        for i in range(len(zs) - 1):
            pairs = bands[i] if i < len(bands) else None
            event = pairs is not None and any(len(pa) != len(pb) for _name, pa, pb in pairs)
            plain = pairs is not None and not event and _plain_ruled(pairs)
            if plain:
                emit_ruled(i, pairs)
                continue
            bay_log: list[dict] = []
            faces = loft_mapped_strips(
                smps[i], smps[i + 1], tmpls[i], tmpls[i + 1],
                on_section=project, strip_log=bay_log, pairs=pairs,
            )
            record_strips(i, bay_log)
            emit_faces(i, faces, {row["region"]: row for row in bay_log})

    n_source = len(nodes)
    started = time.perf_counter()
    for si, (z, tmpl, smp) in enumerate(zip(zs, tmpls, smps)):
        members: list = [m for m in geo.transverse if abs(float(m.z) - z) <= 1e-8]
        rings: list[tuple[str, str, np.ndarray, object]] = []
        for m in members:
            for name, ring in rib_cell_polygons_from_samples(tmpl, smp, m):
                rings.append((name, "rib", np.asarray(ring, float), m))
        if abs(z - tip_closure_at) <= 1e-8 and getattr(geo, "tip_closure", None) is not None:
            tip = geo.tip_closure
            for name, ring in rib_cell_polygons_from_samples(tmpl, smp, str(tip.name)):
                rings.append((name, "tip_closure", np.asarray(ring, float), tip))
        station: list[tuple[tuple, np.ndarray]] = []
        for name, surface, ring, member in rings:
            for p in ring:
                key = _point_key(p)
                k = index.get(key)
                if k is None or k >= n_source:
                    raise ValueError(f"partition ring vertex at z = {z} is not a skin or web node")
            if surface == "rib":
                row = (name, surface, str(member.name), float(member.thickness), 0.0,
                       (str(member.name),), (float(member.thickness),))
            else:
                t_tip = float(geo.tip_closure_thickness())
                row = (name, surface, str(member.partition), t_tip, 0.0,
                       (str(member.partition),), (t_tip,))
            faces = triangulate_ring(
                ring, size=polygon_size, pattern="unstructured", split_boundary=False,
            )
            if len(faces) == 0:
                continue
            station.append((row, np.asarray(faces, float)))
        if not station:
            continue
        welded = _weld_triangles(
            np.concatenate([faces for _row, faces in station], axis=0),
            index, nodes, origins,
        )
        offset = 0
        for row, faces in station:
            block = welded[offset:offset + len(faces)]
            offset += len(faces)
            for tri in block:
                if int(tri[0]) < 0:
                    continue
                ids = [int(k) for k in tri]
                kind = "interior" if any(k >= n_source for k in ids) else "inherited"
                add(ids, row, (si, float(z), float(z), kind, 0, 0))
    _mark("shell_ribs", started, stations=len(zs))

    with stage("implicit_geometry_kernel", "shell_quality"):
        started = time.perf_counter()
        _split_crossed_triangles(elements, rows, meta, nodes)
        _mark("shell_split", started)
        started = time.perf_counter()
        pts = np.asarray(nodes, float).reshape(-1, 3)
        te_at = {float(z): np.asarray(smp["te_su"][-1], float) for z, smp in zip(zs, smps) if len(smp.get("te_su", ()))}
        elements, flip_flags = _oriented(pts, elements, [r[1] for r in rows], zs, te_at)
        _mark("shell_orient", started)
        started = time.perf_counter()
        on_skin = {k for e, s in zip(elements, (r[1] for r in rows)) if s == "skin" for k in e}
        distance = np.zeros(len(pts))
        by_station: dict[float, list[int]] = {}
        for k, origin in enumerate(origins):
            if origin == "station" or k not in on_skin:
                continue
            by_station.setdefault(round(float(pts[k, 2]), 9), []).append(k)
        for ks in by_station.values():
            idx = np.asarray(ks, np.int64)
            loop = loop_at(float(pts[idx[0], 2]))
            if loop is None or len(loop) < 2:
                continue
            distance[idx] = np.linalg.norm(pts[idx] - _closest_on_loop_many(loop, pts[idx]), axis=1)
        _mark("shell_distance", started)
        started = time.perf_counter()
        # Leading and trailing rulings come from the airfoil stack. The offset
        # skin stays out of the mesh until it matches the section midline.
        le, te, le_sag, te_sag, _budget = _span_edges(geo, zs)
        edge_sagitta = _edge_sagitta(
            geo, pts, origins, zs, elements, [r[1] for r in rows],
            span=(le, te, le_sag, te_sag),
        )
        _mark("shell_sagitta", started)
        started = time.perf_counter()
        partition = [r[2] for r in rows]
        thickness = np.array([r[3] for r in rows], float)
        table, jmax = _measure(pts, elements, thickness, with_ratio=False)
        _mark("shell_measure", started)
    started = time.perf_counter()
    edge = tuple(_target_edge(geo, z, h, fraction) for z in zs) if (h is not None or fraction is not None) else ()
    recipe = MeshRecipe(
        tip_h_over_c=fraction,
        h=None if h is None else float(h),
        h_over_c=None if h_over_c is None else float(h_over_c),
        n_span=None if n_span is None else int(n_span),
        region_h=tuple(sorted((str(k), float(v)) for k, v in (region_h or {}).items())),
        refine=tuple((float(a), float(b), float(dz)) for a, b, dz in refine),
        edge=edge,
    )
    mesh = ShellMesh(
        nodes=pts,
        elements=elements,
        name=[r[0] for r in rows],
        surface=[r[1] for r in rows],
        partition=partition,
        thickness=thickness,
        eccentricity=np.array([r[4] for r in rows], float),
        layers=[r[5] for r in rows],
        layer_thicknesses=[r[6] for r in rows],
        stations=tuple(float(z) for z in zs),
        quality=(),
        recipe=recipe,
        element_quality=table,
        n_source=n_source,
        bay=[m[0] for m in meta],
        z0=[m[1] for m in meta],
        z1=[m[2] for m in meta],
        constructor=[m[3] for m in meta],
        n0=[m[4] for m in meta],
        n1=[m[5] for m in meta],
        flipped=flip_flags,
        node_origin=origins,
        section_distance=distance,
        counts=count_log,
        station_log=station_log,
        strips=strips,
        edge_sagitta=edge_sagitta,
    )
    facts = assemble(mesh, float(zs[0]))
    table["volume_ratio"] = volume_ratio_from_edges(facts["edges"], table["shell_volume"])
    _mark("shell_assemble", started, elements=len(elements))
    return replace(
        mesh, element_quality=table, quality=quality_rollup(partition, table, jmax), **facts,
    )


__all__ = [
    "SURFACES",
    "TIP_H_OVER_C",
    "MeshRecipe",
    "MeshRunQuality",
    "ShellMesh",
    "mesh_stations",
    "orient_shell_elements",
    "shell_mesh",
]
