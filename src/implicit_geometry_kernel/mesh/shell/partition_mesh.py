"""The (t, z) partition mesher: a structured quad grid on every partition (decision 0015).

:func:`partition_mesh` meshes each skin and web partition of the shell graph as
its own (t, z) patch: ``t`` across the partition by normalised arc length
(:func:`~implicit_geometry_kernel.blade.skeleton.run_t`), ``z`` along the span
at the mesh stations. Each patch is a structured grid of quads, so element
edges follow the partition's own coordinates.

Nodes are shared by skeleton identity, not welded by coordinate: a partition
boundary is one node per station, used by both partitions, and a web foot is
the skin node at that boundary or vertex. Mesh lines run through every
skeleton point on a partition (web feet) and through every sub-partition band
edge (decision 0016), so a thickness or eccentricity step lies on an element
edge.

Nodes lie on the section mid-line at each station (the cut's run polylines).
The skin loop is meshed as one closed ring from the partition boundary nearest
the trailing edge; it must be closed and have no branches (stringers).
Terminations ride persistent mesh lines (decision 0018): every web foot, every
partition boundary and every band edge seen anywhere along the span is a line
round the loop at every station, held between its neighbouring lines where its
feature is absent. So the loop's lines never change in number: a web that ends
is a plain all-quad strip that stops on its end station, and caps or skin bays
that appear, vanish or renumber move between existing lines. Boundaries are
followed by partition kind and position, not by name. Each skin quad takes the
partition its middle lies in at the lower station of its bay. Ribs and the tip
closure are triangle fills on the station's skin and web nodes. Branches and
Newton junction nodes between stations are later steps of plan row L6.
"""

from __future__ import annotations

import math
import re
from dataclasses import replace
from typing import Any, Sequence

import numpy as np

from implicit_geometry_kernel.blade.skeleton import run_point_at, run_t
from implicit_geometry_kernel.blade.topology.sizing import EPS_OVER_C, KAPPA_FLOOR
from implicit_geometry_kernel.mesh.fill import triangulate_ring
from implicit_geometry_kernel.mesh.shell.orient import _oriented
from implicit_geometry_kernel.mesh.shell.spanwise import _le_te_points
from implicit_geometry_kernel.mesh.quality import _measure, quality_rollup, volume_ratio_from_edges
from implicit_geometry_kernel.mesh.tables import assemble
from implicit_geometry_kernel.mesh.shell.types import MeshRecipe, ShellMesh
from implicit_geometry_kernel.mesh.shell.stations import mesh_stations

_T_TOL = 1e-9


class PartitionMeshError(ValueError):
    """The blade is outside what the partition mesher handles yet."""


def _graph_signature(graph) -> tuple:
    chain = graph.closed_loop_chain()
    if chain is None:
        raise PartitionMeshError(
            f"station z={float(graph.z):g}: the skin loop is not closed; the partition mesher needs a closed loop"
        )
    if graph.branches:
        raise PartitionMeshError(
            f"station z={float(graph.z):g}: branches ({', '.join(b.name for b in graph.branches)}) are not meshed yet"
        )
    return tuple(name for name, _, _ in chain), tuple(sorted(w.name for w in graph.webs))


def _flipped(graph) -> frozenset:
    """Loop partitions whose points run against the loop direction at this station."""
    chain = graph.closed_loop_chain() or ()
    return frozenset(name for name, _, rev in chain if rev)


def _along(t: float, flip: bool) -> float:
    """A run's own ``t`` as ``t`` along the loop direction."""
    return 1.0 - float(t) if flip else float(t)


def _closest_t(points: np.ndarray, xy: np.ndarray) -> tuple[float, float]:
    """``(t, distance)`` of the point on a run polyline closest to ``xy``."""
    pts = np.asarray(points, float)[:, :2]
    tt = run_t(pts)
    a, b = pts[:-1], pts[1:]
    d = b - a
    dd = np.einsum("ij,ij->i", d, d)
    u = np.where(dd > 0.0, np.einsum("ij,ij->i", np.asarray(xy, float)[:2] - a, d) / np.where(dd > 0.0, dd, 1.0), 0.0)
    u = np.clip(u, 0.0, 1.0)
    q = a + u[:, None] * d
    dist = np.linalg.norm(q - np.asarray(xy, float)[:2], axis=1)
    k = int(np.argmin(dist))
    return float(tt[k] + u[k] * (tt[k + 1] - tt[k])), float(dist[k])


def _web_feet(skeleton, web: str, loop_names: Sequence[str], flip: frozenset) -> dict:
    """``{end: (host, t)}`` of a web's feet on the loop, ``t`` along the loop direction."""
    out: dict = {}
    by_key = skeleton.by_key()
    for end in (0, 1):
        node = by_key.get(("junction", web, end))
        if node is None:
            continue
        for host in loop_names:
            t = node.runs.get(host)
            if t is not None:
                out[end] = (host, _along(t, host in flip))
                break
    return out


def _breakpoints(skeleton, run, bands: Sequence[Any], flip: bool = False) -> dict:
    """``{key: t}`` of interior points on ``run`` that a mesh line must pass through.

    ``t`` runs along the loop direction: a run stored against it (``flip``) has
    its own ``t`` reversed, so the mesh walks every loop partition the same way.
    Web junctions are left out: :func:`partition_mesh` places those itself, so
    a web that starts or stops inside the span keeps its mesh line.
    """
    out: dict = {}
    for node in skeleton.nodes:
        if node.key[0] == "junction":
            continue
        t = node.runs.get(run.name)
        if t is not None and _T_TOL < t < 1.0 - _T_TOL:
            out[node.key] = _along(t, flip)
    for s in bands:
        if s.host != run.name:
            continue
        for k, t in enumerate(s.t_range):
            if _T_TOL < t < 1.0 - _T_TOL:
                out[("band", s.name, k)] = _along(t, flip)
    return out


def _segment_breaks(t_by_key: dict) -> list:
    """Sorted ``[(key, t)]`` with the run ends; coincident breakpoints merged."""
    items = sorted(t_by_key.items(), key=lambda kv: kv[1])
    merged: list = [(("end", 0), 0.0)]
    for key, t in items:
        if abs(t - merged[-1][1]) <= _T_TOL:
            continue
        merged.append((key, t))
    if abs(1.0 - merged[-1][1]) <= _T_TOL:
        merged.pop()
    merged.append((("end", 1), 1.0))
    return merged


def _cells(ids: dict, j: int, web_names: Sequence[str]) -> list[list[int]]:
    """Node-id rings of the cells the webs cut the skin loop into at station ``j``."""
    cells = [list(ids[(j, "loop")][:-1])]
    for w in web_names:
        row = ids[(j, w)]
        a, b = row[0], row[-1]
        inner = row[1:-1]
        for c, cell in enumerate(cells):
            if a in cell and b in cell:
                ia, ib = cell.index(a), cell.index(b)
                n = len(cell)
                path_ab = [cell[(ia + k) % n] for k in range((ib - ia) % n + 1)]
                path_ba = [cell[(ib + k) % n] for k in range((ia - ib) % n + 1)]
                cells[c:c + 1] = [path_ab + inner[::-1], path_ba + inner]
                break
        else:
            raise PartitionMeshError(f"web {w!r} does not join two nodes of one cell at station {j}")
    return cells


#: Two mesh lines round the loop closer than this (in normalised loop arc
#: length) are one point: they share a node.
_U_TOL = 1e-6


class _StationLoop:
    """The skin loop at one station as one closed ring along the loop direction.

    ``points`` run from the start of the first partition round to it again;
    ``u`` is their normalised arc length and ``w`` the graded edge measure,
    each partition's own measure laid end to end. ``boundaries`` holds the
    ``u`` of each partition boundary but the ring start.
    """

    def __init__(self, graph, flip: frozenset, loop_names: Sequence[str], h: float, chord: float):
        by_name = {r.name: r for r in graph.loop}
        self.runs = [by_name[n] for n in loop_names]
        pts_all: list[np.ndarray] = []
        u_all: list[np.ndarray] = []
        w_all: list[np.ndarray] = []
        lengths = []
        oriented = []
        for r in self.runs:
            p = np.asarray(r.points, float)
            p = p[::-1].copy() if r.name in flip else p
            oriented.append(p)
            lengths.append(float(np.sum(np.linalg.norm(np.diff(p[:, :2], axis=0), axis=1))) if len(p) > 1 else 0.0)
        self.length = float(sum(lengths)) or 1.0
        self.s0 = np.concatenate([[0.0], np.cumsum(lengths)[:-1]]) / self.length
        self.s1 = np.cumsum(lengths) / self.length
        w_off = 0.0
        for i, (p, n) in enumerate(zip(oriented, lengths)):
            tt, w = _edge_measure(p, h, chord)
            if len(tt) != len(p):  # a run of zero length
                continue
            keep = slice(1, None) if pts_all else slice(None)
            pts_all.append(p[keep])
            u_all.append(self.s0[i] + tt[keep] * n / self.length)
            w_all.append(w[keep] + w_off)
            w_off += float(w[-1])
        self.points = np.vstack(pts_all)
        self.u = np.concatenate(u_all)
        self.w = np.concatenate(w_all)
        self.index = {r.name: i for i, r in enumerate(self.runs)}
        names = list(loop_names)
        #: ``[(kind pair, u)]``: each boundary but the ring start, with the kinds
        #: of the partitions before and after it (see :func:`_kind`).
        self.boundaries = [((_kind(names[i - 1]), _kind(names[i])), float(self.s0[i])) for i in range(1, len(names))]

    def u_of(self, host: str, t: float) -> float:
        """Ring ``u`` of the point at ``t`` (along the loop direction) on partition ``host``."""
        i = self.index[host]
        return float(self.s0[i] + float(t) * (self.s1[i] - self.s0[i]))

    def run_at(self, u: float):
        """The loop partition holding ring parameter ``u``."""
        i = int(np.searchsorted(self.s1, float(u), side="right"))
        return self.runs[min(i, len(self.runs) - 1)]


def _kind(name: str) -> str:
    """A partition's kind without its place: ``skin_suction_surface:bay_3`` and
    ``skin_suction_surface:trailing_edge`` are both ``skin_suction_surface``, and
    ``spar_cap_suction_surface_2`` is ``spar_cap_suction_surface``. Partition
    names renumber as webs and caps come and go; kinds do not."""
    return re.sub(r"_\d+$", "", name.split(":")[0])


#: How far, in normalised loop arc length, a partition boundary may move from
#: one station to the next and still be the same mesh line.
_TRACK_U = 0.03


def _hold_absent(real: list[dict]) -> list[dict]:
    """Every key at every station: where a key is absent it is held by :func:`_held_u`
    from the nearest station that has it."""
    n = len(real)
    lines = [dict(at) for at in real]
    for key in sorted({k for at in real for k in at}, key=repr):
        have = [j for j in range(n) if key in real[j]]
        for j in range(n):
            if key not in real[j]:
                r = min(have, key=lambda i: abs(i - j))
                lines[j][key] = _held_u(real[r], real[j], key)
    return lines


def _slot(foot_lines: dict, u: float):
    """The web-foot line just before ``u`` round the ring (``None`` before the first)."""
    before = [(v, k) for k, v in foot_lines.items() if v < u - _U_TOL]
    return max(before)[1] if before else None


def _track_boundaries(loops: Sequence[Any], feet: dict, foot_lines: list[dict]) -> list[dict]:
    """``{key: u}`` per station of the partition boundaries no web foot sits on.

    Boundaries are followed along the span by kind pair, by the web-foot line
    they follow, and by position, not by partition name: each takes the
    nearest line of the same kind pair and foot seen at an earlier station,
    within :data:`_TRACK_U`, or starts a new line. A cap edge that crosses a
    web foot therefore becomes a new line on the far side, and the old one is
    held on its own side, so lines never cross.
    """
    out: list[dict] = []
    last: dict = {}  # key -> u where last seen
    counter: dict = {}
    for j, lp in enumerate(loops):
        foot_u = [u for (jj, _w, _e), u in feet.items() if jj == j]
        found = [((pair, _slot(foot_lines[j], u)), u) for pair, u in lp.boundaries
                 if not any(min(abs(u - f), 1.0 - abs(u - f)) <= _U_TOL for f in foot_u)]
        cands = sorted(
            (abs(u - lu), i, key) for i, (pair, u) in enumerate(found) for key, lu in last.items()
            if key[1] == pair and abs(u - lu) <= _TRACK_U
        )
        taken: dict = {}
        used: set = set()
        for _d, i, key in cands:
            if i in taken or key in used:
                continue
            taken[i] = key
            used.add(key)
        at: dict = {}
        for i, (pair, u) in enumerate(found):
            key = taken.get(i)
            if key is None:
                counter[pair] = counter.get(pair, -1) + 1
                key = ("partition_boundary", pair, counter[pair])
            at[key] = u
            last[key] = u
        out.append(at)
    return out


def _station_loop(graph, flip: frozenset, loop_names: Sequence[str], h: float, chord: float) -> _StationLoop:
    return _StationLoop(graph, flip, loop_names, h, chord)


def _start_at_te(graph, flip: frozenset, loop_names: Sequence[str], te: np.ndarray) -> tuple:
    """``loop_names`` turned to start at the partition whose start is nearest ``te``."""
    by_name = {r.name: r for r in graph.loop}
    starts = []
    for n in loop_names:
        p = np.asarray(by_name[n].points, float)
        starts.append(p[-1] if n in flip else p[0])
    i0 = int(np.argmin([np.linalg.norm(q[:2] - te[:2]) for q in starts]))
    names = list(loop_names)
    return tuple(names[i0:] + names[:i0])


def _wrap(u: float) -> float:
    """``u`` in ``[0, 1)``, with the ring end taken as the ring start."""
    u = float(u) % 1.0
    return 0.0 if u >= 1.0 - _U_TOL else u


def _held_u(ref: dict, here: dict, key) -> float:
    """``u`` here of a line known only at ``ref``: same fraction between its neighbours there."""
    u = ref[key]
    shared = {k: v for k, v in ref.items() if k != key and k in here}
    below = [(v, k) for k, v in shared.items() if v <= u]
    above = [(v, k) for k, v in shared.items() if v > u]
    lo_u, lo_k = max(below) if below else (0.0, None)
    hi_u, hi_k = min(above) if above else (1.0, None)
    frac = (u - lo_u) / (hi_u - lo_u) if hi_u > lo_u else 0.0
    lo_here = here[lo_k] if lo_k is not None else 0.0
    hi_here = here[hi_k] if hi_k is not None else 1.0
    if hi_here < lo_here:  # the upper neighbour is the ring start, reached going round
        hi_here += 1.0
    return _wrap(lo_here + frac * (hi_here - lo_here))


def _ring_breaks(at: dict) -> list:
    """``[(key, u)]`` round the ring from ``("end", 0)`` to ``("end", 1)``.

    A line within :data:`_U_TOL` of the ring start sits on it; ties sort by key
    so the order is the same at every station.
    """
    items = []
    for key, u in at.items():
        u = float(u)
        if u <= _U_TOL or u >= 1.0 - _U_TOL:
            u = 0.0
        items.append((key, u))
    items.sort(key=lambda kv: (kv[1], repr(kv[0])))
    return [(("end", 0), 0.0)] + items + [(("end", 1), 1.0)]


#: How fast the local edge may grow along a run, per metre of arc (graded sizing).
_GRADING = 0.25


#: Arc length, as a fraction of chord, over which the sizing curvature is measured.
_KAPPA_WINDOW = 0.02


def _window_kappa(pts: np.ndarray, window: float) -> np.ndarray:
    """Curvature at each point as the tangent turn over ``window`` metres of arc, divided by it.

    A nose of radius well above ``window`` gives its true curvature. A single
    kinked vertex (a midline offset step, a polyline glitch) gives only its
    angle over ``window``, not over one tiny segment, so it does not shrink the
    mesh there.
    """
    xy = np.asarray(pts, float)[:, :2]
    d = np.diff(xy, axis=0)
    seg = np.linalg.norm(d, axis=1)
    keep = seg > 1e-15
    if keep.sum() < 2:
        return np.zeros(len(xy))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    mid = 0.5 * (s[:-1] + s[1:])[keep]
    theta = np.unwrap(np.arctan2(d[keep, 1], d[keep, 0]))
    half = 0.5 * float(window)
    lo = np.clip(s - half, mid[0], mid[-1])
    hi = np.clip(s + half, mid[0], mid[-1])
    span = np.maximum(hi - lo, 1e-15)
    turn = np.abs(np.interp(hi, mid, theta) - np.interp(lo, mid, theta))
    return np.where(hi - lo > 1e-12, turn / span, 0.0)


def _edge_measure(pts: np.ndarray, h: float, chord: float | None) -> tuple[np.ndarray, np.ndarray]:
    """``(t, W)`` at the run's points: ``W(t)`` is how many edges fit from the start to ``t``.

    The local edge is ``h``, or on the skin (``chord`` given) the sagitta edge
    ``sqrt(8 EPS_OVER_C c / kappa)`` where that is shorter, the same rule as the
    ruled-strip mesher's region counts. It is graded so it grows by at most
    :data:`_GRADING` per metre of arc, which spreads the leading-edge refinement
    smoothly into the flatter skin. Equal steps in ``W`` give the node spacing.
    """
    pts = np.asarray(pts, float)
    tt = run_t(pts)
    seg = np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)
    length = float(seg.sum())
    if length <= 0.0 or len(pts) < 2:
        return np.array([0.0, 1.0]), np.array([0.0, 0.0])
    edge = np.full(len(pts), float(h))
    if chord is not None and chord > 0.0 and len(pts) >= 3:
        kap = _window_kappa(pts, _KAPPA_WINDOW * float(chord))
        with np.errstate(divide="ignore"):
            sag = np.where(kap > KAPPA_FLOOR, np.sqrt(8.0 * EPS_OVER_C * float(chord) / np.maximum(kap, KAPPA_FLOOR)), np.inf)
        edge = np.minimum(edge, sag)
        for i in range(1, len(edge)):  # grade forwards, then backwards
            edge[i] = min(edge[i], edge[i - 1] + _GRADING * seg[i - 1])
        for i in range(len(edge) - 2, -1, -1):
            edge[i] = min(edge[i], edge[i + 1] + _GRADING * seg[i])
    density = 1.0 / edge  # edges per metre
    w = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * seg)])
    return tt, w


def _ts_between(tt: np.ndarray, w: np.ndarray, t0: float, t1: float, n: int) -> np.ndarray:
    """``n + 1`` run parameters from ``t0`` to ``t1`` at equal steps of the edge measure."""
    w0, w1 = float(np.interp(t0, tt, w)), float(np.interp(t1, tt, w))
    if not w1 > w0:
        return np.linspace(t0, t1, n + 1)
    ts = np.interp(np.linspace(w0, w1, n + 1), w, tt)
    ts[0], ts[-1] = t0, t1
    return ts


#: A span bay may be this many times the station's mean skin edge (as the ruled mesher).
_SPAN_OVER_EDGE = 1.25


def _mean_loop_edge(graph, h: float, chord: float) -> float:
    """Mean graded edge round the skin loop at one station: loop length over its edge measure."""
    length = 0.0
    edges = 0.0
    for run in graph.loop:
        pts = np.asarray(run.points, float)
        if len(pts) < 2:
            continue
        _, w = _edge_measure(pts, h, chord)
        length += float(np.sum(np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)))
        edges += float(w[-1])
    return length / edges if edges > 0.0 else float(h)


def _grade_span(geo, z_list: list[float], graphs: list, h: float) -> tuple[list[float], list]:
    """Split span bays longer than :data:`_SPAN_OVER_EDGE` times the skin edge at either end.

    Curvature sizing shortens the skin edges near the nose; without this the
    quads there would be long in span and short across it. One pass, evenly,
    and never into steps under ``h / 2``.
    """
    mean = [_mean_loop_edge(g, h, float(geo.interpolate(z)[0])) for z, g in zip(z_list, graphs)]
    out_z, out_g = [z_list[0]], [graphs[0]]
    for i in range(len(z_list) - 1):
        za, zb = z_list[i], z_list[i + 1]
        limit = _SPAN_OVER_EDGE * min(mean[i], mean[i + 1])
        n = int(math.ceil((zb - za) / limit - 1e-9)) if limit > 0.0 else 1
        # As the ruled mesher: a step stays at least h / 2. A skin edge already
        # shorter than that is the geometry (a tight nose), and the span step
        # does not chase it.
        while n > 1 and (zb - za) / n < 0.5 * h * (1.0 - 1e-9):
            n -= 1
        for k in range(1, max(n, 1)):
            z = za + (zb - za) * k / n
            out_z.append(z)
            out_g.append(geo.section_at(z).shell_graph())
        out_z.append(zb)
        out_g.append(graphs[i + 1])
    return out_z, out_g


#: Fewest edges round the skin loop, whatever ``h`` is, so the rib and tip faces
#: have a ring to fill.
_MIN_LOOP_EDGES = 8


def partition_mesh(geo, *, h: float, zs: Sequence[float] | None = None,
                   sub_partitions: Sequence[Any] = (), polygon_size: float | None = None) -> ShellMesh:
    """A structured (t, z) quad mesh of the blade's skin and webs, with rib and tip triangles.

    Parameters
    ----------
    geo:
        The blade (anything with ``section_at(z).shell_graph()`` and stations).
    h:
        Target element edge in metres, across each partition. Each segment
        between mesh lines gets the largest count it needs at any station, so
        the edge is at most ``h`` everywhere.
    zs:
        Span stations; default :func:`mesh_stations` with ``h``.
    sub_partitions:
        :class:`~implicit_geometry_kernel.structure.SubPartition` bands whose
        ``t`` edges get mesh lines. Their span windows are not mesh lines.
    polygon_size:
        Triangle spacing on ribs and the tip closure (default ``h``).

    Rib stations are added to ``zs``. Each rib face and the tip closure is
    split into cells by the webs; each cell is a triangle fill whose boundary
    nodes are the skin and web nodes of that station, so only interior nodes
    are new.

    Raises
    ------
    PartitionMeshError
        If the loop is open or branched, a web stops and starts again, or mesh
        lines round the loop would cross or meet at some stations only.
    """
    if not float(h) > 0.0:
        raise ValueError("element size h must be positive")
    z_list = [float(z) for z in (zs if zs is not None else mesh_stations(geo, h=float(h)))]
    for m in getattr(geo, "transverse", ()) or ():
        if min(z_list) - 1e-9 <= float(m.z) <= max(z_list) + 1e-9:
            z_list.append(float(m.z))
    z_list = sorted({round(z, 9): z for z in z_list}.values())
    if len(z_list) < 2:
        raise ValueError("the partition mesher needs at least two stations")
    graphs = [geo.section_at(z).shell_graph() for z in z_list]
    if zs is None:
        z_list, graphs = _grade_span(geo, z_list, graphs, float(h))
    chains = [_graph_signature(g)[0] for g in graphs]
    n_st = len(z_list)
    flips = [_flipped(g) for g in graphs]
    skeletons = [g.skeleton() for g in graphs]
    chords = [float(geo.interpolate(z)[0]) for z in z_list]

    # The skin loop at each station as one closed ring, walked along the loop
    # direction from the start of the first partition. ``u`` is the ring's
    # normalised arc length; each partition holds a ``u`` interval. The edge
    # measure is each partition's own graded measure, laid end to end.
    # The ring starts at the partition boundary nearest the trailing edge, so
    # its start does not move with a cap edge.
    _le, te_pts, _chord = _le_te_points(geo, z_list)
    # Skin partitions may appear or vanish along the span (decision 0018): each
    # station has its own chain, and their boundaries are held lines where absent.
    names_at = [_start_at_te(g, flips[j], chains[j], np.asarray(te_pts[j], float)) for j, g in enumerate(graphs)]
    loop_names = sorted({n for c in chains for n in c})
    loops = [_station_loop(g, flips[j], names_at[j], float(h), chords[j]) for j, g in enumerate(graphs)]

    # Webs may start or stop inside the span. Each web's stations must be one
    # unbroken run. A web that ends on a station (a birth or death anchor) is
    # taken from the section a hair inside its span there, so its strip reaches
    # that station.
    web_names = sorted({w.name for g in graphs for w in g.webs})
    web_src: dict = {}  # (station, web) -> (run, skeleton, or None for a run from just inside the span)
    for j, g in enumerate(graphs):
        for w in g.webs:
            web_src[(j, w.name)] = (w, skeletons[j])
    for name in web_names:
        js = [j for j in range(n_st) if (j, name) in web_src]
        if js != list(range(js[0], js[-1] + 1)):
            raise PartitionMeshError(
                f"web {name!r} stops and starts again inside the span; the partition mesher needs one unbroken run"
            )
        for j, side in ((js[0] - 1, 1.0), (js[-1] + 1, -1.0)):
            if not 0 <= j < n_st:
                continue
            # A found web end is good to about 1e-5 of the span, so look a
            # little further in, up to a thousandth of the bay.
            bay = abs(z_list[j + int(side)] - z_list[j])
            for k in range(7, 2, -1):
                try:
                    gd = geo.section_at(z_list[j] + side * bay * 10.0 ** -k).shell_graph()
                except Exception:  # noqa: BLE001 - a cut that fails there (BG-062): try nearer, else stop
                    continue
                run = next((w for w in gd.webs if w.name == name), None)
                if run is not None:
                    web_src[(j, name)] = (run, None)
                    break
    webs_at = [[w for w in web_names if (j, w) in web_src] for j in range(n_st)]
    hosts = {s.host for s in sub_partitions}
    unknown = hosts - set(loop_names) - set(web_names)
    if unknown:
        raise ValueError(f"sub-partition hosts not in the shell graph: {sorted(unknown)}")

    # Mesh lines round the loop, as ``{key: u}`` per station. Every web foot is
    # a line at every station: where the web is, at its foot (an end station's
    # foot is the web end projected on the ring); where it is not, held at the
    # ``u`` of the nearest station that has it. So the loop's lines never change
    # in number, and a web that starts or stops is an all-quad patch on them.
    # A partition boundary is a line of its own unless a web foot is on it: a
    # cap edge that moves from one web foot to another at a web end moves
    # between two lines, never adds one.
    feet: dict = {}  # (station, web, end) -> u
    for (j, name), (run, sk) in web_src.items():
        lp = loops[j]
        for end in (0, 1):
            if sk is not None:
                ft = _web_feet(sk, name, names_at[j], flips[j]).get(end)
                if ft is not None:
                    feet[(j, name, end)] = lp.u_of(ft[0], ft[1])
            else:
                xy = np.asarray(run.points, float)[0 if end == 0 else -1]
                feet[(j, name, end)] = _closest_t(lp.points, xy)[0]
    # Lines that are really there at each station: web feet, partition
    # boundaries no foot sits on, and sub-partition band edges.
    # On a foot a step rides the foot's line; a boundary's own line is held
    # apart there like any absent feature.
    feet = {k: _wrap(u) for k, u in feet.items()}
    # Web feet first, held among themselves where a web is absent; every other
    # line is then placed relative to them.
    foot_real = [{("junction", w, end): u for (jj, w, end), u in feet.items() if jj == j} for j in range(n_st)]
    foot_lines = _hold_absent(foot_real)
    bounds = _track_boundaries(loops, feet, foot_lines)
    real: list[dict] = []
    for j, lp in enumerate(loops):
        at: dict = dict(foot_real[j])
        at.update(bounds[j])
        for s in sub_partitions:
            if s.host in names_at[j]:
                for k, t in enumerate(s.t_range):
                    if _T_TOL < t < 1.0 - _T_TOL:
                        at[("band", s.name, k)] = lp.u_of(s.host, _along(t, s.host in flips[j]))
        real.append(at)
    # A line whose feature is absent at a station (a web that has ended, a
    # partition that has vanished) is held between the same two neighbouring
    # lines, at the same fraction of the gap, as at the nearest station with
    # the feature; so it moves with the loop and never crosses its neighbours.
    lines = _hold_absent([{**foot_lines[j], **real[j]} for j in range(n_st)])
    ring_breaks = [_ring_breaks(at) for at in lines]
    keys0 = [k for k, _ in ring_breaks[0]]
    for j in range(1, n_st):
        if [k for k, _ in ring_breaks[j]] != keys0:
            raise PartitionMeshError(
                f"the skin loop's mesh lines change order or number at z={z_list[j]:g}"
            )
    # Lines that meet (two feet on one point, or a foot on the ring start) share
    # a node; they must meet at every station or at none.
    n_ring = len(keys0) - 1
    shut = [[ring_breaks[j][k + 1][1] - ring_breaks[j][k][1] <= _U_TOL for k in range(n_ring)] for j in range(n_st)]
    for j in range(1, n_st):
        if shut[j] != shut[0]:
            raise PartitionMeshError(
                f"two of the skin loop's mesh lines meet at some stations but not at z={z_list[j]:g}"
            )
    ring_counts = [0 if shut[0][k] else 1 for k in range(n_ring)]
    for j, lp in enumerate(loops):
        b = ring_breaks[j]
        for k in range(n_ring):
            if not shut[0][k]:
                span = float(np.interp(b[k + 1][1], lp.u, lp.w) - np.interp(b[k][1], lp.u, lp.w))
                ring_counts[k] = max(ring_counts[k], int(math.ceil(span - 1e-9)))
    # A coarse h must still leave a closed loop that a rib or the tip closure
    # can fill (decision 0007): at least _MIN_LOOP_EDGES edges round the loop,
    # added one at a time to the longest-per-edge segment.
    while sum(ring_counts) < _MIN_LOOP_EDGES:
        b, length = ring_breaks[0], loops[0].length
        k_best = max((k for k in range(n_ring) if ring_counts[k] > 0),
                     key=lambda k: (b[k + 1][1] - b[k][1]) * length / ring_counts[k])
        ring_counts[k_best] += 1

    # Mesh lines on each web: the same keys, in the same order, at every
    # station that has the web.
    breaks: dict = {}
    counts: dict = {}
    measure: dict = {}  # (station, web) -> (t, W), the graded edge measure
    for name in web_names:
        stations = [j for j in range(n_st) if (j, name) in web_src]
        per_station: dict = {}
        for j in stations:
            run, sk = web_src[(j, name)]
            if sk is not None:
                bp = _breakpoints(sk, run, sub_partitions)
            else:
                bp = {("band", s.name, k): t for s in sub_partitions if s.host == name
                      for k, t in enumerate(s.t_range) if _T_TOL < t < 1.0 - _T_TOL}
            per_station[j] = _segment_breaks(bp)
        wkeys = [k for k, _ in per_station[stations[0]]]
        for j in stations[1:]:
            if [k for k, _ in per_station[j]] != wkeys:
                raise PartitionMeshError(f"web {name!r}: its mesh lines change order or number at z={z_list[j]:g}")
        need = [1] * (len(wkeys) - 1)
        for j in stations:
            b = per_station[j]
            tt, w = measure[(j, name)] = _edge_measure(np.asarray(web_src[(j, name)][0].points, float), float(h), None)
            for k in range(len(need)):
                span = float(np.interp(b[k + 1][1], tt, w) - np.interp(b[k][1], tt, w))
                need[k] = max(need[k], int(math.ceil(span - 1e-9)))
        breaks[name] = per_station
        counts[name] = need

    nodes: list[np.ndarray] = []

    def new_node(p: np.ndarray) -> int:
        nodes.append(np.asarray(p, float))
        return len(nodes) - 1

    ids: dict = {}  # (station, "loop" or web) -> node ids; the loop ring is closed (last id is the first)
    ring_u: dict = {}  # station -> u of each ring node
    for j, lp in enumerate(loops):
        z = z_list[j]
        b = ring_breaks[j]
        line_ids: dict = {}
        row: list[int] = []
        us: list[float] = []
        for k in range(n_ring):
            n = ring_counts[k]
            if n == 0:
                line_ids[b[k][0]] = None  # filled below by the next line's node
                continue
            ts = _ts_between(lp.u, lp.w, b[k][1], b[k + 1][1], n)
            for m, u in enumerate(ts[:-1]):
                p = run_point_at(lp.points, float(u))
                nid = new_node(np.array([p[0], p[1], z]))
                row.append(nid)
                us.append(float(u))
                if m == 0:
                    line_ids[b[k][0]] = nid
        # A line on a shut segment shares the next line's node.
        for k in range(n_ring - 1, -1, -1):
            if line_ids.get(b[k][0], 0) is None:
                nxt = b[k + 1][0]
                line_ids[b[k][0]] = line_ids[nxt] if nxt in line_ids else row[0]
        row.append(row[0])
        us.append(1.0)
        ids[(j, "loop")] = row
        ring_u[j] = np.asarray(us)
        for name in webs_at[j]:
            run = web_src[(j, name)][0]
            pts = np.asarray(run.points, float)
            wb = breaks[name][j]
            ends = {end: line_ids[("junction", name, end)] for end in (0, 1)
                    if line_ids.get(("junction", name, end)) is not None}
            wrow: list[int] = []
            for k in range(len(wb) - 1):
                tt_w = measure[(j, name)]
                ts = _ts_between(tt_w[0], tt_w[1], wb[k][1], wb[k + 1][1], counts[name][k])
                for m, t in enumerate(ts[:-1] if k < len(wb) - 2 else ts):
                    first = k == 0 and m == 0
                    last_pt = k == len(wb) - 2 and m == len(ts) - 1
                    if first and 0 in ends:
                        wrow.append(ends[0])
                    elif last_pt and 1 in ends:
                        wrow.append(ends[1])
                    else:
                        p = run_point_at(pts, float(t))
                        wrow.append(new_node(np.array([p[0], p[1], z])))
            ids[(j, name)] = wrow

    elements: list[tuple[int, ...]] = []
    rows: list[tuple] = []
    for j in range(n_st - 1):
        # Skin: each quad takes the partition its middle lies in at the lower
        # station of the bay, so a cap that ends at a station ends on that row.
        a, b = ids[(j, "loop")], ids[(j + 1, "loop")]
        ua = ring_u[j]
        for i in range(len(a) - 1):
            um = 0.5 * (ua[i] + ua[i + 1])
            r0 = loops[j].run_at(um)
            r1 = loops[j + 1].run_at(um)
            pair = r1.name == r0.name
            th = 0.5 * (float(r0.thickness) + float(r1.thickness)) if pair else float(r0.thickness)
            ecc = 0.5 * (float(r0.eccentricity_m) + float(r1.eccentricity_m)) if pair else float(r0.eccentricity_m)
            elements.append((a[i], a[i + 1], b[i + 1], b[i]))
            rows.append((r0.name, "skin", th, ecc, tuple(r0.layers), tuple(r0.layer_thicknesses_m), j,
                         z_list[j], z_list[j + 1], r0.name, "partition"))
        for name in web_names:
            if (j, name) not in ids or (j + 1, name) not in ids:
                continue  # the web is absent at one end of the bay
            a, b = ids[(j, name)], ids[(j + 1, name)]
            r0, r1 = web_src[(j, name)][0], web_src[(j + 1, name)][0]
            th = 0.5 * (float(r0.thickness) + float(r1.thickness))
            ecc = 0.5 * (float(r0.eccentricity_m) + float(r1.eccentricity_m))
            for i in range(len(a) - 1):
                elements.append((a[i], a[i + 1], b[i + 1], b[i]))
                rows.append((name, "web", th, ecc, tuple(r0.layers), tuple(r0.layer_thicknesses_m), j,
                             z_list[j], z_list[j + 1], name, "partition"))

    n_source = len(nodes)
    origins = ["station"] * n_source
    if polygon_size is None:
        polygon_size = float(h)
    tip = getattr(geo, "tip_closure", None)
    tip_z = float(getattr(geo, "tip_closure_z", z_list[-1]))
    for j, z in enumerate(z_list):
        members = [m for m in getattr(geo, "transverse", ()) or () if abs(float(m.z) - z) <= 1e-8]
        caps: list = []
        for m in members:
            caps.append(("rib", str(m.name), str(m.name), float(m.thickness), (str(m.name),)))
        if tip is not None and abs(z - tip_z) <= 1e-8 and j == len(z_list) - 1:
            t_tip = float(geo.tip_closure_thickness())
            caps.append(("tip_closure", str(tip.name), str(tip.partition), t_tip, (str(tip.partition),)))
        if not caps:
            continue
        cells = _cells(ids, j, webs_at[j])
        for surface, label, part, th, layers in caps:
            for k, cell in enumerate(cells):
                ring = np.asarray([nodes[i] for i in cell], float)
                faces = triangulate_ring(ring, size=polygon_size, pattern="unstructured", split_boundary=False)
                key_of = {tuple(np.round(nodes[i], 12)): i for i in cell}
                for tri in np.asarray(faces, float):
                    tri_ids = []
                    for p in tri:
                        key = tuple(np.round(p, 12))
                        hit = key_of.get(key)
                        if hit is None:
                            hit = new_node(np.array([p[0], p[1], z]))
                            origins.append("interior")
                            key_of[key] = hit
                        tri_ids.append(hit)
                    if len(set(tri_ids)) < 3:
                        continue
                    a3, b3, c3 = (nodes[i] for i in tri_ids)
                    if np.linalg.norm(np.cross(b3 - a3, c3 - a3)) <= 1e-12 * max(1.0, chords[j]) ** 2:
                        continue  # a sliver on three collinear ring nodes (a straight web's row)
                    elements.append(tuple(tri_ids))
                    rows.append((f"{label}:cell_{k}", surface, th, 0.0, layers, (th,), j, z, z, part, "interior"
                                 if any(i >= n_source for i in tri_ids) else "inherited"))

    pts = np.asarray(nodes, float)
    # The placed trailing edge from the loft (as the ruled mesher), not the loop
    # point farthest from the loop's mean: on a thick section dense near the
    # nose that point can be the leading edge, which flipped web normals.
    te_at = {z: np.asarray(te, float) for z, te in zip(z_list, te_pts)}
    elements, flipped = _oriented(pts, elements, [r[1] for r in rows], z_list, te_at)
    mesh = ShellMesh(
        nodes=pts,
        elements=tuple(tuple(int(i) for i in e) for e in elements),
        name=tuple(r[0] for r in rows),
        surface=tuple(r[1] for r in rows),
        partition=tuple(r[9] for r in rows),
        thickness=np.array([r[2] for r in rows], float),
        eccentricity=np.array([r[3] for r in rows], float),
        layers=tuple(r[4] for r in rows),
        layer_thicknesses=tuple(r[5] for r in rows),
        stations=tuple(z_list),
        n_source=n_source,
        bay=tuple(r[6] for r in rows),
        z0=tuple(r[7] for r in rows),
        z1=tuple(r[8] for r in rows),
        constructor=tuple(r[10] for r in rows),
        flipped=tuple(int(f) for f in flipped),
        node_origin=tuple(origins),
        section_distance=np.zeros(len(pts)),
        recipe=MeshRecipe(tip_h_over_c=None, h=float(h), h_over_c=None, n_span=None, region_h=(),
                          refine=(), edge=tuple(float(h) for _ in z_list)),
    )
    # The same quality table, rollup and assembled tables as the ruled-strip mesh.
    thickness = np.asarray(mesh.thickness, float)
    table, jmax = _measure(np.asarray(pts, float), mesh.elements, thickness, with_ratio=False)
    facts = assemble(mesh, float(z_list[0]))
    table["volume_ratio"] = volume_ratio_from_edges(facts["edges"], table["shell_volume"])
    return replace(mesh, element_quality=table, quality=quality_rollup(mesh.partition, table, jmax), **facts)


__all__ = ["PartitionMeshError", "partition_mesh"]
