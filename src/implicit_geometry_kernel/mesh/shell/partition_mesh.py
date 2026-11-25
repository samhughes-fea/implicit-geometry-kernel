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
This step meshes blades whose shell topology is the same at every station,
with a closed skin loop and no branches (stringers). Ribs and the tip closure
are triangle fills on the station's skin and web nodes. Topology changes
(caps or webs that start or stop inside the span, as transition bands) and
Newton junction nodes between stations are later steps of plan row L6. It
sits beside :func:`shell_mesh`, which stays the default.
"""

from __future__ import annotations

import math
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


def _breakpoints(skeleton, run, bands: Sequence[Any], flip: bool = False) -> dict:
    """``{key: t}`` of interior points on ``run`` that a mesh line must pass through.

    ``t`` runs along the loop direction: a run stored against it (``flip``) has
    its own ``t`` reversed, so the mesh walks every loop partition the same way.
    """
    out: dict = {}
    for node in skeleton.nodes:
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


def _cells(ids: dict, j: int, loop_names: Sequence[str], web_names: Sequence[str]) -> list[list[int]]:
    """Node-id rings of the cells the webs cut the skin loop into at station ``j``."""
    ring: list[int] = []
    for name in loop_names:
        ring.extend(ids[(j, name)][:-1])
    cells = [ring]
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
        If the shell topology changes along the span or the skin loop is open.
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
    signature = _graph_signature(graphs[0])
    for g in graphs[1:]:
        if _graph_signature(g) != signature:
            raise PartitionMeshError(
                f"shell topology changes between z={z_list[0]:g} and z={float(g.z):g}; "
                "the partition mesher handles constant topology for now (use shell_mesh)"
            )
    loop_names, web_names = signature
    hosts = {s.host for s in sub_partitions}
    unknown = hosts - set(loop_names) - set(web_names)
    if unknown:
        raise ValueError(f"sub-partition hosts not in the shell graph: {sorted(unknown)}")

    flips = [_flipped(g) for g in graphs]

    def runs_of(g, flip: frozenset) -> dict:
        # A loop run stored against the loop direction is walked reversed, so
        # node rows always follow the loop; orientation fixes normals later.
        return {r.name: (replace(r, points=np.asarray(r.points)[::-1].copy()) if r.name in flip else r)
                for r in g.loop + g.webs}

    by_station = [runs_of(g, f) for g, f in zip(graphs, flips)]
    skeletons = [g.skeleton() for g in graphs]
    order = list(loop_names) + list(web_names)

    # Mesh lines on each run: the same keys, in the same order, at every station.
    breaks: dict = {}
    counts: dict = {}
    measure: dict = {}  # (station, run) -> (t, W), the graded edge measure
    chords = [float(geo.interpolate(z)[0]) for z in z_list]
    for name in order:
        per_station = [_segment_breaks(_breakpoints(sk, rs[name], sub_partitions, name in fl))
                       for sk, rs, fl in zip(skeletons, by_station, flips)]
        keys0 = [k for k, _ in per_station[0]]
        for g, b in zip(graphs[1:], per_station[1:]):
            if [k for k, _ in b] != keys0:
                raise PartitionMeshError(
                    f"partition {name!r}: its mesh lines change order or number at z={float(g.z):g}"
                )
        n_seg = len(keys0) - 1
        need = [1] * n_seg
        for j, (rs, b) in enumerate(zip(by_station, per_station)):
            tt, w = measure[(j, name)] = _edge_measure(
                np.asarray(rs[name].points, float), float(h), chords[j] if name in loop_names else None)
            for k in range(n_seg):
                span = float(np.interp(b[k + 1][1], tt, w) - np.interp(b[k][1], tt, w))
                need[k] = max(need[k], int(math.ceil(span - 1e-9)))
        breaks[name] = per_station
        counts[name] = need
    # A coarse h must still leave a closed loop that a rib or the tip closure
    # can fill (decision 0007): at least _MIN_LOOP_EDGES edges round the loop,
    # added one at a time to the longest-per-edge segment.
    loop_total = sum(sum(counts[name]) for name in loop_names)
    while loop_total < _MIN_LOOP_EDGES:
        best, best_len = None, -1.0
        for name in loop_names:
            pts = np.asarray(by_station[0][name].points, float)
            length = float(np.sum(np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)))
            b = breaks[name][0]
            for k, n in enumerate(counts[name]):
                per_edge = (b[k + 1][1] - b[k][1]) * length / n
                if per_edge > best_len:
                    best, best_len = (name, k), per_edge
        assert best is not None
        counts[best[0]][best[1]] += 1
        loop_total += 1

    nodes: list[np.ndarray] = []

    def new_node(p: np.ndarray) -> int:
        nodes.append(np.asarray(p, float))
        return len(nodes) - 1

    ids: dict = {}  # (station, run) -> node ids across the run
    for j, (g, rs) in enumerate(zip(graphs, by_station)):
        z = z_list[j]
        line_ids: dict = {}  # skeleton key -> node id at this station
        prev_end: int | None = None
        first_start: int | None = None
        for name in loop_names:
            run = rs[name]
            pts = np.asarray(run.points, float)
            b = breaks[name][j]
            row: list[int] = []
            closing = name == loop_names[-1] and first_start is not None
            for k in range(len(b) - 1):
                tt_w = measure[(j, name)]
                ts = _ts_between(tt_w[0], tt_w[1], b[k][1], b[k + 1][1], counts[name][k])
                for m, t in enumerate(ts[:-1] if k < len(b) - 2 else ts):
                    if k == 0 and m == 0 and prev_end is not None:
                        row.append(prev_end)
                        continue
                    if closing and first_start is not None and k == len(b) - 2 and m == len(ts) - 1:
                        row.append(first_start)  # the ring closes on the first run's start
                        continue
                    key = b[k][0] if m == 0 else None
                    p = run_point_at(pts, float(t))
                    nid = new_node(np.array([p[0], p[1], z]))
                    row.append(nid)
                    if key is not None:
                        line_ids[key] = nid
            if first_start is None:
                first_start = row[0]
            prev_end = row[-1]
            ids[(j, name)] = row
        for name in web_names:
            run = rs[name]
            pts = np.asarray(run.points, float)
            b = breaks[name][j]
            row = []
            ends = {}
            sk = skeletons[j].by_key()
            for end in (0, 1):
                node = sk.get(("junction", name, end))
                if node is None:
                    continue
                for host in loop_names:
                    t = node.runs.get(host)
                    if t is None:
                        continue
                    t = _along(t, host in flips[j])
                    hrow = ids[(j, host)]
                    if t <= _T_TOL:
                        ends[end] = hrow[0]
                    elif t >= 1.0 - _T_TOL:
                        ends[end] = hrow[-1]
                    elif node.key in line_ids:
                        ends[end] = line_ids[node.key]
                    if end in ends:
                        break
            for k in range(len(b) - 1):
                tt_w = measure[(j, name)]
                ts = _ts_between(tt_w[0], tt_w[1], b[k][1], b[k + 1][1], counts[name][k])
                for m, t in enumerate(ts[:-1] if k < len(b) - 2 else ts):
                    first = k == 0 and m == 0
                    last_pt = k == len(b) - 2 and m == len(ts) - 1
                    if first and 0 in ends:
                        row.append(ends[0])
                    elif last_pt and 1 in ends:
                        row.append(ends[1])
                    else:
                        p = run_point_at(pts, float(t))
                        row.append(new_node(np.array([p[0], p[1], z])))
            ids[(j, name)] = row

    elements: list[tuple[int, ...]] = []
    rows: list[tuple] = []
    for j in range(len(z_list) - 1):
        for name in order:
            a, b = ids[(j, name)], ids[(j + 1, name)]
            r0, r1 = by_station[j][name], by_station[j + 1][name]
            th = 0.5 * (float(r0.thickness) + float(r1.thickness))
            ecc = 0.5 * (float(r0.eccentricity_m) + float(r1.eccentricity_m))
            surface = "web" if name in web_names else "skin"
            for i in range(len(a) - 1):
                elements.append((a[i], a[i + 1], b[i + 1], b[i]))
                rows.append((name, surface, th, ecc, tuple(r0.layers), tuple(r0.layer_thicknesses_m), j,
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
        cells = _cells(ids, j, loop_names, web_names)
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
                    elements.append(tuple(tri_ids))
                    rows.append((f"{label}:cell_{k}", surface, th, 0.0, layers, (th,), j, z, z, part, "interior"
                                 if any(i >= n_source for i in tri_ids) else "inherited"))

    pts = np.asarray(nodes, float)
    # The placed trailing edge from the loft (as the ruled mesher), not the loop
    # point farthest from the loop's mean: on a thick section dense near the
    # nose that point can be the leading edge, which flipped web normals.
    _le, te_pts, _chord = _le_te_points(geo, z_list)
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
