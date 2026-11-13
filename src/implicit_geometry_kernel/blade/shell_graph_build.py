"""Build mid-surface runs, rib faces and solid runs for ShellGraph."""
from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING, Dict, Tuple

import numpy as np

from implicit_geometry_kernel.blade.cut import SectionWall

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.cut import SectionCut

from implicit_geometry_kernel.blade.shell_graph import (
    MidSegment,
    MidSurface,
    ShellGraph,
    ShellRun,
    _chain,
    _side_from_name,
)

def _mid_surface(cut: SectionCut) -> MidSurface:
    tol = 1e-9 * max(1.0, max(float(np.abs(w.xy).max()) for w in cut.walls))
    pts = []  # type: ignore
    lookup: Dict[Tuple[int, int], list] = {}

    def index(p):
        key = (int(round(p[0] / (100 * tol))), int(round(p[1] / (100 * tol))))
        for dk in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)):
            for k in lookup.get((key[0] + dk[0], key[1] + dk[1]), ()):
                if np.linalg.norm(pts[k] - p) <= tol:
                    return k
        pts.append(np.array(p, float))
        lookup.setdefault(key, []).append(len(pts) - 1)
        return len(pts) - 1

    skins, webs, caps = cut.walls_of("skin"), cut.walls_of("web"), cut.walls_of("cap")
    # G9: cores / TE inserts are not SectionWalls; CoreWall.to_midline is not implemented.
    # skins: orient every skin so that z x s points out of the section
    flip = {w.label: False for w in skins}
    loop = _chain([(w.label, w.xy) for w in skins], tol)
    if loop is not None:
        ring = np.vstack([xy[::-1] if rev else xy for _, xy, rev in loop])
        area = 0.5 * np.sum(ring[:, 0] * np.roll(ring[:, 1], -1) - np.roll(ring[:, 0], -1) * ring[:, 1])
        for label, _, rev in loop:
            flip[label] = rev != (area > 0.0)  # counter-clockwise: z x s points in
    # caps: which skin vertices carry a stacked footprint (xy stays the skin mid)
    capped: Dict[Tuple[str, int], str] = {}
    skin_of = {w.label: w for w in skins}
    for cap in caps:
        for host_label, (ic, ih) in cap.bonds.items():
            host = skin_of.get(host_label)
            if host is None:
                continue
            if len(ih) < 2:
                continue
            bonded_len = float(np.sum(np.linalg.norm(np.diff(host.xy[ih], axis=0), axis=1)))
            host_len = float(host.length) if host.length else 0.0
            # Spill onto the next panel past a web: skip a sliver shorter than
            # t_cap. Do not skip the primary host (most of the wall is the
            # footprint), or cap colour (and the run name) flickers along span.
            if bonded_len < cap.thickness and bonded_len < 0.5 * host_len:
                continue
            for _a, b in zip(ic, ih):
                capped[(host_label, int(b))] = cap.label
            arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(host.xy, axis=0), axis=1))])
            lo, hi = int(np.min(ih)), int(np.max(ih))
            for edge, rest, gap in ((lo, range(0, lo), arc[lo]), (hi, range(hi + 1, len(host.xy)), arc[-1] - arc[hi])):
                if len(rest) and gap < cap.thickness:
                    for v in rest:
                        capped[(host_label, v)] = cap.label
    segments = []  # type: ignore
    skin_order = [lab for lab, _, _ in loop] if loop is not None else [w.label for w in skins]
    for label in skin_order:
        w = skin_of[label]
        ids = [index(p) for p in w.xy]
        run = []
        for v in range(len(ids) - 1):
            i, j = ids[v], ids[v + 1]
            if i == j:
                continue
            ca, cb = capped.get((w.label, v)), capped.get((w.label, v + 1))
            if ca is not None and ca == cb:
                seg = (i, j, ca, (ca, w.label))
            else:
                seg = (i, j, w.label, (w.label,))  # type: ignore
            run.append(seg)
        if flip[w.label]:
            run = [(j, i, name, layers) for i, j, name, layers in reversed(run)]
        segments.extend(MidSegment(i, j, name, "skin", layers) for i, j, name, layers in run)
    t = np.radians(cut.twist)
    chord_dir = np.array([np.cos(t), -np.sin(t)])
    for w in webs:
        ids = [index(p) for p in w.xy]
        run = [(i, j) for i, j in zip(ids[:-1], ids[1:]) if i != j]  # type: ignore
        if run:
            d = pts[run[-1][1]] - pts[run[0][0]]
            if np.array([-d[1], d[0]]) @ chord_dir < 0.0:
                run = [(j, i) for i, j in reversed(run)]  # type: ignore
        segments.extend(MidSegment(i, j, w.label, "web", (w.label,)) for i, j in run)  # type: ignore
    P = np.array(pts)
    return MidSurface(P, tuple(segments), cut.z)


def _host_panel_suffix(layers: Tuple[str, ...]) -> str | None:
    """Web-bounded host panel (``bay_0``, ``leading_edge``, …) from run layers."""
    for lb in layers:
        if ":" in lb and not str(lb).startswith("spar_cap"):
            return str(lb).split(":", 1)[-1]
    return None


def _offset_to_outer_surface(points: np.ndarray, distance: float) -> np.ndarray:
    """Move a skin run from its mid-surface onto the outer mould line.

    The run is oriented so ``z × tangent`` points out of the section.
    """
    pts = np.array(points, float, copy=True)
    for i in range(len(pts)):
        if i == 0:
            direction = pts[1, :2] - pts[0, :2]
        elif i == len(pts) - 1:
            direction = pts[-1, :2] - pts[-2, :2]
        else:
            direction = pts[i + 1, :2] - pts[i - 1, :2]
        length = float(np.linalg.norm(direction))
        if length <= 1e-15:
            continue
        outward = np.array([-direction[1], direction[0]]) / length
        pts[i, 0] += outward[0] * distance
        pts[i, 1] += outward[1] * distance
    return pts


def _endpoint_gap(a: np.ndarray, b: np.ndarray) -> float:
    pa = np.asarray(a, float).reshape(-1, a.shape[-1])[:, :2]
    pb = np.asarray(b, float).reshape(-1, b.shape[-1])[:, :2]
    return float(min(
        np.linalg.norm(pa[-1] - pb[0]),
        np.linalg.norm(pa[-1] - pb[-1]),
        np.linalg.norm(pa[0] - pb[0]),
        np.linalg.norm(pa[0] - pb[-1]),
    ))


def _join_polylines(a: np.ndarray, b: np.ndarray, tol: float) -> np.ndarray:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    variants = ((a, b), (a, b[::-1]), (a[::-1], b), (a[::-1], b[::-1]))
    dist = [float(np.linalg.norm(x[-1, :2] - y[0, :2])) for x, y in variants]
    x, y = variants[int(np.argmin(dist))]
    if dist[int(np.argmin(dist))] <= tol:
        return np.vstack([x, y[1:]])
    return np.vstack([x, y])


def _colour_kind(surface: str, layers: Tuple[str, ...], walls: Dict[str, SectionWall]) -> str:
    if surface == "web":
        return "web"
    for label in layers:
        w = walls.get(label)
        if w is not None and w.kind == "cap":
            return "cap"
    return "skin"


def _runs_of(ms: MidSurface, segs, cut: SectionCut, surface: str) -> Tuple[ShellRun, ...]:
    if not segs:
        return ()
    by_label = {w.label: w for w in cut.walls}
    z = float(cut.z)

    def xyz(i: int) -> np.ndarray:
        p = ms.points[int(i)]
        return np.array([p[0], p[1], z], dtype=float)

    def pack(name, layers, pts) -> ShellRun:
        layer_t = tuple(float(by_label[lb].thickness) for lb in layers if lb in by_label)
        t_cap = 0.0
        t_skin = 0.0
        for lb in layers:
            wall = by_label.get(lb)
            if wall is None:
                continue
            if wall.kind == "cap":
                t_cap = float(wall.thickness)
            elif wall.kind == "skin":
                t_skin = float(wall.thickness)
        if surface == "skin":
            mesh_t = t_skin if t_skin else (layer_t[-1] if layer_t else 0.0)
        else:
            mesh_t = float(sum(layer_t))
        ecc = 0.5 * t_cap if t_cap else 0.0
        points = np.vstack(pts)
        if (surface == "skin" and getattr(cut, "_wall_reference", "mid") == "oml"
                and mesh_t > 0.0 and len(points) >= 2):
            points = _offset_to_outer_surface(points, 0.5 * float(mesh_t))
        return ShellRun(
            name, surface, _colour_kind(surface, layers, by_label), layers,
            points, float(mesh_t),
            requires_eccentricity=t_cap > 0.0 and surface == "skin",
            eccentricity_m=ecc,
            layer_thicknesses_m=layer_t,
            skin_side=_side_from_name(name),
        )

    runs = []
    name, layers = segs[0].name, segs[0].layers
    pts = [xyz(segs[0].start), xyz(segs[0].end)]
    for s in segs[1:]:
        # Same partition name is one loft run even when a cap footprint
        # spans two host skin walls (layers differ). Duplicate names
        # cannot be lofted.
        if s.name == name:
            nxt = xyz(s.start)
            if np.linalg.norm(np.asarray(pts[-1][:2], float) - nxt[:2]) <= 1e-6 * max(
                1.0, float(np.max(np.abs(nxt[:2])))
            ):
                pts.append(xyz(s.end))
                continue
        runs.append(pack(name, layers, pts))
        name, layers = s.name, s.layers
        pts = [xyz(s.start), xyz(s.end)]
    runs.append(pack(name, layers, pts))
    if surface != "skin":
        return tuple(runs)
    # Join same-name fragments only when endpoints meet (one contiguous
    # arc). A cap that paints two web-bounded panels stays two runs.
    scale = 1.0
    for run in runs:
        scale = max(scale, float(np.max(np.abs(run.points[:, :2]))))
    tol = 1e-6 * scale
    contiguous: list[ShellRun] = []
    for run in runs:
        if contiguous and run.name == contiguous[-1].name:
            prev = contiguous[-1]
            if _endpoint_gap(prev.points, run.points) <= tol:
                joined = _join_polylines(prev.points, run.points, tol)
                contiguous[-1] = replace(prev, points=joined)
                continue
        contiguous.append(run)
    counts: Dict[str, int] = {}
    for run in contiguous:
        counts[run.name] = counts.get(run.name, 0) + 1
    used: set[str] = set()
    unique: list[ShellRun] = []
    for run in contiguous:
        name = run.name
        if counts.get(name, 1) > 1:
            suffix = _host_panel_suffix(run.layers)
            name = f"{run.name}:{suffix}" if suffix else name
            base = name
            i = 0
            while name in used:
                i += 1
                name = f"{base}_{i}"
        used.add(name)
        unique.append(run if name == run.name else replace(run, name=name))
    return tuple(unique)


def _closed_polyline(points: np.ndarray) -> np.ndarray:
    pts = np.asarray(points, float)
    if len(pts) == 0:
        return pts.reshape(0, pts.shape[1] if pts.ndim == 2 else 3)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    return pts


def _skin_loop_polyline(graph: ShellGraph) -> np.ndarray:
    if not graph.loop:
        return np.zeros((0, 3))
    chunks = [np.asarray(graph.loop[0].points, float)]
    for run in graph.loop[1:]:
        pts = np.asarray(run.points, float)
        if np.allclose(chunks[-1][-1], pts[0]):
            chunks.append(pts[1:])
        else:
            chunks.append(pts)
    return _closed_polyline(np.vstack(chunks))


def _rib_faces(cut: SectionCut, graph: ShellGraph) -> Tuple[ShellRun, ...]:
    """Rib faces at this station, then the tip closure on the tip station.

    Each fills the cut's cells (or the skin loop when there are none), so a face
    meets every web that reaches the station. Rib faces are ``surface="rib"``. The
    tip closure is part of the skin: ``surface="skin"``, ``kind="tip_closure"`` and
    the skin partition in ``layers``.
    """
    members = tuple(getattr(cut, "transverse", ()) or ())
    tip = getattr(cut, "tip_closure", None)
    if not members and tip is None:
        return ()
    loops: list[tuple[str | None, np.ndarray]] = []
    if cut.cells:
        for cell in cut.cells:
            loops.append((str(cell.id), _closed_polyline(cell.points)))
    else:
        pts = _skin_loop_polyline(graph)
        if len(pts):
            loops.append((None, pts))
    faces = []
    for member in members:
        for cell_id, pts in loops:
            name = member.name if len(loops) == 1 else f"{member.name}:cell_{cell_id}"
            faces.append(ShellRun(name, "rib", member.name, (name,), pts, float(member.thickness)))
    if tip is not None:
        for cell_id, pts in loops:
            name = tip.name if len(loops) == 1 else f"{tip.name}:cell_{cell_id}"
            faces.append(ShellRun(name, "skin", "tip_closure", (tip.partition,), pts, float(tip.thickness)))
    return tuple(faces)


def _chord_fractions(cut: SectionCut, xy: np.ndarray) -> np.ndarray:
    t = np.radians(float(cut.twist))
    c, s = np.cos(t), np.sin(t)
    x = xy[:, 0] - float(cut.offset[0])
    y = xy[:, 1] - float(cut.offset[1])
    local_x = c * x - s * y
    return (local_x + float(cut.pitch_axis) * float(cut.chord)) / float(cut.chord)


def _stringer_branches(cut: SectionCut) -> Tuple[ShellRun, ...]:
    members = tuple(getattr(cut, "stringers", ()) or ())
    if not members:
        return ()
    out = []
    for member in members:
        walls = [w for w in cut.walls if w.kind == "skin" and (w.side == member.side or w.side is None)]
        if not walls:
            continue
        pts = np.vstack([w.points for w in walls])
        frac = _chord_fractions(cut, pts[:, :2])
        i = int(np.argmin(np.abs(frac - member.position)))
        foot = np.asarray(pts[i], float)
        if member.free_end:
            j0, j1 = max(0, i - 1), min(len(pts) - 1, i + 1)
            tangent = pts[j1, :2] - pts[j0, :2]
            nrm = float(np.linalg.norm(tangent))
            if nrm < 1e-15:
                continue
            tangent = tangent / nrm
            normal = np.array([-tangent[1], tangent[0]])
            if cut.cells:
                centroid = np.mean(np.vstack([c.points[:, :2] for c in cut.cells]), axis=0)
                if float(np.dot(centroid - foot[:2], normal)) < 0.0:
                    normal = -normal
            tip_xy = foot[:2] + normal * float(member.depth)
            tip = np.array([tip_xy[0], tip_xy[1], foot[2]])
        else:
            opp = [w for w in cut.walls if w.kind == "skin" and w.side and w.side != member.side]
            if not opp:
                continue
            opts = np.vstack([w.points for w in opp])
            k = int(np.argmin(np.abs(_chord_fractions(cut, opts[:, :2]) - member.position)))
            tip = np.asarray(opts[k], float)
        out.append(ShellRun(member.name, "stringer", member.name, (member.name,),
                            np.vstack([foot, tip]), float(member.thickness)))
    return tuple(out)


def _exterior_xy(geom) -> np.ndarray:
    if geom is None or getattr(geom, "is_empty", True):
        return np.zeros((0, 2))
    poly = geom
    if getattr(geom, "geom_type", "") == "MultiPolygon":
        geoms = [g for g in geom.geoms if not g.is_empty]
        if not geoms:
            return np.zeros((0, 2))
        poly = max(geoms, key=lambda g: g.area)
    return np.asarray(poly.exterior.coords, float)


def _solid_runs(cut: SectionCut) -> Tuple[ShellRun, ...]:
    """Core and insert loops from constructive polygons, never a fake midline."""
    host = getattr(cut, "host", None)
    if host is None:
        return ()
    from implicit_geometry_kernel.section.boundary.cores import emit_constructive_polygons
    from implicit_geometry_kernel.section.boundary.labels import is_core_label, is_insert_label
    from implicit_geometry_kernel.section.boundary.walls import (
        _component_map,
        _insert_poly_via_recipe,
        _label_key,
    )

    place = getattr(cut, "place_host_xy", None)
    out: list[ShellRun] = []
    try:
        polys = emit_constructive_polygons(host)
    except Exception:
        polys = {}
    for key, geom in polys.items():
        if not is_core_label(key):
            continue
        xy = _exterior_xy(geom)
        if len(xy) < 3:
            continue
        pts = place(xy) if callable(place) else np.column_stack([xy, np.full(len(xy), cut.z)])
        out.append(ShellRun(str(key), "core", "core", (str(key),), pts, 0.0))
    try:
        comps = _component_map(host)
    except Exception:
        comps = {}
    for label, comp in comps.items():
        key = _label_key(label)
        if not is_insert_label(key):
            continue
        try:
            geom = _insert_poly_via_recipe(host, key, comp)
        except Exception:
            continue
        xy = _exterior_xy(geom)
        if len(xy) < 3:
            continue
        pts = place(xy) if callable(place) else np.column_stack([xy, np.full(len(xy), cut.z)])
        out.append(ShellRun(str(key), "insert", "insert", (str(key),), pts, 0.0))
    return tuple(out)

