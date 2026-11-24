"""Element orientation, crossed-triangle repair and welding for the shell mesh.

Part of :mod:`implicit_geometry_kernel.mesh.shell.shell_mesh`, split out of it in 0.9.4.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np



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
            i0, i1, i2 = (int(v) for v in elements[ei])
            for u, v, w in ((i0, i1, i2), (i1, i2, i0), (i2, i0, i1)):
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
                hit_edge = edge_of.get((ia, ib))
                if hit_edge is None:
                    continue
                edge = hit_edge
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

    # Near-coincident points weld by distance, not by rounding: two points a
    # hair apart can round to different keys when they straddle a rounding
    # boundary, which left fill nodes lying on boundary nodes (an unshared seam).
    weld = 1e-8

    def _cell(p) -> tuple[int, int, int]:
        return (int(np.floor(p[0] / weld)), int(np.floor(p[1] / weld)), int(np.floor(p[2] / weld)))

    grid: dict[tuple[int, int, int], list[int]] = {}

    def _add(k: int) -> None:
        grid.setdefault(_cell(nodes[k]), []).append(k)

    def _near(p) -> int | None:
        cx, cy, cz = _cell(p)
        best, best_d = None, weld
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for k in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        d = float(np.linalg.norm(np.asarray(nodes[k], float) - p))
                        if d < best_d:
                            best, best_d = k, d
        return best

    for k in range(len(nodes)):
        _add(k)
    ids = np.empty(len(first), np.int64)
    for i, p in enumerate(flat[first]):
        key = _point_key(p)
        found = index.get(key)
        if found is None:
            found = _near(p)
        if found is not None:
            k = found
        else:
            k = len(nodes)
            index[key] = k
            nodes.append(np.asarray(key, float))
            origins.append("interior")
            _add(k)
        ids[i] = k
    corners = ids[inverse].reshape(-1, 3)
    collapsed = (
        (corners[:, 0] == corners[:, 1])
        | (corners[:, 1] == corners[:, 2])
        | (corners[:, 0] == corners[:, 2])
    )
    corners[collapsed] = -1
    return corners
