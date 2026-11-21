"""Section cuts between sites: placing, blending and sampling the cut at any span station.

Part of :mod:`implicit_geometry_kernel.mesh.shell.shell_mesh`, split out of it in 0.9.4.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology import (
    place_station_keys,
    section_cut_near,
)
from implicit_geometry_kernel.blade.topology.map_walk import (
    _skin_mid_loop, te_face_open, te_wrap_for_mesh,
)
from implicit_geometry_kernel.blade.topology.ribs import _template_from_cut
from implicit_geometry_kernel.blade.topology.sizing import EPS_OVER_C
from implicit_geometry_kernel.mesh.shell.stations import (
    _anchor_stations,
    _target_edge,
    _unique_z,
)
from implicit_geometry_kernel.mesh.shell.spanwise import (
    _closest_on_loop,
    _closest_on_loop_many,
    _span_edges,
    _web_fraction_pairs,
)


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
            pending = [int(i) for i in chosen]
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
        parts: list[float] = []
        for side in sorted(raw):
            window = raw[side]
            if window is None:
                continue
            parts.extend(float(v) for v in window)
        return tuple(round(v, 6) for v in parts)
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
    if webs and width and raw is not None:
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


def _ensure_site_cuts(geo) -> list[tuple[float, Any]]:
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


def _run_units(runs, cut) -> list[tuple[Any, np.ndarray | None]]:
    """Chord-normal polyline of each run, or ``None`` when the run has no points."""
    from implicit_geometry_kernel.blade.loft.span_place import to_chord_normal

    out: list[tuple[Any, np.ndarray | None]] = []
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
