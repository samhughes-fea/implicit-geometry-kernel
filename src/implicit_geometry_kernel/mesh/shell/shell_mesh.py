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
from typing import Any, Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology import (
    rib_cell_polygons_from_samples,
    sample_region_polylines,
)
from implicit_geometry_kernel.blade.topology.region_counts import region_counts
from implicit_geometry_kernel.mesh.fill import triangulate_ring
from implicit_geometry_kernel.mesh.quality import (
    QUALITY_DTYPE, MeshRunQuality, _measure, quality_rollup, volume_ratio_from_edges,
)
from implicit_geometry_kernel.mesh.tables import (
    CORNER_DTYPE, EDGE_DTYPE, StationDecision, StripDecision, assemble,
)
from implicit_geometry_kernel.profiling import report_stage, stage
from implicit_geometry_kernel.mesh.shell.stations import (
    _chord_at,
    _grade_stations,
    _size_fraction,
    _target_edge,
    _unique_z,
    mesh_stations,
)
from implicit_geometry_kernel.mesh.shell.orient import (
    _nearest_runs,
    _oriented,
    _plain_ruled,
    _point_key,
    _ruled_station_queries,
    _split_crossed_triangles,
    _weld_triangles,
)
from implicit_geometry_kernel.mesh.shell.spanwise import (
    _closest_on_loop_many,
    _freeze_sagitta,
    _sagitta_pending,
    _span_edges,
)
from implicit_geometry_kernel.mesh.shell.section_blend import (
    _bay_weight,
    _branch_names,
    _edge_sagitta,
    _ensure_site_cuts,
    _lerp_polylines,
    _pin_station_skin_ends,
    _place_polylines,
    _placed_cuts,
    _resample_counts,
    _section_projector,
    _shared_count_gap,
    _station_kind,
    _template_at,
)
from implicit_geometry_kernel.mesh.shell.stations import (  # noqa: F401 (re-exported helpers)
    TIP_H_OVER_C,
    _anchor_stations,
    _drop_stations,
    _fill_gaps,
    _median_skin_edge,
    _n_bays,
    _refine_stations,
)
from implicit_geometry_kernel.mesh.shell.orient import (  # noqa: F401 (re-exported helpers)
    _area_vector,
    _areas_means,
    _nearest_run,
    _run_distances,
    orient_shell_elements,
)
from implicit_geometry_kernel.mesh.shell.spanwise import (  # noqa: F401 (re-exported helpers)
    _closest_on_loop,
    _fill_le_te_cache,
    _le_te_point,
    _le_te_points,
    _span_sagitta,
    _web_fraction_pairs,
)
from implicit_geometry_kernel.mesh.shell.section_blend import (  # noqa: F401 (re-exported helpers)
    _bay_of,
    _blend_frames,
    _blend_open,
    _blend_open_many,
    _blend_wall_points,
    _cache_graph,
    _cut_at_z,
    _frame_at,
    _gauge_key,
    _lerp_cut,
    _lerp_cuts_bay,
    _lerp_runs,
    _lerp_runs_many,
    _live_fractions,
    _loft_gap,
    _loop_wind,
    _named_window,
    _nearest_at,
    _nearest_key,
    _nearest_stations,
    _open_unit_loop,
    _place_graph,
    _place_runs,
    _placed_unit,
    _refine_shape_sites,
    _run_units,
    _section_sites,
    _shape_keys,
    _sibling_cut,
    _sibling_each,
    _stacked_shell,
    _station_wrap,
    _wall_units,
    _window_key,
)


SURFACES = ("skin", "web", "rib", "tip_closure")


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
    if zs is not None:
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
        run_of: dict[tuple[int, int], Any] = {}
        for (at_tip, surface), items in grouped.items():
            station = i + 1 if at_tip else i
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
        rings: list[tuple[str, str, np.ndarray, Any]] = []
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
                pkey = _point_key(p)
                hit = index.get(pkey)
                if hit is None or hit >= n_source:
                    raise ValueError(f"partition ring vertex at z = {z} is not a skin or web node")
            if surface == "rib":
                ring_row = (name, surface, str(member.name), float(member.thickness), 0.0,
                       (str(member.name),), (float(member.thickness),))
            else:
                t_tip = float(geo.tip_closure_thickness())
                ring_row = (name, surface, str(member.partition), t_tip, 0.0,
                       (str(member.partition),), (t_tip,))
            ring_faces = triangulate_ring(
                ring, size=polygon_size, pattern="unstructured", split_boundary=False,
            )
            if len(ring_faces) == 0:
                continue
            station.append((ring_row, np.asarray(ring_faces, float)))
        if not station:
            continue
        welded = _weld_triangles(
            np.concatenate([ring_faces for _ring_row, ring_faces in station], axis=0),
            index, nodes, origins,
        )
        offset = 0
        for ring_row, ring_faces in station:
            block = welded[offset:offset + len(ring_faces)]
            offset += len(ring_faces)
            for tri in block:
                if int(tri[0]) < 0:
                    continue
                ids = [int(k) for k in tri]
                kind = "interior" if any(k >= n_source for k in ids) else "inherited"
                add(ids, ring_row, (si, float(z), float(z), kind, 0, 0))
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
        elements=tuple(elements),
        name=tuple(r[0] for r in rows),
        surface=tuple(r[1] for r in rows),
        partition=tuple(partition),
        thickness=thickness,
        eccentricity=np.array([r[4] for r in rows], float),
        layers=tuple(r[5] for r in rows),
        layer_thicknesses=tuple(r[6] for r in rows),
        stations=tuple(float(z) for z in zs),
        quality=(),
        recipe=recipe,
        element_quality=table,
        n_source=n_source,
        bay=tuple(m[0] for m in meta),
        z0=tuple(m[1] for m in meta),
        z1=tuple(m[2] for m in meta),
        constructor=tuple(m[3] for m in meta),
        n0=tuple(m[4] for m in meta),
        n1=tuple(m[5] for m in meta),
        flipped=tuple(flip_flags),
        node_origin=tuple(origins),
        section_distance=distance,
        counts=tuple(count_log),
        station_log=station_log,
        strips=tuple(strips),
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
