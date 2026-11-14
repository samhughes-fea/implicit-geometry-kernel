"""The blade's conforming shell mesh: quads on the skin and webs, polygons on ribs and the tip.

:func:`shell_mesh` is the kernel's mesh product. Skin and web mid-surfaces are
all-quad mapped strips between span stations; each rib face and the tip
closure is a polygon mesh (every cell has at least five nodes) whose perimeter
nodes are the skin and web nodes of its station. Every element carries the
partition it lies on and that partition's shell data from the cut
(:class:`~implicit_geometry_kernel.blade.shell_graph.ShellRun`): mesh
thickness, layer names and thicknesses, and the cap eccentricity. Nodes stay on
the skin mid-surface; a consumer applies the eccentricity in its constitutive
law.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology import (
    _template_from_cut,
    loft_mapped_strips,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
    section_cut_near,
)
from implicit_geometry_kernel.blade.topology.crossings import key_crossing_stations
from implicit_geometry_kernel.blade.topology.region_counts import region_counts
from implicit_geometry_kernel.mesh.polygon_voronoi import polygon_mesh_ring
from implicit_geometry_kernel.stage_log import stage

SURFACES = ("skin", "web", "rib", "tip_closure")


@dataclass(frozen=True)
class ShellMesh:
    """Nodes and elements of the blade shell, with per-element shell data.

    ``elements[i]`` lists node indices (4 for a quad, 5 or more for a polygon)
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

    @property
    def quads(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) == 4)

    @property
    def polygons(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) > 4)

    def edge_use(self) -> dict[tuple[int, int], int]:
        """How many elements use each edge (2 on a smooth sheet, more at a junction)."""
        out: dict[tuple[int, int], int] = {}
        for e in self.elements:
            for a, b in zip(e, e[1:] + e[:1]):
                key = (a, b) if a < b else (b, a)
                out[key] = out.get(key, 0) + 1
        return out


def mesh_stations(geo, n_span: int | None = None, *, h: float | None = None,
                  refine: Sequence[tuple[float, float, float]] = ()) -> tuple[float, ...]:
    """Span stations of the mesh: evenly spaced, every forced station, and key crossings.

    ``n_span`` evenly spaced stations, or about one per ``h`` of span when only
    ``h`` is given. ``refine`` holds ``(z0, z1, dz)`` windows that get stations
    at most ``dz`` apart. Forced stations (ribs, ply drops, stringer ends), the
    tip closure and the stations where two skin keys meet are always added.
    """
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    if n_span is None:
        if h is None:
            raise ValueError("give n_span or h")
        n_span = max(int(np.ceil((z1 - z0) / float(h))) + 1, 2)
    zs = {round(float(z), 9) for z in np.linspace(z0, z1, int(n_span))}
    for a, b, dz in refine:
        a, b = max(float(min(a, b)), z0), min(float(max(a, b)), z1)
        if float(dz) <= 0.0:
            raise ValueError("refinement spacing dz must be positive")
        if b > a:
            zs |= {round(float(z), 9) for z in np.linspace(a, b, int(np.ceil((b - a) / float(dz) - 1e-9)) + 1)}
    zs |= {round(float(z), 9) for z in geo.forced_stations()}
    zs |= {round(float(m.z), 9) for m in geo.transverse}
    zs.add(round(float(geo.tip_closure_z), 9))
    base = sorted(z for z in zs if z0 - 1e-9 <= z <= z1 + 1e-9)
    return tuple(sorted(set(base) | set(key_crossing_stations(geo, base))))


def _nearest_run(runs, p: np.ndarray):
    best, best_d = None, np.inf
    for run in runs:
        pts = np.asarray(run.points, float)[:, :2]
        if len(pts) == 1:
            d = float(np.linalg.norm(pts[0] - p))
        else:
            a, ab = pts[:-1], np.diff(pts, axis=0)
            t = np.clip(((p - a) * ab).sum(axis=1) / np.maximum((ab * ab).sum(axis=1), 1e-30), 0.0, 1.0)
            d = float(np.min(np.linalg.norm(a + t[:, None] * ab - p, axis=1)))
        if d < best_d:
            best, best_d = run, d
    return best


def _area_vector(p: np.ndarray) -> np.ndarray:
    """Newell normal of a polygon: its area times its unit normal (right-hand rule on the node order)."""
    return 0.5 * np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)


def _oriented(nodes: np.ndarray, elements: list[tuple[int, ...]], surface: Sequence[str],
              zs: Sequence[float], te_at: Mapping[float, np.ndarray]) -> list[tuple[int, ...]]:
    """Elements reordered so that their normals follow one rule per surface.

    Skin: out of the section. The skin is made consistent across every edge two
    skin elements share, then each connected sheet takes the sign for which
    ``sum(A_e . (x_e - c))`` is positive (``A_e`` area vector, ``c`` any point
    of the station): for a closed section that sum is twice its area with
    outward normals, whatever the section's concavity. Webs: towards the
    trailing edge of their station. Ribs and the tip closure: along +z.
    """
    out = list(elements)
    z_st = np.asarray(sorted(zs), float)

    def flip(i: int) -> None:
        e = out[i]
        out[i] = (e[0],) + tuple(reversed(e[1:]))

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
        flux = 0.0
        for i in sheet:
            p = nodes[list(out[i])]
            x = p.mean(axis=0)
            c = centre[station(x)]
            flux += float(_area_vector(p)[:2] @ (x - c)[:2])
        if flux < 0.0:
            for i in sheet:
                flip(i)
    for i, s in enumerate(surface):
        p = nodes[list(out[i])]
        area = _area_vector(p)
        if s == "web":
            te = te_at.get(station(p.mean(axis=0)))
            if te is not None and float(area[:2] @ (te - p.mean(axis=0))[:2]) < 0.0:
                flip(i)
        elif s in ("rib", "tip_closure") and area[2] < 0.0:
            flip(i)
    return out


def shell_mesh(geo, n_span: int | None = None, *, h: float | None = None, n_chord: int = 8,
               zs: Sequence[float] | None = None, polygon_size: float | None = None,
               region_h: Mapping[str, float] | None = None,
               refine: Sequence[tuple[float, float, float]] = ()) -> ShellMesh:
    """The blade's conforming shell mesh.

    Stations are ``zs`` or :func:`mesh_stations`. With ``h`` (target element
    edge, m) each region's node count follows its length
    (:func:`~implicit_geometry_kernel.blade.topology.region_counts.region_counts`);
    without it every region has ``n_chord`` nodes (odd). ``region_h`` sets the
    element size of named regions (``{"le_su": h / 3, "le_pr": h / 3}`` refines
    the leading edge) and ``refine`` adds span stations in ``(z0, z1, dz)``
    windows (see :func:`mesh_stations`). ``polygon_size`` is the rib and tip
    cell size (default: from the perimeter spacing).
    """
    with stage("implicit_geometry_kernel", "shell_mesh"):
        return _shell_mesh_body(
            geo, n_span, h=h, n_chord=n_chord, zs=zs,
            polygon_size=polygon_size, region_h=region_h, refine=refine)


def _shell_mesh_body(geo, n_span: int | None = None, *, h: float | None = None, n_chord: int = 8,
                     zs: Sequence[float] | None = None, polygon_size: float | None = None,
                     region_h: Mapping[str, float] | None = None,
                     refine: Sequence[tuple[float, float, float]] = ()) -> ShellMesh:
    if region_h and h is None:
        raise ValueError("region_h needs h")
    zs = tuple(sorted(float(z) for z in (
        zs if zs is not None else mesh_stations(geo, n_span, h=h, refine=refine))))
    cuts, tmpls, smps = [], [], []
    for z in zs:
        cut = section_cut_near(geo, z)
        if cut is None:
            raise ValueError(f"no section at z = {z}")
        tmpl = _template_from_cut(geo.region_graph, cut, z, n_chord)
        cuts.append(cut)
        tmpls.append(tmpl)
        smps.append(sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph))
    if h is not None:
        tmpls = region_counts(tmpls, smps, float(h), region_h=region_h)
        smps = [sample_region_polylines(c, c.shell_graph(), t, geo.region_graph) for c, t in zip(cuts, tmpls)]
    graphs = [c.shell_graph() for c in cuts]

    index: dict[tuple[float, float, float], int] = {}
    nodes: list[np.ndarray] = []

    def node(p: np.ndarray) -> int:
        key = (float(p[0]), float(p[1]), float(p[2]))
        k = index.get(key)
        if k is None:
            k = index[key] = len(nodes)
            nodes.append(np.array(key))
        return k

    elements: list[tuple[int, ...]] = []
    rows: list[tuple[str, str, str, float, float, tuple[str, ...], tuple[float, ...]]] = []

    def add(face: np.ndarray, row) -> None:
        ids: list[int] = []
        for p in np.asarray(face, float):
            k = node(p)
            if not ids or ids[-1] != k:
                ids.append(k)
        if len(ids) > 1 and ids[0] == ids[-1]:
            ids.pop()
        if len(ids) < 4:
            raise ValueError(f"element of {row[0]!r} has {len(ids)} nodes; the mesh is all quads and polygons")
        elements.append(tuple(ids))
        rows.append(row)

    web_names = {w.name for t in tmpls for w in t.webs}
    for i in range(len(zs) - 1):
        for name, face, _collapsing in loft_mapped_strips(smps[i], smps[i + 1], tmpls[i], tmpls[i + 1]):
            face = np.asarray(face, float)
            zc = float(np.mean(face[:, 2]))
            j = i if abs(zc - zs[i]) <= abs(zs[i + 1] - zc) else i + 1
            on = face[np.abs(face[:, 2] - zs[j]) <= 1e-12]
            p = (on if len(on) else face)[:, :2].mean(axis=0)
            surface = "web" if name in web_names else "skin"
            run = _nearest_run(graphs[j].webs if surface == "web" else graphs[j].loop, p)
            if run is None:
                raise ValueError(f"no {surface} partition at z = {zs[j]} for {name!r}")
            add(face, (name, surface, run.name, float(run.thickness), float(run.eccentricity_m),
                       tuple(run.layers), tuple(float(t) for t in run.layer_thicknesses_m)))

    tip_z = float(geo.tip_closure_z)
    for z, tmpl, smp in zip(zs, tmpls, smps):
        members: list = [m for m in geo.transverse if abs(float(m.z) - z) <= 1e-8]
        for m in members:
            for name, ring in rib_cell_polygons_from_samples(tmpl, smp, m):
                for f in polygon_mesh_ring(np.asarray(ring, float), size=polygon_size):
                    add(f, (name, "rib", str(m.name), float(m.thickness), 0.0, (str(m.name),),
                            (float(m.thickness),)))
        if abs(z - tip_z) <= 1e-8 and getattr(geo, "tip_closure", None) is not None:
            tip = geo.tip_closure
            t_tip = float(geo.tip_closure_thickness())
            for name, ring in rib_cell_polygons_from_samples(tmpl, smp, str(tip.name)):
                for f in polygon_mesh_ring(np.asarray(ring, float), size=polygon_size):
                    add(f, (name, "tip_closure", str(tip.partition), t_tip, 0.0, (str(tip.partition),),
                            (t_tip,)))

    pts = np.asarray(nodes, float).reshape(-1, 3)
    te_at = {float(z): np.asarray(smp["te_su"][-1], float) for z, smp in zip(zs, smps) if len(smp.get("te_su", ()))}
    elements = _oriented(pts, elements, [r[1] for r in rows], zs, te_at)
    return ShellMesh(
        nodes=pts,
        elements=tuple(elements),
        name=tuple(r[0] for r in rows),
        surface=tuple(r[1] for r in rows),
        partition=tuple(r[2] for r in rows),
        thickness=np.array([r[3] for r in rows], float),
        eccentricity=np.array([r[4] for r in rows], float),
        layers=tuple(r[5] for r in rows),
        layer_thicknesses=tuple(r[6] for r in rows),
        stations=zs,
    )


__all__ = ["SURFACES", "ShellMesh", "mesh_stations", "shell_mesh"]
