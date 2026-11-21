"""Mesh health registry: classify and flag issues for engineer judgement.

Detectors emit :class:`HealthHit` records under registered ``code`` values.
Nothing here decides that a mesh is good or bad; dumps stay inspectable.
"""
from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

HEALTH_LAYERS = ("topology", "quality", "geometry", "construction")
HEALTH_SEVERITIES = ("info", "warn", "fail")
HEALTH_SCOPES = ("edge", "node", "element", "band", "region_pair")

MISSING_ID = -1


@dataclass(frozen=True)
class HealthClass:
    """One registered flag kind. Severity is a hint, not a publish gate."""

    code: str
    layer: str
    severity: str
    scope: str
    detector: str
    doc: str
    inspect: str = ""


@dataclass(frozen=True)
class HealthHit:
    """One classified observation on a mesh entity."""

    code: str
    severity: str
    element: int = MISSING_ID
    node: int = MISSING_ID
    edge: int = MISSING_ID
    region: str = ""
    bay: int = -1
    constructor: str = ""
    z0: float = 0.0
    z1: float = 0.0
    measure: float = float("nan")
    threshold: float = float("nan")
    detail: str = ""


@dataclass(frozen=True)
class HealthThresholds:
    """Cut-offs used only to classify hits; not a verdict."""

    angle: float = 20.0
    skew: float | None = None
    volume_ratio: float | None = None
    sagitta: float | None = None
    distance: float | None = None
    trailing_diagonal: float = 1e-3
    trailing_edge: float = 5e-3
    area_ratio: float = 2.0
    near_miss: float = 1e-8
    t_junction: float = 1e-5
    dihedral_fold: float = 45.0
    edge_over_target: float = 1.5
    sheet_gap: float = 1e-6


HEALTH_CLASSES: tuple[HealthClass, ...] = (
    HealthClass("open_loop_non_root", "topology", "fail", "band", "run_mesh_health",
                "Open edge loop whose z is not the root station.",
                "inspect loops.csv; confirm root ring vs leak z"),
    HealthClass("open_edge_off_root", "topology", "fail", "edge", "run_mesh_health",
                "use==1 edge with a node off the root station.",
                "inspect failures.csv edge id → edges.csv; region/kind/z"),
    HealthClass("near_miss", "topology", "warn", "node", "run_mesh_health",
                "Distinct nodes closer than the near-miss tolerance.",
                "inspect junctions.csv kind=near_miss; nodes.csv"),
    HealthClass("t_junction", "topology", "fail", "node", "run_mesh_health",
                "A node lies on a foreign edge (T-junction).",
                "inspect junctions.csv kind=t_junction"),
    HealthClass("ruled_valence_defect", "topology", "warn", "node", "run_mesh_health",
                "Interior ruled-skin node whose valence is not 4.",
                "inspect nodes.csv valence; VTK valence scalar"),
    HealthClass("quad_angle", "quality", "warn", "element", "run_mesh_health",
                "Quad corner angle under the angle threshold (or reflex).",
                "inspect elements.csv angle_min; corners.csv"),
    HealthClass("quad_skew", "quality", "warn", "element", "run_mesh_health",
                "Quad skew above the skew threshold.",
                "inspect elements.csv skew"),
    HealthClass("quad_volume_ratio", "quality", "warn", "element", "run_mesh_health",
                "Quad volume_ratio above the threshold.",
                "inspect elements.csv volume_ratio"),
    HealthClass("quad_inverted", "quality", "fail", "element", "run_mesh_health",
                "Element jacobian_min is negative.",
                "inspect elements.csv jacobian_min; VTK jacobian_min"),
    HealthClass("section_distance", "quality", "warn", "element", "run_mesh_health",
                "A node of the element is farther from the section than allowed.",
                "inspect nodes.csv section_distance"),
    HealthClass("edge_sagitta", "quality", "warn", "element", "run_mesh_health",
                "An element edge sagitta exceeds the threshold.",
                "inspect edges.csv sagitta"),
    HealthClass("trailing_face_angle", "quality", "warn", "element", "run_mesh_health",
                "Trailing-wrap (region te) quad is shallow.",
                "inspect elements.csv region=te angle_min"),
    HealthClass("trailing_face_inverted", "quality", "fail", "element", "run_mesh_health",
                "Trailing-wrap quad has negative Jacobian.",
                "inspect elements.csv region=te jacobian_min"),
    HealthClass("trailing_face_diagonal", "quality", "warn", "element", "run_mesh_health",
                "Trailing-wrap quad has a short diagonal across a long edge.",
                "inspect elements.csv diagonal_min for region te"),
    HealthClass("transition_band_angle", "quality", "warn", "element", "run_mesh_health",
                "Skin/web quad in a bay is shallow.",
                "inspect elements.csv angle_min by region/bay"),
    HealthClass("transition_band_area", "quality", "warn", "element", "run_mesh_health",
                "Skin/web quad area jumps vs the bay median.",
                "inspect elements.csv area by region/bay"),
    HealthClass("sheet_intersection", "geometry", "fail", "region_pair", "detect_sheet_intersection",
                "Two non-adjacent region sheets intersect off a shared edge.",
                "inspect failures.csv region_pair; VTK both regions"),
    HealthClass("ridge_mismatch", "geometry", "fail", "band", "detect_ridge_mismatch",
                "te_su and te_pr aft corners should coincide but do not.",
                "inspect station TE nodes on te_su vs te_pr"),
    HealthClass("mid_surface_fold", "geometry", "warn", "edge", "detect_mid_surface_fold",
                "Dihedral across a smooth skin seam exceeds the fold threshold.",
                "inspect edges.csv dihedral_deg meet=skin"),
    HealthClass("edge_over_target", "quality", "warn", "edge", "detect_edge_over_target",
                "Edge length exceeds the local recipe sizing target.",
                "inspect edges.csv length vs recipe.edge at station z"),
)

HEALTH_CODES: tuple[str, ...] = tuple(entry.code for entry in HEALTH_CLASSES)
_BY_CODE: dict[str, HealthClass] = {entry.code: entry for entry in HEALTH_CLASSES}


def health_class(code: str) -> HealthClass:
    """Return the registry entry for ``code``."""
    return _BY_CODE[code]


def _hit(code: str, **kwargs) -> HealthHit:
    entry = _BY_CODE[code]
    return HealthHit(code=code, severity=entry.severity, **kwargs)


def _attr_len(mesh, name: str, element: int, default):
    values = getattr(mesh, name, ())
    if values is None or len(values) <= element:
        return default
    return values[element]


def _meta(mesh, element: int) -> dict:
    region = _attr_len(mesh, "name", element, "")
    bay = int(_attr_len(mesh, "bay", element, -1))
    constructor = _attr_len(mesh, "constructor", element, "")
    z0 = float(_attr_len(mesh, "z0", element, 0.0))
    z1 = float(_attr_len(mesh, "z1", element, 0.0))
    return {"region": str(region), "bay": bay, "constructor": str(constructor), "z0": z0, "z1": z1}


def _edge_index(mesh) -> dict[tuple[int, int], int]:
    edges = np.asarray(getattr(mesh, "edges", ()))
    if len(edges) == 0:
        return {}
    return {(int(edges["a"][i]), int(edges["b"][i])): i for i in range(len(edges))}


def run_mesh_health(mesh, thresholds: HealthThresholds | None = None, *,
                    include_geometry: bool = True, **overrides) -> tuple[HealthHit, ...]:
    """Classify the mesh into registered hits. Thresholds classify; they do not judge.

    ``include_geometry`` runs sheet/ridge/fold/sizing detectors. Legacy failure
    adapters turn it off so they only reshape the catalogue classes.
    """
    thr = thresholds or HealthThresholds()
    if overrides:
        thr = HealthThresholds(**{**thr.__dict__, **overrides})
    hits: list[HealthHit] = []
    hits.extend(_topology_hits(mesh, thr))
    hits.extend(_quality_hits(mesh, thr))
    hits.extend(_trailing_hits(mesh, thr))
    hits.extend(_transition_hits(mesh, thr))
    if include_geometry:
        hits.extend(detect_ridge_mismatch(mesh))
        hits.extend(detect_mid_surface_fold(mesh, fold=thr.dihedral_fold))
        hits.extend(detect_edge_over_target(mesh, ratio=thr.edge_over_target))
        hits.extend(detect_sheet_intersection(mesh, gap=thr.sheet_gap))
    return tuple(hits)


def _topology_hits(mesh, thr: HealthThresholds) -> list[HealthHit]:
    from implicit_geometry_kernel.mesh.tables import build_edges, near_misses, t_junctions

    hits: list[HealthHit] = []
    for loop in getattr(mesh, "loops", ()):
        if not loop.is_root:
            hits.append(_hit(
                "open_loop_non_root", z0=float(loop.z), z1=float(loop.z),
                detail=f"n_edges={loop.n_edges}",
            ))
    edges = getattr(mesh, "edges", None)
    if edges is None or len(edges) == 0:
        edges = build_edges(mesh)
    nodes = np.asarray(mesh.nodes, float)
    root_z = float(np.min(nodes[:, 2])) if len(nodes) else 0.0
    names = getattr(mesh, "name", ())
    if names is None or len(names) != len(mesh.elements):
        names = ("",) * len(mesh.elements)
    for index, edge in enumerate(np.asarray(edges)):
        if int(edge["use"]) != 1:
            continue
        a, b = int(edge["a"]), int(edge["b"])
        za, zb = float(nodes[a, 2]), float(nodes[b, 2])
        if abs(za - root_z) <= 1e-9 and abs(zb - root_z) <= 1e-9:
            continue
        el = int(edge["e0"])
        region = str(names[el]) if 0 <= el < len(names) else ""
        hits.append(_hit(
            "open_edge_off_root", edge=index, element=el if el >= 0 else MISSING_ID,
            region=region, z0=min(za, zb), z1=max(za, zb),
            detail=str(edge["kind"]),
        ))
    for seam in getattr(mesh, "near_miss", None) or near_misses(nodes, thr.near_miss):
        hits.append(_hit(
            "near_miss", node=int(seam.node),
            measure=float(seam.distance), threshold=float(thr.near_miss),
            detail=f"a={seam.a} b={seam.b}",
        ))
    for seam in getattr(mesh, "junctions", None) or t_junctions(nodes, edges, thr.t_junction):
        hits.append(_hit(
            "t_junction", node=int(seam.node),
            measure=float(seam.distance), threshold=float(thr.t_junction),
            detail=f"edge={seam.a}-{seam.b}",
        ))
    for node in getattr(mesh, "valence_defects", ()):
        hits.append(_hit("ruled_valence_defect", node=int(node)))
    return hits


def _quality_hits(mesh, thr: HealthThresholds) -> list[HealthHit]:
    table = np.asarray(mesh.element_quality)
    if len(table) != len(mesh.elements):
        return []
    hits: list[HealthHit] = []
    dist = np.asarray(mesh.section_distance, float) if len(mesh.section_distance) else np.zeros(len(mesh.nodes))
    edge_rows = getattr(mesh, "edges", None)
    if edge_rows is not None and getattr(getattr(edge_rows, "dtype", None), "names", None):
        edge_sag = {(int(edge["a"]), int(edge["b"])): float(edge["sagitta"]) for edge in edge_rows}
    else:
        edge_sag = {}
    for i, element in enumerate(mesh.elements):
        meta = _meta(mesh, i)
        quad = len(element) == 4
        if quad and (float(table["angle_min"][i]) < thr.angle or float(table["angle_max"][i]) > 180.0 - thr.angle):
            hits.append(_hit(
                "quad_angle", element=i, measure=float(table["angle_min"][i]),
                threshold=float(thr.angle), **meta,
            ))
        if quad and thr.skew is not None and float(table["skew"][i]) > thr.skew:
            hits.append(_hit(
                "quad_skew", element=i, measure=float(table["skew"][i]),
                threshold=float(thr.skew), **meta,
            ))
        if quad and thr.volume_ratio is not None and float(table["volume_ratio"][i]) > thr.volume_ratio:
            hits.append(_hit(
                "quad_volume_ratio", element=i, measure=float(table["volume_ratio"][i]),
                threshold=float(thr.volume_ratio), **meta,
            ))
        if float(table["jacobian_min"][i]) < 0.0:
            hits.append(_hit(
                "quad_inverted", element=i, measure=float(table["jacobian_min"][i]),
                threshold=0.0, **meta,
            ))
        if thr.distance is not None and any(float(dist[k]) > thr.distance for k in element):
            hits.append(_hit(
                "section_distance", element=i, threshold=float(thr.distance), **meta,
            ))
        if thr.sagitta is not None:
            for a, b in zip(element, element[1:] + element[:1]):
                key = (a, b) if a < b else (b, a)
                sag = edge_sag.get(key, 0.0)
                if sag > thr.sagitta:
                    hits.append(_hit(
                        "edge_sagitta", element=i, measure=float(sag),
                        threshold=float(thr.sagitta), **meta,
                    ))
                    break
    return hits


def _trailing_hits(mesh, thr: HealthThresholds) -> list[HealthHit]:
    table = np.asarray(getattr(mesh, "element_quality", ()))
    nodes = np.asarray(mesh.nodes, float)
    diags = getattr(mesh, "diagonal_min", None)
    if diags is None:
        diags = ()
    names = getattr(mesh, "name", ())
    if names is None or len(names) == 0:
        names = ("",) * len(mesh.elements)
    hits: list[HealthHit] = []
    for i, element in enumerate(mesh.elements):
        if names[i] != "te" or len(element) != 4:
            continue
        meta = _meta(mesh, i)
        if len(table) == len(mesh.elements) and float(table["angle_min"][i]) < thr.angle:
            hits.append(_hit(
                "trailing_face_angle", element=i, measure=float(table["angle_min"][i]),
                threshold=float(thr.angle), **meta,
            ))
        if len(table) == len(mesh.elements) and float(table["jacobian_min"][i]) < 0.0:
            hits.append(_hit(
                "trailing_face_inverted", element=i, measure=float(table["jacobian_min"][i]),
                threshold=0.0, **meta,
            ))
        if i < len(diags) and float(diags[i]) < thr.trailing_diagonal:
            pts = nodes[[int(k) for k in element]]
            lengths = np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1)
            if float(np.max(lengths)) > thr.trailing_edge:
                hits.append(_hit(
                    "trailing_face_diagonal", element=i, measure=float(diags[i]),
                    threshold=float(thr.trailing_diagonal), **meta,
                ))
    return hits


def _transition_hits(mesh, thr: HealthThresholds) -> list[HealthHit]:
    table = np.asarray(getattr(mesh, "element_quality", ()))
    surfaces = getattr(mesh, "surface", None)
    if surfaces is None or len(surfaces) == 0:
        surfaces = ("",) * len(mesh.elements)
    names = getattr(mesh, "name", None)
    if names is None or len(names) == 0:
        names = ("",) * len(mesh.elements)
    bays = getattr(mesh, "bay", None)
    if bays is None or len(bays) == 0:
        bays = (-1,) * len(mesh.elements)
    groups: dict[tuple, list[int]] = {}
    for i, element in enumerate(mesh.elements):
        if surfaces[i] not in ("skin", "web") or len(element) != 4:
            continue
        groups.setdefault((names[i], int(bays[i])), []).append(i)
    hits: list[HealthHit] = []
    for (_name, _bay), members in groups.items():
        areas = [float(table["area"][i]) for i in members] if len(table) == len(mesh.elements) else []
        for index, i in enumerate(members):
            meta = _meta(mesh, i)
            if len(table) == len(mesh.elements) and float(table["angle_min"][i]) < thr.angle:
                hits.append(_hit(
                    "transition_band_angle", element=i, measure=float(table["angle_min"][i]),
                    threshold=float(thr.angle), **meta,
                ))
            others = areas[:index] + areas[index + 1:]
            if others:
                median = float(np.median(others))
                if median > 0.0:
                    area = areas[index]
                    if area < median / thr.area_ratio or area > median * thr.area_ratio:
                        hits.append(_hit(
                            "transition_band_area", element=i, measure=float(area),
                            threshold=float(median), **meta,
                        ))
    return hits


def detect_ridge_mismatch(mesh, tol: float = 1e-8) -> tuple[HealthHit, ...]:
    """Flag stations where te_su and te_pr aft corners do not share a node (no wrap)."""
    names = list(mesh.name) if len(mesh.name) else []
    if "te" in names:
        return ()
    nodes = np.asarray(mesh.nodes, float)
    origins = list(mesh.node_origin) if len(mesh.node_origin) else ["station"] * len(nodes)
    su_aft: dict[float, list[int]] = {}
    pr_aft: dict[float, list[int]] = {}
    for element, name in zip(mesh.elements, names):
        if name not in ("te_su", "te_pr") or len(element) != 4:
            continue
        station = [int(k) for k in element if origins[int(k)] == "station"]
        if len(station) < 2:
            continue
        # aft-most by max chord-ish x as a stable TE proxy on station nodes
        aft = max(station, key=lambda k: float(nodes[k, 0]))
        z = round(float(nodes[aft, 2]), 9)
        (su_aft if name == "te_su" else pr_aft).setdefault(z, []).append(aft)
    hits = []
    for z in sorted(set(su_aft) & set(pr_aft)):
        su = set(su_aft[z])
        pr = set(pr_aft[z])
        if su & pr:
            continue
        # nearest pair distance
        best = float("inf")
        for a in su:
            for b in pr:
                best = min(best, float(np.linalg.norm(nodes[a] - nodes[b])))
        if best > tol:
            hits.append(_hit(
                "ridge_mismatch", z0=z, z1=z, measure=best, threshold=float(tol),
                region="te_su|te_pr", detail=f"su={sorted(su)} pr={sorted(pr)}",
            ))
    return tuple(hits)


def detect_mid_surface_fold(mesh, fold: float = 45.0) -> tuple[HealthHit, ...]:
    """Flag skin–skin edges whose absolute dihedral exceeds ``fold`` degrees."""
    edges = np.asarray(getattr(mesh, "edges", ()))
    if len(edges) == 0:
        return ()
    hits = []
    for index, edge in enumerate(edges):
        if str(edge["meet"]) != "skin":
            continue
        angle = abs(float(edge["dihedral_deg"]))
        if angle <= fold:
            continue
        el = int(edge["e0"])
        meta = _meta(mesh, el) if 0 <= el < len(mesh.elements) else {}
        hits.append(_hit(
            "mid_surface_fold", edge=index, element=el if el >= 0 else MISSING_ID,
            measure=angle, threshold=float(fold), **meta,
        ))
    return tuple(hits)


def detect_edge_over_target(mesh, ratio: float = 1.5) -> tuple[HealthHit, ...]:
    """Flag edges longer than ``ratio`` times the recipe edge at that station."""
    recipe = getattr(mesh, "recipe", None)
    edge_targets = () if recipe is None else tuple(getattr(recipe, "edge", ()) or ())
    if not edge_targets:
        return ()
    stations = np.asarray(mesh.stations, float)
    target = np.asarray(edge_targets, float)
    if len(stations) == 0 or len(target) == 0:
        return ()
    edges = np.asarray(getattr(mesh, "edges", ()))
    nodes = np.asarray(mesh.nodes, float)
    if len(edges) == 0:
        return ()
    hits = []
    for index, edge in enumerate(edges):
        a, b = int(edge["a"]), int(edge["b"])
        mid_z = 0.5 * (float(nodes[a, 2]) + float(nodes[b, 2]))
        j = int(np.argmin(np.abs(stations - mid_z)))
        if float(target[j]) <= 0.0:
            continue
        length = float(edge["length"])
        over = length / float(target[j])
        if over <= ratio:
            continue
        el = int(edge["e0"])
        meta = _meta(mesh, el) if 0 <= el < len(mesh.elements) else {}
        meta = {**meta, "z0": mid_z, "z1": mid_z}
        hits.append(_hit(
            "edge_over_target", edge=index, element=el if el >= 0 else MISSING_ID,
            measure=over, threshold=float(ratio), **meta,
        ))
    return tuple(hits)


def detect_sheet_intersection(mesh, gap: float = 1e-6,
                              region_pairs: Sequence[tuple[str, str]] = (("te_su", "te_pr"),)
                              ) -> tuple[HealthHit, ...]:
    """Flag approximate face intersections between named region pairs off shared nodes.

    A hit means at least one triangle pair from the two regions has crossing
    half-spaces and no shared element edge. Coarse and conservative; engineer
    confirms in VTK.
    """
    nodes = np.asarray(mesh.nodes, float)
    names = list(mesh.name) if len(mesh.name) else []
    by_region: dict[str, list[tuple[int, np.ndarray]]] = {}
    for i, (element, name) in enumerate(zip(mesh.elements, names)):
        if len(element) != 4:
            continue
        by_region.setdefault(str(name), []).append((i, nodes[list(element)]))
    # shared node pairs that are mesh edges (exclude ridge)
    edge_keys = set(_edge_index(mesh))

    def _tris(quad: np.ndarray):
        return ((quad[0], quad[1], quad[2]), (quad[0], quad[2], quad[3]))

    def _orient(p, q, r, s):
        return float(np.dot(np.cross(q - p, r - p), s - p))

    def _cross(t1, t2) -> bool:
        a, b, c = t1
        d, e, f = t2
        o1, o2, o3 = _orient(d, e, f, a), _orient(d, e, f, b), _orient(d, e, f, c)
        if (o1 > gap and o2 > gap and o3 > gap) or (o1 < -gap and o2 < -gap and o3 < -gap):
            return False
        p1, p2, p3 = _orient(a, b, c, d), _orient(a, b, c, e), _orient(a, b, c, f)
        if (p1 > gap and p2 > gap and p3 > gap) or (p1 < -gap and p2 < -gap and p3 < -gap):
            return False
        # require mixed signs on both planes (proper pierce), not coplanar overlap
        if o1 * o2 > 0 and o1 * o3 > 0:
            return False
        if p1 * p2 > 0 and p1 * p3 > 0:
            return False
        return abs(o1) + abs(o2) + abs(o3) > gap and abs(p1) + abs(p2) + abs(p3) > gap

    hits = []
    for left_name, right_name in region_pairs:
        left = by_region.get(left_name, [])
        right = by_region.get(right_name, [])
        if not left or not right:
            continue
        bins: dict[int, list[tuple[int, np.ndarray]]] = {}
        for item in right:
            zc = 0.5 * (float(item[1][:, 2].min()) + float(item[1][:, 2].max()))
            bins.setdefault(int(zc * 4), []).append(item)
        seen: set[tuple[int, int]] = set()
        for ei, quad in left:
            zc = 0.5 * (float(quad[:, 2].min()) + float(quad[:, 2].max()))
            cand = []
            for key in (int(zc * 4) - 1, int(zc * 4), int(zc * 4) + 1):
                cand.extend(bins.get(key, []))
            for ej, other in cand:
                # skip if quads share an edge (ridge / seam)
                ids_a = set(int(k) for k in mesh.elements[ei])
                ids_b = set(int(k) for k in mesh.elements[ej])
                shared = ids_a & ids_b
                if len(shared) >= 2:
                    pair = tuple(sorted(shared))
                    if len(pair) == 2 and ((pair[0], pair[1]) in edge_keys):
                        continue
                if (quad[:, 0].max() < other[:, 0].min() - gap
                        or other[:, 0].max() < quad[:, 0].min() - gap
                        or quad[:, 1].max() < other[:, 1].min() - gap
                        or other[:, 1].max() < quad[:, 1].min() - gap):
                    continue
                crossed = False
                for t1 in _tris(quad):
                    for t2 in _tris(other):
                        if _cross(t1, t2):
                            crossed = True
                            break
                    if crossed:
                        break
                if not crossed:
                    continue
                key = (min(ei, ej), max(ei, ej))
                if key in seen:
                    continue
                seen.add(key)
                meta = _meta(mesh, ei)
                hits.append(_hit(
                    "sheet_intersection", element=ei,
                    region=f"{left_name}|{right_name}",
                    bay=meta["bay"], constructor=meta["constructor"],
                    z0=meta["z0"], z1=meta["z1"],
                    detail=f"elements={ei},{ej}",
                ))
    return tuple(hits)


def health_counts(hits: Sequence[HealthHit]) -> dict[str, int]:
    """Count hits by registry code (missing codes omitted)."""
    return dict(Counter(hit.code for hit in hits))


def health_counts_table(hits: Sequence[HealthHit]) -> list[dict]:
    """Rows for ``health_counts.csv``: code, severity, layer, count, inspect."""
    counts = health_counts(hits)
    rows = []
    for entry in HEALTH_CLASSES:
        n = int(counts.get(entry.code, 0))
        if n == 0:
            continue
        rows.append({
            "code": entry.code,
            "severity": entry.severity,
            "layer": entry.layer,
            "count": n,
            "inspect": entry.inspect,
        })
    return rows


def hits_to_rows(hits: Sequence[HealthHit]) -> list[dict]:
    """Dict rows for ``failures.csv`` (joinable ids)."""
    return [{
        "code": hit.code,
        "severity": hit.severity,
        "element": "" if hit.element < 0 else hit.element,
        "node": "" if hit.node < 0 else hit.node,
        "edge": "" if hit.edge < 0 else hit.edge,
        "region": hit.region,
        "bay": hit.bay,
        "constructor": hit.constructor,
        "z0": hit.z0,
        "z1": hit.z1,
        "measure": hit.measure,
        "threshold": hit.threshold,
        "detail": hit.detail,
    } for hit in hits]


def element_health_codes(hits: Sequence[HealthHit], n_elements: int) -> np.ndarray:
    """Per-element registry index (first hit), or ``len(HEALTH_CODES)`` when clean."""
    out = np.full(n_elements, len(HEALTH_CODES), np.int64)
    lookup = {code: i for i, code in enumerate(HEALTH_CODES)}
    for hit in hits:
        if hit.element < 0 or hit.element >= n_elements:
            continue
        if out[hit.element] == len(HEALTH_CODES):
            out[hit.element] = lookup.get(hit.code, len(HEALTH_CODES))
    return out


def edge_health_codes(hits: Sequence[HealthHit], n_edges: int) -> np.ndarray:
    """Per-edge registry index (first hit), or ``len(HEALTH_CODES)`` when clean."""
    out = np.full(n_edges, len(HEALTH_CODES), np.int64)
    lookup = {code: i for i, code in enumerate(HEALTH_CODES)}
    for hit in hits:
        if hit.edge < 0 or hit.edge >= n_edges:
            continue
        if out[hit.edge] == len(HEALTH_CODES):
            out[hit.edge] = lookup.get(hit.code, len(HEALTH_CODES))
    return out


def top_hits(hits: Sequence[HealthHit], n: int = 5) -> list[HealthHit]:
    """Worst-first sample: fail, then warn, then by |measure| when present."""
    rank = {"fail": 0, "warn": 1, "info": 2}

    def key(hit: HealthHit):
        measure = abs(hit.measure) if hit.measure == hit.measure else -1.0
        return (rank.get(hit.severity, 9), -measure, hit.code, hit.element)

    return sorted(hits, key=key)[:n]


def diff_health_counts(before: Mapping[str, int], after: Mapping[str, int]) -> list[dict]:
    """Per-code count deltas between two studies (before → after)."""
    codes = sorted(set(before) | set(after) | set(HEALTH_CODES))
    rows = []
    for code in codes:
        b, a = int(before.get(code, 0)), int(after.get(code, 0))
        if b == 0 and a == 0:
            continue
        rows.append({"code": code, "before": b, "after": a, "delta": a - b})
    return rows


def diff_failure_elements(before_rows: Iterable[dict], after_rows: Iterable[dict]
                          ) -> dict[str, list[int]]:
    """Element ids that appeared or disappeared between two ``failures.csv`` dumps."""
    def elements(rows):
        out: set[int] = set()
        for row in rows:
            value = row.get("element", "")
            if value == "" or value is None:
                continue
            out.add(int(value))
        return out

    b, a = elements(before_rows), elements(after_rows)
    return {"new": sorted(a - b), "gone": sorted(b - a)}


def read_health_counts_csv(path: Path | str) -> dict[str, int]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["code"]: int(row["count"]) for row in csv.DictReader(handle)}


def read_failures_csv(path: Path | str) -> list[dict]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def diff_study_folders(before: Path | str, after: Path | str) -> dict:
    """Compare two example-15 style dump folders (health_counts + failures)."""
    before, after = Path(before), Path(after)
    counts = diff_health_counts(
        read_health_counts_csv(before / "health_counts.csv"),
        read_health_counts_csv(after / "health_counts.csv"),
    )
    elements = diff_failure_elements(
        read_failures_csv(before / "failures.csv"),
        read_failures_csv(after / "failures.csv"),
    )
    return {"counts": counts, "elements": elements}


def format_health_summary(hits: Sequence[HealthHit], *, top: int = 5) -> list[str]:
    """Extra summary lines: counts by code and top-N location pointers."""
    lines = []
    table = health_counts_table(hits)
    total = sum(int(row["count"]) for row in table)
    lines.append(f"health hits: {total}  (failures.csv / health_counts.csv)")
    for row in table:
        lines.append(
            f"  {row['code']} [{row['severity']}]: {row['count']}  — {row['inspect']}"
        )
    for hit in top_hits(hits, top):
        where = []
        if hit.element >= 0:
            where.append(f"element {hit.element}")
        if hit.edge >= 0:
            where.append(f"edge {hit.edge}")
        if hit.node >= 0:
            where.append(f"node {hit.node}")
        if hit.region:
            where.append(hit.region)
        if hit.z0 or hit.z1:
            where.append(f"z {hit.z0:.4g}–{hit.z1:.4g}")
        lines.append(f"top: {hit.code} " + " ".join(where))
    return lines


__all__ = [
    "HEALTH_CLASSES", "HEALTH_CODES", "HEALTH_LAYERS", "HEALTH_SCOPES", "HEALTH_SEVERITIES",
    "MISSING_ID", "HealthClass", "HealthHit", "HealthThresholds",
    "detect_edge_over_target", "detect_mid_surface_fold", "detect_ridge_mismatch",
    "detect_sheet_intersection", "diff_failure_elements", "diff_health_counts",
    "diff_study_folders", "edge_health_codes", "element_health_codes",
    "format_health_summary", "health_class", "health_counts", "health_counts_table",
    "hits_to_rows", "read_failures_csv", "read_health_counts_csv", "run_mesh_health",
    "top_hits",
]
