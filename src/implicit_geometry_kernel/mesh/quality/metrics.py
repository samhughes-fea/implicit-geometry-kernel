"""Aspect ratio and Jacobian of a shell mesh, one row per partition run."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class MeshRunQuality:
    """Quality of the elements that lie on one partition (a shell-graph run).

    ``aspect`` is longest edge over shortest edge (1 is equal-sided).
    ``jacobian`` is the scaled isoparametric Jacobian at the corners: 1 is a
    rectangle, values in (0, 1] are valid, and a negative value is inverted.
    """

    partition: str
    n: int
    aspect_min: float
    aspect_max: float
    jacobian_min: float
    jacobian_max: float
    inverted: int
    skew_min: float = 0.0
    skew_max: float = 0.0
    angle_min: float = 0.0
    angle_max: float = 0.0
    volume_ratio_max: float = 1.0
    warp_max: float = 0.0


def _plane_xy(p: np.ndarray) -> np.ndarray:
    origin = p[0]
    normal = np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)
    length = float(np.linalg.norm(normal))
    if length <= 1e-18:
        return np.zeros((len(p), 2))
    normal = normal / length
    helper = np.array([0.0, 0.0, 1.0]) if abs(normal[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    axis = np.cross(helper, normal)
    axis_len = float(np.linalg.norm(axis))
    if axis_len <= 1e-18:
        return np.zeros((len(p), 2))
    axis = axis / axis_len
    binormal = np.cross(normal, axis)
    return np.column_stack([(p - origin) @ axis, (p - origin) @ binormal])


def _quad_jacobians(p: np.ndarray) -> np.ndarray:
    """Signed scaled Jacobian at the four corners of a bilinear quad."""
    xy = _plane_xy(p)
    corners = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))
    dets = np.empty(4, float)
    for i, (xi, eta) in enumerate(corners):
        d_xi = 0.25 * np.array([-(1.0 - eta), (1.0 - eta), (1.0 + eta), -(1.0 + eta)])
        d_eta = 0.25 * np.array([-(1.0 - xi), -(1.0 + xi), (1.0 + xi), (1.0 - xi)])
        a, b = d_xi @ xy, d_eta @ xy
        dets[i] = float(a[0] * b[1] - a[1] * b[0])
    peak = float(np.max(np.abs(dets)))
    return dets / peak if peak > 0.0 else dets


def _polygon_jacobians(p: np.ndarray) -> np.ndarray:
    """Overall signed area of a polygon, scaled to ±1 (concave fans are not inversions)."""
    xy = _plane_xy(p)
    nxt = np.roll(xy, -1, axis=0)
    area = 0.5 * float(np.sum(xy[:, 0] * nxt[:, 1] - xy[:, 1] * nxt[:, 0]))
    if abs(area) <= 1e-18:
        return np.zeros(1)
    return np.array([1.0 if area > 0.0 else -1.0])


def _element_metrics(p: np.ndarray) -> tuple[float, np.ndarray]:
    lengths = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    short = float(np.min(lengths))
    aspect = float(np.max(lengths) / short) if short > 0.0 else float("inf")
    jac = _quad_jacobians(p) if len(p) == 4 else _polygon_jacobians(p)
    return aspect, jac


QUALITY_DTYPE = np.dtype([
    ("aspect", float),
    ("skew", float),
    ("angle_min", float),
    ("angle_max", float),
    ("jacobian_min", float),
    ("area", float),
    ("shell_volume", float),
    ("volume_ratio", float),
    ("warp", float),
])
"""One row per element. ``shell_volume`` is area times thickness, not a solid volume."""


def _corner_angles(p: np.ndarray) -> np.ndarray:
    n = len(p)
    angles = np.empty(n, float)
    for i in range(n):
        a = p[(i - 1) % n] - p[i]
        b = p[(i + 1) % n] - p[i]
        la, lb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
        if la <= 1e-15 or lb <= 1e-15:
            angles[i] = 0.0
            continue
        angles[i] = np.degrees(np.arccos(np.clip(float(a @ b) / (la * lb), -1.0, 1.0)))
    return angles


def _equiangle_skew(angles: np.ndarray) -> float:
    """Equiangle skew against the ideal corner of an n-gon. 0 is ideal, 1 is collapsed."""
    n = len(angles)
    if n < 3:
        return 1.0
    ideal = ((n - 2) / n) * 180.0
    room = 180.0 - ideal
    if ideal <= 1e-9 or room <= 1e-9:
        return 1.0
    worst = 0.0
    for theta in angles:
        if theta <= ideal:
            skew = (ideal - float(theta)) / ideal
        else:
            skew = (float(theta) - ideal) / room
        if skew > worst:
            worst = skew
    return float(min(max(worst, 0.0), 1.0))


def _cell_area(p: np.ndarray) -> float:
    normal = np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)
    return 0.5 * float(np.linalg.norm(normal))


def _warp(p: np.ndarray) -> float:
    """Distance between the diagonal midpoints of a quad, over its longest edge."""
    if len(p) != 4:
        return 0.0
    gap = 0.5 * (p[0] + p[2]) - 0.5 * (p[1] + p[3])
    longest = float(np.max(np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)))
    if longest <= 1e-15:
        return 0.0
    return float(np.linalg.norm(gap) / longest)


def _fill_block(table: np.ndarray, jmax: np.ndarray, area_out: np.ndarray,
                idx: np.ndarray, p: np.ndarray, quads: bool) -> None:
    """Corner angles, area, skew, Jacobian, and warp for one block of equal-sized cells."""
    prev = np.roll(p, 1, axis=1) - p
    nxt = np.roll(p, -1, axis=1) - p
    la = np.linalg.norm(prev, axis=2)
    lb = np.linalg.norm(nxt, axis=2)
    cos = np.sum(prev * nxt, axis=2) / np.maximum(la * lb, 1e-30)
    angles = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    angles = np.where((la <= 1e-15) | (lb <= 1e-15), 0.0, angles)
    n = p.shape[1]
    ideal = ((n - 2) / n) * 180.0
    room = 180.0 - ideal
    skew = np.clip(np.where(angles <= ideal, (ideal - angles) / ideal, (angles - ideal) / room), 0.0, 1.0)
    lengths = np.linalg.norm(np.roll(p, -1, axis=1) - p, axis=2)
    short = np.min(lengths, axis=1)
    aspect = np.divide(np.max(lengths, axis=1), short, out=np.full(len(p), np.inf), where=short > 0.0)
    cell = 0.5 * np.linalg.norm(np.sum(np.cross(p, np.roll(p, -1, axis=1)), axis=1), axis=1)
    origin = p[:, 0]
    normal = np.sum(np.cross(p, np.roll(p, -1, axis=1)), axis=1)
    length = np.linalg.norm(normal, axis=1)
    ok = length > 1e-18
    normal = np.divide(normal, length[:, None], out=np.zeros_like(normal), where=ok[:, None])
    helper = np.zeros_like(normal)
    helper[:, 2] = 1.0
    helper[np.abs(normal[:, 2]) >= 0.9] = (1.0, 0.0, 0.0)
    axis = np.cross(helper, normal)
    axis_len = np.linalg.norm(axis, axis=1)
    axis = np.divide(axis, axis_len[:, None], out=np.zeros_like(axis), where=(axis_len > 1e-18)[:, None])
    rel = p - origin[:, None, :]
    binormal = np.cross(normal, axis)
    xy = np.stack([
        np.sum(rel * axis[:, None, :], axis=2),
        np.sum(rel * binormal[:, None, :], axis=2),
    ], axis=2)
    xy[~(ok & (axis_len > 1e-18))] = 0.0
    if quads:
        xi = np.array([-1.0, 1.0, 1.0, -1.0])
        eta = np.array([-1.0, -1.0, 1.0, 1.0])
        d_xi = 0.25 * np.stack([-(1.0 - eta), (1.0 - eta), (1.0 + eta), -(1.0 + eta)], axis=1)
        d_eta = 0.25 * np.stack([-(1.0 - xi), -(1.0 + xi), (1.0 + xi), (1.0 - xi)], axis=1)
        a = np.einsum("cj,njd->ncd", d_xi, xy)
        b = np.einsum("cj,njd->ncd", d_eta, xy)
        jac = a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
        peak = np.max(np.abs(jac), axis=1)
        jac = np.divide(jac, peak[:, None], out=jac.copy(), where=(peak > 0.0)[:, None])
        longest = np.max(lengths, axis=1)
        gap = 0.5 * (p[:, 0] + p[:, 2]) - 0.5 * (p[:, 1] + p[:, 3])
        warp = np.where(longest <= 1e-15, 0.0, np.linalg.norm(gap, axis=1) / np.maximum(longest, 1e-30))
    else:
        rolled = np.roll(xy, -1, axis=1)
        signed = 0.5 * np.sum(xy[..., 0] * rolled[..., 1] - xy[..., 1] * rolled[..., 0], axis=1)
        jac = np.where(np.abs(signed) <= 1e-18, 0.0, np.where(signed > 0.0, 1.0, -1.0))[:, None]
        warp = np.zeros(len(p), float)
    area_out[idx] = cell
    table["aspect"][idx] = aspect
    table["skew"][idx] = skew.max(axis=1)
    table["angle_min"][idx] = np.min(angles, axis=1)
    table["angle_max"][idx] = np.max(angles, axis=1)
    table["jacobian_min"][idx] = np.min(jac, axis=1)
    table["area"][idx] = cell
    table["warp"][idx] = warp
    jmax[idx] = np.max(jac, axis=1)


def _volume_ratio(elements: Sequence[Sequence[int]], volume: np.ndarray) -> np.ndarray:
    """Largest neighbour jump of ``volume``. A cell with no neighbour stays at 1."""
    ratio = np.ones(len(elements), float)
    by_edge: dict[tuple[int, int], list[int]] = {}
    for i, element in enumerate(elements):
        ring = list(element)
        for a, b in zip(ring, ring[1:] + ring[:1]):
            by_edge.setdefault((a, b) if a < b else (b, a), []).append(i)
    for users in by_edge.values():
        if len(users) < 2:
            continue
        for i in range(len(users)):
            for j in range(i + 1, len(users)):
                a, b = users[i], users[j]
                lo, hi = float(volume[a]), float(volume[b])
                if hi < lo:
                    lo, hi = hi, lo
                jump = 1.0 if lo <= 1e-30 else hi / lo
                if jump > ratio[a]:
                    ratio[a] = jump
                if jump > ratio[b]:
                    ratio[b] = jump
    return ratio


def volume_ratio_from_edges(edges: np.ndarray, volume: np.ndarray) -> np.ndarray:
    """Largest neighbour jump of ``volume``, read from an edge table.

    An edge stores up to three element ids. A cell with no neighbour stays at 1.
    """
    ratio = np.ones(len(volume), float)
    if len(edges) == 0:
        return ratio
    users = np.column_stack([edges["e0"], edges["e1"], edges["e2"]])
    use = np.asarray(edges["use"], np.int64)
    for width in (2, 3):
        chosen = users[use >= width][:, :width]
        if len(chosen) == 0:
            continue
        values = volume[chosen]
        for left in range(width):
            for right in range(left + 1, width):
                lo = np.minimum(values[:, left], values[:, right])
                hi = np.maximum(values[:, left], values[:, right])
                jump = np.where(lo <= 1e-30, 1.0, hi / lo)
                np.maximum.at(ratio, chosen[:, left], jump)
                np.maximum.at(ratio, chosen[:, right], jump)
    return ratio


def _measure(nodes: np.ndarray, elements: Sequence[Sequence[int]], thickness: np.ndarray,
             *, with_ratio: bool = True):
    """Per-element quality table and the max scaled Jacobian of each element."""
    nodes = np.asarray(nodes, float).reshape(-1, 3)
    thickness = np.asarray(thickness, float)
    table = np.zeros(len(elements), dtype=QUALITY_DTYPE)
    table["volume_ratio"] = 1.0
    jmax = np.zeros(len(elements), float)
    area = np.zeros(len(elements), float)
    grouped: dict[int, list[int]] = {}
    for i, element in enumerate(elements):
        grouped.setdefault(len(element), []).append(i)
    for n, idx_list in grouped.items():
        idx = np.asarray(idx_list, np.int64)
        if n in (3, 4) and n > 0:
            corners = np.asarray([list(elements[i]) for i in idx_list], np.int64)
            _fill_block(table, jmax, area, idx, nodes[corners], n == 4)
            continue
        for i in idx_list:
            p = nodes[list(elements[i])]
            angles = _corner_angles(p)
            aspect, jac = _element_metrics(p)
            cell_area = _cell_area(p)
            area[i] = cell_area
            table[i] = (
                aspect,
                _equiangle_skew(angles),
                float(np.min(angles)) if len(angles) else 0.0,
                float(np.max(angles)) if len(angles) else 0.0,
                float(np.min(jac)) if len(jac) else 0.0,
                cell_area,
                0.0,
                1.0,
                _warp(p),
            )
            jmax[i] = float(np.max(jac)) if len(jac) else 0.0
    volume = area * thickness
    table["shell_volume"] = volume
    if with_ratio:
        table["volume_ratio"] = _volume_ratio(elements, volume)
    return table, jmax


def element_quality(nodes: np.ndarray, elements: Sequence[Sequence[int]],
                    thickness: Sequence[float]) -> np.ndarray:
    """One :data:`QUALITY_DTYPE` row per element."""
    table, _jmax = _measure(nodes, elements, np.asarray(thickness, float))
    return table


def quality_rollup(partition: Sequence[str], table: np.ndarray, jacobian_max: np.ndarray
                   ) -> tuple[MeshRunQuality, ...]:
    """One :class:`MeshRunQuality` per partition, read from the element table."""
    groups: dict[str, list[int]] = {}
    for i, name in enumerate(partition):
        groups.setdefault(str(name), []).append(i)
    rows = []
    for name in sorted(groups):
        idx = groups[name]
        jac = table["jacobian_min"][idx]
        rows.append(MeshRunQuality(
            name,
            len(idx),
            float(np.min(table["aspect"][idx])),
            float(np.max(table["aspect"][idx])),
            float(np.min(jac)),
            float(np.max(jacobian_max[idx])),
            int(np.sum(jac < 0.0)),
            float(np.min(table["skew"][idx])),
            float(np.max(table["skew"][idx])),
            float(np.min(table["angle_min"][idx])),
            float(np.max(table["angle_max"][idx])),
            float(np.max(table["volume_ratio"][idx])),
            float(np.max(table["warp"][idx])),
        ))
    return tuple(rows)


def mesh_quality(
    nodes: np.ndarray,
    elements: Sequence[Sequence[int]],
    partition: Sequence[str],
    thickness: Sequence[float] | None = None,
) -> tuple[MeshRunQuality, ...]:
    """One :class:`MeshRunQuality` per partition that has elements."""
    thick = np.ones(len(elements), float) if thickness is None else np.asarray(thickness, float)
    table, jmax = _measure(nodes, elements, thick)
    return quality_rollup(partition, table, jmax)


@dataclass(frozen=True)
class MeshCheck:
    """Transparent rollup a study reads off a finished shell mesh.

    Counts and extrema only — no pass/fail gate. ``hit_counts`` maps registry
    health codes to hit totals (see :mod:`implicit_geometry_kernel.mesh.health`).
    ``edge_over_s_max`` is the longest element edge divided by the recipe edge
    at that edge's station, or ``None`` when there is no target edge.
    ``boundary_in_source`` reports whether rib/tip sheet nodes sit in the source
    block; it is a fact, not a verdict.
    """

    edge_over_s_max: float | None
    quad_angle_min: float
    open_edges: int
    volume_ratio_max: float
    skew_max: float
    boundary_in_source: bool
    hit_counts: dict[str, int] = field(default_factory=dict)


def mesh_check(mesh) -> MeshCheck:
    """Read :class:`MeshCheck` from a shell mesh, quality table, and health hits."""
    from implicit_geometry_kernel.mesh.health import health_counts, run_mesh_health

    table = np.asarray(mesh.element_quality)
    quads = [i for i, e in enumerate(mesh.elements) if len(e) == 4]
    quad_angle_min = float(np.min(table["angle_min"][quads])) if quads and len(table) else float("nan")
    open_edges = sum(1 for n in mesh.edge_use().values() if n == 1)
    volume_ratio_max = float(np.max(table["volume_ratio"])) if len(table) else 1.0
    skew_max = float(np.max(table["skew"])) if len(table) else 0.0
    recipe = getattr(mesh, "recipe", None)
    edge = () if recipe is None else tuple(recipe.edge)
    edge_over_s_max = None
    if edge:
        stations = np.asarray(mesh.stations, float)
        target = np.asarray(edge, float)
        ratios = []
        nodes = np.asarray(mesh.nodes, float)
        for element in mesh.elements:
            p = nodes[list(element)]
            for a, b in zip(p, np.roll(p, -1, axis=0)):
                j = int(np.argmin(np.abs(stations - 0.5 * (a[2] + b[2]))))
                if target[j] > 0.0:
                    ratios.append(float(np.linalg.norm(a - b)) / float(target[j]))
        if ratios:
            edge_over_s_max = float(max(ratios))
    n_source = int(mesh.n_source)
    sheet: set[int] = set()
    for element, surface in zip(mesh.elements, mesh.surface):
        if surface in ("skin", "web"):
            sheet.update(int(k) for k in element)
    boundary_in_source = True
    for element, surface in zip(mesh.elements, mesh.surface):
        if surface not in ("rib", "tip_closure"):
            continue
        for k in element:
            k = int(k)
            if (k in sheet) != (k < n_source):
                boundary_in_source = False
                break
        if not boundary_in_source:
            break
    if any(k >= n_source for k in sheet):
        boundary_in_source = False
    counts = health_counts(run_mesh_health(mesh))
    return MeshCheck(edge_over_s_max, quad_angle_min, open_edges, volume_ratio_max, skew_max,
                     boundary_in_source, counts)


__all__ = [
    "MeshCheck", "MeshRunQuality", "QUALITY_DTYPE", "element_quality", "mesh_check",
    "mesh_quality", "quality_rollup", "volume_ratio_from_edges",
]
