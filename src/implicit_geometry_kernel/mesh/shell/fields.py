"""Per-element (t, z) coordinates, thickness fields and per-node thickness on a shell mesh.

Decision 0014 makes thickness a field ``t = f(s, z)`` on each partition, with
``t`` the normalised arc length across the partition. These helpers place every
skin and web element on its partition's (t, z) patch, sample
:class:`~implicit_geometry_kernel.structure.ThicknessField` objects there, and
give each node a thickness. Properties only: node positions never change
(decision 0016).
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from implicit_geometry_kernel.blade.skeleton import run_t


def _project_t(points: np.ndarray, xy: np.ndarray) -> float:
    """Normalised arc length of the closest point on a polyline to ``xy``."""
    pts = np.asarray(points, float)[:, :2]
    if len(pts) < 2:
        return 0.0
    tt = run_t(pts)
    a, b = pts[:-1], pts[1:]
    ab = b - a
    denom = np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)
    w = np.clip(np.einsum("ij,ij->i", xy[None, :] - a, ab) / denom, 0.0, 1.0)
    foot = a + w[:, None] * ab
    k = int(np.argmin(np.linalg.norm(foot - xy[None, :], axis=1)))
    return float(tt[k] + w[k] * (tt[k + 1] - tt[k]))


def element_tz(mesh, geo) -> tuple[np.ndarray, np.ndarray]:
    """``(t, z)`` of each element's centroid on its partition's patch.

    ``t`` is NaN for rib and tip-closure elements, which have no (t, z) patch.
    ``z`` is the centroid's span position for every element.
    """
    nodes = np.asarray(mesh.nodes, float)
    n = len(mesh.elements)
    t = np.full(n, np.nan)
    z = np.empty(n)
    cache: dict = {}
    for i, conn in enumerate(mesh.elements):
        c = nodes[list(conn)].mean(axis=0)
        z[i] = float(c[2])
        if mesh.surface[i] not in ("skin", "web"):
            continue
        key = round(float(c[2]), 9)
        graph = cache.get(key)
        if graph is None:
            graph = geo.section_at(float(c[2])).shell_graph()
            cache[key] = graph
        runs = {r.name: r for r in graph.loop + graph.webs}
        run = runs.get(mesh.partition[i])
        if run is not None:
            t[i] = _project_t(np.asarray(run.points, float), c[:2])
    return t, z


def sample_thickness(mesh, geo, fields: Mapping[str, Any]) -> np.ndarray:
    """Element thickness with ``fields`` (``{partition name: ThicknessField}``) applied.

    Elements on a partition in ``fields`` take the field at their centroid's
    (t, z); every other element keeps :attr:`ShellMesh.thickness`.
    """
    out = np.array(mesh.thickness, float, copy=True)
    if not fields:
        return out
    t, z = element_tz(mesh, geo)
    for i, part in enumerate(mesh.partition):
        field = fields.get(part)
        if field is not None and np.isfinite(t[i]):
            out[i] = float(np.asarray(field(t[i], z[i])))
    if np.any(out < 0.0):
        raise ValueError("a thickness field returned a negative thickness")
    return out


def node_thickness(mesh, thickness: np.ndarray | None = None, *, partition: str | None = None) -> np.ndarray:
    """Area-weighted mean thickness of the elements around each node, metres.

    ``thickness`` defaults to :attr:`ShellMesh.thickness`. With ``partition``,
    only that partition's elements count, so a node on a partition boundary gets
    that side's value (no averaging across the step); nodes off the partition
    are NaN.
    """
    th = np.asarray(mesh.thickness if thickness is None else thickness, float)
    nodes = np.asarray(mesh.nodes, float)
    total = np.zeros(len(nodes))
    weight = np.zeros(len(nodes))
    for i, conn in enumerate(mesh.elements):
        if partition is not None and mesh.partition[i] != partition:
            continue
        pts = nodes[list(conn)]
        c = pts.mean(axis=0)
        area = 0.0
        for a, b in zip(pts, np.roll(pts, -1, axis=0)):
            area += 0.5 * float(np.linalg.norm(np.cross(a - c, b - c)))
        for j in conn:
            total[j] += area * th[i]
            weight[j] += area
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(weight > 0.0, total / weight, np.nan)


def element_curvature(mesh, geo, *, n_t: int = 33, curvature_tol: float = 1e-6) -> dict[str, np.ndarray]:
    """Principal curvatures, direction and kind at each element centroid (decision 0017).

    Skin and web elements are evaluated on a :class:`~implicit_geometry_kernel.blade.patch.PartitionPatch`
    of their partition, built over the mesh stations where the partition exists
    and interpolated along the span like the blade (``span_interp``). Rib and
    tip-closure elements lie in a plane of constant z, so they are ``FLAT`` with
    exactly zero curvature by construction.

    Returns
    -------
    dict
        ``k1``, ``k2`` (1/m, ``k1 <= k2``), ``theta`` (radians, from the
        partition's t direction to the first principal direction) and ``kind``
        (``FLAT``, ``SINGLY``, ``DOUBLY``, ``UMBILIC``; ``UNSET`` only where a
        patch could not be built).
    """
    from implicit_geometry_kernel.blade.patch import PartitionPatch

    n = len(mesh.elements)
    k1 = np.zeros(n)
    k2 = np.zeros(n)
    theta = np.zeros(n)
    kind = np.array(["FLAT"] * n, dtype=object)
    t, z = element_tz(mesh, geo)
    stations = sorted(float(v) for v in mesh.stations)
    patches: dict = {}
    for i, part in enumerate(mesh.partition):
        if mesh.surface[i] not in ("skin", "web"):
            continue
        patch = patches.get(part, False)
        if patch is False:
            present = []
            for zs in stations:
                graph = geo.section_at(zs).shell_graph()
                if any(r.name == part for r in graph.loop + graph.webs):
                    present.append(zs)
            try:
                patch = PartitionPatch(geo, part, present, n_t=n_t, sampling="flag",
                                       curvature_tol=curvature_tol)
            except ValueError:
                patch = None
            patches[part] = patch
        if patch is None or not np.isfinite(t[i]):
            k1[i] = k2[i] = theta[i] = np.nan
            kind[i] = "UNSET"
            continue
        sp = patch.at(float(t[i]), float(z[i]))
        k1[i], k2[i], theta[i], kind[i] = sp.k1, sp.k2, sp.theta, sp.kind
    return {"k1": k1, "k2": k2, "theta": theta, "kind": kind}


__all__ = ["element_curvature", "element_tz", "node_thickness", "sample_thickness"]
