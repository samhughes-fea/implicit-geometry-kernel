"""Marching squares on an ``SDFNode`` at fixed ``z``, then Newton onto ``phi = 0``."""

from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.sdf.nodes.base import SDFNode
from implicit_geometry_kernel.sdf.solve import GeometrySolver

_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0))
_CORNERS = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])


def _interp(p0: np.ndarray, p1: np.ndarray, v0: float, v1: float) -> np.ndarray:
    t = v0 / (v0 - v1 + 1e-30)
    return p0 + t * (p1 - p0)


def extract_zero_loops(field: SDFNode, z: float, bounds: tuple[float, float, float, float],
                       n: int = 48, solver: GeometrySolver | None = None) -> list[np.ndarray]:
    """Closed polylines ``(M, 3)`` of ``phi = 0`` in ``(xmin, xmax, ymin, ymax)`` at ``z``."""
    xmin, xmax, ymin, ymax = bounds
    xs = np.linspace(xmin, xmax, n)
    ys = np.linspace(ymin, ymax, n)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, float(z))])
    phi = field.evaluate(pts).phi.reshape(n, n)
    segs = []
    for j in range(n - 1):
        for i in range(n - 1):
            corners = np.array([[xs[i], ys[j]], [xs[i + 1], ys[j]],
                                [xs[i + 1], ys[j + 1]], [xs[i], ys[j + 1]]])
            vals = np.array([phi[j, i], phi[j, i + 1], phi[j + 1, i + 1], phi[j + 1, i]])
            inside = vals < 0.0
            if inside.all() or not inside.any():
                continue
            hits = {}
            for k, (a, b) in enumerate(_EDGES):
                if inside[a] != inside[b]:
                    hits[k] = _interp(corners[a], corners[b], vals[a], vals[b])
            if len(hits) == 2:
                first, second = hits.values()
                segs.append((first, second))
            elif len(hits) == 4:
                # Saddle cell: two segments. The asymptotic decider (R4.31) gives the
                # sign of the bilinear interpolant at its saddle point; if it has the
                # sign of corners 0 and 2, those corners connect and the segments cut
                # off corners 1 and 3, otherwise they cut off corners 0 and 2.
                v0, v1, v2, v3 = vals
                denom = v0 + v2 - v1 - v3
                saddle = (v0 * v2 - v1 * v3) / denom if denom != 0.0 else 0.25 * float(np.sum(vals))
                if (saddle < 0.0) == inside[0]:
                    pairs = ((0, 1), (2, 3))
                else:
                    pairs = ((3, 0), (1, 2))
                for p, q in pairs:
                    segs.append((hits[p], hits[q]))
    if not segs:
        return []
    loops = _chain(segs)
    solver = solver or GeometrySolver(tolerance=1e-8, max_step=min(xmax - xmin, ymax - ymin) / n)
    out = []
    for loop in loops:
        xyz = np.column_stack([loop, np.full(len(loop), float(z))])
        snapped = solver.project(field, xyz, plane=True).points
        if np.linalg.norm(snapped[0, :2] - snapped[-1, :2]) > 1e-9:
            snapped = np.vstack([snapped, snapped[0]])
        out.append(snapped)
    return out


def _chain(segs: list) -> list[np.ndarray]:
    unused = [ (np.asarray(a, float), np.asarray(b, float)) for a, b in segs ]
    loops = []
    while unused:
        a, b = unused.pop()
        chain = [a, b]
        moved = True
        while moved:
            moved = False
            tip = chain[-1]
            for k, (p, q) in enumerate(unused):
                if np.linalg.norm(p - tip) < 1e-9:
                    chain.append(q)
                    unused.pop(k)
                    moved = True
                    break
                if np.linalg.norm(q - tip) < 1e-9:
                    chain.append(p)
                    unused.pop(k)
                    moved = True
                    break
        loops.append(np.asarray(chain))
    loops.sort(key=len, reverse=True)
    return loops
