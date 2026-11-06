"""Spanwise feature-curve locators and trailing-edge classification."""

from __future__ import annotations

import numpy as np

from implicit_geometry_kernel.sdf.nodes.base import SDFNode
from implicit_geometry_kernel.sdf.solve import GeometrySolver


def feature_curve(a: SDFNode, b: SDFNode, zs, guess_xy, solver: GeometrySolver | None = None) -> np.ndarray:
    """Trace ``phi_a = phi_b = 0`` along ``zs`` (in-plane intersect at each station).

    Named partition edges (cap/skin keys) are ``ShellGraph.partition_edges`` /
    ``partition_edge_curve`` — graph vertices, not this Newton intersect.
    """
    solver = solver or GeometrySolver(tolerance=1e-8, max_step=0.05)
    g = np.asarray(guess_xy, float).reshape(2)
    out = []
    for z in np.asarray(zs, float):
        start = np.array([[g[0], g[1], float(z)]])
        hit = solver.intersect(a, b, start)
        p = hit.points[0]
        out.append(p)
        g = p[:2]
    return np.asarray(out)


def classify_trailing_edge(ring_xy: np.ndarray, blunt_tol: float = 1e-4) -> str:
    """``sharp`` if the TE is a point; ``blunt`` if a finite chordwise edge remains."""
    pts = np.asarray(ring_xy, float)[:, :2]
    if np.linalg.norm(pts[0] - pts[-1]) < 1e-12:
        pts = pts[:-1]
    i_te = int(np.argmax(pts[:, 0]))
    x_te = pts[i_te, 0]
    band = pts[np.abs(pts[:, 0] - x_te) <= blunt_tol]
    span = float(np.ptp(band[:, 1])) if len(band) else 0.0
    return "blunt" if span > blunt_tol else "sharp"
