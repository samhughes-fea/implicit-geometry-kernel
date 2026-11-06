"""Marching-squares zero contour on quadtree leaves, and mesh summaries."""

from __future__ import annotations



import numpy as np

from implicit_geometry_kernel.integration.contour import ensure_closed, smooth_contour
from implicit_geometry_kernel.integration.quadtree.tree import QuadCell, QuadTreeMesh


def _interpolate_zero_crossing(
    p1: np.ndarray,
    p2: np.ndarray,
    v1: float,
    v2: float,
) -> np.ndarray:
    if abs(v2 - v1) < 1e-15:
        return 0.5 * (p1 + p2)
    t = -v1 / (v2 - v1)
    return p1 + t * (p2 - p1)


def _connect_segments_to_contour(
    segments: list[tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    if not segments:
        raise ValueError("No zero contour found in quadtree leaf cells.")
    chains: list[list[np.ndarray]] = []
    tol = 1e-9
    for a, b in segments:
        placed = False
        for chain in chains:
            if np.allclose(chain[-1], a, atol=tol, rtol=0.0):
                chain.append(b)
                placed = True
                break
            if np.allclose(chain[-1], b, atol=tol, rtol=0.0):
                chain.append(a)
                placed = True
                break
            if np.allclose(chain[0], b, atol=tol, rtol=0.0):
                chain.insert(0, a)
                placed = True
                break
            if np.allclose(chain[0], a, atol=tol, rtol=0.0):
                chain.insert(0, b)
                placed = True
                break
        if not placed:
            chains.append([a, b])
    merged = True
    while merged and len(chains) > 1:
        merged = False
        for i in range(len(chains)):
            if merged:
                break
            for j in range(i + 1, len(chains)):
                c1, c2 = chains[i], chains[j]
                if np.allclose(c1[-1], c2[0], atol=tol, rtol=0.0):
                    chains[i] = c1 + c2[1:]
                elif np.allclose(c1[-1], c2[-1], atol=tol, rtol=0.0):
                    chains[i] = c1 + c2[-2::-1]
                elif np.allclose(c1[0], c2[-1], atol=tol, rtol=0.0):
                    chains[i] = c2[:-1] + c1
                elif np.allclose(c1[0], c2[0], atol=tol, rtol=0.0):
                    chains[i] = c2[::-1][:-1] + c1
                else:
                    continue
                del chains[j]
                merged = True
                break
    longest = max(chains, key=len)
    return ensure_closed(np.asarray(longest, dtype=float))


def contour_from_quadtree_marching_squares(
    leaves: list[QuadCell],
    passes: int = 2,
) -> np.ndarray:
    """Extract the zero-level contour from quadtree leaf cells."""
    edge_pairs = [(0, 1), (1, 2), (2, 3), (3, 0)]
    leaves_with_phi = [c for c in leaves if c.phi_corners is not None]
    if not leaves_with_phi:
        raise ValueError("No leaf cells with phi_corners found.")
    phi_array = np.array([c.phi_corners for c in leaves_with_phi], dtype=float)
    sign_change_mask = (phi_array.min(axis=1) <= 0.0) & (phi_array.max(axis=1) >= 0.0)
    active = [c for c, active in zip(leaves_with_phi, sign_change_mask) if active]
    segments = []
    for cell in active:
        phi = cell.phi_corners
        corners = cell.corners()
        crossings = []
        for i, j in edge_pairs:
            if (phi[i] <= 0.0) != (phi[j] <= 0.0):  # type: ignore
                crossings.append(
                    _interpolate_zero_crossing(corners[i], corners[j], phi[i], phi[j])  # type: ignore
                )
        if len(crossings) == 2:
            segments.append((crossings[0], crossings[1]))
        elif len(crossings) == 4:
            segments.append((crossings[0], crossings[1]))
            segments.append((crossings[2], crossings[3]))
    contour = _connect_segments_to_contour(segments)
    return smooth_contour(contour, passes=passes)


def summarize_quadtree(mesh: QuadTreeMesh) -> dict[str, object]:
    """Summary statistics for a refined quadtree mesh."""
    if getattr(mesh, "_is_soa", False):
        bounds = mesh._leaf_bounds  # type: ignore[attr-defined]
        levels_arr = mesh._leaf_levels  # type: ignore[attr-defined]
        widths = bounds[:, 1] - bounds[:, 0] if bounds.size else np.empty(0)
        uniq, counts = np.unique(levels_arr, return_counts=True)
        level_counts = {int(k): int(v) for k, v in zip(uniq, counts)}
        return {
            "leaf_count": int(bounds.shape[0]),
            "level_counts": dict(sorted(level_counts.items())),
            "max_level": int(np.max(levels_arr)) if levels_arr.size else 0,
            "min_cell_width": float(np.min(widths)) if widths.size else 0.0,
            "max_cell_width": float(np.max(widths)) if widths.size else 0.0,
            "refine_stats": dict(mesh.refine_stats),
        }
    leaves = mesh.leaf_cells()
    levels: dict[int, int] = {}
    for cell in leaves:
        levels[cell.level] = levels.get(cell.level, 0) + 1
    return {
        "leaf_count": len(leaves),
        "level_counts": dict(sorted(levels.items())),
        "max_level": max(levels) if levels else 0,
        "min_cell_width": min((c.width for c in leaves), default=0.0),
        "max_cell_width": max((c.width for c in leaves), default=0.0),
        "refine_stats": dict(mesh.refine_stats),
    }
