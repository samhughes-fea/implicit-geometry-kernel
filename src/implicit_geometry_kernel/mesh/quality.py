"""Aspect ratio and Jacobian of a shell mesh, one row per partition run."""

from __future__ import annotations

from dataclasses import dataclass
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


def mesh_quality(
    nodes: np.ndarray,
    elements: Sequence[Sequence[int]],
    partition: Sequence[str],
) -> tuple[MeshRunQuality, ...]:
    """One :class:`MeshRunQuality` per partition that has elements."""
    nodes = np.asarray(nodes, float).reshape(-1, 3)
    groups: dict[str, list[int]] = {}
    for i, name in enumerate(partition):
        groups.setdefault(str(name), []).append(i)
    rows = []
    for name in sorted(groups):
        aspects, jmins, jmaxs = [], [], []
        inverted = 0
        for i in groups[name]:
            p = nodes[list(elements[i])]
            aspect, jac = _element_metrics(p)
            aspects.append(aspect)
            lo, hi = float(np.min(jac)), float(np.max(jac))
            jmins.append(lo)
            jmaxs.append(hi)
            if lo < 0.0:
                inverted += 1
        rows.append(
            MeshRunQuality(
                name,
                len(groups[name]),
                float(np.min(aspects)),
                float(np.max(aspects)),
                float(np.min(jmins)),
                float(np.max(jmaxs)),
                inverted,
            )
        )
    return tuple(rows)


__all__ = ["MeshRunQuality", "mesh_quality"]
