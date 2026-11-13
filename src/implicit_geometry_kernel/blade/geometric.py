"""Geometric section integrals at one station.

Area, centroid and second moments of the wall mid-lines. A spar cap is
integrated on its own mid-line, so its offset from the skin is in the
moments. This is geometry: no elastic centre, shear centre, principal
axes, mass centre or warping.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from implicit_geometry_kernel.blade.cut import SectionCut, SectionWall


@dataclass(frozen=True)
class GeometricIntegrals:
    """Wall area and second moments in the section frame (metres).

    ``ixx``, ``iyy`` and ``ixy`` are about the centroid. The ``*_origin``
    moments are about the pitch-axis origin, so a cap offset of ``e``
    contributes ``area * e**2`` there.
    """

    area: float
    centroid: tuple[float, float]
    ixx: float
    iyy: float
    ixy: float
    ixx_origin: float
    iyy_origin: float
    ixy_origin: float


def geometric_integrals(source: SectionCut | list[SectionWall]) -> GeometricIntegrals:
    """Integrate every wall of a cut, caps on their own mid-line."""
    walls = source.walls if isinstance(source, SectionCut) else list(source)
    area = 0.0
    sx = sy = 0.0
    ixx_o = iyy_o = ixy_o = 0.0
    for wall in walls:
        xy = np.asarray(wall.xy, dtype=float)
        if xy.ndim != 2 or len(xy) < 2:
            continue
        t = float(wall.thickness)
        if t == 0.0:
            continue
        p0 = xy[:-1, :2]
        p1 = xy[1:, :2]
        d = p1 - p0
        length = np.hypot(d[:, 0], d[:, 1])
        live = length > 0.0
        if not np.any(live):
            continue
        p0, p1, length = p0[live], p1[live], length[live]
        da = t * length
        area += float(da.sum())
        sx += float((0.5 * (p0[:, 0] + p1[:, 0]) * da).sum())
        sy += float((0.5 * (p0[:, 1] + p1[:, 1]) * da).sum())
        iyy_o += float((t * length * (p0[:, 0] ** 2 + p0[:, 0] * p1[:, 0] + p1[:, 0] ** 2) / 3.0).sum())
        ixx_o += float((t * length * (p0[:, 1] ** 2 + p0[:, 1] * p1[:, 1] + p1[:, 1] ** 2) / 3.0).sum())
        ixy_o += float((t * length * (
            2.0 * p0[:, 0] * p0[:, 1] + p0[:, 0] * p1[:, 1] + p1[:, 0] * p0[:, 1] + 2.0 * p1[:, 0] * p1[:, 1]
        ) / 6.0).sum())
    if area == 0.0:
        return GeometricIntegrals(0.0, (0.0, 0.0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    xc, yc = sx / area, sy / area
    return GeometricIntegrals(
        area=area,
        centroid=(xc, yc),
        ixx=ixx_o - area * yc * yc,
        iyy=iyy_o - area * xc * xc,
        ixy=ixy_o - area * xc * yc,
        ixx_origin=ixx_o,
        iyy_origin=iyy_o,
        ixy_origin=ixy_o,
    )
