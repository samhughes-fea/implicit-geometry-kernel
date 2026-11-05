"""
CellWallLoop
============
FE-agnostic DTO for one void-cell **wall mid-surface** circuit (skin arc + web
closures) produced by ``section_geometry``.

This is constructive geometry only — not an SDF medial axis of the cavity, not
a sandwich-core bay clip, and not a thin-wall shear graph.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class CellWallLoop:
    """Closed wall mid-surface loop for one multicell void bay.

    Parameters
    ----------
    id
        Cell index LE → TE (``0 … n_webs``).
    polyline_twisted_pitch_axis_blade_frame
        ``(N, 2)`` closed polyline in ``twisted_pitch_axis_blade_frame`` (**B**)
        (first point repeated at the end).
    area_m2
        Absolute shoelace area of the polyline [m²].
    web_labels
        Labels of webs that bound this cell (empty when ``n_webs == 0``).
    piece_labels
        Ordered role tags along the loop construction (e.g. ``skin``,
        ``shear_web_0``).
    wall_labels
        The midline strip labels that bound this cell, in loop order (e.g.
        ``skin_suction_surface:bay_1``, ``shear_web_2``,
        ``skin_pressure_surface:bay_1``, ``shear_web_1``). This **is** the
        wall-to-cell membership: a strip label appearing in two loops is a wall
        shared by those cells. Consumers read membership from it and never
        re-derive it from coordinates (SG-014).
    """

    id: int
    polyline_twisted_pitch_axis_blade_frame: NDArray[np.float64]
    area_m2: float
    web_labels: tuple[str, ...] = ()
    piece_labels: tuple[str, ...] = ()
    wall_labels: tuple[str, ...] = ()

    def is_closed(self, *, tol: float = 1e-9) -> bool:
        pts = np.asarray(self.polyline_twisted_pitch_axis_blade_frame, dtype=float)
        if pts.shape[0] < 3:
            return False
        return float(np.hypot(*(pts[0] - pts[-1]))) <= tol


__all__ = ["CellWallLoop"]
