"""JSON-serializable debug payload for void-cell wall loop demos."""

from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from blade_structure_utilities.geometry_utilities.section_contract.cell_wall_loop import CellWallLoop
from blade_structure_utilities.geometry_utilities.section_contract.shell_midline_strip import ShellMidlineStrip
from blade_structure_utilities.geometry_utilities.section_geometry.cell_wall_loops.cell_wall_loop_builders import (
    build_cell_wall_loops,
)


def _polyline_to_list(pts: NDArray[np.float64]) -> list[list[float]]:
    arr = np.asarray(pts, dtype=float)
    return [[float(x), float(y)] for x, y in arr]


def cell_wall_loops_case_payload(
    strips: Sequence[ShellMidlineStrip],
    *,
    name: str = "",
    loops: Sequence[CellWallLoop] | None = None,
) -> dict:
    """JSON-serializable payload for one midline â†’ cell-loop case (debug / examples)."""
    strip_list = list(strips)
    cell_loops = list(loops) if loops is not None else list(build_cell_wall_loops(strip_list))
    n_webs = sum(1 for s in strip_list if s.kind == "web")
    return {
        "name": str(name) if name else f"n_webs_{n_webs}",
        "n_webs": int(n_webs),
        "n_loops": len(cell_loops),
        "strips": [
            {
                "label": str(s.label),
                "kind": str(s.kind),
                "thickness_m": float(s.thickness_m),
                "closed": bool(s.closed),
                "midline_twisted_pitch_axis_blade_frame": _polyline_to_list(s.midline_twisted_pitch_axis_blade_frame),
            }
            for s in strip_list
            if s.kind in ("skin", "web")  # caps omitted from debug payload
        ],
        "loops": [
            {
                "id": int(lp.id),
                "area_m2": float(lp.area_m2),
                "web_labels": list(lp.web_labels),
                "piece_labels": list(lp.piece_labels),
                "wall_labels": list(lp.wall_labels),
                "polyline_twisted_pitch_axis_blade_frame": _polyline_to_list(lp.polyline_twisted_pitch_axis_blade_frame),
            }
            for lp in cell_loops
        ],
    }


__all__ = ["cell_wall_loops_case_payload"]

