"""Plot helpers for shell midline strip topology overlays."""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np

from blade_geometry_model.section_contract.shell_midline_strip import ShellMidlineStrip
from blade_geometry_model.coordinate_frames import (
    section_plot_axes,
)


def plot_shell_midline_strip_topology(
    strips: Sequence[ShellMidlineStrip],
    *,
    ax: Any | None = None,
    title: str = "Shell midline strip topology",
    kind_colour: dict[str, str] | None = None,
    figsize: tuple[float, float] = (8.0, 6.0),
) -> tuple[Any, Any]:
    """Overlay ``ShellMidlineStrip`` polylines coloured by strip kind."""
    import matplotlib.pyplot as plt

    if kind_colour is None:
        kind_colour = {"skin": "#2196F3", "web": "#4CAF50", "cap": "#F44336"}

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    label_seen: set[str] = set()
    for strip in strips:
        pts = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=np.float64)
        if pts.shape[0] < 2:
            continue
        colour = kind_colour.get(str(strip.kind), "#607D8B")
        legend = str(strip.label) if str(strip.label) not in label_seen else None
        label_seen.add(str(strip.label))
        # midline_twisted_pitch_axis_blade_frame is canonical B (X_B, Y_B);
        # plots are horizontal = Y_B.
        h, v = section_plot_axes(pts)
        ax.plot(h, v, color=colour, lw=1.5, label=legend, alpha=0.9)
        if strip.closed:
            ax.plot(
                [pts[-1, 1], pts[0, 1]],
                [pts[-1, 0], pts[0, 0]],
                color=colour,
                lw=1.5,
                alpha=0.9,
            )

    ax.set_title(title)
    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
    )

    apply_structural_section_axes(ax)
    if label_seen:
        ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig, ax


__all__ = ["plot_shell_midline_strip_topology"]
