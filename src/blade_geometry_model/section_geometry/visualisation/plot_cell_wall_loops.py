"""Cell-wall loop overlays for section geometry plots.

Promoted from ``visualisation/demos/cell_wall_loops_demo.py`` so campaign
consumers can draw closed cell loops without the demo case payload.
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np
from matplotlib.lines import Line2D

LOOP_STYLES: tuple[tuple[str, str], ...] = (
    ("#e41a1c", "#fbb4ae"),
    ("#377eb8", "#b3cde3"),
    ("#4daf4a", "#ccebc5"),
    ("#ff7f00", "#fed9a6"),
    ("#984ea3", "#decbe4"),
    ("#a65628", "#e5d8bd"),
)


def draw_component_midlines(ax: Any, strips: Sequence[Any]) -> list[Line2D]:
    """Faint skin/web midlines over a section plot."""
    from blade_geometry_model.coordinate_frames import (
        section_plot_axes,
    )

    handles: list[Line2D] = []
    web_drawn = False
    for s in strips:
        kind = str(getattr(s, "kind", "") or "")
        xy = np.asarray(
            getattr(s, "midline_twisted_pitch_axis_blade_frame", None), dtype=float
        )
        if xy.ndim != 2 or xy.shape[0] < 2:
            continue
        h, v = section_plot_axes(xy)
        if kind == "skin":
            (ln,) = ax.plot(
                h, v, color="#444444", lw=1.4, ls=":", zorder=1, label="skin midline"
            )
            if not any(handle.get_label() == "skin midline" for handle in handles):
                handles.append(ln)
        elif kind == "web":
            (ln,) = ax.plot(
                h, v, color="#111111", lw=1.8, ls="--", zorder=5, label="web midline"
            )
            if not web_drawn:
                handles.append(ln)
                web_drawn = True
    return handles


def draw_cell_wall_loops_dashed(
    ax: Any,
    loops: Sequence[Any],
    *,
    lw: float = 1.35,
) -> list[Line2D]:
    """Draw each cell wall loop as a thin dashed centreline."""
    from blade_geometry_model.coordinate_frames import (
        section_plot_axes,
    )

    handles: list[Line2D] = []
    for i, lp in enumerate(loops):
        xy = np.asarray(lp.polyline_twisted_pitch_axis_blade_frame, dtype=float)
        edge, _fill = LOOP_STYLES[i % len(LOOP_STYLES)]
        h, v = section_plot_axes(xy)
        area = float(getattr(lp, "area_m2", 0.0) or 0.0)
        cid = getattr(lp, "id", i)
        (ln,) = ax.plot(
            h,
            v,
            color=edge,
            lw=lw,
            ls=(0, (2.8, 1.6)),
            solid_capstyle="round",
            zorder=20,
            label=f"cell_{cid}  A={area:.4g} m²",
        )
        handles.append(ln)
    return handles


def draw_cell_wall_loops_filled(
    ax: Any,
    strips: Sequence[Any],
    loops: Sequence[Any],
) -> None:
    """Filled cell loops over midlines (standalone midline view)."""
    from blade_geometry_model.coordinate_frames import (
        section_plot_axes,
    )

    draw_component_midlines(ax, strips)
    for i, lp in enumerate(loops):
        xy = np.asarray(lp.polyline_twisted_pitch_axis_blade_frame, dtype=float)
        edge, fill = LOOP_STYLES[i % len(LOOP_STYLES)]
        h, v = section_plot_axes(xy)
        ax.fill(
            h, v, facecolor=fill, edgecolor=edge, linewidth=2.2, alpha=0.55, zorder=3
        )
        ax.plot(h, v, color=edge, lw=2.2, solid_capstyle="round", zorder=4)


def plot_section_with_cell_wall_loops(
    section: Any,
    grid: Any,
    strips: Sequence[Any],
    loops: Sequence[Any],
    ax: Any = None,
    *,
    title: str = "",
    **plot_section_kwargs: Any,
) -> Any:
    """``plot_section`` plus midlines and dashed cell-wall loops."""
    import matplotlib.pyplot as plt
    from blade_geometry_model.section_geometry import plot_section

    if ax is None:
        _fig, ax = plt.subplots(figsize=(8.0, 6.0))

    kwargs = {
        "alpha": 0.42,
        "show_airfoil": False,
        "show_legend": True,
        "omit_cores": True,
        "auto_limits": True,
        "limit_pad_frac": 0.06,
        "linewidths": 0.85,
        **plot_section_kwargs,
    }
    if title:
        kwargs["title"] = title
    plot_section(section, grid, ax=ax, **kwargs)
    mid_handles = draw_component_midlines(ax, strips)
    loop_handles = draw_cell_wall_loops_dashed(ax, loops)
    existing = ax.get_legend()
    if existing is not None:
        ax.add_artist(existing)
    ax.legend(
        mid_handles + loop_handles,
        [h.get_label() for h in mid_handles + loop_handles],
        fontsize=9,
        loc="center left",
        bbox_to_anchor=(1.01, 0.35),
        frameon=True,
        fancybox=False,
        edgecolor="0.55",
        framealpha=0.97,
        title="centrelines / cell loops",
        title_fontsize=9,
    )
    return ax


__all__ = [
    "LOOP_STYLES",
    "draw_component_midlines",
    "draw_cell_wall_loops_dashed",
    "draw_cell_wall_loops_filled",
    "plot_section_with_cell_wall_loops",
]
