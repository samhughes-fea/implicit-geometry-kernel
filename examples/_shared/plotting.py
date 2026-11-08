"""Shared helpers for the worked examples: plot style, wall drawing and output paths.

Each example writes its figures to a ``figures/`` folder next to its script. Walls
keep the same colour in every figure: skins blue, spar caps red, webs green.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from implicit_geometry_kernel.plotting.palette import NAVY, RED, TEAL  # noqa: E402
from implicit_geometry_kernel.section.boundary import LineString, OPS  # noqa: E402

TEXT_SECONDARY = "#52514e"
REFERENCE = "#0b0b0b"  # closed-form references are drawn in ink, dashed
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
WALL_COLOURS = {"skin": NAVY, "cap": RED, "web": TEAL}
# rust, mustard, olive, teal, burgundy, ochre, slate, terracotta,
# forest, gold, brick, plum, cadet, sienna, pine, wine
RETRO_PANELS = (
    "#c45c26", "#d4a017", "#6b7f3a", "#2e6b6a",
    "#8b2942", "#c4893c", "#4a6080", "#b85c38",
    "#3d5c45", "#c9a227", "#9c3b2a", "#6e4a7a",
    "#3d6e7a", "#a0673b", "#2f4f3e", "#7a3a4a",
)


def panel_colour(role: str) -> str:
    """Stable loft colour for a region / web / rib role (not draw order)."""
    key = str(role)
    h = 0
    for ch in key:
        h = (h * 16777619) ^ ord(ch)
    return RETRO_PANELS[h % len(RETRO_PANELS)]

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.grid": False})


def _thickness_band(xy, thickness):
    band = OPS.buffer(LineString(xy), 0.5 * float(thickness), cap_style=2, join_style=2)
    if getattr(band, "geom_type", "") == "MultiPolygon":
        return band.geoms
    return (band,)


def draw_cut(ax, cut, midlines: bool = True, alpha: float = 0.35) -> None:
    """Draw each wall of a ``SectionCut`` as a band of its thickness, with its midline."""
    for wall in cut.walls:
        colour = WALL_COLOURS[wall.kind]
        for poly in _thickness_band(wall.xy, wall.thickness):
            if poly.is_empty:
                continue
            ax.fill(*poly.exterior.xy, color=colour, alpha=alpha, lw=0)
        if midlines:
            ax.plot(wall.xy[:, 0], wall.xy[:, 1], color=colour, lw=0.8)
    ax.set_aspect("equal")


_CAP_FILL_THICKNESS_OVER_CHORD = 0.04


def draw_cap_solids(ax, cut, alpha: float = 0.35) -> None:
    """Stacked cap walls as thickness bands (solid midline, not mesh).

    Fill is skipped when ``t/chord`` is large (tiny tip sections) so OFFSET
    midlines stay readable next to ``draw_shell``.
    """
    chord = max(float(getattr(cut, "chord", 0.0)), 1e-12)
    for wall in cut.walls:
        if wall.kind != "cap":
            continue
        colour = WALL_COLOURS["cap"]
        if float(wall.thickness) / chord < _CAP_FILL_THICKNESS_OVER_CHORD:
            for poly in _thickness_band(wall.xy, wall.thickness):
                if poly.is_empty:
                    continue
                ax.fill(*poly.exterior.xy, color=colour, alpha=alpha, lw=0)
        ax.plot(wall.xy[:, 0], wall.xy[:, 1], color=colour, lw=0.8)
    ax.set_aspect("equal")


def draw_shell(ax, cut, alpha: float = 0.35) -> None:
    """Conformal skin-mid loop: offset regions recoloured, no second cap polyline."""
    graph = cut.shell_graph()
    offset = {r.name for r in graph.offset_regions()}
    for run in graph.loop:
        xy = np.asarray(run.points[:, :2], float)
        colour = WALL_COLOURS["cap"] if run.name in offset else WALL_COLOURS["skin"]
        for poly in _thickness_band(xy, run.thickness):
            if poly.is_empty:
                continue
            ax.fill(*poly.exterior.xy, color=colour, alpha=alpha, lw=0)
        ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=0.8)
    for run in graph.webs:
        xy = np.asarray(run.points[:, :2], float)
        colour = WALL_COLOURS["web"]
        for poly in _thickness_band(xy, run.thickness):
            if poly.is_empty:
                continue
            ax.fill(*poly.exterior.xy, color=colour, alpha=alpha, lw=0)
        ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=0.8)
    ax.set_aspect("equal")


def note(fig, text: str) -> None:
    """A grey note under the figure: the model and what the figure checks against."""
    fig.text(0.01, 0.005, text, fontsize=7.5, color=TEXT_SECONDARY, ha="left", va="bottom")


def save(fig, script: str, name: str) -> Path:
    """Write ``figures/<name>.png`` next to ``script`` and close the figure."""
    out = Path(script).resolve().parent / "figures"
    out.mkdir(exist_ok=True)
    path = out / f"{name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path.relative_to(Path(script).resolve().parents[2])}")
    return path
