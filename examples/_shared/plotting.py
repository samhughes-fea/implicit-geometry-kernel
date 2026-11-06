"""Shared helpers for the worked examples: plot style, wall drawing and output paths.

Each example writes its figures to a ``figures/`` folder next to its script. Walls
keep the same colour in every figure: skins blue, spar caps red, webs green.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from shapely.geometry import LineString  # noqa: E402

from plotting.palette import NAVY, RED, TEAL  # noqa: E402

TEXT_SECONDARY = "#52514e"
REFERENCE = "#0b0b0b"  # closed-form references are drawn in ink, dashed
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
WALL_COLOURS = {"skin": NAVY, "cap": RED, "web": TEAL}

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.grid": False})


def draw_cut(ax, cut, midlines: bool = True, alpha: float = 0.35) -> None:
    """Draw each wall of a ``SectionCut`` as a band of its thickness, with its midline."""
    for wall in cut.walls:
        colour = WALL_COLOURS[wall.kind]
        band = LineString(wall.xy).buffer(0.5 * wall.thickness, cap_style=2, join_style=2)
        for poly in getattr(band, "geoms", [band]):
            ax.fill(*poly.exterior.xy, color=colour, alpha=alpha, lw=0)
        if midlines:
            ax.plot(wall.xy[:, 0], wall.xy[:, 1], color=colour, lw=0.8)
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
