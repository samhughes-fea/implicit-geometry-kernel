"""Shared helpers for the worked examples: plot style, wall drawing and output paths.

Each example writes its figures to a ``figures/`` folder next to its script. Walls
keep the same colour in every figure: skins navy, spar caps red, webs teal.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, to_rgba  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from implicit_geometry_kernel.plotting.palette import NAVY, RED, RETRO_PANELS, TEAL, panel_colour  # noqa: E402
from implicit_geometry_kernel.section.boundary import LineString, OPS  # noqa: E402

TEXT_SECONDARY = "#52514e"
REFERENCE = "#0b0b0b"  # closed-form references are drawn in ink, dashed
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
WALL_COLOURS = {"skin": NAVY, "cap": RED, "web": TEAL}
MESH_EDGE = "#1f2937"
MESH_EDGE_ALPHA = 0.82
MESH_LW_3D = 0.40
MESH_LW_FLAT = 0.28
MESH_LW_RIB = MESH_LW_3D


def mesh_edge_rgba(n: int):
    """Broadcast loft-quad stroke for ``n`` faces."""
    return np.broadcast_to(to_rgba(MESH_EDGE, MESH_EDGE_ALPHA), (max(int(n), 0), 4)).copy()


def add_faces_3d(ax, faces, *, scale: bool = True, **kwargs):
    """Add polygons to a 3D axes and set the data limits from their real vertices.

    Matplotlib pads faces of unequal vertex count (loft quads next to rib polygons) into
    one array with ``np.empty`` and autoscales on the whole array, padding included, so
    the limits pick up whatever was in that memory: wrong limits, or "Axis limits cannot
    be NaN or Inf". The collection is therefore added without autoscaling and the limits
    come from the vertices themselves.

    A ``(n, k, 3)`` array is one equal-length sheet (the quad mesh). ``scale=False``
    leaves limits to the caller, which already holds that array.
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    if isinstance(faces, np.ndarray) and faces.ndim == 3:
        verts = np.asarray(faces, float)
        pts = verts.reshape(-1, 3)
    else:
        verts = [np.asarray(f, float) for f in faces]
        pts = np.concatenate(verts) if verts else np.zeros((0, 3))
    coll = Poly3DCollection(verts, **kwargs)
    ax.add_collection3d(coll, autolim=False)
    if scale and len(pts):
        ax.auto_scale_xyz(pts[:, 0], pts[:, 1], pts[:, 2], had_data=ax.has_data())
    return coll


def is_rib_role(role: str) -> bool:
    key = str(role).lower()
    return "rib" in key or key.startswith("tip")


def mesh_linewidths(kinds, *, base: float = MESH_LW_3D, rib: float | None = None) -> list[float]:
    """Loft-quad stroke on every role (ribs match webs and skins)."""
    weight = base if rib is None else rib
    return [weight if is_rib_role(k) else base for k in kinds]

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


def legend_outside(ax, **kwargs):
    """Legend to the right of the axes, not over the geometry."""
    loc = kwargs.pop("loc", "center left")
    bbox = kwargs.pop("bbox_to_anchor", (1.02, 0.5))
    return ax.legend(loc=loc, bbox_to_anchor=bbox, fontsize=kwargs.pop("fontsize", 7), **kwargs)


def meta_cmap(n: int, colours=None):
    """Discrete colormap for integer feature ids 0 .. n-1 (not tab10 on raw ids)."""
    n = max(int(n), 1)
    palette = list(colours) if colours is not None else list(SERIES) + list(RETRO_PANELS)
    while len(palette) < n:
        palette = palette + palette
    cmap = ListedColormap(palette[:n])
    bounds = np.arange(n + 1) - 0.5
    return cmap, BoundaryNorm(bounds, n)


def meta_legend(ax, names: dict[int, str], colours=None):
    """Patches for integer id → name."""
    cmap, _ = meta_cmap(max(names) + 1 if names else 1, colours)
    handles = [Patch(facecolor=cmap(i), edgecolor="none", label=f"{i}:{names[i]}") for i in sorted(names)]
    if handles:
        legend_outside(ax, handles=handles)
    return handles


def note(fig, text: str) -> None:
    """A grey note under the figure: the model and what the figure checks against."""
    fig.text(0.01, 0.005, text, fontsize=7.5, color=TEXT_SECONDARY, ha="left", va="bottom")


def save(fig, script: str, name: str, *, dpi: int = 150) -> Path:
    """Write ``figures/<name>.png`` next to ``script`` and close the figure."""
    out = Path(script).resolve().parent / "figures"
    out.mkdir(exist_ok=True)
    path = out / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path.relative_to(Path(script).resolve().parents[2])}")
    return path
