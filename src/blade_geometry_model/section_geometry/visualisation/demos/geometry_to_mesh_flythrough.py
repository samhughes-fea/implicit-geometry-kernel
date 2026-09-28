"""
visualisation/demos/geometry_to_mesh_flythrough.py
=======================================
The handoff that makes the implicit representation useful: solids -> shell mesh.

Animates the four stages that carry an implicitly-defined section into a
finite-element shell mesh, with the camera diving into a web-skin junction
where the difference between the stages is actually visible:

1. **CSG solids** - named regions, each a compiled ``Expr`` tree evaluated as a
   signed-distance field. Nothing here is a curve or a node yet.
2. **Raw midlines** - ``midlines.build_shell_midline_strips``. One centreline
   per component, taken from the *recipe* (aerofoil offset, web capsule axis,
   IML offset), **not** from a post-boolean zero contour. At this stage the
   strips are independent: a web centreline runs to its own end and the skin
   centreline passes straight through it. They do not share endpoints.
3. **Junction-prepared strips** - ``topology.prepare_shell_midline_junctions``.
   Skin arcs are split at web and cap feet, endpoints are merged onto shared
   junction points, and cap-web T-junctions are inserted. This is the step that
   turns a pile of curves into a connected 1-D topology.
4. **Q4 shell mesh** - ``section_mesh.SectionMeshRunner``. Each prepared panel
   is discretised into shell elements; endpoint clusters become shared nodes,
   so the mesh is conforming across every junction.

Note on ordering: stages 2 and 3 are both "midlines" and are easy to conflate.
The distinction that matters is **raw** (per-component, endpoints free) versus
**prepared** (split and snapped). Meshing consumes the prepared form only;
``SectionGeometryResult.midline_strips`` is already prepared.

Writes under ``visualisation/output/geometry_to_mesh/``:

* ``geometry_to_mesh_flythrough.gif`` - the animation.
* ``geometry_to_mesh_stages.png`` - the same four stages as a static strip.
* ``geometry_to_mesh.json`` - per-stage counts and junction metrics.

Requires ImageMagick (``magick``) on PATH for GIF encoding; ``--no-gif`` skips it.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/geometry_to_mesh_flythrough.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

_demos_dir = Path(__file__).resolve().parent


def _load_sibling(name: str, filename: str):
    path = _demos_dir / filename
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


_ensure_mod = _load_sibling("_ensure_repo_path", "_ensure_repo_path.py")
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    PRESENTATION_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)

# The render/specs helpers use plain sibling imports and dataclasses, so they
# have to be imported normally rather than exec'd out of a spec - a module
# loaded that way is absent from sys.modules and @dataclass cannot resolve it.
if str(_demos_dir) not in sys.path:
    sys.path.insert(0, str(_demos_dir))

import feature_flythrough_render as _render  # noqa: E402
import feature_flythrough_specs as _specs  # noqa: E402

from blade_structure_utilities.geometry_utilities.coordinate_frames import (
    section_plot_axes,
)
from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SDFGrid,
    SectionGeometryRunner,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.ir import (
    eval_expr_batch,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
    _component_color,
    rotated_naca_grid_to_section_plot_axes,
)
from blade_structure_utilities.mesh_utilities.section_mesh import SectionMeshRunner

_STAGE_NAMES = ("csg", "raw", "prepared", "mesh")

_STAGE_TITLES = {
    "csg": "1 - CSG solids",
    "raw": "2 - Raw midlines",
    "prepared": "3 - Junction-prepared strips",
    "mesh": "4 - Q4 shell mesh",
}

# Kept short: these are rendered inside fixed-width sidebar boxes.
_STAGE_SUBTITLES = {
    "csg": "named regions as distance fields",
    "raw": "one per component; ends free",
    "prepared": "skin split at feet; ends merged",
    "mesh": "shell elements; shared nodes",
}

_STAGE_MODULE = {
    "csg": "sections / sdf.ir",
    "raw": "midlines.build_shell_midline_strips",
    "prepared": "topology.prepare_shell_midline_junctions",
    "mesh": "section_mesh.SectionMeshRunner",
}

# Callout text is built per-context in PipelineContext.stage_callouts so the
# numbers quoted on screen are the ones actually measured for this section.

_KIND_COLOR = {"skin": "#2f6f9f", "web": "#c0392b", "cap": "#d9a021"}
_INK = "#1d262d"
_MUTED = "#7d8894"
_ACCENT = "#7a4fa3"


# --------------------------------------------------------------------------
# data assembly
# --------------------------------------------------------------------------


@dataclass
class PipelineContext:
    """Everything the animation draws, built once."""

    key: str
    chord_m: float
    twist_rad: float
    section: Any
    grid: SDFGrid
    x_plot: np.ndarray
    y_plot: np.ndarray
    phi_by_label: dict[str, np.ndarray]
    raw_strips: tuple
    prepared_strips: tuple
    mesh: Any
    junction_point: tuple[float, float]
    wide_center: tuple[float, float]
    wide_half_width: float
    zoom_half_width: float
    frame_aspect: float

    def windows(
        self, aspect: float | None = None
    ) -> tuple[tuple[float, float, float, float], tuple[float, float, float, float]]:
        """``(wide, zoom)`` windows sized for a panel of the given aspect.

        Axes are drawn with ``aspect='equal'``, so a window whose aspect does not
        match its panel gets padded with dead space. The GIF view and the static
        strip have different panel shapes, so each asks for its own.
        """
        use = float(self.frame_aspect if aspect is None else aspect)
        return (
            _window_around(self.wide_center, self.wide_half_width, use),
            _window_around(self.junction_point, self.zoom_half_width, use),
        )

    @property
    def wide_window(self) -> tuple[float, float, float, float]:
        return self.windows()[0]

    @property
    def zoom_window(self) -> tuple[float, float, float, float]:
        return self.windows()[1]

    @property
    def web_foot_gap_max_mm(self) -> tuple[float, float]:
        """Worst web-foot-to-skin-node gap, raw then prepared, in millimetres."""
        raw = endpoint_to_skin_vertex_gaps(self.raw_strips, "web")
        prepared = endpoint_to_skin_vertex_gaps(self.prepared_strips, "web")
        as_mm = lambda gaps: (float(gaps.max()) * 1000.0 if gaps.size else float("nan"))
        return as_mm(raw), as_mm(prepared)

    @property
    def stage_metrics(self) -> dict[str, str]:
        summary = self.mesh.summary()
        raw_gap_mm, prepared_gap_mm = self.web_foot_gap_max_mm
        return {
            "csg": f"{len(self.phi_by_label)} implicit regions",
            "raw": (
                f"{len(self.raw_strips)} centrelines, "
                f"feet {raw_gap_mm:.1f} mm off-node"
            ),
            "prepared": (
                f"{len(self.prepared_strips)} panels, "
                f"{int(summary['n_clusters'])} junctions, "
                f"gap {prepared_gap_mm:.0f} mm"
            ),
            "mesh": (
                f"{int(summary['n_total_elements'])} elements, "
                f"{int(summary['n_total_nodes'])} nodes"
            ),
        }

    @property
    def stage_callouts(self) -> dict[str, str]:
        raw_gap_mm, prepared_gap_mm = self.web_foot_gap_max_mm
        return {
            "csg": "the junction exists only as overlapping fields",
            "raw": (
                f"web foot lands {raw_gap_mm:.1f} mm along the skin: no node to share"
            ),
            "prepared": (
                f"skin split at the foot; web endpoint merged ({prepared_gap_mm:.0f} mm)"
            ),
            "mesh": "elements on both sides terminate on that shared node",
        }


def _to_key(label) -> str:
    return str(getattr(label, "value", label))


def _plot_points(points_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Canonical-B ``(X_B, Y_B)`` -> plot ``(horizontal, vertical)``."""
    horizontal, vertical = section_plot_axes(
        np.asarray(points_b, dtype=float).reshape(-1, 2)
    )
    return np.asarray(horizontal, dtype=float), np.asarray(vertical, dtype=float)


def _window_around(
    center: tuple[float, float], half_width: float, aspect: float
) -> tuple[float, float, float, float]:
    half_height = half_width / max(aspect, 1e-6)
    return (
        center[0] - half_width,
        center[0] + half_width,
        center[1] - half_height,
        center[1] + half_height,
    )


def build_context(args: argparse.Namespace) -> PipelineContext:
    layout = resolve_system_type(str(args.key))
    twist_rad = float(np.deg2rad(float(args.twist_deg)))
    airfoil = AirfoilSDF.from_naca(
        str(args.naca), chord=float(args.chord), n_points=200
    )

    result = (
        SectionGeometryRunner(
            layout=layout,
            airfoil=airfoil,
            skin_thickness=float(args.skin_thickness),
            spar_cap_thickness=float(args.spar_cap_thickness),
            shear_web_thickness=float(args.shear_web_thickness),
            twist_rad=twist_rad,
            outputs=("midlines",),
        )
        .execute()
        .get_results()
    )

    # Raw emit, before topology preparation. SectionGeometryResult only ever
    # exposes the prepared form, so the raw stage is rebuilt from the host.
    raw_strips = build_shell_midline_strips(
        result.section, twist_rad=twist_rad, chord_m=float(args.chord)
    )
    prepared_strips = tuple(result.midline_strips or ())

    mesh = (
        SectionMeshRunner(
            geometry_result=result,
            layout_key=str(args.key),
            target_element_length_frac=float(args.element_length_frac),
        )
        .execute()
        .get_results()
    )

    # Fields are evaluated once and reused for every frame; re-evaluating per
    # frame would dominate render time.
    airfoil_for_grid = (
        airfoil.rotate(twist_rad) if abs(twist_rad) > 1e-12 else airfoil
    )
    grid = SDFGrid.from_airfoil(
        airfoil_for_grid, nx=int(args.nx), ny=int(args.ny), padding=0.05
    )
    x_plot, y_plot = rotated_naca_grid_to_section_plot_axes(
        grid.X, grid.Y, twist_rad=twist_rad, chord_m=float(args.chord)
    )
    phi_by_label: dict[str, np.ndarray] = {}
    for label in result.section:
        expression = getattr(result.section[label], "expr", None)
        if expression is None:
            continue
        phi_by_label[_to_key(label)] = np.asarray(
            eval_expr_batch(expression, grid.X, grid.Y, {}), dtype=float
        )

    junction_point = _select_junction(mesh, prefer_size=4)

    aspect = float(args.frame_aspect)
    all_horizontal, all_vertical = _plot_points(
        np.concatenate(
            [
                np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
                for s in prepared_strips
            ]
        )
    )
    span = max(
        float(all_horizontal.max() - all_horizontal.min()),
        float(all_vertical.max() - all_vertical.min()),
    )
    wide_center = (
        0.5 * float(all_horizontal.max() + all_horizontal.min()),
        0.5 * float(all_vertical.max() + all_vertical.min()),
    )

    return PipelineContext(
        key=str(args.key),
        chord_m=float(args.chord),
        twist_rad=twist_rad,
        section=result.section,
        grid=grid,
        x_plot=x_plot,
        y_plot=y_plot,
        phi_by_label=phi_by_label,
        raw_strips=tuple(raw_strips),
        prepared_strips=prepared_strips,
        mesh=mesh,
        junction_point=junction_point,
        wide_center=wide_center,
        wide_half_width=0.56 * span,
        zoom_half_width=float(args.zoom_half_width_frac) * float(args.chord),
        frame_aspect=aspect,
    )


def _select_junction(mesh: Any, *, prefer_size: int) -> tuple[float, float]:
    """Pick a junction to dive into: a web foot if there is one.

    Web feet are the interesting case - three or four panels meet there, and the
    raw-vs-prepared difference is largest. Falls back to the busiest cluster.
    """
    clusters = list(mesh.clusters)
    if not clusters:
        raise RuntimeError("mesh has no junction clusters to zoom on")
    candidates = [c for c in clusters if len(c.members) >= prefer_size]
    chosen = max(candidates or clusters, key=lambda c: len(c.members))
    horizontal, vertical = _plot_points(np.asarray(chosen.point_yz, dtype=float))
    return (float(horizontal[0]), float(vertical[0]))


# --------------------------------------------------------------------------
# stage drawing
# --------------------------------------------------------------------------


def draw_junction_locator(axis, context: PipelineContext, alpha: float) -> None:
    """Ring on the junction the camera is heading for.

    Without it the CSG stage is an undifferentiated field of colour once the
    camera is in close, and the viewer has nothing to fix on across the fade.
    """
    if alpha <= 0.01:
        return
    axis.scatter(
        [context.junction_point[0]],
        [context.junction_point[1]],
        s=340,
        facecolor="none",
        edgecolor=_ACCENT,
        linewidths=2.0,
        alpha=0.75 * alpha,
        zorder=11,
    )


def draw_csg(axis, context: PipelineContext, alpha: float) -> None:
    if alpha <= 0.01:
        return
    for label, phi in context.phi_by_label.items():
        axis.contourf(
            context.x_plot,
            context.y_plot,
            phi,
            levels=[float(phi.min()) - 1.0, 0.0],
            colors=[_component_color(label)],
            alpha=0.55 * alpha,
            zorder=1,
        )
        axis.contour(
            context.x_plot,
            context.y_plot,
            phi,
            levels=[0.0],
            colors=[_component_color(label)],
            linewidths=1.1,
            alpha=0.95 * alpha,
            zorder=2,
        )


def draw_strips(
    axis,
    strips,
    alpha: float,
    *,
    show_endpoints: bool,
    linewidth: float,
    zorder: int,
) -> None:
    if alpha <= 0.01:
        return
    for strip in strips:
        points = np.asarray(
            strip.midline_twisted_pitch_axis_blade_frame, dtype=float
        )
        if points.ndim != 2 or points.shape[0] < 2:
            continue
        horizontal, vertical = _plot_points(points)
        color = _KIND_COLOR.get(str(strip.kind), _MUTED)
        axis.plot(
            horizontal,
            vertical,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            solid_capstyle="round",
            zorder=zorder,
        )
        if show_endpoints:
            axis.scatter(
                [horizontal[0], horizontal[-1]],
                [vertical[0], vertical[-1]],
                s=44,
                facecolor="white",
                edgecolor=color,
                linewidths=1.6,
                alpha=alpha,
                zorder=zorder + 1,
            )


def draw_junction_markers(axis, context: PipelineContext, alpha: float) -> None:
    if alpha <= 0.01:
        return
    for cluster in context.mesh.clusters:
        horizontal, vertical = _plot_points(
            np.asarray(cluster.point_yz, dtype=float)
        )
        axis.scatter(
            horizontal,
            vertical,
            s=150,
            facecolor="none",
            edgecolor=_ACCENT,
            linewidths=2.2,
            alpha=alpha,
            zorder=9,
        )


def draw_mesh(axis, context: PipelineContext, alpha: float) -> None:
    if alpha <= 0.01:
        return
    mesh = context.mesh
    for panel_mesh in mesh.panel_meshes:
        nodes = np.asarray(panel_mesh.yz_nodes, dtype=float)
        if nodes.size == 0:
            continue
        horizontal, vertical = _plot_points(nodes)
        color = _KIND_COLOR.get(str(panel_mesh.kind), _MUTED)
        axis.plot(
            horizontal,
            vertical,
            color=color,
            linewidth=1.6,
            alpha=0.85 * alpha,
            zorder=5,
        )
        # Element boundaries: a tick across the midline at every node.
        axis.scatter(
            horizontal,
            vertical,
            s=13,
            color=color,
            alpha=alpha,
            zorder=6,
            linewidths=0,
        )
    for cluster in mesh.clusters:
        horizontal, vertical = _plot_points(
            np.asarray(cluster.point_yz, dtype=float)
        )
        axis.scatter(
            horizontal,
            vertical,
            s=170,
            marker="o",
            facecolor="white",
            edgecolor=_ACCENT,
            linewidths=2.4,
            alpha=alpha,
            zorder=10,
        )


_STAGE_DRAW = {
    "csg": lambda ax, ctx, a: (draw_csg(ax, ctx, a), draw_junction_locator(ax, ctx, a)),
    "raw": lambda ax, ctx, a: draw_strips(
        ax, ctx.raw_strips, a, show_endpoints=True, linewidth=2.4, zorder=4
    ),
    "prepared": lambda ax, ctx, a: (
        draw_strips(
            ax, ctx.prepared_strips, a, show_endpoints=True, linewidth=2.4, zorder=4
        ),
        draw_junction_markers(ax, ctx, a),
    ),
    "mesh": lambda ax, ctx, a: draw_mesh(ax, ctx, a),
}


# --------------------------------------------------------------------------
# timeline
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Beat:
    """One timeline segment: a camera move and/or a stage cross-fade."""

    frames: int
    stage_from: str
    stage_to: str
    zoom_from: float  # 0 = wide, 1 = junction
    zoom_to: float


def default_timeline(hold_scale: float = 1.0) -> list[Beat]:
    hold = max(int(round(6 * hold_scale)), 2)
    return [
        Beat(hold, "csg", "csg", 0.0, 0.0),
        Beat(16, "csg", "csg", 0.0, 1.0),
        Beat(4, "csg", "csg", 1.0, 1.0),
        Beat(12, "csg", "raw", 1.0, 1.0),
        Beat(hold, "raw", "raw", 1.0, 1.0),
        Beat(12, "raw", "prepared", 1.0, 1.0),
        Beat(hold, "prepared", "prepared", 1.0, 1.0),
        Beat(12, "prepared", "mesh", 1.0, 1.0),
        Beat(hold, "mesh", "mesh", 1.0, 1.0),
        Beat(16, "mesh", "mesh", 1.0, 0.0),
        Beat(hold + 2, "mesh", "mesh", 0.0, 0.0),
    ]


def smoothstep(t: float) -> float:
    t = min(max(float(t), 0.0), 1.0)
    return t * t * (3.0 - 2.0 * t)


def expand_timeline(beats: list[Beat]) -> list[tuple[dict[str, float], float, str]]:
    """Flatten beats into per-frame ``(stage weights, zoom, headline stage)``."""
    frames: list[tuple[dict[str, float], float, str]] = []
    for beat in beats:
        for index in range(int(beat.frames)):
            progress = (
                0.0 if beat.frames <= 1 else index / float(beat.frames - 1)
            )
            eased = smoothstep(progress)
            weights = {name: 0.0 for name in _STAGE_NAMES}
            if beat.stage_from == beat.stage_to:
                weights[beat.stage_to] = 1.0
            else:
                weights[beat.stage_from] = 1.0 - eased
                weights[beat.stage_to] = eased
            zoom = beat.zoom_from + (beat.zoom_to - beat.zoom_from) * eased
            headline = beat.stage_to if eased >= 0.5 else beat.stage_from
            frames.append((weights, zoom, headline))
    return frames


def interpolate_window(
    wide: tuple[float, float, float, float],
    zoom: tuple[float, float, float, float],
    t: float,
) -> tuple[float, float, float, float]:
    """Blend the two windows geometrically so the zoom feels linear in scale."""
    t = min(max(float(t), 0.0), 1.0)
    wide_center = (0.5 * (wide[0] + wide[1]), 0.5 * (wide[2] + wide[3]))
    zoom_center = (0.5 * (zoom[0] + zoom[1]), 0.5 * (zoom[2] + zoom[3]))
    wide_half = 0.5 * (wide[1] - wide[0])
    zoom_half = 0.5 * (zoom[1] - zoom[0])
    half = float(np.exp(np.log(wide_half) * (1.0 - t) + np.log(zoom_half) * t))
    center_x = wide_center[0] + (zoom_center[0] - wide_center[0]) * t
    center_y = wide_center[1] + (zoom_center[1] - wide_center[1]) * t
    aspect = (wide[1] - wide[0]) / max(wide[3] - wide[2], 1e-9)
    return _window_around((center_x, center_y), half, aspect)


# --------------------------------------------------------------------------
# frame rendering
# --------------------------------------------------------------------------


def draw_sidebar(
    axis, context: PipelineContext, headline: str, weights: dict[str, float]
) -> None:
    """Pipeline step list with the active step lit."""
    axis.set_axis_off()
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    metrics = context.stage_metrics

    axis.text(
        0.02,
        0.965,
        "section_geometry -> section_mesh",
        fontsize=11.5,
        fontweight="bold",
        color=_INK,
        va="top",
    )
    axis.text(
        0.02,
        0.918,
        f"layout {context.key}   chord {context.chord_m:g} m",
        fontsize=9,
        color=_MUTED,
        va="top",
    )

    top = 0.855
    row_height = 0.185
    for index, name in enumerate(_STAGE_NAMES):
        weight = float(weights.get(name, 0.0))
        active = name == headline
        y_top = top - index * row_height
        emphasis = 0.28 + 0.72 * weight
        box = FancyBboxPatch(
            (0.02, y_top - row_height + 0.028),
            0.96,
            row_height - 0.038,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            facecolor="#eef2f6" if active else "#f7f9fb",
            edgecolor=_ACCENT if active else "#dde3e9",
            linewidth=2.0 if active else 1.0,
            zorder=1,
        )
        axis.add_patch(box)
        axis.text(
            0.055,
            y_top - 0.012,
            _STAGE_TITLES[name],
            fontsize=10.2,
            fontweight="bold" if active else "normal",
            color=_INK if active else _MUTED,
            alpha=emphasis,
            va="top",
            zorder=2,
        )
        axis.text(
            0.055,
            y_top - 0.056,
            _STAGE_SUBTITLES[name],
            fontsize=8.2,
            color=_MUTED,
            alpha=emphasis,
            va="top",
            zorder=2,
        )
        axis.text(
            0.055,
            y_top - 0.094,
            metrics[name],
            fontsize=8.6,
            color=_KIND_COLOR["web"] if active else _MUTED,
            alpha=emphasis,
            va="top",
            zorder=2,
        )
        axis.text(
            0.055,
            y_top - 0.130,
            _STAGE_MODULE[name],
            fontsize=7.4,
            family="monospace",
            color=_MUTED,
            alpha=0.55 * emphasis + 0.15,
            va="top",
            zorder=2,
        )


def render_frame(
    context: PipelineContext,
    weights: dict[str, float],
    zoom: float,
    headline: str,
    profile,
    out_path: Path,
) -> None:
    figure = plt.figure(figsize=(13.0, 6.0))
    grid_spec = figure.add_gridspec(1, 2, width_ratios=[2.55, 1.0], wspace=0.04)
    axis_view = figure.add_subplot(grid_spec[0, 0])
    axis_side = figure.add_subplot(grid_spec[0, 1])

    for name in _STAGE_NAMES:
        weight = float(weights.get(name, 0.0))
        if weight > 0.01:
            _STAGE_DRAW[name](axis_view, context, weight)

    window = interpolate_window(context.wide_window, context.zoom_window, zoom)
    axis_view.set_xlim(window[0], window[1])
    axis_view.set_ylim(window[2], window[3])
    if not axis_view.xaxis_inverted():
        axis_view.invert_xaxis()  # LE on the right, matching the package figures
    axis_view.set_aspect("equal")
    axis_view.set_xticks([])
    axis_view.set_yticks([])
    for spine in axis_view.spines.values():
        spine.set_edgecolor("#c9d1d9")

    axis_view.set_title(
        f"{_STAGE_TITLES[headline]} - {_STAGE_SUBTITLES[headline]}",
        fontsize=12.5,
        color=_INK,
        pad=10,
    )
    scale_note = "whole section" if zoom < 0.35 else "junction zoom"
    axis_view.text(
        0.015,
        0.02,
        scale_note,
        transform=axis_view.transAxes,
        fontsize=8.5,
        color=_MUTED,
    )
    # Only the dominant stage's callout is drawn, and it fades through zero at
    # the midpoint of a cross-fade. Drawing both at partial alpha overlaps two
    # different sentences in the same place and neither is readable.
    if zoom > 0.55:
        headline_weight = float(weights.get(headline, 0.0))
        text_alpha = min(max(2.0 * headline_weight - 1.0, 0.0), 1.0)
        if text_alpha > 0.02:
            axis_view.text(
                0.5,
                0.955,
                context.stage_callouts[headline],
                transform=axis_view.transAxes,
                ha="center",
                va="top",
                fontsize=10.5,
                color=_INK,
                alpha=text_alpha,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.80 * text_alpha,
                    "edgecolor": "none",
                    "boxstyle": "round,pad=0.35",
                },
            )

    draw_sidebar(axis_side, context, headline, weights)

    save_figure(
        figure,
        out_path,
        dpi=int(profile.dpi),
        pad_inches=0.18,
    )
    plt.close(figure)


def render_static_strip(
    context: PipelineContext, out_path: Path, *, presentation: bool = False
) -> None:
    """The same four stages side by side, for slides that cannot play a GIF."""
    figure, axes = plt.subplots(2, 4, figsize=(19.0, 8.4))
    metrics = context.stage_metrics
    # Panels here are near-square, unlike the GIF's wide view axis.
    wide_window, zoom_window = context.windows(aspect=1.55)
    for column, name in enumerate(_STAGE_NAMES):
        for row, window in enumerate((wide_window, zoom_window)):
            axis = axes[row, column]
            _STAGE_DRAW[name](axis, context, 1.0)
            axis.set_xlim(window[0], window[1])
            axis.set_ylim(window[2], window[3])
            if not axis.xaxis_inverted():
                axis.invert_xaxis()
            axis.set_aspect("equal")
            axis.set_xticks([])
            axis.set_yticks([])
            if row == 0:
                axis.set_title(
                    f"{_STAGE_TITLES[name]}\n{_STAGE_SUBTITLES[name]}",
                    fontsize=11,
                    color=_INK,
                )
                axis.text(
                    0.5,
                    -0.06,
                    metrics[name],
                    transform=axis.transAxes,
                    ha="center",
                    va="top",
                    fontsize=9,
                    color=_KIND_COLOR["web"],
                )
            else:
                axis.text(
                    0.5,
                    -0.05,
                    _STAGE_MODULE[name],
                    transform=axis.transAxes,
                    ha="center",
                    va="top",
                    fontsize=8,
                    family="monospace",
                    color=_MUTED,
                )
    axes[0, 0].set_ylabel("whole section", fontsize=10, color=_MUTED)
    axes[1, 0].set_ylabel("junction zoom", fontsize=10, color=_MUTED)
    figure.suptitle(
        "Implicit solids to a conforming shell mesh: the four stages, "
        "whole section (top) and at a web-skin junction (bottom)",
        fontsize=14,
    )
    figure.tight_layout(rect=(0, 0.01, 1, 0.93))
    save_figure(figure, out_path, presentation=presentation)
    plt.close(figure)
    print(f"Wrote {out_path}")


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------


def endpoint_to_skin_vertex_gaps(strips, kind: str) -> np.ndarray:
    """Distance from each ``kind`` endpoint to the nearest skin **vertex**.

    This is what topology preparation fixes for web feet, and a pairwise minimum
    over all endpoints does not measure it: raw skin strips already meet each
    other exactly at the LE and TE, so that statistic reads 0 either way.

    A web foot landing part-way along a skin segment is not a mesh junction -
    there is no node there to share. Preparation splits the skin arc at the foot,
    so the web gap goes to exactly zero.

    Cap endpoints are reported separately and are *expected* to stay offset:
    a spar cap runs inboard of the skin centreline by roughly
    ``(skin_thickness + cap_thickness) / 2``. That residual is geometry, not a
    junction defect, and calling it one would misrepresent the pipeline.
    """
    skin_arrays = [
        np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for strip in strips
        if str(strip.kind) == "skin"
    ]
    if not skin_arrays:
        return np.asarray([], dtype=float)
    skin_points = np.concatenate(skin_arrays, axis=0)

    gaps: list[float] = []
    for strip in strips:
        if str(strip.kind) != kind:
            continue
        points = np.asarray(
            strip.midline_twisted_pitch_axis_blade_frame, dtype=float
        )
        if points.ndim != 2 or points.shape[0] < 2:
            continue
        for endpoint in (points[0], points[-1]):
            gaps.append(
                float(np.linalg.norm(skin_points - endpoint[None, :], axis=1).min())
            )
    return np.asarray(gaps, dtype=float)


def _describe_gaps(gaps: np.ndarray) -> dict[str, float]:
    if gaps.size == 0:
        return {"count": 0, "max_m": float("nan"), "mean_m": float("nan")}
    return {
        "count": int(gaps.size),
        "max_m": float(gaps.max()),
        "mean_m": float(gaps.mean()),
    }


def junction_report(context: PipelineContext) -> dict[str, Any]:
    """Web-foot and cap-end offsets from the nearest skin node, before and after."""
    return {
        "metric": "distance from each endpoint to the nearest skin vertex [m]",
        "raw_strip_count": int(len(context.raw_strips)),
        "prepared_strip_count": int(len(context.prepared_strips)),
        "web_feet": {
            "raw": _describe_gaps(
                endpoint_to_skin_vertex_gaps(context.raw_strips, "web")
            ),
            "prepared": _describe_gaps(
                endpoint_to_skin_vertex_gaps(context.prepared_strips, "web")
            ),
            "note": "preparation splits the skin at the foot; prepared max is 0",
        },
        "cap_ends": {
            "raw": _describe_gaps(
                endpoint_to_skin_vertex_gaps(context.raw_strips, "cap")
            ),
            "prepared": _describe_gaps(
                endpoint_to_skin_vertex_gaps(context.prepared_strips, "cap")
            ),
            "note": (
                "expected to stay offset by ~(skin_thickness + cap_thickness)/2; "
                "a cap runs inboard of the skin centreline, not on it"
            ),
        },
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--key", default="2D-F")
    parser.add_argument("--naca", default="2412")
    parser.add_argument("--chord", type=float, default=1.0)
    parser.add_argument("--twist-deg", type=float, default=0.0)
    parser.add_argument("--nx", type=int, default=760)
    parser.add_argument("--ny", type=int, default=340)
    parser.add_argument("--skin-thickness", type=float, default=0.006)
    parser.add_argument("--spar-cap-thickness", type=float, default=0.014)
    parser.add_argument("--shear-web-thickness", type=float, default=0.008)
    parser.add_argument("--element-length-frac", type=float, default=0.012)
    parser.add_argument(
        "--zoom-half-width-frac",
        type=float,
        default=0.035,
        help="Junction zoom half-width as a fraction of chord.",
    )
    parser.add_argument("--frame-aspect", type=float, default=1.60)
    parser.add_argument(
        "--profile", default="balanced", choices=sorted(_specs.PROFILES)
    )
    parser.add_argument("--hold-scale", type=float, default=1.0)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--no-gif", action="store_true", help="Render frames only.")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    apply_presentation_style(bool(args.presentation))
    out_dir = args.out_dir or demo_output_dir("geometry_to_mesh")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    profile = _specs.PROFILES[str(args.profile)]
    if bool(args.presentation) and int(profile.dpi) < PRESENTATION_DPI:
        from dataclasses import replace

        profile = replace(profile, dpi=PRESENTATION_DPI)
    context = build_context(args)
    print(
        f"Built {context.key}: {len(context.phi_by_label)} regions, "
        f"{len(context.raw_strips)} raw strips, "
        f"{len(context.prepared_strips)} prepared panels, "
        f"{context.mesh.n_total_elements} elements"
    )

    render_static_strip(
        context, out_dir / "geometry_to_mesh_stages.png", presentation=bool(args.presentation)
    )

    paths = _render.prepare_scenario_paths(
        out_dir, "geometry_to_mesh", "geometry_to_mesh_flythrough.gif"
    )
    _render.reset_frame_dir(paths.frame_dir)

    timeline = expand_timeline(default_timeline(float(args.hold_scale)))
    selected = _render.select_frame_indices(len(timeline), args.max_frames)
    frame_paths: list[Path] = []
    for output_index, timeline_index in enumerate(selected):
        weights, zoom, headline = timeline[timeline_index]
        frame_file = _render.frame_path(paths.frame_dir, output_index)
        render_frame(context, weights, zoom, headline, profile, frame_file)
        frame_paths.append(frame_file)
    print(f"Rendered {len(frame_paths)} frames -> {paths.frame_dir}")

    payload: dict[str, Any] = {
        "key": context.key,
        "naca": str(args.naca),
        "chord_m": context.chord_m,
        "twist_deg": float(args.twist_deg),
        "profile": str(args.profile),
        "frames": len(frame_paths),
        "stages": [
            {
                "name": name,
                "title": _STAGE_TITLES[name],
                "module": _STAGE_MODULE[name],
                "metric": context.stage_metrics[name],
            }
            for name in _STAGE_NAMES
        ],
        "junction": junction_report(context),
        "mesh_summary": context.mesh.summary(),
        "zoom_junction_plot_xy": list(context.junction_point),
    }

    if not args.no_gif:
        held = _render.apply_holds(
            frame_paths, head=profile.hold_head, tail=profile.hold_tail
        )
        _render.encode_gif(held, paths.gif_path, profile)
        payload["gif"] = _render.validate_gif(paths.gif_path)
        print(f"Wrote {paths.gif_path}")

    _render.write_metadata(paths.metadata_path, payload)
    print(f"Wrote {paths.metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
