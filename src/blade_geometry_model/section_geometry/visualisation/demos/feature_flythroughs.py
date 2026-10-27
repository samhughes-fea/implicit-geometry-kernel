"""Generate high-quality feature fly-through GIFs for section geometry."""

from __future__ import annotations

import argparse
import importlib.util
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

_examples_dir = Path(__file__).resolve().parent


def _load_sibling(name: str, filename: str):
    path = _examples_dir / filename
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


_ensure_mod = _load_sibling("_ensure_repo_path", "_ensure_repo_path.py")
repo_root = _ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    PRESENTATION_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)

from blade_structure_precompute.module.section_selection.processing.adapters.scale_adapter import (
    airfoil_sdf_from_station,
)
from blade_geometry_model.coordinate_frames import section_plot_axes
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_geometry_model.section_geometry.props_emit.area_mesh import (
    build_shared_section_quadtree,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.junction_preparation.shell_midline_junction_preparation import (
    prepare_shell_midline_junctions,
)
from blade_geometry_model.section_geometry.visualisation.plot import (
    apply_structural_section_axes,
    plot_section,
    plot_section_hybrid,
    plot_section_with_quadtree,
)
from blade_geometry_model.section_geometry.visualisation.plot_section import (
    rotated_naca_grid_to_section_plot_axes,
)

import spanwise_web_anchors as spanwise_probe
from feature_flythrough_render import (
    ScenarioPaths,
    apply_holds,
    encode_gif,
    frame_path,
    prepare_scenario_paths,
    reset_frame_dir,
    select_frame_indices,
    validate_gif,
    write_metadata,
)
from feature_flythrough_specs import (
    DEFAULT_PROFILE,
    DEFAULT_ROTOR,
    PROFILES,
    SCENARIOS,
    SCENARIO_ORDER,
    RenderProfile,
    ScenarioSpec,
)
from spanwise_web_anchors_shared import DEFAULT_MATERIALS, build_stations

DEFAULT_T_SKIN = spanwise_probe.T_SKIN
DEFAULT_T_SPAR = spanwise_probe.T_SPAR
DEFAULT_T_WEB = spanwise_probe.T_WEB
SPANWISE_KEYS = ("5B-F", "5C-CN", "5D-F")


@dataclass
class FlythroughContext:
    repo_root: Path
    rotor: str
    stations: list[dict[str, Any]] = field(default_factory=list)
    probe_by_key: dict[str, list[dict[str, Any]]] = field(default_factory=dict)

    def ensure_spanwise(self, keys: tuple[str, ...] = SPANWISE_KEYS) -> None:
        if not self.stations:
            self.stations = build_stations(
                repo_root=self.repo_root,
                rotor=self.rotor,
                materials=DEFAULT_MATERIALS,
                structural_system_key="5D-F",
            )
        for key in keys:
            if key not in self.probe_by_key:
                self.probe_by_key[key] = spanwise_probe.probe_layout(key, self.stations)


BuilderFn = Callable[[ScenarioSpec, RenderProfile, argparse.Namespace, ScenarioPaths, FlythroughContext], dict[str, Any]]


def _plot_strip_family(ax, strips: list[object], title: str) -> None:
    colors = {"skin": "#1E88E5", "web": "#43A047", "cap": "#E53935"}
    for strip in strips:
        points = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if points.shape[0] < 2:
            continue
        horizontal, vertical = section_plot_axes(points)
        ax.plot(horizontal, vertical, color=colors.get(str(strip.kind), "#607D8B"), lw=1.9)
    apply_structural_section_axes(ax)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.2)
    ax.set_title(title, fontsize=10)


def _max_web_skin_gap(strips: list[object]) -> float:
    skins = [
        np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for strip in strips
        if str(strip.kind) == "skin"
    ]
    if not skins:
        return float("nan")
    skin_points = np.vstack(skins)
    max_gap = 0.0
    for strip in strips:
        if str(strip.kind) != "web":
            continue
        points = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if points.shape[0] < 2:
            continue
        for row in (0, -1):
            point = points[row]
            distance = float(np.min(np.linalg.norm(skin_points - point.reshape(1, 2), axis=1)))
            max_gap = max(max_gap, distance)
    return max_gap


def _add_footer(fig: plt.Figure, text: str) -> None:
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=9)


def _plot_root_anchor_overlay(
    ax,
    airfoil,
    root_all: list[float],
    root_survivors: list[float],
) -> None:
    vertices = np.asarray(airfoil.vertices, dtype=float)
    closed = np.vstack([vertices, vertices[:1]])
    horizontal, vertical = rotated_naca_grid_to_section_plot_axes(
        closed[:, 0], closed[:, 1], twist_rad=0.0, chord_m=float(airfoil.chord)
    )
    ax.plot(horizontal, vertical, color="black", lw=1.4)

    y_low = float(vertices[:, 1].min()) * 1.25
    y_high = float(vertices[:, 1].max()) * 1.25

    for x_pos in root_all:
        x_values = np.array([x_pos, x_pos], dtype=float)
        y_values = np.array([y_low, y_high], dtype=float)
        h_line, v_line = rotated_naca_grid_to_section_plot_axes(
            x_values,
            y_values,
            twist_rad=0.0,
            chord_m=float(airfoil.chord),
        )
        ax.plot(h_line, v_line, color="#9E9E9E", lw=1.0, linestyle="--", alpha=0.8)

    for x_pos in root_survivors:
        x_values = np.array([x_pos, x_pos], dtype=float)
        y_values = np.array([y_low, y_high], dtype=float)
        h_line, v_line = rotated_naca_grid_to_section_plot_axes(
            x_values,
            y_values,
            twist_rad=0.0,
            chord_m=float(airfoil.chord),
        )
        ax.plot(h_line, v_line, color="#D32F2F", lw=1.5)

    apply_structural_section_axes(ax)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.2)


def _build_root_section(layout, airfoil, root_survivors: list[float]):
    if not root_survivors:
        return None
    try:
        return MultiCellSection(
            airfoil_sdf=airfoil,
            web_x_positions=sorted(float(value) for value in root_survivors),
            web_thickness=float(DEFAULT_T_WEB),
            web_alignment=str(layout.web_orientation),
            cap_height=float(DEFAULT_T_SPAR),
            skin_thickness=float(DEFAULT_T_SKIN),
            twist_angle=0.0,
            structural_family=str(layout.structural_family),
            core_enabled=True,
        )
    except Exception:
        return None


def _build_station_map(stations: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    return {int(station["idx"]): station for station in stations}


def _records_by_station(records: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    return {int(record["idx"]): record for record in records if "idx" in record}


def _render_anchoring_cascade(
    spec: ScenarioSpec,
    profile: RenderProfile,
    args: argparse.Namespace,
    paths: ScenarioPaths,
    context: FlythroughContext,
) -> dict[str, Any]:
    context.ensure_spanwise(keys=("5D-F",))
    records = context.probe_by_key["5D-F"]
    station_map = _build_station_map(context.stations)
    layout = resolve_system_type("5D-F")

    indices = select_frame_indices(len(records), args.max_frames)
    frame_files: list[Path] = []
    source_indices: list[int] = []

    for out_idx, source_idx in enumerate(indices):
        record = records[source_idx]
        if "error" in record:
            continue
        station = station_map.get(int(record["idx"]))
        if station is None:
            continue
        airfoil = airfoil_sdf_from_station(station["coords"], station["metadata"])
        grid = SDFGrid.from_airfoil(airfoil, nx=profile.grid_nx, ny=profile.grid_ny)

        section_station = build_section_view(
            airfoil,
            layout,
            twist_rad=0.0,
            skin_thickness=float(DEFAULT_T_SKIN),
            spar_cap_thickness=float(DEFAULT_T_SPAR),
            shear_web_thickness=float(DEFAULT_T_WEB),
        )

        root_state = record.get("root_anchored") or {}
        root_survivors = [float(value) for value in root_state.get("xs_survivors_m", [])]
        root_all = [float(value) for value in root_state.get("xs_root_m", [])]
        root_section = _build_root_section(layout, airfoil, root_survivors)

        fig, axes = plt.subplots(1, 2, figsize=(14.6, 5.5))
        plot_section(
            section_station,
            grid,
            ax=axes[0],
            alpha=0.52,
            show_airfoil=False,
            show_legend=False,
            auto_limits=True,
            title="Station-relative geometry",
        )

        if root_section is not None:
            plot_section(
                root_section,
                grid,
                ax=axes[1],
                alpha=0.52,
                show_airfoil=False,
                show_legend=False,
                auto_limits=True,
                title="Root-anchored counterfactual",
            )
        else:
            _plot_root_anchor_overlay(axes[1], airfoil, root_all, root_survivors)
            axes[1].set_title("Root anchors (line overlay)")
            axes[1].text(
                0.02,
                0.96,
                "Full section infeasible\nfor this survivor count",
                transform=axes[1].transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
            )

        figure_title = (
            f"{spec.title} | station {record['idx']} (z/L={record['zL']:.3f}) "
            f"chord={record['chord'] * 1e3:.1f} mm"
        )
        fig.suptitle(figure_title, fontsize=12)
        footer = (
            f"station-relative survivors={record['n_webs_resolved']}/{record['n_webs_nominal']} "
            f"(dropped={record['n_dropped']}); root-anchored survivors={root_state.get('n_survivors', 0)}"
        )
        _add_footer(fig, footer)

        frame_file = frame_path(paths.frame_dir, out_idx)
        fig.tight_layout(rect=[0.0, 0.04, 1.0, 0.95])
        save_figure(fig, frame_file, dpi=profile.dpi)
        plt.close(fig)
        frame_files.append(frame_file)
        source_indices.append(source_idx)

    held_frames = apply_holds(
        frame_files,
        head=profile.hold_head,
        tail=profile.hold_tail,
    )
    encode_gif(held_frames, paths.gif_path, profile)

    return {
        "scenario": spec.name,
        "layout_key": "5D-F",
        "source_frame_indices": source_indices,
        "raw_frame_count": len(frame_files),
        "encoded_frame_count": len(held_frames),
        "notes": "Left panel is station-relative implicit geometry; right panel is root-anchored counterfactual.",
    }


def _render_key_cascade_spanwise(
    spec: ScenarioSpec,
    profile: RenderProfile,
    args: argparse.Namespace,
    paths: ScenarioPaths,
    context: FlythroughContext,
) -> dict[str, Any]:
    context.ensure_spanwise(keys=SPANWISE_KEYS)
    station_map = _build_station_map(context.stations)
    station_ids = sorted(station_map)
    selected = select_frame_indices(len(station_ids), args.max_frames)

    probe_maps = {key: _records_by_station(context.probe_by_key[key]) for key in SPANWISE_KEYS}
    frame_files: list[Path] = []
    source_indices: list[int] = []

    for out_idx, pos in enumerate(selected):
        station_id = int(station_ids[pos])
        station = station_map[station_id]
        airfoil = airfoil_sdf_from_station(station["coords"], station["metadata"])
        grid = SDFGrid.from_airfoil(airfoil, nx=profile.grid_nx, ny=profile.grid_ny)

        fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.4))
        for axis, key in zip(axes, SPANWISE_KEYS):
            record = probe_maps[key].get(station_id, {})
            layout = resolve_system_type(key)
            try:
                section = build_section_view(
                    airfoil,
                    layout,
                    twist_rad=0.0,
                    skin_thickness=float(DEFAULT_T_SKIN),
                    spar_cap_thickness=float(DEFAULT_T_SPAR),
                    shear_web_thickness=float(DEFAULT_T_WEB),
                )
                plot_section(
                    section,
                    grid,
                    ax=axis,
                    alpha=0.52,
                    show_airfoil=False,
                    show_legend=False,
                    auto_limits=True,
                    title=key,
                )
            except Exception as exc:
                _plot_root_anchor_overlay(axis, airfoil, [], [])
                axis.set_title(key)
                axis.text(
                    0.5,
                    0.5,
                    f"Unavailable\n{type(exc).__name__}",
                    transform=axis.transAxes,
                    ha="center",
                    va="center",
                    fontsize=9,
                )

            root_state = record.get("root_anchored") or {}
            cap_info = record.get("cap") or {}
            axis.text(
                0.02,
                0.98,
                (
                    f"sr={record.get('n_webs_resolved', 0)}/{record.get('n_webs_nominal', 0)}\n"
                    f"root={root_state.get('n_survivors', 0)}\n"
                    f"cap={cap_info.get('cap_h_allowed_m', float('nan')) * 1e3:.2f} mm"
                ),
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
            )

        fig.suptitle(
            f"{spec.title} | station {station_id} (z/L={station['zL']:.3f})",
            fontsize=12,
        )
        _add_footer(
            fig,
            "Panels show station-relative geometry with root-survivor and cap-allowance diagnostics.",
        )

        frame_file = frame_path(paths.frame_dir, out_idx)
        fig.tight_layout(rect=[0.0, 0.04, 1.0, 0.95])
        save_figure(fig, frame_file, dpi=profile.dpi)
        plt.close(fig)
        frame_files.append(frame_file)
        source_indices.append(pos)

    held_frames = apply_holds(frame_files, head=profile.hold_head, tail=profile.hold_tail)
    encode_gif(held_frames, paths.gif_path, profile)

    return {
        "scenario": spec.name,
        "keys": list(SPANWISE_KEYS),
        "source_frame_indices": source_indices,
        "raw_frame_count": len(frame_files),
        "encoded_frame_count": len(held_frames),
    }


def _render_topology_twist(
    spec: ScenarioSpec,
    profile: RenderProfile,
    args: argparse.Namespace,
    paths: ScenarioPaths,
    _context: FlythroughContext,
) -> dict[str, Any]:
    twists_deg = np.linspace(0.0, 30.0, 31)
    selected = select_frame_indices(len(twists_deg), args.max_frames)
    airfoil = AirfoilSDF.from_naca("2412", chord=1.0)

    frame_files: list[Path] = []
    source_indices: list[int] = []

    for out_idx, source_idx in enumerate(selected):
        twist_deg = float(twists_deg[source_idx])
        twist_rad = float(np.deg2rad(twist_deg))
        section = MultiCellSection.torsion_box(
            airfoil,
            x_fore=0.15,
            x_mid=0.35,
            x_aft=0.55,
            skin_thickness=0.005,
            web_thickness=0.006,
            cap_height=(0.014, 0.012),
            web_alignment="chord_normal",
            twist_angle=twist_rad,
            structural_family="D",
            core_enabled=False,
        )
        raw = build_shell_midline_strips(section, twist_rad=twist_rad)
        snapped = prepare_shell_midline_junctions(raw)

        raw_gap = _max_web_skin_gap(raw)
        snapped_gap = _max_web_skin_gap(snapped)

        fig, axes = plt.subplots(1, 2, figsize=(14.4, 5.4))
        _plot_strip_family(axes[0], raw, f"Raw strips at {twist_deg:.1f} deg")
        _plot_strip_family(axes[1], snapped, f"Topology snapped at {twist_deg:.1f} deg")
        fig.suptitle(spec.title, fontsize=12)
        _add_footer(fig, f"max gap raw={raw_gap:.3e} m, snapped={snapped_gap:.3e} m")

        frame_file = frame_path(paths.frame_dir, out_idx)
        fig.tight_layout(rect=[0.0, 0.04, 1.0, 0.95])
        save_figure(fig, frame_file, dpi=profile.dpi)
        plt.close(fig)
        frame_files.append(frame_file)
        source_indices.append(source_idx)

    held_frames = apply_holds(frame_files, head=profile.hold_head, tail=profile.hold_tail)
    encode_gif(held_frames, paths.gif_path, profile)

    return {
        "scenario": spec.name,
        "source_frame_indices": source_indices,
        "raw_frame_count": len(frame_files),
        "encoded_frame_count": len(held_frames),
        "twist_range_deg": [float(twists_deg[0]), float(twists_deg[-1])],
    }


def _render_orientation_comparator(
    spec: ScenarioSpec,
    profile: RenderProfile,
    args: argparse.Namespace,
    paths: ScenarioPaths,
    context: FlythroughContext,
) -> dict[str, Any]:
    context.ensure_spanwise(keys=("5D-F",))
    station_map = _build_station_map(context.stations)
    station_ids = sorted(station_map)
    selected = select_frame_indices(len(station_ids), args.max_frames)

    layout_cn = resolve_system_type("3D-CN")
    layout_flapwise = resolve_system_type("3D-F")

    frame_files: list[Path] = []
    source_indices: list[int] = []

    for out_idx, source_idx in enumerate(selected):
        station_id = int(station_ids[source_idx])
        station = station_map[station_id]
        airfoil = airfoil_sdf_from_station(station["coords"], station["metadata"])
        grid = SDFGrid.from_airfoil(airfoil, nx=profile.grid_nx, ny=profile.grid_ny)

        section_cn = build_section_view(
            airfoil,
            layout_cn,
            twist_rad=0.0,
            skin_thickness=float(DEFAULT_T_SKIN),
            spar_cap_thickness=float(DEFAULT_T_SPAR),
            shear_web_thickness=float(DEFAULT_T_WEB),
        )
        section_flapwise = build_section_view(
            airfoil,
            layout_flapwise,
            twist_rad=0.0,
            skin_thickness=float(DEFAULT_T_SKIN),
            spar_cap_thickness=float(DEFAULT_T_SPAR),
            shear_web_thickness=float(DEFAULT_T_WEB),
        )

        fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.4))
        plot_section(
            section_cn,
            grid,
            ax=axes[0],
            alpha=0.52,
            show_airfoil=False,
            show_legend=False,
            auto_limits=True,
            title="3D-CN (chord-normal)",
        )
        plot_section(
            section_flapwise,
            grid,
            ax=axes[1],
            alpha=0.52,
            show_airfoil=False,
            show_legend=False,
            auto_limits=True,
            title="3D-F (flapwise)",
        )

        fig.suptitle(
            f"{spec.title} | station {station_id} (z/L={station['zL']:.3f})",
            fontsize=12,
        )
        _add_footer(fig, "Matched station and thicknesses; only web orientation changes.")

        frame_file = frame_path(paths.frame_dir, out_idx)
        fig.tight_layout(rect=[0.0, 0.04, 1.0, 0.95])
        save_figure(fig, frame_file, dpi=profile.dpi)
        plt.close(fig)
        frame_files.append(frame_file)
        source_indices.append(source_idx)

    held_frames = apply_holds(frame_files, head=profile.hold_head, tail=profile.hold_tail)
    encode_gif(held_frames, paths.gif_path, profile)

    return {
        "scenario": spec.name,
        "source_frame_indices": source_indices,
        "raw_frame_count": len(frame_files),
        "encoded_frame_count": len(held_frames),
        "keys": ["3D-CN", "3D-F"],
    }


def _render_insert_dynamics(
    spec: ScenarioSpec,
    profile: RenderProfile,
    args: argparse.Namespace,
    paths: ScenarioPaths,
    _context: FlythroughContext,
) -> dict[str, Any]:
    airfoil = AirfoilSDF.from_naca("2412", chord=1.0)
    grid = SDFGrid.from_airfoil(airfoil.rotate(float(np.deg2rad(15.0))), nx=profile.grid_nx, ny=profile.grid_ny)
    quadtree_spec = quadtree_spec_for_airfoil(
        airfoil,
        chord=1.0,
        max_depth=7,
        max_cells=140000,
    )

    forward_le = np.linspace(0.07, 0.15, 10)
    forward_te = np.linspace(0.72, 0.84, 10)
    sequence: list[tuple[float, float]] = []
    for le_value, te_value in zip(forward_le, forward_te):
        sequence.append((float(le_value), float(te_value)))
    for le_value, te_value in zip(forward_le[-2:0:-1], forward_te[-2:0:-1]):
        sequence.append((float(le_value), float(te_value)))

    selected = select_frame_indices(len(sequence), args.max_frames)
    frame_files: list[Path] = []
    source_indices: list[int] = []

    for out_idx, source_idx in enumerate(selected):
        le_x, te_x = sequence[source_idx]
        section = MultiCellSection.twin_web(
            airfoil,
            x_fore=0.20,
            x_aft=0.50,
            skin_thickness=float(DEFAULT_T_SKIN),
            cap_height=float(DEFAULT_T_SPAR),
            web_thickness=float(DEFAULT_T_WEB),
            web_alignment="chord_normal",
            trailing_edge_insert_x_m=te_x,
            leading_edge_insert_x_m=le_x,
            structural_family="D",
            twist_angle=float(np.deg2rad(15.0)),
            core_enabled=True,
        )
        _props, metadata, mesh = compute_section_properties(
            section,
            backend="hybrid",
            quadtree_spec=quadtree_spec,
            with_mesh=True,
            chord=1.0,
            twist_rad=float(np.deg2rad(15.0)),
            airfoil_for_spec=airfoil,
        )

        fig, axes = plt.subplots(1, 2, figsize=(14.6, 5.4))
        plot_section(
            section,
            grid,
            ax=axes[0],
            alpha=0.52,
            show_airfoil=False,
            show_legend=False,
            auto_limits=True,
            title="Geometry with moving LE/TE inserts",
        )
        plot_section_hybrid(
            section,
            residual_mesh=mesh,
            ax=axes[1],
            alpha=0.52,
            show_legend=False,
            title="Hybrid constructive + residual path",
        )

        residual_labels = metadata.get("residual_labels", [])
        fig.suptitle(spec.title, fontsize=12)
        _add_footer(
            fig,
            f"LE={le_x:.3f} m, TE={te_x:.3f} m, residual labels={','.join(str(v) for v in residual_labels)}",
        )

        frame_file = frame_path(paths.frame_dir, out_idx)
        fig.tight_layout(rect=[0.0, 0.04, 1.0, 0.95])
        save_figure(fig, frame_file, dpi=profile.dpi)
        plt.close(fig)
        frame_files.append(frame_file)
        source_indices.append(source_idx)

    held_frames = apply_holds(frame_files, head=profile.hold_head, tail=profile.hold_tail)
    encode_gif(held_frames, paths.gif_path, profile)

    return {
        "scenario": spec.name,
        "source_frame_indices": source_indices,
        "raw_frame_count": len(frame_files),
        "encoded_frame_count": len(held_frames),
    }


def _render_quadtree_refinement(
    spec: ScenarioSpec,
    profile: RenderProfile,
    args: argparse.Namespace,
    paths: ScenarioPaths,
    _context: FlythroughContext,
) -> dict[str, Any]:
    depth_sequence = [4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4]
    selected = select_frame_indices(len(depth_sequence), args.max_frames)

    airfoil = AirfoilSDF.from_naca("2412", chord=1.0)
    layout = resolve_system_type("2D-F")
    section = build_section_view(
        airfoil,
        layout,
        twist_rad=float(np.deg2rad(15.0)),
        skin_thickness=float(DEFAULT_T_SKIN),
        spar_cap_thickness=float(DEFAULT_T_SPAR),
        shear_web_thickness=float(DEFAULT_T_WEB),
    )

    frame_files: list[Path] = []
    source_indices: list[int] = []

    for out_idx, source_idx in enumerate(selected):
        depth = int(depth_sequence[source_idx])
        quadtree_spec = quadtree_spec_for_airfoil(
            airfoil,
            chord=1.0,
            max_depth=depth,
            max_cells=180000,
        )
        mesh = build_shared_section_quadtree(
            section,
            quadtree_spec,
            twist_rad=float(np.deg2rad(15.0)),
        )
        leaves = mesh.leaf_cells()

        fig, _axes = plot_section_with_quadtree(
            mesh,
            section=section,
            title=f"Geometry view (depth <= {depth})",
            mesh_title=f"Union quadtree mesh (depth <= {depth})",
            show_legend=False,
            alpha=0.50,
            figsize=(15.2, 5.5),
            twist_rad=float(np.deg2rad(15.0)),
            chord_m=1.0,
        )
        fig.suptitle(spec.title, fontsize=12)
        _add_footer(fig, f"depth={depth}, leaves={len(leaves)}")

        frame_file = frame_path(paths.frame_dir, out_idx)
        fig.tight_layout(rect=[0.0, 0.04, 1.0, 0.95])
        save_figure(fig, frame_file, dpi=profile.dpi)
        plt.close(fig)
        frame_files.append(frame_file)
        source_indices.append(source_idx)

    held_frames = apply_holds(frame_files, head=profile.hold_head, tail=profile.hold_tail)
    encode_gif(held_frames, paths.gif_path, profile)

    return {
        "scenario": spec.name,
        "source_frame_indices": source_indices,
        "raw_frame_count": len(frame_files),
        "encoded_frame_count": len(held_frames),
        "depth_sequence": [int(depth_sequence[index]) for index in source_indices],
    }


BUILDERS: dict[str, BuilderFn] = {
    "anchoring_cascade": _render_anchoring_cascade,
    "key_cascade_spanwise": _render_key_cascade_spanwise,
    "topology_twist": _render_topology_twist,
    "orientation_comparator": _render_orientation_comparator,
    "insert_dynamics": _render_insert_dynamics,
    "quadtree_refinement": _render_quadtree_refinement,
}


def _selected_scenarios(requested: str) -> list[str]:
    if requested == "all":
        return list(SCENARIO_ORDER)
    if requested not in SCENARIOS:
        raise SystemExit(f"Unknown scenario {requested!r}; options: {list(SCENARIOS)}")
    return [requested]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", default="all", help="Scenario name or 'all'.")
    parser.add_argument(
        "--profile",
        choices=sorted(PROFILES),
        default=DEFAULT_PROFILE,
        help="Render profile preset.",
    )
    parser.add_argument(
        "--rotor",
        default=DEFAULT_ROTOR,
        help="Rotor directory key for spanwise scenarios.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "output" / "flythroughs",
        help="Output directory for GIFs and scenario metadata.",
    )
    parser.add_argument(
        "--keep-frames",
        action="store_true",
        help="Keep per-scenario PNG frame folders.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional cap on frame count after deterministic down-sampling.",
    )
    parser.add_argument(
        "--skip-validate",
        action="store_true",
        help="Skip post-render GIF metadata validation.",
    )
    parser.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    from dataclasses import replace

    args = parse_args(argv)
    apply_presentation_style(bool(args.presentation))
    profile = PROFILES[str(args.profile)]
    if bool(args.presentation) and int(profile.dpi) < PRESENTATION_DPI:
        profile = replace(profile, dpi=PRESENTATION_DPI)
    scenario_names = _selected_scenarios(str(args.scenario))
    context = FlythroughContext(repo_root=Path(repo_root), rotor=str(args.rotor))

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    validations: list[dict[str, Any]] = []

    for scenario_name in scenario_names:
        spec = SCENARIOS[scenario_name]
        paths = prepare_scenario_paths(out_dir, scenario_name, spec.output_name)
        reset_frame_dir(paths.frame_dir)

        print(f"Rendering scenario: {scenario_name}")
        metadata = BUILDERS[scenario_name](spec, profile, args, paths, context)

        scenario_payload = {
            "scenario": scenario_name,
            "title": spec.title,
            "description": spec.description,
            "frame_source": spec.frame_source,
            "panel_layout": spec.panel_layout,
            "target_frame_count": spec.frame_count,
            "sweep_params": spec.sweep_params,
            "metric_overlays": list(spec.metric_overlays),
            "profile": profile.name,
            "width_px": profile.width_px,
            "delay_cs": profile.frame_delay_cs,
            "hold_head": profile.hold_head,
            "hold_tail": profile.hold_tail,
            "rotor": str(args.rotor),
            "max_frames": args.max_frames,
            "output_gif": paths.gif_path.name,
            **metadata,
        }

        if not args.skip_validate:
            validation = validate_gif(paths.gif_path)
            validations.append(validation)
            scenario_payload["validation"] = validation
            print(
                f"  done: {paths.gif_path.name} ({validation['frame_count']} frames, "
                f"{validation['width_px']}x{validation['height_px']})"
            )
        else:
            print(f"  done: {paths.gif_path.name}")

        write_metadata(paths.metadata_path, scenario_payload)

        if not args.keep_frames:
            reset_frame_dir(paths.frame_dir)
            paths.frame_dir.rmdir()

    if not args.keep_frames:
        frames_root = out_dir / "frames"
        if frames_root.exists() and not any(frames_root.iterdir()):
            frames_root.rmdir()

    print("Completed scenarios:")
    for scenario_name in scenario_names:
        print(f"  - {scenario_name}: {SCENARIOS[scenario_name].output_name}")

    if validations:
        print("Validation summary:")
        for row in validations:
            print(
                f"  - {Path(row['path']).name}: frames={row['frame_count']}, "
                f"width={row['width_px']}, loop={row['infinite_loop']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
