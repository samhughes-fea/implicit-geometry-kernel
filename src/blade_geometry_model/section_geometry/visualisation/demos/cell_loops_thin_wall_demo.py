"""
visualisation/demos/cell_loops_thin_wall_demo.py
================================================
Cell-wall loops, for any layout, become the thin-wall equations.

``section_geometry`` emits a closed loop around every empty cell of every
layout, with the walls that bound it named. That is exactly the input the
multi-cell Bredt problem and the shear-centre problem need, so neither has to
be set up by hand per architecture:

* ``section_properties.reference_axes.wall_network.network_from_midlines_props``
  turns midlines + loops into a wall network, one cell per loop;
* ``section_properties.torsion.multi_cell.torsion_from_wall_network`` solves the
  n-cell compatibility system for the cell shear flows and the torsion constant;
* ``section_properties.reference_axes.shear_centre.shear_centre_from_midlines``
  solves for the shear centre on the same loops.

Top row: four layouts with one to six cells, each wall coloured by the net
shear flow it carries under unit twist (shared webs carry the *difference* of
their two cells' flows), shear centre and elastic centroid marked. Bottom row:
one layout with its aft web swept across the chord. Every point is a fresh
geometry build, fresh loops and a fresh solve, so the curves are smooth only
because the loops are regenerated correctly at every parameter value.

Walls carry ``G = 1``, so the torsion result is the geometric constant ``J``.
Every solve is checked for torque equilibrium, ``2 sum(A_i q_i) = J``; a
mismatch raises.

Writes under ``visualisation/output/cell_loops_thin_wall/``:

* ``cell_loops_thin_wall.png``
* ``cell_loops_thin_wall.json``

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/cell_loops_thin_wall_demo.py --presentation
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.coordinate_frames import section_plot_axes
from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SectionGeometryRunner,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    NAVY,
    RED,
    TEAL,
    strip_axes_chrome,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
    apply_structural_section_axes,
)
from blade_structure_utilities.geometry_utilities.section_properties.reference_axes.shear_centre import (
    shear_centre_from_midlines,
)
from blade_structure_utilities.geometry_utilities.section_properties.reference_axes.wall_network import (
    network_from_midlines_props,
)
from blade_structure_utilities.geometry_utilities.section_properties.torsion.multi_cell import (
    torsion_from_wall_network,
)

_GALLERY_KEYS = ("0A", "2D-F", "3B-F", "5A-F")
_SKIN_M = 0.004
_CAP_M = 0.012
_WEB_M = 0.005
_EQUILIBRIUM_RTOL = 1e-9


def _solve(airfoil: AirfoilSDF, key: str, *, web_x_positions_m=None) -> dict:
    """Geometry → loops → wall network → Bredt + shear centre, for one section."""
    result = (
        SectionGeometryRunner(
            layout=resolve_system_type(key),
            airfoil=airfoil,
            skin_thickness=_SKIN_M,
            spar_cap_thickness=_CAP_M,
            shear_web_thickness=_WEB_M,
            twist_rad=0.0,
            web_x_positions_m=web_x_positions_m,
            outputs=("properties", "midlines", "cell_wall_loops"),
        )
        .execute()
        .get_results()
    )
    midlines = result.midline_strips
    loops = result.cell_wall_loops
    if not loops:
        raise RuntimeError(f"{key}: no cell wall loops emitted")

    network = network_from_midlines_props(midlines, cell_wall_loops=loops)
    torsion = torsion_from_wall_network(network)
    torque = 2.0 * sum(area * flow for area, flow in zip(torsion.cell_areas, torsion.cell_shear_flows))
    if not np.isclose(torque, torsion.GJ, rtol=_EQUILIBRIUM_RTOL, atol=0.0):
        raise RuntimeError(f"{key}: torque equilibrium fails, 2ΣAq = {torque:.6e} vs J = {torsion.GJ:.6e}")
    centre = shear_centre_from_midlines(midlines, cell_wall_loops=loops)

    # torsion_from_wall_network keeps, in segment order, only the segments some
    # cell owns; rebuild that ownership to put each wall's flow back on its segment.
    owned = set()
    for cell in network.cells:
        owned.update(int(segment_id) for segment_id in cell.segment_ids)
    owned_segments = [segment for index, segment in enumerate(network.segments) if index in owned]
    if len(owned_segments) != len(torsion.walls):
        raise RuntimeError(f"{key}: {len(owned_segments)} owned segments but {len(torsion.walls)} walls")
    segment_flows = [
        ((segment.x1, segment.y1), (segment.x2, segment.y2), wall.shear_flow)
        for segment, wall in zip(owned_segments, torsion.walls)
    ]

    return {
        "key": key,
        "midlines": midlines,
        "loops": loops,
        "segment_flows": segment_flows,
        "n_cells": int(torsion.n_cells),
        "J_m4": float(torsion.GJ),
        "cell_areas_m2": [float(area) for area in torsion.cell_areas],
        "cell_shear_flows": [float(flow) for flow in torsion.cell_shear_flows],
        "shear_centre": (float(centre.xs), float(centre.ys)),
        "elastic_centroid": tuple(float(value) for value in centre.elastic_centroid),
    }


def _draw_section(axis, solved: dict, norm: Normalize, colour_map) -> None:
    for loop in solved["loops"]:
        polyline = np.asarray(loop.polyline_twisted_pitch_axis_blade_frame, dtype=float)
        horizontal, vertical = section_plot_axes(polyline)
        axis.fill(horizontal, vertical, color="#eef1f5", zorder=1)

    for strip in solved["midlines"]:
        points = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        horizontal, vertical = section_plot_axes(points)
        axis.plot(horizontal, vertical, color="#d5d9df", lw=1.0, zorder=2)

    segments = []
    flows = []
    for start, end, flow in solved["segment_flows"]:
        start_h, start_v = section_plot_axes(np.asarray([start], dtype=float))
        end_h, end_v = section_plot_axes(np.asarray([end], dtype=float))
        segments.append([(start_h[0], start_v[0]), (end_h[0], end_v[0])])
        flows.append(abs(flow))
    collection = LineCollection(segments, cmap=colour_map, norm=norm, linewidths=3.0, zorder=3,
                                capstyle="round")
    collection.set_array(np.asarray(flows))
    axis.add_collection(collection)

    (centroid_h,), (centroid_v,) = section_plot_axes(np.asarray([solved["elastic_centroid"]]))
    (shear_h,), (shear_v,) = section_plot_axes(np.asarray([solved["shear_centre"]]))
    axis.scatter([centroid_h], [centroid_v], s=60, color=INK, zorder=5)
    axis.scatter([shear_h], [shear_v], s=110, marker="X", color=RED, edgecolors="white",
                 linewidths=0.8, zorder=6)

    apply_structural_section_axes(axis)
    axis.autoscale_view()
    axis.set_aspect("equal", adjustable="box")
    strip_axes_chrome(axis)
    n = solved["n_cells"]
    axis.set_title(
        f"{solved['key']}  ·  {n} cell{'s' if n != 1 else ''} → {n}×{n} system",
        fontsize=13, color=INK, loc="left",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="cell_loops_thin_wall")
    parser.add_argument("--sweep-key", default="2D-F")
    parser.add_argument("--sweep-min", type=float, default=0.30, help="Aft web start, fraction of chord.")
    parser.add_argument("--sweep-max", type=float, default=0.72, help="Aft web end, fraction of chord.")
    parser.add_argument("--sweep-count", type=int, default=22)
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    airfoil = AirfoilSDF.from_naca("2412", chord=1.0)
    gallery = [_solve(airfoil, key) for key in _GALLERY_KEYS]

    sweep_layout = resolve_system_type(str(args.sweep_key))
    fixed_fractions = list(sweep_layout.web_chord_fracs[:-1])
    aft_fractions = np.linspace(float(args.sweep_min), float(args.sweep_max), int(args.sweep_count))
    sweep = []
    for fraction in aft_fractions:
        positions = tuple(float(value) for value in (*fixed_fractions, fraction))
        solved = _solve(airfoil, str(args.sweep_key), web_x_positions_m=positions)
        sweep.append({"aft_web_fraction": float(fraction), **solved})

    peak_flow = max(abs(flow) for solved in gallery for _, _, flow in solved["segment_flows"])
    norm = Normalize(vmin=0.0, vmax=peak_flow)
    colour_map = plt.get_cmap("magma_r")

    fig = plt.figure(figsize=(15.0, 6.2))
    outer = fig.add_gridspec(3, 1, height_ratios=(0.42, 0.05, 1.0), hspace=0.85)
    top = outer[0].subgridspec(1, 4, wspace=0.06)
    middle = outer[1].subgridspec(1, 4, wspace=0.06)
    bottom = outer[2].subgridspec(1, 2, wspace=0.2)
    for column, solved in enumerate(gallery):
        _draw_section(fig.add_subplot(top[0, column]), solved, norm, colour_map)

    colourbar_axis = fig.add_subplot(middle[0, 1:3])
    colourbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=Normalize(0.0, peak_flow * 1e4), cmap=colour_map),
        cax=colourbar_axis, orientation="horizontal",
    )
    colourbar.set_label("|net shear flow| per unit twist rate  [×10⁻⁴, G = 1]", fontsize=11, color=INK)
    colourbar.ax.xaxis.set_label_position("top")
    colourbar.outline.set_visible(False)
    legend_axis = fig.add_subplot(middle[0, 3])
    legend_axis.axis("off")
    legend_axis.scatter([0.05], [0.5], s=110, marker="X", color=RED)
    legend_axis.text(0.1, 0.5, "shear centre", va="center", fontsize=11, color=INK)
    legend_axis.scatter([0.55], [0.5], s=60, color=INK)
    legend_axis.text(0.6, 0.5, "elastic centroid", va="center", fontsize=11, color=INK)
    legend_axis.set_xlim(0, 1)
    legend_axis.set_ylim(0, 1)

    fractions = np.array([row["aft_web_fraction"] for row in sweep])
    stiffness = np.array([row["J_m4"] for row in sweep])
    shear_edgewise = np.array([row["shear_centre"][1] for row in sweep])
    centroid_edgewise = np.array([row["elastic_centroid"][1] for row in sweep])

    ax_stiffness = fig.add_subplot(bottom[0, 0])
    ax_stiffness.plot(fractions, stiffness * 1e6, "o-", color=NAVY, lw=2.4, ms=5)
    ax_stiffness.axvline(sweep_layout.web_chord_fracs[-1], color=MUTED, ls=":", lw=1.2)
    ax_stiffness.set_xlabel(f"{args.sweep_key} aft web position [fraction of chord]")
    ax_stiffness.set_ylabel("torsion constant J  [×10⁻⁶ m⁴]")
    ax_stiffness.set_title(
        f"{len(sweep)} rebuilds → {len(sweep)} {sweep[0]['n_cells']}-cell Bredt solves", loc="left",
    )
    ax_stiffness.grid(True, alpha=0.4)

    ax_centre = fig.add_subplot(bottom[0, 1])
    ax_centre.plot(fractions, shear_edgewise * 1e3, "o-", color=RED, lw=2.4, ms=5, label="shear centre")
    ax_centre.plot(fractions, centroid_edgewise * 1e3, "o-", color=INK, lw=1.6, ms=4, label="elastic centroid")
    ax_centre.axvline(sweep_layout.web_chord_fracs[-1], color=MUTED, ls=":", lw=1.2)
    ax_centre.set_xlabel(f"{args.sweep_key} aft web position [fraction of chord]")
    ax_centre.set_ylabel("edgewise position Y_B [mm]")
    ax_centre.set_title("The shear centre follows the box", loc="left")
    ax_centre.legend(frameon=False, fontsize=11)
    ax_centre.grid(True, alpha=0.4)

    png_path = out_dir / "cell_loops_thin_wall.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    def summary(row: dict) -> dict:
        return {name: value for name, value in row.items() if name not in ("midlines", "loops", "segment_flows")}

    payload = {
        "airfoil": "NACA 2412, chord 1 m",
        "gauges_m": {"skin": _SKIN_M, "spar_cap": _CAP_M, "shear_web": _WEB_M},
        "shear_modulus": "G = 1 on every wall; J is geometric",
        "equilibrium_rtol": _EQUILIBRIUM_RTOL,
        "gallery": [summary(row) for row in gallery],
        "sweep_key": str(args.sweep_key),
        "sweep": [summary(row) for row in sweep],
    }
    json_path = out_dir / "cell_loops_thin_wall.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
