"""
visualisation/demos/downstream_consumers_demo.py
================================================
One station, three downstream consumers of the same geometry result.

Everything here is the production chain, not a stand-in: catalogue composites,
realised gauges, ``SectionGeometryRunner`` for midlines + cell-wall loops +
per-region properties, then

1. **mesh → FE.** ``section_mesh`` strips → ``SectionBindRunner`` frozen zones →
   ``SectionLinearStaticRunner`` (MITC4, thin-wall load channel) → Hashin FI.
   Drawn: the peak ply FI of every element, on the element.
2. **loads.** ``SectionLoadsRunner`` shapes the beam resultants onto the
   midlines. It needs the cell-wall loops for the closed-cell shear flow.
   Drawn: the shear flow ``q`` round the section.
3. **properties.** ``SectionPropertiesRunner`` aggregates the per-region
   integrals with the materials and solves Bredt torsion on the loops.
   Drawn: the section with its elastic centroid and shear centre, and the
   stiffnesses in the title — the beam the section becomes.

Load case: flapwise bending moment and shear at one station, stated on the figure.

Writes under ``visualisation/output/downstream_consumers/``:

* ``downstream_consumers.png``
* ``downstream_consumers.json`` — the numbers the slide cards quote

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/downstream_consumers_demo.py --presentation
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import time
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

from blade_structure_inputs.module.materials.materials import Materials
from blade_structure_utilities.fem_utilities.section_shell_fea import SectionLinearStaticRunner
from blade_structure_utilities.fem_utilities.section_shell_fea.pre_processing.mesh.api.bound_section_shell_inputs import (
    bind_section_shell_inputs,
)
from blade_structure_utilities.geometry_utilities.coordinate_frames import section_plot_axes
from blade_structure_utilities.geometry_utilities.section_contract import BeamResultantVector
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
from blade_structure_utilities.geometry_utilities.section_properties.api.runner import (
    SectionPropertiesRunner,
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
from blade_structure_utilities.load_utilities.section_loads.api.runner import SectionLoadsRunner
from blade_structure_utilities.material_utilities.section_bind import (
    SectionBindRunner,
    realized_thicknesses_for_fi,
)
from blade_structure_utilities.mesh_utilities.section_mesh import shell_mesh_inputs_from_geometry_result
from blade_structure_utilities.mesh_utilities.section_mesh.api import (
    MIN_ELEMENTS_PER_PANEL,
    TARGET_ELEMENT_LENGTH_FRAC,
    build_section_strip_mesh,
)

_LAYOUT_KEY = "2D-F"
_SKIN_MATERIAL_ID = 13
_SPAR_MATERIAL_ID = 17
_WEB_MATERIAL_ID = 13


def _run_chain(*, moment_flapwise: float, shear_flapwise: float) -> dict:
    catalog = Materials().load()
    materials = {
        "skin": catalog.get_material_by_id(_SKIN_MATERIAL_ID),
        "spar": catalog.get_material_by_id(_SPAR_MATERIAL_ID),
        "web": catalog.get_material_by_id(_WEB_MATERIAL_ID),
    }
    skin_t, spar_t, web_t = realized_thicknesses_for_fi(
        materials, skin_thickness=0.004, spar_cap_thickness=0.012, shear_web_thickness=0.005
    )

    started = time.perf_counter()
    geometry = (
        SectionGeometryRunner(
            layout=resolve_system_type(_LAYOUT_KEY),
            airfoil=AirfoilSDF.from_naca("2412", chord=1.0),
            skin_thickness=skin_t,
            spar_cap_thickness=spar_t,
            shear_web_thickness=web_t,
            outputs=("properties", "midlines", "cell_wall_loops"),
        )
        .execute()
        .get_results()
    )
    geometry_seconds = time.perf_counter() - started

    shell = shell_mesh_inputs_from_geometry_result(geometry, layout_key=_LAYOUT_KEY)
    binding = (
        SectionBindRunner(
            shell_inputs=shell,
            component_materials=materials,
            skin_thickness=skin_t,
            spar_cap_thickness=spar_t,
            shear_web_thickness=web_t,
        )
        .execute()
        .get_results()
    )
    topology = bind_section_shell_inputs(shell, zones=binding.zones, n_cells=binding.n_cells)
    mesh = build_section_strip_mesh(
        topology.shell_inputs,
        target_element_length_frac=TARGET_ELEMENT_LENGTH_FRAC,
        min_elements_per_panel=MIN_ELEMENTS_PER_PANEL,
    )
    resultants = BeamResultantVector.from_flap_edge(M_flapwise=moment_flapwise, V_flapwise=shear_flapwise)
    section_forces = (
        SectionLoadsRunner(
            resultants=resultants,
            midlines=topology.shell_inputs.midlines,
            cell_wall_loops=geometry.cell_wall_loops,
        )
        .execute()
        .get_results()
    )
    solve = (
        SectionLinearStaticRunner(
            panel_topology=topology,
            resultants=resultants,
            section_force_vector=section_forces,
            section_strip_mesh=mesh,
            load_channel="thinwall",
            compute_gaussian=True,
        )
        .execute()
        .get_results()
    )
    # EA / EI: aggregation stays in the geometry worker's own (NACA) frame,
    # where the per-region integrals live, so no B-frame midlines go in.
    properties = (
        SectionPropertiesRunner(
            properties_map=geometry.properties_map,
            materials=materials,
            jt="skipped",
            attach_axes=False,
        )
        .execute()
        .get_results()
    )
    # GJ: multi-cell Bredt on the loops, each wall carrying its bound laminate's
    # in-plane shear stiffness A66. Zone labels are "<kind>:<wall label>".
    shear_stiffness_by_label = {
        zone.label.split(":", 1)[1]: float(np.asarray(zone.lam.build_ABD())[2, 2]) for zone in binding.zones
    }
    network = network_from_midlines_props(geometry.midline_strips, cell_wall_loops=geometry.cell_wall_loops)
    torsion = torsion_from_wall_network(network, shear_stiffness_by_label=shear_stiffness_by_label)
    if torsion.meta.get("gj_shear_source") != "shear_stiffness":
        raise RuntimeError(f"GJ not wholly from laminate A66: {torsion.meta.get('gj_shear_source')}")
    # Centroid and shear centre in the B frame the midlines use.
    centres = shear_centre_from_midlines(
        geometry.midline_strips, materials=materials, cell_wall_loops=geometry.cell_wall_loops
    )
    return {
        "materials": materials,
        "gauges_m": {"skin": skin_t, "spar_cap": spar_t, "shear_web": web_t},
        "geometry": geometry,
        "geometry_seconds": geometry_seconds,
        "mesh": mesh,
        "resultants": resultants,
        "section_forces": section_forces,
        "solve": solve,
        "properties": properties,
        "torsion": torsion,
        "centres": centres,
    }


def _element_segments_and_fi(mesh, solve) -> tuple[list, np.ndarray]:
    """One plot segment and one peak-ply FI per element, in the solver's row order."""
    fi_rows = np.asarray(solve.quaternary_results.hashin_fi.in_plane, dtype=float)
    fi_rows = fi_rows.reshape(fi_rows.shape[0], -1)
    rows = list(solve.secondary_results.elemental_results.resultant_rows)
    if len(rows) != fi_rows.shape[0] or len(rows) != int(mesh.n_total_elements):
        raise RuntimeError(
            f"{len(rows)} result rows, {fi_rows.shape[0]} FI rows, {mesh.n_total_elements} elements"
        )
    segments: list = []
    row_index = 0
    for panel in mesh.panel_meshes:
        nodes = np.asarray(panel.yz_nodes, dtype=float)
        for element in range(int(panel.n_elements)):
            if int(rows[row_index].panel_index) != int(panel.panel_index):
                raise RuntimeError(f"row {row_index} belongs to panel {rows[row_index].panel_index}")
            horizontal, vertical = section_plot_axes(nodes[element : element + 2])
            segments.append(np.column_stack([horizontal, vertical]))
            row_index += 1
    return segments, np.nanmax(fi_rows, axis=1)


def _finish_section_axis(axis, title: str) -> None:
    apply_structural_section_axes(axis)
    axis.autoscale_view()
    axis.set_aspect("equal", adjustable="box")
    strip_axes_chrome(axis)
    axis.set_title(title, loc="left", fontsize=14, color=INK, pad=4)


def _outline(axis, geometry) -> None:
    for strip in geometry.midline_strips:
        points = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        horizontal, vertical = section_plot_axes(points)
        axis.plot(horizontal, vertical, color="#d5d9df", lw=1.0, zorder=1)


def _row_colourbar(fig, axis, mappable, label: str) -> None:
    colourbar = fig.colorbar(mappable, ax=axis, orientation="vertical", fraction=0.02, pad=0.01)
    colourbar.set_label(label, fontsize=11, color=INK)
    colourbar.outline.set_visible(False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="downstream_consumers")
    parser.add_argument("--moment-flapwise", type=float, default=8.0e4, help="M_flapwise [N m].")
    parser.add_argument("--shear-flapwise", type=float, default=1.0e4, help="V_flapwise [N].")
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    chain = _run_chain(moment_flapwise=float(args.moment_flapwise), shear_flapwise=float(args.shear_flapwise))
    geometry = chain["geometry"]
    mesh = chain["mesh"]
    solve = chain["solve"]
    section_forces = chain["section_forces"]
    properties = chain["properties"]
    centres = chain["centres"]
    torsion = chain["torsion"]

    fig, axes = plt.subplots(3, 1, figsize=(13.0, 7.6), gridspec_kw={"hspace": 0.32})
    ax_fe, ax_loads, ax_properties = axes

    # 1 · FE: peak ply Hashin FI on every element.
    segments, element_fi = _element_segments_and_fi(mesh, solve)
    fi_norm = Normalize(vmin=0.0, vmax=max(float(np.nanmax(element_fi)), 1e-12))
    fi_lines = LineCollection(segments, cmap="magma_r", norm=fi_norm, linewidths=4.0, capstyle="butt")
    fi_lines.set_array(element_fi)
    ax_fe.add_collection(fi_lines)
    governing = dict(solve.quaternary_results.governing or {})
    _finish_section_axis(
        ax_fe,
        # BL-044 / FEA-EQ-002: the shell model currently carries no axial
        # membrane force, so this FI is the shear-carried part only.
        f"section_mesh → section_shell_fea · {int(mesh.n_total_elements)} MITC4 elements · "
        f"Hashin FI {float(governing.get('fi', np.nanmax(element_fi))):.2f}, shear only (BL-044)",
    )
    _row_colourbar(fig, ax_fe, fi_lines, "peak ply FI")

    # 2 · Loads: shear flow round the closed cells.
    samples = list(section_forces.samples)
    load_points = np.array([[sample.y, sample.z] for sample in samples], dtype=float)
    load_horizontal, load_vertical = section_plot_axes(load_points)
    shear_flow = np.array([abs(sample.q) for sample in samples], dtype=float) * 1e-3
    _outline(ax_loads, geometry)
    loads_scatter = ax_loads.scatter(
        load_horizontal, load_vertical, c=shear_flow, cmap="viridis", s=9, linewidths=0, zorder=3,
    )
    reconstruction_error = float(
        np.max(np.abs(section_forces.reconstructed.R - section_forces.resultants.R))
        / max(float(np.max(np.abs(section_forces.resultants.R))), 1e-30)
    )
    _finish_section_axis(
        ax_loads,
        f"section_loads · closed-cell shear flow from the cell-wall loops · "
        f"resultants recovered to {reconstruction_error:.0e}",
    )
    _row_colourbar(fig, ax_loads, loads_scatter, "|q| [kN/m]")

    # 3 · Properties: the beam the section becomes.
    _outline(ax_properties, geometry)
    for strip in geometry.midline_strips:
        points = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        horizontal, vertical = section_plot_axes(points)
        ax_properties.plot(horizontal, vertical, color=NAVY, lw=2.2, zorder=2)
    (centroid_h,), (centroid_v,) = section_plot_axes(np.array([centres.elastic_centroid], dtype=float))
    (shear_h,), (shear_v,) = section_plot_axes(np.array([[centres.xs, centres.ys]], dtype=float))
    ax_properties.scatter([centroid_h], [centroid_v], s=80, color=INK, zorder=5, label="elastic centroid")
    ax_properties.scatter([shear_h], [shear_v], s=130, marker="X", color=RED, edgecolors="white",
                          linewidths=0.8, zorder=6, label="shear centre")
    _finish_section_axis(
        ax_properties,
        f"section_properties · EI_flap {properties.EI_flap / 1e6:.2f} MN·m², "
        f"EI_edge {properties.EI_edge / 1e6:.1f} MN·m², GJ {torsion.GJ / 1e6:.2f} MN·m² ({torsion.n_cells}-cell Bredt on the loops)",
    )
    # Reserve the colourbar slot the other rows use, so all three sections share
    # one width and line up; the legend sits where the bar would be.
    spacer = fig.colorbar(plt.cm.ScalarMappable(), ax=ax_properties, fraction=0.02, pad=0.01)
    spacer.ax.set_visible(False)
    ax_properties.legend(loc="center left", bbox_to_anchor=(0.97, 0.5), frameon=False, fontsize=11)

    fig.text(
        0.5, 0.005,
        f"{_LAYOUT_KEY}, NACA 2412, 1 m chord · catalogue composites · "
        f"M_flap = {args.moment_flapwise / 1e3:g} kN·m, V_flap = {args.shear_flapwise / 1e3:g} kN",
        ha="center", va="bottom", fontsize=11, color=MUTED,
    )

    png_path = out_dir / "downstream_consumers.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    payload = {
        "layout": _LAYOUT_KEY,
        "load_case": {"M_flapwise_Nm": float(args.moment_flapwise), "V_flapwise_N": float(args.shear_flapwise)},
        "gauges_m": chain["gauges_m"],
        "geometry_seconds": float(chain["geometry_seconds"]),
        "cell_wall_loops": len(geometry.cell_wall_loops),
        "fe": {
            "elements": int(mesh.n_total_elements),
            "nodes": int(mesh.n_total_nodes),
            "governing": {name: (value if isinstance(value, (int, float, str)) else str(value))
                          for name, value in governing.items()},
        },
        "loads": {
            "channel": str(section_forces.meta.channel),
            "resultant_relative_error": reconstruction_error,
            "peak_shear_flow_N_per_m": float(np.max(shear_flow) * 1e3),
        },
        "properties": {
            "EA_N": float(properties.EA),
            "EI_flap_Nm2": float(properties.EI_flap),
            "EI_edge_Nm2": float(properties.EI_edge),
            "GJ_Nm2": float(torsion.GJ),
            "GJ_source": "multi-cell Bredt on cell-wall loops, wall A66 from bound laminates",
            "n_cells": int(torsion.n_cells),
            "elastic_centroid_B": [float(value) for value in centres.elastic_centroid],
            "shear_centre_B": [float(centres.xs), float(centres.ys)],
        },
    }
    json_path = out_dir / "downstream_consumers.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
