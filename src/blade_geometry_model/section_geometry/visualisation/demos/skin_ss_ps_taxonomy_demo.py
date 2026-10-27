"""
visualisation/demos/skin_ss_ps_taxonomy_demo.py
====================================
Demonstrate nested outer-skin taxonomy: suction/pressure solids and
leading_edge / bay / trailing_edge arcs.

Writes:

* ``skin_ss_ps_section.png`` — hybrid section with distinct skin colours
* ``skin_ss_ps_midline_arcs.png`` — topology midlines labelled by arc
* ``skin_ss_ps_taxonomy.json`` — strip inventory

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/skin_ss_ps_taxonomy_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    demo_output_dir,
    save_figure,
)

from blade_geometry_model.coordinate_frames import (
    section_plot_axes,
)
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
)
from blade_geometry_model.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_geometry_model.section_geometry.junction_preparation.shell_midline_junction_preparation import (
    prepare_shell_midline_junctions,
)
from blade_geometry_model.section_geometry.visualisation.plot import (
    plot_section_hybrid,
    plot_shell_midline_strip_topology,
)

_ARC_COLOUR = {
    ("suction_surface", "leading_edge"): "#0D47A1",
    ("suction_surface", "bay"): "#1976D2",
    ("suction_surface", "trailing_edge"): "#64B5F6",
    ("pressure_surface", "leading_edge"): "#4A148C",
    ("pressure_surface", "bay"): "#7B1FA2",
    ("pressure_surface", "trailing_edge"): "#CE93D8",
}


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "skin_ss_ps_taxonomy"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    args = parser.parse_args(argv)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", chord=1.0)
    section = MultiCellSection.twin_web(
        af,
        x_fore=0.20,
        x_aft=0.50,
        skin_thickness=0.005,
        web_thickness=0.006,
        cap_height=(0.014, 0.012),
        structural_family="D",
        core_enabled=False,
    )

    fig, ax = plt.subplots(figsize=(12.0, 5.0))
    plot_section_hybrid(
        section,
        residual_mesh=None,
        ax=ax,
        alpha=0.55,
        title=(
            "Suction/pressure outer-skin solids "
            "(camber-upper = suction_surface, camber-lower = pressure_surface)"
        ),
        show_legend=True,
    )
    fig.tight_layout()
    path1 = out_dir / "skin_ss_ps_section.png"
    save_figure(fig, path1)
    plt.close(fig)
    print(f"Wrote {path1}")

    raw = build_shell_midline_strips(section, twist_rad=0.0)
    topo = prepare_shell_midline_junctions(raw)
    skins = [s for s in topo if s.kind == "skin"]

    fig2, ax2 = plt.subplots(figsize=(12.0, 5.2))
    non_skin = [s for s in topo if s.kind != "skin"]
    plot_shell_midline_strip_topology(
        non_skin,
        ax=ax2,
        title=(
            "Nested skin arcs: "
            "skin_{suction|pressure}_surface:{leading_edge|bay_i|trailing_edge}"
        ),
    )
    for s in skins:
        pts = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        col = _ARC_COLOUR.get((s.skin_side or "", s.skin_arc or ""), "#607D8B")
        _h, _v = section_plot_axes(pts)
        ax2.plot(_h, _v, color=col, lw=2.4, label=s.label)
    ax2.legend(fontsize=7, loc="best", ncol=2)
    fig2.tight_layout()
    path2 = out_dir / "skin_ss_ps_midline_arcs.png"
    save_figure(fig2, path2)
    plt.close(fig2)
    print(f"Wrote {path2}")

    payload = {
        "airfoil": "NACA 2412",
        "n_webs": 2,
        "skin_strips": [
            {
                "label": s.label,
                "skin_side": s.skin_side,
                "skin_arc": s.skin_arc,
                "bay_index": s.bay_index,
                "n_points": int(
                    np.asarray(s.midline_twisted_pitch_axis_blade_frame).shape[0]
                ),
                "arc_length_m": s.arc_length_m(),
            }
            for s in skins
        ],
        "notes": (
            "SDF solids are skin_suction_surface / skin_pressure_surface only; "
            "leading_edge / bay / trailing_edge are panelisation of the open "
            "midlines after web-foot topology split."
        ),
    }
    json_path = out_dir / "skin_ss_ps_taxonomy.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return out_dir


if __name__ == "__main__":
    main()
