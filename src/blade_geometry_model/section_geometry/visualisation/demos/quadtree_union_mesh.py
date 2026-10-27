"""
visualisation/demos/quadtree_union_mesh.py
===============================
Hybrid demo: build section geometry from the orchestration ``SYSTEM_TYPE_KEYS``
registry, then build ONE adaptive quadtree refined to the UNION of all
subcomponent boundaries (skin, webs, caps, cores) and render:

  (left)  the built geometry (filled subcomponents)
  (right) the adaptive quadtree mesh + every subcomponent phi=0 contour
          (the union of all internal interfaces)

By default renders **all 40** registered layout keys (``0A``, ``0B``, and every
multicell ``{X}{Y}-{Z}`` token) at **15 deg** section twist (``--twist-deg``),
including skin-only layouts, chord-normal webs, and flapwise webs.

Outputs under ``visualisation/output/quadtree_union/``:

  * ``quadtree_union_{key}.png`` per layout (e.g. ``quadtree_union_2D-CN.png``)
  * ``quadtree_union_summary.json`` — coverage metadata and leaf counts

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/quadtree_union_mesh.py
    python src/blade_geometry_model/section_geometry/visualisation/demos/quadtree_union_mesh.py 2D-CN 2C-F
    python src/blade_geometry_model/section_geometry/visualisation/demos/quadtree_union_mesh.py --list-keys
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

_EXAMPLES_DIR = Path(__file__).resolve().parent


def _load_sibling_module(name: str, filename: str):
    path = _EXAMPLES_DIR / filename
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


_ensure_mod = _load_sibling_module("_ensure_repo_path", "_ensure_repo_path.py")
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)
_taxonomy = _load_sibling_module("_taxonomy_cases", "_taxonomy_cases.py")
airfoil_for_twist = _taxonomy.airfoil_for_twist
system_type_display_name = _taxonomy.system_type_display_name
twist_deg_for_layout = _taxonomy.twist_deg_for_layout
twist_rad_for_layout = _taxonomy.twist_rad_for_layout

from blade_geometry_model.section_geometry.api.factory import (
    SYSTEM_TYPE_KEYS,
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.props_emit.area_mesh import build_shared_section_quadtree
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SDFGrid,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.section_geometry.visualisation.plot import (
    _component_color,
    plot_section,
)


def render_case(
    airfoil_base: Any,
    key: str,
    *,
    chord: float,
    max_depth: int,
    out_dir: Path,
    twist_deg: float = 15.0,
    write_png: bool = True,
    presentation: bool = False,
) -> dict[str, Any]:
    """Build and render one taxonomy case; return summary metadata."""
    layout = resolve_system_type(key)
    twist_deg = twist_deg_for_layout(layout, twist_deg=twist_deg)
    twist_rad = twist_rad_for_layout(layout, twist_deg=twist_deg)
    section = build_section_view(airfoil_base, layout, twist_angle_rad=twist_rad)
    airfoil_b = airfoil_for_twist(airfoil_base, twist_rad)
    grid = SDFGrid.from_airfoil(airfoil_b, nx=600, ny=240)
    labels = list(section)
    spec = quadtree_spec_for_airfoil(
        airfoil_b, chord=float(chord), max_depth=int(max_depth), max_cells=120000
    )
    mesh = build_shared_section_quadtree(section, spec, twist_rad=twist_rad)
    leaves = mesh.leaf_cells()
    display_name = system_type_display_name(key)
    geom_title = f"{display_name} (twist={twist_deg:g} deg)"

    if presentation:
        # Slide variant: one panel, cropped to the section — the refinement
        # hugging the walls is the whole message, and the square root cell
        # around it is dead space.
        fig, axq = plt.subplots(figsize=(13.0, 4.6))
    else:
        fig, (axg, axq) = plt.subplots(1, 2, figsize=(20, 6))
        plot_section(
            section,
            grid,
            ax=axg,
            alpha=0.5,
            show_airfoil=False,
            title=f"Built geometry — {geom_title}",
            show_legend=True,
        )
    for cell in leaves:
        axq.add_patch(
            Rectangle(
                (cell.x_min, cell.y_min),
                cell.width,
                cell.height,
                fill=False,
                edgecolor="0.6",
                linewidth=0.25 if presentation else 0.3,
            )
        )
    for label in labels:
        phi = grid.eval(section[label])
        axq.contour(
            grid.X, grid.Y, phi, levels=[0.0], colors=[_component_color(label)],
            linewidths=2.0 if presentation else 1.4,
        )
    if presentation:
        outline = np.asarray(airfoil_b.vertices, dtype=float)
        pad = 0.03 * float(chord)
        axq.set_xlim(outline[:, 0].min() - pad, outline[:, 0].max() + pad)
        axq.set_ylim(outline[:, 1].min() - pad, outline[:, 1].max() + pad)
        axq.set_aspect("equal")
        axq.set_title(
            f"{key}, {twist_deg:g}° twist · one quadtree shared by every region · "
            f"{len(leaves):,} leaves, depth ≤ {max_depth}",
            loc="left", fontsize=13,
        )
    else:
        axq.set_xlim(spec.x_min, spec.x_max)
        axq.set_ylim(spec.y_min, spec.y_max)
        axq.set_title(
            f"Square union-refined quadtree — {display_name} "
            f"(depth<={max_depth}, leaves={len(leaves)})"
        )
    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
    )

    apply_structural_section_axes(axq)
    if presentation:
        from blade_geometry_model.section_geometry.visualisation.palette import (
            strip_axes_chrome,
        )

        strip_axes_chrome(axq)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    png_name = f"quadtree_union_{key}.png"
    out_png = out_dir / png_name
    if write_png:
        save_figure(fig, out_png, presentation=presentation)
        print(f"Wrote {out_png}  (leaves={len(leaves)})")
    plt.close(fig)

    z_token = "CN" if layout.web_orientation == "chord_normal" else "F"
    if layout.n_webs == 0:
        z_token = None

    return {
        "system_type_key": key,
        "structural_family": layout.structural_family,
        "web_orientation": layout.web_orientation,
        "z_token": z_token,
        "n_webs": layout.n_webs,
        "n_cells": layout.n_cells,
        "twist_deg": float(twist_deg),
        "labels": labels,
        "leaf_count": int(len(leaves)),
        "png": png_name,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "keys",
        nargs="*",
        metavar="KEY",
        help="Layout keys to render (default: all SYSTEM_TYPE_KEYS).",
    )
    p.add_argument("--list-keys", action="store_true", help="Print registered keys and exit.")
    p.add_argument("--naca", default="2412")
    p.add_argument("--chord", type=float, default=1.0)
    p.add_argument("--max-depth", type=int, default=8)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument(
        "--twist-deg",
        type=float,
        default=15.0,
        help="Section twist for every layout [deg] (default 15).",
    )
    p.add_argument("--presentation", action="store_true", help="Single cropped slide panel.")
    args = p.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    if args.list_keys:
        for key in SYSTEM_TYPE_KEYS:
            print(key)
        return 0

    keys = tuple(args.keys) if args.keys else SYSTEM_TYPE_KEYS
    for key in keys:
        resolve_system_type(key)

    out_dir = args.out_dir or (Path(__file__).resolve().parents[1] / "output" / "quadtree_union")
    out_dir = Path(out_dir).resolve()
    chord = float(args.chord)
    airfoil = AirfoilSDF.from_naca(str(args.naca).zfill(4), chord=chord)

    cases: list[dict[str, Any]] = []
    for key in keys:
        meta = render_case(
            airfoil,
            key,
            chord=chord,
            max_depth=args.max_depth,
            out_dir=out_dir,
            twist_deg=args.twist_deg,
            presentation=bool(args.presentation),
        )
        cases.append(meta)

    summary = {
        "naca": str(args.naca).zfill(4),
        "chord_m": chord,
        "max_depth": int(args.max_depth),
        "twist_deg": float(args.twist_deg),
        "expected_keys": list(SYSTEM_TYPE_KEYS),
        "rendered_keys": list(keys),
        "cases": cases,
    }
    summary_path = out_dir / "quadtree_union_summary.json"
    with summary_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
