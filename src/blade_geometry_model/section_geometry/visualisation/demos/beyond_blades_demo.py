"""
visualisation/demos/beyond_blades_demo.py
=========================================
The kernel knows thin-walled multi-cell sections, not blades.

The same ``SectionGeometryRunner`` call, on outlines that are not aerofoils —
built with ``AirfoilSDF.from_array`` — then the same cell-wall loops into the
same multi-cell Bredt and shear-centre solves (``cell_loops_thin_wall_demo._solve``,
which checks torque equilibrium on every solve).

One constraint is real and stated rather than hidden: each surface must be
single-valued along the chord, because the offset branches are traced as
graphs ``y(x)``. A rounded rectangle with straight vertical sides is therefore
rejected with ``InfeasibleSectionGeometry`` for capped layouts; the demo tries
one and records the refusal in its JSON.

Writes under ``visualisation/output/beyond_blades/``:

* ``beyond_blades.png``
* ``beyond_blades.json``

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/beyond_blades_demo.py --presentation
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_contract import InfeasibleSectionGeometry
from blade_structure_utilities.geometry_utilities.section_geometry import AirfoilSDF
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import INK

_loops_spec = importlib.util.spec_from_file_location(
    "cell_loops_thin_wall_demo", Path(__file__).resolve().parent / "cell_loops_thin_wall_demo.py"
)
_loops_mod = importlib.util.module_from_spec(_loops_spec)
assert _loops_spec.loader is not None
_loops_spec.loader.exec_module(_loops_mod)


def _ellipse(chord: float, depth: float, n: int = 240) -> np.ndarray:
    angle = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.column_stack([0.5 * chord * (1.0 + np.cos(angle)), 0.5 * depth * np.sin(angle)])


def _superellipse(chord: float, depth: float, power: float = 4.0, n: int = 320) -> np.ndarray:
    """|x|^p + |y|^p = 1, scaled: flat-sided like a box girder, but smooth everywhere."""
    angle = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    cosine, sine = np.cos(angle), np.sin(angle)
    x = np.sign(cosine) * np.abs(cosine) ** (2.0 / power)
    y = np.sign(sine) * np.abs(sine) ** (2.0 / power)
    return np.column_stack([0.5 * chord * (1.0 + x), 0.5 * depth * y])


def _rounded_rectangle(chord: float, depth: float, corner: float, n: int = 20) -> np.ndarray:
    points = []
    half = 0.5 * depth
    for centre_x, centre_y, start in (
        (chord - corner, half - corner, 0.0),
        (corner, half - corner, 0.5 * np.pi),
        (corner, -half + corner, np.pi),
        (chord - corner, -half + corner, 1.5 * np.pi),
    ):
        for angle in np.linspace(start, start + 0.5 * np.pi, n):
            points.append((centre_x + corner * np.cos(angle), centre_y + corner * np.sin(angle)))
    return np.asarray(points)


#: (title, outline, layout key)
_CASES = (
    ("Tidal-turbine blade · NACA 2412", AirfoilSDF.from_naca("2412", chord=1.0), "2D-F"),
    ("Elliptical strut", AirfoilSDF.from_array(_ellipse(1.0, 0.22)), "3C-CN"),
    ("Box girder (superellipse, p = 4)", AirfoilSDF.from_array(_superellipse(1.0, 0.30)), "2D-F"),
    ("Symmetric wing · NACA 0018", AirfoilSDF.from_naca("0018", chord=1.0), "3B-F"),
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="beyond_blades")
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    solved = [(title, _loops_mod._solve(outline, key)) for title, outline, key in _CASES]

    refusal: str | None = None
    try:
        _loops_mod._solve(AirfoilSDF.from_array(_rounded_rectangle(1.0, 0.30, 0.05)), "2D-F")
    except InfeasibleSectionGeometry as exc:
        refusal = str(exc).split(" (first")[0]

    peak = max(abs(flow) for _, row in solved for _, _, flow in row["segment_flows"])
    norm = Normalize(vmin=0.0, vmax=peak)
    colour_map = plt.get_cmap("magma_r")

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 5.6), gridspec_kw={"hspace": 0.55, "wspace": 0.08})
    for axis, (title, row) in zip(axes.ravel(), solved):
        _loops_mod._draw_section(axis, row, norm, colour_map)
        axis.set_title(
            f"{title}\n{row['key']}, {row['n_cells']} cells  ·  J = {row['J_m4'] * 1e6:.0f}×10⁻⁶ m⁴",
            loc="left", fontsize=12, color=INK,
        )

    png_path = out_dir / "beyond_blades.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    payload = {
        "cases": [
            {"title": title, "layout": row["key"], "n_cells": row["n_cells"], "J_m4": row["J_m4"],
             "shear_centre": row["shear_centre"]}
            for title, row in solved
        ],
        "rejected_rounded_rectangle": refusal,
    }
    (out_dir / "beyond_blades.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
