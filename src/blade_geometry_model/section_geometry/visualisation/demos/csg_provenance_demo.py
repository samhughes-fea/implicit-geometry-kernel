"""
visualisation/demos/csg_provenance_demo.py
===============================
Demonstrate implicit-CSG component interactions via distance-field diagnostics.

Writes under ``visualisation/output/csg_provenance/``:

* ``csg_provenance_map.png`` - winner-label map, nearest-interface map and
  junction map (distance to the third-nearest region), stacked, sharing one
  region legend.
* ``csg_provenance_map.json`` - run metadata and winner-pixel counts.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from scipy.ndimage import center_of_mass, distance_transform_edt
from scipy.ndimage import label as label_connected

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)
from blade_geometry_model.section_geometry.visualisation.palette import (
    NAVY,
    RED,
    strip_axes_chrome,
    tighten_to_section,
)

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SDFGrid,
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.visualisation.plot import (
    _component_color,
)
from blade_geometry_model.section_geometry.visualisation.plot_section import (
    apply_structural_section_axes,
    rotated_naca_grid_to_section_plot_axes,
)

_DPI = DEFAULT_DPI


def _to_key(label) -> str:
    return str(getattr(label, "value", label))


def _save(fig, path: Path, *, presentation: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--key", default="2D-F")
    p.add_argument("--naca", default="2412")
    p.add_argument("--chord", type=float, default=1.0)
    p.add_argument("--twist-deg", type=float, default=15.0)
    p.add_argument("--nx", type=int, default=1600)
    p.add_argument("--ny", type=int, default=700)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument(
        "--junction-clip-mm",
        type=float,
        default=1.5,
        help="A pixel within this distance of a third resolved region belongs to a junction.",
    )
    p.add_argument("--presentation", action="store_true", help="Deck styling.")
    args = p.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = args.out_dir or (Path(__file__).resolve().parents[1] / "output" / "csg_provenance")
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = resolve_system_type(str(args.key))
    twist_rad = float(np.deg2rad(float(args.twist_deg)))
    airfoil = AirfoilSDF.from_naca(str(args.naca), chord=float(args.chord), n_points=180)
    section = build_section_view(
        airfoil,
        layout,
        twist_rad=twist_rad,
        skin_thickness=0.0035,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.004,
    )

    airfoil_for_grid = airfoil.rotate(twist_rad) if abs(twist_rad) > 1e-12 else airfoil
    grid = SDFGrid.from_airfoil(airfoil_for_grid, nx=int(args.nx), ny=int(args.ny), padding=0.08)

    labels = list(section)
    phi_stack: list[np.ndarray] = []
    label_keys: list[str] = []
    for label in labels:
        key = _to_key(label)
        try:
            phi_i = np.asarray(grid.eval(section[label]), dtype=float)
        except Exception:
            continue
        if phi_i.shape == grid.X.shape:
            phi_stack.append(phi_i)
            label_keys.append(key)
    if not phi_stack:
        raise RuntimeError("No evaluable label SDFs found for provenance map")

    phi = np.stack(phi_stack, axis=0)
    n_labels = int(phi.shape[0])
    inside = phi < 0.0
    any_inside = np.any(inside, axis=0)
    winner_inside = np.argmin(phi, axis=0)
    nearest = np.argmin(np.abs(phi), axis=0)
    winner = np.where(any_inside, winner_inside, nearest)

    # Where does a mesher need a shared node? Where three materials meet. The
    # raw fields cannot answer that: implicit regions overlap (cores run under
    # the walls, adjacent cores overlap each other), so "k-th nearest |phi|"
    # measures overlaps, not junctions. Answer it on the resolved partition
    # instead: the winner map, then each label's distance transform. The third
    # smallest of those distances is ~0 exactly where a third label arrives.
    if n_labels < 3:
        raise RuntimeError("A junction map needs at least three regions.")
    cell_spacing = (
        float(abs(grid.Y[1, 0] - grid.Y[0, 0])),
        float(abs(grid.X[0, 1] - grid.X[0, 0])),
    )
    resolved = np.where(any_inside, winner_inside, -1)
    label_distances = np.stack(
        [distance_transform_edt(resolved != index, sampling=cell_spacing) for index in range(n_labels)],
        axis=0,
    )
    third_nearest = np.partition(label_distances, kth=2, axis=0)[2]

    x_plot, y_plot = rotated_naca_grid_to_section_plot_axes(
        grid.X, grid.Y, twist_rad=twist_rad, chord_m=float(args.chord)
    )

    outline = np.asarray(airfoil_for_grid.vertices, dtype=float)
    outline_x, outline_y = rotated_naca_grid_to_section_plot_axes(
        np.r_[outline[:, 0], outline[0, 0]],
        np.r_[outline[:, 1], outline[0, 1]],
        twist_rad=twist_rad,
        chord_m=float(args.chord),
    )

    colors = [_component_color(k) for k in label_keys]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, n_labels + 0.5, 1.0), cmap.N)

    fig, axes = plt.subplots(3, 1, figsize=(11.0, 8.4), gridspec_kw={"hspace": 0.32})
    ax_winner, ax_nearest, ax_margin = axes

    # Winner only means something inside a region: outside, leave the page blank
    # rather than painting the bounding box with whichever label is closest.
    ax_winner.pcolormesh(
        x_plot, y_plot, np.ma.masked_where(~any_inside, winner_inside),
        cmap=cmap, norm=norm, shading="auto",
    )
    ax_winner.set_title("Winner: most negative φ — which laminate is at this point", loc="left")

    ax_nearest.pcolormesh(
        x_plot, y_plot, np.ma.masked_where(~any_inside, nearest),
        cmap=cmap, norm=norm, shading="auto",
    )
    ax_nearest.set_title("Nearest interface — what the mesher asks for", loc="left")

    # Each junction is a handful of pixels; cluster them and ring the centroid
    # so they read at slide scale. The faded partition underneath says which
    # materials are meeting.
    junction_clip = float(args.junction_clip_mm) * 1e-3
    junction_mask = any_inside & (third_nearest < junction_clip)
    junction_labels, junction_count = label_connected(junction_mask, structure=np.ones((3, 3)))
    junction_centres = center_of_mass(junction_mask, junction_labels, range(1, junction_count + 1))

    ax_margin.pcolormesh(
        x_plot, y_plot, np.ma.masked_where(~any_inside, winner_inside),
        cmap=cmap, norm=norm, shading="auto", alpha=0.3,
    )
    if junction_count:
        rows = np.array([centre[0] for centre in junction_centres])
        columns = np.array([centre[1] for centre in junction_centres])
        row_index = np.clip(np.round(rows).astype(int), 0, x_plot.shape[0] - 1)
        column_index = np.clip(np.round(columns).astype(int), 0, x_plot.shape[1] - 1)
        ax_margin.scatter(
            x_plot[row_index, column_index], y_plot[row_index, column_index],
            s=150, facecolors="none", edgecolors=RED, linewidths=2.2, zorder=4,
        )
    ax_margin.set_title(
        f"Junctions: {junction_count} points where three materials meet — each needs a shared node",
        loc="left",
    )

    for axis in axes:
        axis.plot(outline_x, outline_y, color=NAVY, lw=1.0)
        apply_structural_section_axes(axis)
        tighten_to_section(axis, outline_x, outline_y, pad_fraction=0.03)
        strip_axes_chrome(axis)

    # One legend for all three panels, by region, in drawing colours.
    legend_handles = [
        Patch(facecolor=colors[i], edgecolor="none", label=label_keys[i].replace("_", " "))
        for i in range(n_labels)
    ]
    fig.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.06),
        ncol=5, frameon=False, fontsize=10.5,
    )

    out_png = out_dir / "csg_provenance_map.png"
    _save(fig, out_png, presentation=bool(args.presentation))

    payload = {
        "key": str(args.key),
        "naca": str(args.naca),
        "chord_m": float(args.chord),
        "twist_deg": float(args.twist_deg),
        "grid": {"nx": int(args.nx), "ny": int(args.ny)},
        "labels": label_keys,
        "winner_pixel_counts": {k: int(np.sum(winner == i)) for i, k in enumerate(label_keys)},
        "nearest_pixel_counts": {k: int(np.sum(nearest == i)) for i, k in enumerate(label_keys)},
        "inside_fraction": float(np.mean(any_inside)),
        "third_nearest_region_inside_p05_m": float(np.quantile(third_nearest[any_inside], 0.05)),
        "junction_count": int(junction_count),
        "png": out_png.name,
    }
    out_json = out_dir / "csg_provenance_map.json"
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

