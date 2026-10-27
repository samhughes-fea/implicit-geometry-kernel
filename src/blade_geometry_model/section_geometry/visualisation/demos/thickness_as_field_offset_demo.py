"""
visualisation/demos/thickness_as_field_offset_demo.py
==========================================
Wall thickness is a field operation, not a mesh operation.

Resizing a wall in an explicit-geometry pipeline means offsetting a boundary
curve and repairing whatever self-intersections that creates. Here it is one
scalar in one ``Offset`` node: the expression tree is byte-for-byte the same
shape, only ``amount`` changes, and the area responds smoothly with no boolean
solver and no repair pass.

The honest half is the fourth panel. After CSG, the assembled field is
signed-distance-*like*, not metric: ``Offset`` is a true geometric offset only
where the Eikonal condition ``|grad phi| = 1`` holds. Along CSG seams (where
``min``/``max`` switch branches) the gradient collapses and the offset under-cuts.
That is why offsets are taken against the aerofoil polygon - which is metric -
rather than against an already-assembled union.

Writes under ``visualisation/output/thickness_as_field_offset/``:

* ``thickness_sweep.png`` - nested interfaces, field extremes, area response.
* ``offset_metric_validity.png`` - |grad phi| and where offsetting is trustworthy.
* ``thickness_as_field_offset.json`` - thickness/area table, Eikonal percentiles.

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/thickness_as_field_offset_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize

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
    INK,
    NAVY,
    RED,
    region_color,
    strip_axes_chrome,
)

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SDFGrid,
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.sdf.ir import (
    Offset,
    eval_expr_batch,
    union_expr,
)
from blade_geometry_model.section_geometry.sections.clip_algebra import (
    airfoil_expr,
)

_DPI = DEFAULT_DPI


def _to_key(label) -> str:
    return str(getattr(label, "value", label))


def _save(fig, path: Path, *, presentation: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def union_field_for(section, grid: SDFGrid, prefixes: tuple[str, ...] | None = None) -> np.ndarray:
    """Union of the section's regions, optionally restricted by label prefix.

    ``prefixes=None`` unions everything. Note that the full union's *outer*
    boundary is always the aerofoil OML regardless of gauge - it is the inner
    boundaries that move - so the interesting sweeps restrict to the skin or to
    the cores.
    """
    expressions = []
    for label in section:
        key = str(getattr(label, "value", label))
        if prefixes is not None and not key.startswith(prefixes):
            continue
        expression = getattr(section[label], "expr", None)
        if expression is not None:
            expressions.append(expression)
    if not expressions:
        return np.full(grid.X.shape, np.inf, dtype=float)
    shared_memo: dict = {}
    return np.asarray(
        eval_expr_batch(union_expr(expressions), grid.X, grid.Y, shared_memo),
        dtype=float,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--key", default="2D-F")
    parser.add_argument("--naca", default="2412")
    parser.add_argument("--chord", type=float, default=1.0)
    parser.add_argument("--twist-deg", type=float, default=0.0)
    parser.add_argument("--nx", type=int, default=900)
    parser.add_argument("--ny", type=int, default=400)
    parser.add_argument(
        "--skin-thicknesses-mm",
        type=float,
        nargs="+",
        default=[2.0, 5.0, 9.0, 14.0, 20.0],
        help="Skin gauges to sweep [mm].",
    )
    parser.add_argument("--spar-cap-thickness", type=float, default=0.012)
    parser.add_argument("--shear-web-thickness", type=float, default=0.006)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--presentation", action="store_true", help="Deck styling.")
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)

    out_dir = args.out_dir or (
        Path(__file__).resolve().parents[1] / "output" / "thickness_as_field_offset"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = resolve_system_type(str(args.key))
    twist_rad = float(np.deg2rad(float(args.twist_deg)))
    airfoil = AirfoilSDF.from_naca(str(args.naca), chord=float(args.chord), n_points=200)
    airfoil_for_grid = airfoil.rotate(twist_rad) if abs(twist_rad) > 1e-12 else airfoil
    grid = SDFGrid.from_airfoil(
        airfoil_for_grid, nx=int(args.nx), ny=int(args.ny), padding=0.06
    )

    thicknesses_m = [float(value) / 1000.0 for value in args.skin_thicknesses_mm]

    # --- the sweep: same tree, one scalar changed --------------------------
    sweep_rows: list[dict] = []
    sweep_fields: list[np.ndarray] = []
    sweep_skin_fields: list[np.ndarray] = []
    for skin_thickness in thicknesses_m:
        section = build_section_view(
            airfoil,
            layout,
            twist_rad=twist_rad,
            skin_thickness=skin_thickness,
            spar_cap_thickness=float(args.spar_cap_thickness),
            shear_web_thickness=float(args.shear_web_thickness),
        )
        phi_union = union_field_for(section, grid)
        phi_skin = union_field_for(section, grid, prefixes=("skin_",))
        phi_core = union_field_for(section, grid, prefixes=("sandwich_core_",))
        sweep_fields.append(phi_union)
        sweep_skin_fields.append(phi_skin)
        sweep_rows.append(
            {
                "skin_thickness_m": skin_thickness,
                "skin_thickness_mm": skin_thickness * 1000.0,
                "skin_area_m2": float(grid.area(phi_skin)),
                "core_area_m2": float(grid.area(phi_core)),
                "all_regions_area_m2": float(grid.area(phi_union)),
                "labels": [_to_key(label) for label in section],
            }
        )

    # --- the pure-offset family on the metric aerofoil field ---------------
    aerofoil_expression = airfoil_expr(airfoil_for_grid)
    phi_aerofoil = np.asarray(
        eval_expr_batch(aerofoil_expression, grid.X, grid.Y, {}), dtype=float
    )

    # --- figure 1 -----------------------------------------------------------
    figure = plt.figure(figsize=(15.0, 6.6))
    outer = figure.add_gridspec(2, 1, height_ratios=[1.0, 0.55], hspace=0.12)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.6, 1.0], wspace=0.14)
    axis_nested = figure.add_subplot(top[0, 0])
    axis_area = figure.add_subplot(top[0, 1])

    # Deck ramp: thin gauges pale, thick gauges in the skin's own navy.
    color_map = matplotlib.colors.LinearSegmentedColormap.from_list("gauge", ["#9ec5e8", NAVY])
    normalise = Normalize(vmin=min(thicknesses_m), vmax=max(thicknesses_m))

    # The exterior skin-shell surface never moves; the interior surface is what the scalar drives.
    axis_nested.contour(
        grid.X, grid.Y, phi_aerofoil, levels=[0.0], colors="#111418", linewidths=2.4
    )
    for thickness, phi_skin in zip(thicknesses_m, sweep_skin_fields):
        axis_nested.contour(
            grid.X,
            grid.Y,
            phi_skin,
            levels=[0.0],
            colors=[color_map(normalise(thickness))],
            linewidths=1.4,
        )
    axis_nested.set_aspect("equal")
    axis_nested.set_title(
        "Exterior skin-shell surface fixed; interior per gauge, 2 → 20 mm", loc="left", fontsize=13, color=INK,
    )
    strip_axes_chrome(axis_nested)
    mappable = cm.ScalarMappable(norm=Normalize(1e3 * min(thicknesses_m), 1e3 * max(thicknesses_m)), cmap=color_map)
    color_bar = figure.colorbar(
        mappable, ax=axis_nested, fraction=0.045, pad=0.02, orientation="horizontal", aspect=45
    )
    color_bar.set_label("skin gauge [mm]")
    color_bar.outline.set_visible(False)

    thickness_mm = [row["skin_thickness_mm"] for row in sweep_rows]
    axis_area.plot(
        thickness_mm,
        [row["skin_area_m2"] * 1e4 for row in sweep_rows],
        marker="o",
        color=region_color("skin_suction_surface"),
        lw=2.4,
        label="skin wall",
    )
    axis_area.plot(
        thickness_mm,
        [row["core_area_m2"] * 1e4 for row in sweep_rows],
        marker="s",
        color=region_color("sandwich_core"),
        lw=2.4,
        label="enclosed cores (cavity)",
    )
    axis_area.set_xlabel("skin thickness [mm]")
    axis_area.set_ylabel("area [cm²]")
    axis_area.set_title("Skin grows, cavity gives way", loc="left", fontsize=13, color=INK)
    axis_area.grid(alpha=0.4)
    axis_area.legend(fontsize=11, frameon=False)

    bottom = outer[1].subgridspec(1, len(thicknesses_m), wspace=0.06)
    for panel_index, (thickness, phi_union) in enumerate(zip(thicknesses_m, sweep_fields)):
        axis_panel = figure.add_subplot(bottom[0, panel_index])
        limit = 0.035
        axis_panel.pcolormesh(
            grid.X,
            grid.Y,
            phi_union,
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            shading="auto",
        )
        axis_panel.contour(
            grid.X, grid.Y, phi_union, levels=[0.0], colors="k", linewidths=1.2
        )
        axis_panel.set_aspect("equal")
        strip_axes_chrome(axis_panel)
        axis_panel.set_title(f"{thickness * 1000.0:.0f} mm", fontsize=12, color=INK)

    _save(figure, out_dir / "thickness_sweep.png", presentation=presentation)

    # --- figure 2: where offsetting is actually valid ----------------------
    reference_thickness = thicknesses_m[len(thicknesses_m) // 2]
    reference_section = build_section_view(
        airfoil,
        layout,
        twist_rad=twist_rad,
        skin_thickness=reference_thickness,
        spar_cap_thickness=float(args.spar_cap_thickness),
        shear_web_thickness=float(args.shear_web_thickness),
    )
    phi_assembled = union_field_for(reference_section, grid)

    gradient_aerofoil = grid.grad_magnitude(phi_aerofoil)
    gradient_assembled = grid.grad_magnitude(phi_assembled)
    residual_aerofoil = np.abs(gradient_aerofoil - 1.0)
    residual_assembled = np.abs(gradient_assembled - 1.0)

    metric_figure, metric_axes = plt.subplots(
        2, 2, figsize=(15.0, 6.4), gridspec_kw={"hspace": 0.32, "wspace": 0.12}
    )
    (axis_grad_aero, axis_grad_asm), (axis_offset_family, axis_residual) = metric_axes

    # Map the residual, not |grad phi| itself: the field is 1.0 almost
    # everywhere, so plotting the magnitude spends the whole colour range on
    # the uninteresting part and leaves the seams as hairlines.
    for axis, residual_field, title in (
        (
            axis_grad_aero,
            residual_aerofoil,
            "Aerofoil polygon - metric almost everywhere",
        ),
        (
            axis_grad_asm,
            residual_assembled,
            "After CSG assembly - the min/max seam skeleton lights up",
        ),
    ):
        mesh = axis.pcolormesh(
            grid.X,
            grid.Y,
            residual_field,
            cmap="inferno",
            vmin=0.0,
            vmax=0.6,
            shading="auto",
        )
        axis.contour(
            grid.X, grid.Y, phi_assembled, levels=[0.0], colors="#7fd4ff", linewidths=0.8
        )
        axis.set_aspect("equal")
        strip_axes_chrome(axis)
        axis.set_title(title, loc="left", fontsize=12.5, color=INK)
        color_bar_residual = metric_figure.colorbar(
            mesh, ax=axis, fraction=0.028, pad=0.015
        )
        color_bar_residual.set_label("| ‖∇φ‖ − 1 |", fontsize=11)
        color_bar_residual.outline.set_visible(False)

    # Offsetting the metric aerofoil field gives a genuine parallel family.
    for offset_amount in (-0.005, -0.015, -0.030, -0.050):
        phi_offset = np.asarray(
            eval_expr_batch(
                Offset(aerofoil_expression, amount=offset_amount), grid.X, grid.Y, {}
            ),
            dtype=float,
        )
        axis_offset_family.contour(
            grid.X,
            grid.Y,
            phi_offset,
            levels=[0.0],
            colors=[NAVY],
            linewidths=1.6,
        )
    axis_offset_family.contour(
        grid.X, grid.Y, phi_aerofoil, levels=[0.0], colors="#111418", linewidths=2.0
    )
    axis_offset_family.set_aspect("equal")
    strip_axes_chrome(axis_offset_family)
    axis_offset_family.set_title(
        "Offsets of the aerofoil: a true parallel family (5, 15, 30, 50 mm)", loc="left", fontsize=12.5,
        color=INK,
    )

    percentiles = (50, 90, 99)
    residual_summary = {
        "aerofoil_polygon": {
            f"p{p}": float(np.percentile(residual_aerofoil, p)) for p in percentiles
        },
        "assembled_section": {
            f"p{p}": float(np.percentile(residual_assembled, p)) for p in percentiles
        },
    }
    axis_residual.hist(
        residual_aerofoil.ravel(),
        bins=80,
        range=(0.0, 1.0),
        histtype="step",
        linewidth=2.2,
        color=NAVY,
        label="aerofoil polygon",
        density=True,
    )
    axis_residual.hist(
        residual_assembled.ravel(),
        bins=80,
        range=(0.0, 1.0),
        histtype="step",
        linewidth=2.2,
        color=RED,
        label="assembled section",
        density=True,
    )
    axis_residual.set_yscale("log")
    axis_residual.set_xlabel("Eikonal residual  | ‖∇φ‖ − 1 |")
    axis_residual.set_ylabel("density")
    # Both medians are ~1e-15, so the medians say nothing; the tail is the story.
    axis_residual.set_title(
        f"99th percentile: aerofoil {residual_summary['aerofoil_polygon']['p99']:.4f} "
        f"vs assembled {residual_summary['assembled_section']['p99']:.3f}",
        loc="left", fontsize=12.5, color=INK,
    )
    axis_residual.grid(alpha=0.4)
    axis_residual.legend(fontsize=11, frameon=False)

    _save(metric_figure, out_dir / "offset_metric_validity.png", presentation=presentation)

    payload = {
        "key": str(args.key),
        "naca": str(args.naca),
        "chord_m": float(args.chord),
        "twist_deg": float(args.twist_deg),
        "grid": {"nx": int(args.nx), "ny": int(args.ny)},
        "spar_cap_thickness_m": float(args.spar_cap_thickness),
        "shear_web_thickness_m": float(args.shear_web_thickness),
        "sweep": sweep_rows,
        "eikonal_residual": residual_summary,
        "reference_skin_thickness_m": float(reference_thickness),
    }
    json_path = out_dir / "thickness_as_field_offset.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
