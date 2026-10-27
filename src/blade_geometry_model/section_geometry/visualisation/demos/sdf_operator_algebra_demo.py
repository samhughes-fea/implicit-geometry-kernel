"""
visualisation/demos/sdf_operator_algebra_demo.py
=====================================
The whole shape language, one line of numpy per operator.

Implicit geometry buys its expressiveness with a very small algebra. Every
region in every one of the 40 catalogue layouts is assembled from the handful
of ``Expr`` operators below, and each one is a single vectorised numpy
expression in ``sdf/ir.py``'s ``eval_expr_batch`` dispatch. There is no mesh,
no boolean solver, and no repair pass.

Figure 1 covers the operators the blade recipes actually use. Figure 2 covers
smooth/blending operators that exist in ``sdf/csg.py`` and are exported from the
kernel but are **not** used by any current blade recipe - shown as available
capability, labelled as such.

Writes under ``visualisation/output/sdf_operator_algebra/``:

* ``sdf_operator_algebra.png`` - one panel per operator in use.
* ``sdf_kernel_only_operators.png`` - smooth union / subtract / blend sweeps.
* ``sdf_operator_algebra.json`` - operator inventory with numpy formulas.

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/sdf_operator_algebra_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    demo_output_dir,
    save_figure,
)
from blade_geometry_model.section_geometry.visualisation.palette import (
    draw_distance_field,
)

from blade_geometry_model.section_geometry import SDFGrid
from blade_geometry_model.section_geometry.sdf.csg import (
    blend,
    smooth_subtract,
    smooth_union,
)
from blade_geometry_model.section_geometry.sdf.ir import (
    Capsule,
    Circle,
    Intersect,
    Offset,
    Rotate,
    Shell,
    Slab,
    Subtract,
    Union,
    eval_expr_batch,
)

_DPI = DEFAULT_DPI

# Operand shapes reused across every panel so the reader compares operators,
# not inputs. Deliberately overlapping.
_LEFT_SHAPE = Circle(cx=-0.22, cy=0.0, r=0.42)
_RIGHT_SHAPE = Circle(cx=0.22, cy=0.0, r=0.42)
_BAR_SHAPE = Capsule(ax=-0.55, ay=-0.30, bx=0.55, by=0.30, r=0.16)

# One colour scale across every panel: the reader is comparing operators, not
# inputs, so a per-panel norm would make identical fields look different. The
# inner limit is the operand radius, so a full circle spends the whole blue
# half; the outer is clipped at twice that, which is what stops the far field
# flooding the panel with saturated colour and hiding the banding.
_FIELD_INNER_LIMIT = 0.42
_FIELD_OUTER_LIMIT = 0.84
_FIELD_BAND_INTERVAL = 0.105


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, path)
    plt.close(fig)
    print(f"Wrote {path}")


def draw_operator_panel(
    axis,
    grid: SDFGrid,
    phi: np.ndarray,
    *,
    title: str,
    formula: str,
    operand_fields: tuple[np.ndarray, ...] = (),
) -> None:
    """Field, interface, and the operand interfaces it was built from."""
    draw_distance_field(
        axis,
        grid.X,
        grid.Y,
        phi,
        norm=TwoSlopeNorm(
            vcenter=0.0, vmin=-_FIELD_INNER_LIMIT, vmax=_FIELD_OUTER_LIMIT
        ),
        band_interval=_FIELD_BAND_INTERVAL,
        zero_line_width=2.0,
    )
    for operand in operand_fields:
        axis.contour(
            grid.X,
            grid.Y,
            operand,
            levels=[0.0],
            colors="#59636e",
            linewidths=0.9,
            linestyles="dashed",
        )
    axis.contour(grid.X, grid.Y, phi, levels=[0.0], colors="k", linewidths=2.0)
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_title(title, fontsize=11.5, pad=6)
    # Inside the panel, not below it: captions placed under the axes collide
    # with the next row's titles once the grid has more than one row.
    axis.text(
        0.5,
        0.028,
        formula,
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.0,
        family="monospace",
        color="#1d262d",
        bbox={
            "facecolor": "white",
            "alpha": 0.86,
            "edgecolor": "none",
            "boxstyle": "round,pad=0.32",
        },
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=420)
    parser.add_argument("--ny", type=int, default=420)
    parser.add_argument("--extent", type=float, default=0.95)
    parser.add_argument(
        "--smooth-k",
        type=float,
        default=0.28,
        help="Blending radius for the kernel-only smooth operators.",
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args(argv)

    out_dir = args.out_dir or (
        Path(__file__).resolve().parents[1] / "output" / "sdf_operator_algebra"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    extent = float(args.extent)
    grid = SDFGrid.from_bbox(
        -extent, extent, -extent, extent, nx=int(args.nx), ny=int(args.ny)
    )

    def evaluate(expression) -> np.ndarray:
        return np.asarray(
            eval_expr_batch(expression, grid.X, grid.Y, {}), dtype=float
        )

    phi_left = evaluate(_LEFT_SHAPE)
    phi_right = evaluate(_RIGHT_SHAPE)
    phi_bar = evaluate(_BAR_SHAPE)

    # --- figure 1: the operators the blade recipes use ---------------------
    operators_in_use = [
        (
            "Union",
            "np.minimum(a, b)",
            evaluate(Union(_LEFT_SHAPE, _RIGHT_SHAPE)),
            (phi_left, phi_right),
        ),
        (
            "Intersect",
            "np.maximum(a, b)",
            evaluate(Intersect(_LEFT_SHAPE, _RIGHT_SHAPE)),
            (phi_left, phi_right),
        ),
        (
            "Subtract",
            "np.maximum(base, -cutter)",
            evaluate(Subtract(_LEFT_SHAPE, _RIGHT_SHAPE)),
            (phi_left, phi_right),
        ),
        (
            "Offset",
            "base - amount",
            evaluate(Offset(_BAR_SHAPE, amount=-0.10)),
            (phi_bar,),
        ),
        (
            "Shell",
            "np.abs(base) - thickness/2",
            evaluate(Shell(_BAR_SHAPE, thickness=0.12)),
            (phi_bar,),
        ),
        (
            "Slab",
            "np.maximum(base, -base - depth)",
            evaluate(Slab(_BAR_SHAPE, depth=0.14)),
            (phi_bar,),
        ),
        (
            "Rotate",
            "base(R(-angle) @ xy)",
            evaluate(Rotate(_BAR_SHAPE, angle=float(np.deg2rad(55.0)))),
            (phi_bar,),
        ),
    ]

    columns = 4
    rows = int(np.ceil((len(operators_in_use) + 1) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(4.05 * columns, 4.75 * rows))
    flat_axes = np.asarray(axes).ravel()

    for panel_index, (name, formula, phi, operands) in enumerate(operators_in_use):
        draw_operator_panel(
            flat_axes[panel_index],
            grid,
            phi,
            title=name,
            formula=formula,
            operand_fields=operands,
        )

    # Closing panel: the same algebra composed two levels deep.
    composed = Subtract(
        Intersect(Offset(_LEFT_SHAPE, amount=0.06), _BAR_SHAPE), _RIGHT_SHAPE
    )
    draw_operator_panel(
        flat_axes[len(operators_in_use)],
        grid,
        evaluate(composed),
        title="Composed",
        formula="Subtract(Intersect(Offset(a, .06), bar), b)",
        operand_fields=(phi_left, phi_right, phi_bar),
    )

    for spare_index in range(len(operators_in_use) + 1, flat_axes.size):
        flat_axes[spare_index].set_axis_off()

    figure.suptitle(
        "The complete shape algebra: every region of every layout is built "
        "from these operators\n"
        "solid line = resulting interface   |   dashed = operand interfaces   |   "
        "red outside, blue inside",
        fontsize=14,
        y=0.985,
    )
    # Equal-aspect panels shrink under tight_layout, so titles need explicit
    # vertical room or the row below collides with the row above.
    figure.tight_layout(rect=(0, 0, 1, 0.94), h_pad=3.0)
    _save(figure, out_dir / "sdf_operator_algebra.png")

    # --- figure 2: kernel-only smooth operators ----------------------------
    smooth_k = float(args.smooth_k)

    def as_callable(expression):
        def _field(x, y):
            return eval_expr_batch(expression, x, y, {})

        return _field

    left_callable = as_callable(_LEFT_SHAPE)
    right_callable = as_callable(_RIGHT_SHAPE)

    kernel_figure, kernel_axes = plt.subplots(1, 4, figsize=(16.6, 4.9))

    draw_operator_panel(
        kernel_axes[0],
        grid,
        np.asarray(smooth_union(left_callable, right_callable, k=smooth_k)(grid.X, grid.Y)),
        title=f"smooth_union (k={smooth_k:g})",
        formula="a*h + b*(1-h) - k*h*(1-h)",
        operand_fields=(phi_left, phi_right),
    )
    draw_operator_panel(
        kernel_axes[1],
        grid,
        np.asarray(
            smooth_subtract(left_callable, right_callable, k=smooth_k)(grid.X, grid.Y)
        ),
        title=f"smooth_subtract (k={smooth_k:g})",
        formula="smooth_intersect(base, -cutter, k)",
        operand_fields=(phi_left, phi_right),
    )
    for panel_index, blend_t in enumerate((0.35, 0.65)):
        draw_operator_panel(
            kernel_axes[2 + panel_index],
            grid,
            np.asarray(blend(left_callable, as_callable(_BAR_SHAPE), blend_t)(grid.X, grid.Y)),
            title=f"blend (t={blend_t:g})",
            formula="(1-t)*a + t*b",
            operand_fields=(phi_left, phi_bar),
        )

    kernel_figure.suptitle(
        "Kernel-only operators: available in sdf/csg.py, NOT used by any current blade recipe",
        fontsize=14,
        y=1.0,
        color="#8a4b12",
        fontweight="bold",
    )
    kernel_figure.text(
        0.5,
        0.935,
        "Sharp min/max give creases; the smooth variants trade the exact "
        "distance property for C-infinity fillets. blend() is a linear morph and "
        "its result is not a true SDF.",
        ha="center",
        fontsize=10,
        color="#41505c",
    )
    kernel_figure.tight_layout(rect=(0, 0, 1, 0.90))
    _save(kernel_figure, out_dir / "sdf_kernel_only_operators.png")

    payload = {
        "grid": {"nx": int(args.nx), "ny": int(args.ny), "extent": extent},
        "smooth_k": smooth_k,
        "operators_used_by_blade_recipes": [
            {"node": name, "numpy": formula} for name, formula, _, _ in operators_in_use
        ],
        "kernel_only_operators": [
            {
                "name": "smooth_union",
                "numpy": "a*h + b*(1-h) - k*h*(1-h)",
                "used_by_blade_recipes": False,
            },
            {
                "name": "smooth_intersect",
                "numpy": "a*h + b*(1-h) + k*h*(1-h)",
                "used_by_blade_recipes": False,
            },
            {
                "name": "smooth_subtract",
                "numpy": "smooth_intersect(base, -cutter, k)",
                "used_by_blade_recipes": False,
            },
            {
                "name": "blend",
                "numpy": "(1-t)*a + t*b",
                "used_by_blade_recipes": False,
                "note": "linear morph; result is not a true SDF",
            },
        ],
    }
    json_path = out_dir / "sdf_operator_algebra.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
