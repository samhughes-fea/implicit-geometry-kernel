"""
visualisation/demos/evaluation_cost_pruning_demo.py
========================================
Why the implicit path is cheap: shared subtrees and analytic support boxes.

Two independent savings, measured rather than asserted:

1. **Shared-subtree memoisation.** ``Expr`` nodes are frozen dataclasses, so
   they are hashable and compare structurally. Passing one ``memo`` dict through
   ``eval_expr_batch`` for the whole section means a node shared by several
   regions - the aerofoil polygon above all - is evaluated once, not once per
   region. Here we count the distinct field evaluations both ways and time them.

2. **Analytic support boxes.** ``expr_bbox`` derives a conservative bounding box
   for ``expr <= 0`` straight from the IR graph, with **zero** SDF evaluations.
   Cells outside it cannot contain any part of the region, so they can be culled
   before any field work happens.

Note on scope: the culling figure reports what ``expr_bbox`` *makes available*.
``eval_expr_batch`` itself does not prune - it evaluates the grid it is handed.
Bbox culling is applied by the adaptive area mesh in ``props_emit/area_mesh.py``.

Writes under ``visualisation/output/evaluation_cost_pruning/``:

* ``support_box_culling.png`` - per-region support boxes and culled fractions.
* ``evaluation_scaling.png`` - naive walk vs per-region memo vs one shared memo,
  and section cost against region count across the catalogue.
* ``evaluation_cost_pruning.json`` - all measured numbers.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/evaluation_cost_pruning_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    NAVY,
    RED,
)

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SDFGrid,
    build_section_view,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.ir import (
    eval_expr_batch,
    expr_bbox,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
    _component_color,
)

_DPI = DEFAULT_DPI

# Node-visit counting lives with the expression-tree demo; reuse it rather than
# re-deriving what "a node" means.
_tree_spec = importlib.util.spec_from_file_location(
    "csg_expression_tree_demo", Path(__file__).resolve().parent / "csg_expression_tree_demo.py"
)
_tree_mod = importlib.util.module_from_spec(_tree_spec)
assert _tree_spec.loader is not None
_tree_spec.loader.exec_module(_tree_mod)
count_nodes = _tree_mod.count_nodes

_OPERATOR_COLOR = "#2f6f9f"
_LEAF_COLOR = "#c2703d"
_SHARED_COLOR = "#7a4fa3"


def _to_key(label) -> str:
    return str(getattr(label, "value", label))


def _save(fig, path: Path, *, presentation: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def region_expressions(section) -> list[tuple[str, object]]:
    """``(label, expr)`` for every region that carries an IR tree."""
    pairs: list[tuple[str, object]] = []
    for label in section:
        expression = getattr(section[label], "expr", None)
        if expression is not None:
            pairs.append((_to_key(label), expression))
    return pairs


def evaluate_section(
    section, grid: SDFGrid, *, share_memo: bool
) -> tuple[float, int]:
    """Evaluate every region; return (seconds, distinct field evaluations).

    ``len(memo)`` is the exact count of nodes that reached a numpy kernel,
    because ``eval_expr_batch`` stores one entry per node it computes.
    """
    pairs = region_expressions(section)
    shared_memo: dict = {}
    started = time.perf_counter()
    distinct = 0
    for _, expression in pairs:
        memo = shared_memo if share_memo else {}
        eval_expr_batch(expression, grid.X, grid.Y, memo)
        if not share_memo:
            distinct += len(memo)
    elapsed = time.perf_counter() - started
    if share_memo:
        distinct = len(shared_memo)
    return elapsed, int(distinct)


def timed_median(function, repeats: int) -> float:
    samples = []
    for _ in range(max(1, repeats)):
        started = time.perf_counter()
        function()
        samples.append(time.perf_counter() - started)
    return float(np.median(samples))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--key", default="2D-F")
    parser.add_argument("--naca", default="2412")
    parser.add_argument("--chord", type=float, default=1.0)
    parser.add_argument("--twist-deg", type=float, default=0.0)
    parser.add_argument("--nx", type=int, default=700)
    parser.add_argument("--ny", type=int, default=300)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--scaling-keys",
        nargs="+",
        default=["0A", "1C-F", "2C-F", "3C-F", "4C-F", "5C-F"],
        help=(
            "Layout keys for the web-count scaling panel. Family C is used "
            "because it spans 1-5 webs; family D has no single-web key."
        ),
    )
    parser.add_argument(
        "--grid-widths",
        type=int,
        nargs="+",
        default=[200, 350, 500, 700, 950],
        help="Grid widths for the resolution scaling panel (height = width * 3/7).",
    )
    parser.add_argument("--skin-thickness", type=float, default=0.0035)
    parser.add_argument("--spar-cap-thickness", type=float, default=0.012)
    parser.add_argument("--shear-web-thickness", type=float, default=0.006)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--presentation", action="store_true", help="Deck styling.")
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)

    out_dir = args.out_dir or (
        Path(__file__).resolve().parents[1] / "output" / "evaluation_cost_pruning"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    twist_rad = float(np.deg2rad(float(args.twist_deg)))
    airfoil = AirfoilSDF.from_naca(str(args.naca), chord=float(args.chord), n_points=180)
    airfoil_for_grid = airfoil.rotate(twist_rad) if abs(twist_rad) > 1e-12 else airfoil

    def build(key: str):
        return build_section_view(
            airfoil,
            resolve_system_type(str(key)),
            twist_rad=twist_rad,
            skin_thickness=float(args.skin_thickness),
            spar_cap_thickness=float(args.spar_cap_thickness),
            shear_web_thickness=float(args.shear_web_thickness),
        )

    grid = SDFGrid.from_airfoil(
        airfoil_for_grid, nx=int(args.nx), ny=int(args.ny), padding=0.06
    )
    section = build(str(args.key))
    pairs = region_expressions(section)

    # --- figure 1: analytic support boxes ----------------------------------
    domain_x = float(grid.X.max() - grid.X.min())
    domain_y = float(grid.Y.max() - grid.Y.min())
    domain_area = domain_x * domain_y

    culling_rows: list[dict] = []
    figure_boxes, axes_boxes = plt.subplots(
        1, 2, figsize=(16.4, 5.4), gridspec_kw={"width_ratios": [1.5, 1.0]}
    )
    axis_boxes, axis_culled = axes_boxes

    for label, expression in pairs:
        bbox = expr_bbox(expression)
        phi = np.asarray(eval_expr_batch(expression, grid.X, grid.Y, {}), dtype=float)
        axis_boxes.contour(
            grid.X, grid.Y, phi, levels=[0.0], colors=[_component_color(label)], linewidths=1.1
        )
        if bbox is None:
            culling_rows.append(
                {"label": label, "bounded": False, "culled_fraction": 0.0}
            )
            continue
        x_min, x_max, y_min, y_max = bbox
        axis_boxes.add_patch(
            Rectangle(
                (x_min, y_min),
                x_max - x_min,
                y_max - y_min,
                facecolor="none",
                edgecolor=_component_color(label),
                linewidth=1.2,
                linestyle="dashed",
            )
        )
        box_area = max(x_max - x_min, 0.0) * max(y_max - y_min, 0.0)
        culled_fraction = 1.0 - min(box_area / domain_area, 1.0)
        culling_rows.append(
            {
                "label": label,
                "bounded": True,
                "bbox": [float(v) for v in bbox],
                "box_area_m2": float(box_area),
                "culled_fraction": float(culled_fraction),
            }
        )

    axis_boxes.set_aspect("equal")
    axis_boxes.set_title(
        f"{args.key}: dashed = each region's analytic bound, zero field evaluations",
        loc="left", fontsize=12.5,
    )
    axis_boxes.set_xticks([])
    axis_boxes.set_yticks([])
    for spine in axis_boxes.spines.values():
        spine.set_visible(False)

    bounded_rows = [row for row in culling_rows if row["bounded"]]
    bounded_rows.sort(key=lambda row: row["culled_fraction"])
    axis_culled.barh(
        np.arange(len(bounded_rows)),
        [row["culled_fraction"] * 100.0 for row in bounded_rows],
        color=[_component_color(row["label"]) for row in bounded_rows],
        height=0.66,
    )
    axis_culled.set_yticks(np.arange(len(bounded_rows)))
    axis_culled.set_yticklabels(
        [row["label"].replace("_surface", "").replace("_", " ") for row in bounded_rows], fontsize=11
    )
    for bar_index, row in enumerate(bounded_rows):
        percentage = row["culled_fraction"] * 100.0
        axis_culled.text(
            percentage - 1.5,
            bar_index,
            f"{percentage:.1f}%",
            ha="right",
            va="center",
            fontsize=10.5,
            color="white",
            fontweight="bold",
        )
    axis_culled.set_xlabel("% of the sampling domain provably outside the region")
    axis_culled.set_xlim(0, 100)
    axis_culled.set_title(
        "Work that never has to happen (applied by the adaptive area mesh)",
        loc="left", fontsize=12.5,
    )
    axis_culled.grid(axis="x", alpha=0.28)

    figure_boxes.tight_layout()
    _save(figure_boxes, out_dir / "support_box_culling.png", presentation=presentation)

    # --- figure 2: memo sharing and scaling --------------------------------
    shared_seconds = timed_median(
        lambda: evaluate_section(section, grid, share_memo=True), int(args.repeats)
    )
    isolated_seconds = timed_median(
        lambda: evaluate_section(section, grid, share_memo=False), int(args.repeats)
    )
    _, shared_evaluations = evaluate_section(section, grid, share_memo=True)
    _, isolated_evaluations = evaluate_section(section, grid, share_memo=False)

    grid_rows: list[dict] = []
    for width in args.grid_widths:
        height = max(int(round(width * 3 / 7)), 8)
        scaling_grid = SDFGrid.from_airfoil(
            airfoil_for_grid, nx=int(width), ny=height, padding=0.06
        )
        seconds = timed_median(
            lambda g=scaling_grid: evaluate_section(section, g, share_memo=True),
            max(2, int(args.repeats) // 2),
        )
        grid_rows.append(
            {
                "nx": int(width),
                "ny": int(height),
                "points": int(width * height),
                "seconds": float(seconds),
            }
        )

    key_rows: list[dict] = []
    for key in args.scaling_keys:
        try:
            key_section = build(str(key))
        except Exception as exc:  # a key may be infeasible at these gauges
            print(f"  skipping {key}: {exc}")
            continue
        seconds = timed_median(
            lambda s=key_section: evaluate_section(s, grid, share_memo=True),
            max(2, int(args.repeats) // 2),
        )
        _, evaluations = evaluate_section(key_section, grid, share_memo=True)
        key_rows.append(
            {
                "key": str(key),
                "regions": len(region_expressions(key_section)),
                "distinct_evaluations": int(evaluations),
                "seconds": float(seconds),
            }
        )

    scaling_figure, scaling_axes = plt.subplots(
        1, 2, figsize=(14.0, 5.0), gridspec_kw={"width_ratios": (1.15, 1.0), "wspace": 0.28}
    )
    axis_memo, axis_key_scaling = scaling_axes

    # Three ways to evaluate the same section. A naive walk visits every node of
    # every region's tree; a memo per region stops repeats inside one tree; one
    # memo for the section also stops repeats *between* trees, which is where
    # the shared aerofoil and cavity subtrees live.
    naive_visits = int(sum(count_nodes(expression) for _, expression in pairs))
    speedup = isolated_seconds / shared_seconds if shared_seconds > 0 else float("nan")
    steps = [
        ("walk every tree", naive_visits, None, MUTED),
        ("memo per region", isolated_evaluations, isolated_seconds, "#e89aae"),
        ("one memo for\nthe section", shared_evaluations, shared_seconds, RED),
    ]
    positions = np.arange(len(steps))
    axis_memo.bar(positions, [step[1] for step in steps], color=[step[3] for step in steps], width=0.6)
    for position, (_, count, seconds, _) in zip(positions, steps):
        axis_memo.text(position, count, f"{count}", ha="center", va="bottom", fontsize=15,
                       fontweight="bold", color=INK)
        if seconds is not None:
            axis_memo.text(position, count / 2.0, f"{seconds * 1e3:.0f} ms", ha="center",
                           va="center", fontsize=12, color="white", fontweight="bold")
    axis_memo.set_xticks(positions)
    axis_memo.set_xticklabels([step[0] for step in steps])
    axis_memo.set_ylim(0, naive_visits * 1.18)
    axis_memo.set_ylabel("field evaluations on the grid")
    axis_memo.set_title(
        f"{args.key}: {naive_visits} → {shared_evaluations} evaluations "
        f"({naive_visits / max(shared_evaluations, 1):.1f}× fewer)",
        loc="left",
    )
    axis_memo.grid(axis="y", alpha=0.4)

    regions = np.array([row["regions"] for row in key_rows], dtype=float)
    milliseconds = np.array([row["seconds"] * 1e3 for row in key_rows], dtype=float)
    axis_key_scaling.plot(regions, milliseconds, marker="o", ms=9, lw=2.4, color=NAVY)
    # The proportional line is what "structure costs per region" would look like
    # if nothing were shared: it makes the flatness of the measured curve legible.
    if regions.size:
        proportional = milliseconds[0] * regions / regions[0]
        axis_key_scaling.plot(regions, proportional, ls="--", lw=1.4, color=MUTED)
        axis_key_scaling.text(
            regions[-1], proportional[-1], "  if cost ∝ regions", va="center",
            fontsize=11, color=MUTED,
        )
    for row in key_rows:
        axis_key_scaling.annotate(
            row["key"], (row["regions"], row["seconds"] * 1e3), textcoords="offset points",
            xytext=(0, -18), ha="center", fontsize=11, color=INK,
        )
    if regions.size >= 2:
        growth = milliseconds[-1] / milliseconds[0] - 1.0
        axis_key_scaling.set_title(
            f"{int(regions[0])} → {int(regions[-1])} regions costs +{growth:.0%}, "
            f"not {regions[-1] / regions[0]:.0f}×",
            loc="left",
        )
    axis_key_scaling.set_xlabel("regions in the section")
    axis_key_scaling.set_ylabel("section evaluation [ms]")
    axis_key_scaling.set_xlim(0, regions.max() * 1.35 if regions.size else 1)
    axis_key_scaling.set_ylim(0, None)
    axis_key_scaling.grid(alpha=0.4)

    _save(scaling_figure, out_dir / "evaluation_scaling.png", presentation=presentation)

    payload = {
        "key": str(args.key),
        "naca": str(args.naca),
        "chord_m": float(args.chord),
        "twist_deg": float(args.twist_deg),
        "grid": {"nx": int(args.nx), "ny": int(args.ny)},
        "repeats": int(args.repeats),
        "memoisation": {
            "distinct_evaluations_shared_memo": int(shared_evaluations),
            "distinct_evaluations_memo_per_region": int(isolated_evaluations),
            "seconds_shared_memo": float(shared_seconds),
            "seconds_memo_per_region": float(isolated_seconds),
            "speedup": float(speedup),
        },
        "support_boxes": culling_rows,
        "grid_scaling": grid_rows,
        "layout_scaling": key_rows,
        "naive_node_visits": naive_visits,
        "note": (
            "expr_bbox culling is reported as available work reduction; "
            "eval_expr_batch does not prune, props_emit/area_mesh.py applies it."
        ),
    }
    json_path = out_dir / "evaluation_cost_pruning.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
