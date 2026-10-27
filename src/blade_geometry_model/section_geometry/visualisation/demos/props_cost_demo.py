"""
visualisation/demos/props_cost_demo.py
===========================
Section-property backends: wall time vs accuracy (Pareto).

For a canonical 2D-CN NACA 2412 section (0° and 15° twist), times:

* hybrid (constructive wall solids + residual QT)
* full-label quadtree at depths 5..8
* dense uniform grid (reference accuracy / cost ceiling)

Writes:

* ``props_cost_pareto.png`` — relative area error vs wall time
* ``props_cost_table.png`` — summary table figure
* ``props_cost.json`` — raw timings and errors

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/props_cost_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

os.environ.setdefault("SECTION_GEOMETRY_PROFILE", "1")

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
    INK,
    MUTED,
    NAVY,
    RED,
    SURFACE,
    TEAL,
)

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
)
from blade_geometry_model.section_geometry.constructive import (
    emit_constructive_polygons,
    section_properties_from_polygons,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.sdf.profiling import (
    get_stats,
    reset_stats,
)

_DPI = DEFAULT_DPI
_DEPTHS = (5, 6, 7, 8)
_TRACK = (
    "skin_suction_surface",
    "skin_pressure_surface",
    "shear_web_0",
    "shear_web_1",
    "spar_cap_suction_surface",
    "spar_cap_pressure_surface",
)


def _build(af: AirfoilSDF, *, twist_rad: float) -> MultiCellSection:
    return MultiCellSection.twin_web(
        af,
        x_fore=0.20,
        x_aft=0.50,
        skin_thickness=0.004,
        web_thickness=0.004,
        cap_height=(0.014, 0.012),
        web_alignment="chord_normal",
        twist_angle=float(twist_rad),
        structural_family="D",
        core_enabled=True,
    )


def _mean_rel_area_err(
    props: dict[str, dict[str, float]],
    ref: dict[str, dict[str, float]],
    labels: tuple[str, ...] = _TRACK,
) -> float:
    errs = []
    for lab in labels:
        if lab not in props or lab not in ref:
            continue
        a = float(props[lab]["area"])
        ar = float(ref[lab]["area"])
        if not np.isfinite(a) or not np.isfinite(ar) or ar <= 0.0:
            continue
        errs.append(abs(a - ar) / ar)
    return float(np.mean(errs)) if errs else float("nan")


def _run_case(
    section: MultiCellSection,
    af: AirfoilSDF,
    *,
    twist_rad: float,
    tag: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    # Constructive-only reference for wall-solid labels (hybrid without residual cost).
    t0 = time.perf_counter()
    polys = emit_constructive_polygons(section)
    constr = section_properties_from_polygons(polys)
    constr_s = time.perf_counter() - t0
    rows.append(
        {
            "tag": tag,
            "backend": "constructive",
            "max_depth": None,
            "wall_time_s": constr_s,
            "phi_eval_count": 0,
            "leaf_count": 0,
            "mean_rel_area_err_vs_constructive": 0.0,
            "n_labels": len(constr),
        }
    )

    # Hybrid at depth 7 (production default class).
    spec7 = quadtree_spec_for_airfoil(af, chord=1.0, max_depth=7, max_cells=150000)
    reset_stats()
    t0 = time.perf_counter()
    hy, hy_meta = compute_section_properties(
        section,
        backend="hybrid",
        quadtree_spec=spec7,
        airfoil_for_spec=af,
        twist_rad=float(twist_rad),
    )
    hy_s = time.perf_counter() - t0
    st = get_stats()
    rows.append(
        {
            "tag": tag,
            "backend": "hybrid",
            "max_depth": 7,
            "wall_time_s": hy_s,
            "phi_eval_count": int(st.phi_eval_count) if st else 0,
            "leaf_count": int(
                (hy_meta.get("quadtree_stats") or {}).get("leaf_count", 0)
            ),
            "mean_rel_area_err_vs_constructive": _mean_rel_area_err(hy, constr),
            "constructive_s": hy_meta.get("constructive_s"),
            "residual_quadtree_s": hy_meta.get("residual_quadtree_s"),
            "n_labels": len(hy),
        }
    )

    for depth in _DEPTHS:
        spec = quadtree_spec_for_airfoil(
            af, chord=1.0, max_depth=int(depth), max_cells=200000
        )
        reset_stats()
        t0 = time.perf_counter()
        qt, qt_meta = compute_section_properties(
            section,
            backend="quadtree",
            quadtree_spec=spec,
            airfoil_for_spec=af,
            twist_rad=float(twist_rad),
        )
        qt_s = time.perf_counter() - t0
        st = get_stats()
        qt_stats = qt_meta.get("quadtree_stats") or {}
        rows.append(
            {
                "tag": tag,
                "backend": "quadtree",
                "max_depth": int(depth),
                "wall_time_s": qt_s,
                "phi_eval_count": int(st.phi_eval_count) if st else 0,
                "leaf_count": int(qt_stats.get("leaf_count", 0)),
                "mean_rel_area_err_vs_constructive": _mean_rel_area_err(qt, constr),
                "n_labels": len(qt),
            }
        )

    # Dense grid — accuracy/cost ceiling (same labels as constructive).
    grid = SDFGrid.from_airfoil(af, padding=0.06, nx=600, ny=240)
    reset_stats()
    t0 = time.perf_counter()
    grid_props, grid_meta = compute_section_properties(
        section,
        backend="grid",
        grid=grid,
        airfoil_for_spec=af,
        twist_rad=float(twist_rad),
    )
    grid_s = time.perf_counter() - t0
    st = get_stats()
    rows.append(
        {
            "tag": tag,
            "backend": "grid",
            "max_depth": None,
            "wall_time_s": grid_s,
            "phi_eval_count": int(st.phi_eval_count) if st else int(grid.X.size),
            "leaf_count": 0,
            "mean_rel_area_err_vs_constructive": _mean_rel_area_err(grid_props, constr),
            "grid_nx": grid.nx,
            "grid_ny": grid.ny,
            "n_labels": len(grid_props),
            "meta_backend": grid_meta.get("backend"),
        }
    )
    return rows


def _save_fig(fig, path: Path) -> None:
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)
    print(f"Wrote {path}")


#: Below this the two exact paths disagree by less than the comparison can
#: resolve, so their reported error is a property of the measurement and not of
#: the backend. Points are clamped here and the floor is drawn, because a
#: clamped value plotted bare reads as a result.
_ERROR_FLOOR = 1e-12

_BACKEND_COLOURS = {
    "constructive": NAVY,
    "hybrid": TEAL,
    "quadtree": "#EF6C00",
    "grid": RED,
}
# Distinct markers so constructive and hybrid stay readable where they land on
# top of each other, which is the whole point of the comparison.
_BACKEND_MARKERS = {"constructive": "o", "hybrid": "s", "quadtree": "^", "grid": "D"}


def _pareto_frontier(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Non-dominated points, ascending in wall time.

    A point survives when nothing else is both faster and more accurate.
    """
    frontier: list[tuple[float, float]] = []
    best_error = float("inf")
    for x, y in sorted(points):
        if y < best_error:
            frontier.append((x, y))
            best_error = y
    return frontier


def _plot_pareto(all_rows: list[dict[str, Any]], path: Path) -> None:
    """Accuracy against cost, on one panel.

    The empty middle of this plot is the message, not a defect: there is nothing
    between the sampled backends and the exact ones, and the gap is annotated
    rather than left as white space. Two twists share the panel — filled markers
    at 0 degrees, hollow at 15 — because splitting them into two panels halved
    each to say the same thing twice.
    """
    figure, axis = plt.subplots(figsize=(11.0, 4.8))

    for row in all_rows:
        backend = str(row["backend"])
        x = max(float(row["wall_time_s"]), 1e-6)
        y = max(float(row["mean_rel_area_err_vs_constructive"]), _ERROR_FLOOR)
        twisted = row["tag"] != "twist_0deg"
        colour = _BACKEND_COLOURS.get(backend, INK)
        axis.scatter(
            x,
            y,
            s=132,
            facecolors="none" if twisted else colour,
            edgecolors=colour,
            marker=_BACKEND_MARKERS.get(backend, "o"),
            linewidths=1.6,
            zorder=4,
        )

    # Depth labels only on the untwisted series, and offset by rank rather than
    # by depth: deep quadtrees converge onto each other, so a fixed offset puts
    # d=7 and d=8 in the same few pixels and they overprint.
    quadtree_rows = sorted(
        (r for r in all_rows if r["backend"] == "quadtree" and r["tag"] == "twist_0deg"),
        key=lambda r: float(r["wall_time_s"]),
    )
    for rank, row in enumerate(quadtree_rows):
        axis.annotate(
            f"d={int(row['max_depth'])}",
            (max(float(row["wall_time_s"]), 1e-6),
             max(float(row["mean_rel_area_err_vs_constructive"]), _ERROR_FLOOR)),
            textcoords="offset points",
            xytext=(12, 9) if rank % 2 == 0 else (12, -18),
            fontsize=8.5,
            color=_BACKEND_COLOURS["quadtree"],
        )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_ylim(3e-14, 3.0)

    # The frontier, drawn. It was called a Pareto plot and had no frontier on it.
    # Here the frontier usually collapses to a single point, because the exact
    # paths are simultaneously the fastest and the most accurate — which is the
    # argument the slide is making, so say it rather than drawing a line of one.
    frontier = _pareto_frontier(
        [
            (
                max(float(r["wall_time_s"]), 1e-6),
                max(float(r["mean_rel_area_err_vs_constructive"]), _ERROR_FLOOR),
            )
            for r in all_rows
            if r["tag"] == "twist_0deg"
        ]
    )
    if len(frontier) > 1:
        axis.plot(
            [p[0] for p in frontier],
            [p[1] for p in frontier],
            color=MUTED,
            linewidth=1.2,
            linestyle="--",
            zorder=2,
            label="Pareto frontier (0°)",
        )
        axis.legend(fontsize=8.5, loc="lower right")
    else:
        # Mid-right, not bottom-right: the floor caption already owns the bottom
        # edge and the two overprinted.
        axis.text(
            0.985,
            0.40,
            "the frontier is one point: the exact paths are\n"
            "both the fastest and the most accurate",
            transform=axis.transAxes,
            ha="right",
            va="center",
            fontsize=9.5,
            color=INK,
        )

    # The measurement floor, declared.
    axis.axhspan(3e-14, _ERROR_FLOOR, color=SURFACE, zorder=1)
    axis.axhline(_ERROR_FLOOR, color=MUTED, linewidth=0.9, linestyle=":", zorder=2)
    axis.text(
        0.985,
        _ERROR_FLOOR * 1.35,
        f"at or below the measurement floor — clamped to {_ERROR_FLOOR:g}",
        transform=axis.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8.5,
        style="italic",
        color=MUTED,
    )

    # The two decisive numbers, computed rather than asserted.
    untwisted = [r for r in all_rows if r["tag"] == "twist_0deg"]
    exact = [r for r in untwisted if r["backend"] in ("constructive", "hybrid")]
    grid_rows = [r for r in untwisted if r["backend"] == "grid"]
    if exact and grid_rows:
        exact_time = min(float(r["wall_time_s"]) for r in exact)
        grid_row = grid_rows[0]
        grid_time = float(grid_row["wall_time_s"])
        grid_error = max(float(grid_row["mean_rel_area_err_vs_constructive"]), _ERROR_FLOOR)
        speedup = grid_time / exact_time if exact_time > 0 else float("nan")
        orders = np.log10(grid_error / _ERROR_FLOOR)

        axis.annotate(
            "",
            xy=(exact_time, _ERROR_FLOOR),
            xytext=(exact_time, grid_error),
            arrowprops={"arrowstyle": "<->", "color": INK, "linewidth": 1.3},
            zorder=5,
        )
        # "At least": the exact paths' measured error is zero, so the ratio
        # against the clamp is a lower bound on the true separation, not an
        # estimate of it.
        axis.text(
            exact_time * 1.22,
            np.sqrt(grid_error * _ERROR_FLOOR),
            f"≥{orders:.0f} orders\nmore accurate",
            fontsize=10.5,
            fontweight="bold",
            color=INK,
            va="center",
        )
        # Well below the quadtree cluster, in the empty band: at the grid's own
        # error this arrow ran straight through the d=5..d=8 points and their
        # labels, and the caption sat on top of a marker.
        speed_arrow_y = 1e-5
        axis.annotate(
            "",
            xy=(exact_time, speed_arrow_y),
            xytext=(grid_time, speed_arrow_y),
            arrowprops={"arrowstyle": "<->", "color": INK, "linewidth": 1.3},
            zorder=5,
        )
        axis.text(
            np.sqrt(exact_time * grid_time),
            speed_arrow_y * 3.0,
            f"~{speedup:.0f}x faster",
            fontsize=10.5,
            fontweight="bold",
            color=INK,
            ha="center",
        )

    # Direct labels instead of a legend box: the backends sit in well-separated
    # clusters, so naming them where they are beats a key in a corner.
    # (offset, which point to hang the label on). The quadtree label hangs off a
    # mid-depth point and sits to its left: anchored on the shallowest point it
    # landed under the 15-degree markers, which cluster at the top of the panel.
    label_anchors = {
        "grid": ((-6, 20), "fastest"),
        "quadtree": ((-52, 2), "middle"),
        "constructive": ((4, -26), "fastest"),
        "hybrid": ((-30, 6), "fastest"),
    }
    for backend, (offset, anchor_rule) in label_anchors.items():
        rows = sorted(
            (r for r in untwisted if r["backend"] == backend),
            key=lambda r: float(r["wall_time_s"]),
        )
        if not rows:
            continue
        row = rows[0] if anchor_rule == "fastest" else rows[len(rows) // 2]
        axis.annotate(
            backend,
            (
                max(float(row["wall_time_s"]), 1e-6),
                max(float(row["mean_rel_area_err_vs_constructive"]), _ERROR_FLOOR),
            ),
            textcoords="offset points",
            xytext=offset,
            fontsize=10.5,
            fontweight="bold",
            color=_BACKEND_COLOURS[backend],
            ha="center",
        )

    axis.set_xlabel("wall time [s]  →  slower")
    axis.set_ylabel("mean rel. area error vs constructive  →  worse")
    axis.set_title(
        "Section property backends (2D-CN NACA 2412) — lower-left is better\n"
        "filled = 0° twist · hollow = 15° twist",
        fontsize=12,
    )
    axis.grid(True, which="major", alpha=0.30)

    _save_fig(figure, path)


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "props_cost"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    args = parser.parse_args(argv)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", n_points=120, chord=1.0)
    all_rows: list[dict[str, Any]] = []
    for tag, twist in (("twist_0deg", 0.0), ("twist_15deg", float(np.deg2rad(15.0)))):
        section = _build(af, twist_rad=twist)
        all_rows.extend(_run_case(section, af, twist_rad=twist, tag=tag))

    _plot_pareto(all_rows, out_dir / "props_cost_pareto.png")

    # Table figure
    fig2, ax2 = plt.subplots(figsize=(11.5, 0.45 * len(all_rows) + 1.8))
    ax2.axis("off")
    col_labels = [
        "case",
        "backend",
        "depth",
        "t [s]",
        "φ-evals",
        "leaves",
        "mean |ΔA|/A",
    ]
    cell = []
    for r in all_rows:
        cell.append(
            [
                r["tag"],
                r["backend"],
                "-" if r["max_depth"] is None else str(r["max_depth"]),
                f"{r['wall_time_s']:.3f}",
                f"{r['phi_eval_count']:,}",
                f"{r['leaf_count']:,}",
                f"{r['mean_rel_area_err_vs_constructive']:.2e}",
            ]
        )
    table = ax2.table(
        cellText=cell,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.25)
    ax2.set_title("Property backend cost / accuracy summary", fontsize=12, pad=12)
    _save_fig(fig2, out_dir / "props_cost_table.png")

    payload = {
        "airfoil": "NACA 2412",
        "layout": "twin_web family D (2D-CN class)",
        "track_labels": list(_TRACK),
        "reference": "constructive wall-solid polygons (shoelace)",
        "rows": all_rows,
        "notes": (
            "Hybrid = constructive wall solids + residual QT for inserts/cores path. "
            "Area error is mean relative |dA|/A over track labels vs constructive. "
            "At non-zero twist, full-label QT can under-resolve vs constructive "
            "(supports/refinement); hybrid remains the production path."
        ),
    }
    json_path = out_dir / "props_cost.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    for r in all_rows:
        print(
            f"  {r['tag']:<12} {r['backend']:<12} d={str(r['max_depth']):<4} "
            f"t={r['wall_time_s']:.3f}s  err={r['mean_rel_area_err_vs_constructive']:.2e}"
        )
    return out_dir


if __name__ == "__main__":
    main()
