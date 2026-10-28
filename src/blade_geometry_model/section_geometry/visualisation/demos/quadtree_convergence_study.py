"""
visualisation/demos/quadtree_convergence_study.py
======================================
Quadtree depth / section-property convergence study.

Sweeps ``max_depth`` for the shared adaptive quadtree backend and reports how
component section properties (area, centroid, second moments) converge toward
a high-resolution uniform-grid reference, alongside cost metrics (leaf count,
phi-eval count, wall time).

Default sweep: depth 4..9 for two cases (untwisted and 15 deg twist). Expect
~5--15 minutes depending on hardware (dominated by depth 9 and the 1200x600
grid baseline).

Run from the repository root:
    python src/blade_geometry_model/section_geometry/visualisation/demos/quadtree_convergence_study.py
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

# Enable phi-eval profiling before engine imports.
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

from blade_geometry_model.section_geometry.props_emit.area_mesh import section_properties_via_quadtree
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import quadtree_spec_for_airfoil
from blade_geometry_model.section_geometry.sdf.profiling import get_stats, reset_stats

PROP_KEYS = ("area", "cx", "cy", "Ixx", "Iyy")
TRACK_LABELS = ("outer_skin", "shear_web_0")


def make_section(airfoil, web_x, twist_rad):
    return MultiCellSection(
        airfoil_sdf=airfoil,
        web_x_positions=list(web_x),
        web_thickness=0.004,
        web_alignment="chord_normal",
        cap_height=0.012,
        skin_thickness=0.003,
        twist_angle=twist_rad,
        core_enabled=True,
    )


def _label_key(label: Any) -> str:
    val = getattr(label, "value", label)
    return str(val)


def _rel_err(value: float, reference: float) -> float:
    if not np.isfinite(value) or not np.isfinite(reference):
        return float("nan")
    denom = max(abs(reference), 1e-30)
    return abs(value - reference) / denom


def _pick_props(props: dict[str, dict[str, float]], label: str) -> dict[str, float] | None:
    if label in props:
        return props[label]
    for key, val in props.items():
        if _label_key(key) == label:
            return val
    return None


def grid_reference_properties(
    section: Any,
    grid: SDFGrid,
    labels: tuple[str, ...] = TRACK_LABELS,
) -> dict[str, dict[str, float]]:
    """High-resolution uniform-grid reference for selected components."""
    out: dict[str, dict[str, float]] = {}
    for label in labels:
        comp = None
        if label in section:
            comp = section[label]
        else:
            for lbl in section:
                if _label_key(lbl) == label:
                    comp = section[lbl]
                    break
        if comp is None:
            continue
        phi = grid.eval(comp)
        out[label] = grid.section_properties(phi)
    return out


@dataclass
class DepthRow:
    max_depth: int
    leaf_count: int
    min_cell_width: float
    phi_eval_count: int
    time_s: float
    properties: dict[str, dict[str, float]] = field(default_factory=dict)
    err_vs_grid: dict[str, dict[str, float]] = field(default_factory=dict)
    err_vs_finest_qt: dict[str, dict[str, float]] = field(default_factory=dict)


def sweep_depths(
    section: Any,
    airfoil_b: Any,
    *,
    chord: float,
    twist_rad: float,
    depths: list[int],
    max_cells: int,
) -> list[DepthRow]:
    rows: list[DepthRow] = []
    for depth in depths:
        spec = quadtree_spec_for_airfoil(
            airfoil_b,
            chord=float(chord),
            max_depth=int(depth),
            max_cells=int(max_cells),
        )
        reset_stats()
        t0 = time.perf_counter()
        props, _mesh, qt_stats = section_properties_via_quadtree(
            section, spec, twist_rad=float(twist_rad)
        )
        elapsed = time.perf_counter() - t0
        prof = get_stats()
        phi_evals = int(prof.phi_eval_count) if prof is not None else 0
        picked: dict[str, dict[str, float]] = {}
        for label in TRACK_LABELS:
            p = _pick_props(props, label)
            if p is not None:
                picked[label] = {k: float(p[k]) for k in PROP_KEYS if k in p}
        rows.append(
            DepthRow(
                max_depth=int(depth),
                leaf_count=int(qt_stats["leaf_count"]),
                min_cell_width=float(qt_stats["min_cell_width"]),
                phi_eval_count=phi_evals,
                time_s=float(elapsed),
                properties=picked,
            )
        )
    return rows


def attach_errors(
    rows: list[DepthRow],
    grid_ref: dict[str, dict[str, float]],
) -> None:
    if not rows:
        return
    finest = rows[-1].properties
    for row in rows:
        for label in TRACK_LABELS:
            props = row.properties.get(label)
            if props is None:
                continue
            gref = grid_ref.get(label, {})
            fref = finest.get(label, {})
            row.err_vs_grid[label] = {
                k: _rel_err(props[k], gref[k]) for k in PROP_KEYS if k in props and k in gref
            }
            row.err_vs_finest_qt[label] = {
                k: _rel_err(props[k], fref[k]) for k in PROP_KEYS if k in props and k in fref
            }


def _series(rows: list[DepthRow], label: str, err_key: str, prop: str) -> np.ndarray:
    vals = []
    for row in rows:
        bucket = row.err_vs_grid if err_key == "grid" else row.err_vs_finest_qt
        v = bucket.get(label, {}).get(prop, float("nan"))
        vals.append(v)
    return np.asarray(vals, dtype=float)


def render_convergence_figure(
    rows: list[DepthRow],
    *,
    tag: str,
    twist_deg: float,
    grid_ref: dict[str, dict[str, float]],
    out_png: Path,
) -> None:
    depths = np.array([r.max_depth for r in rows], dtype=float)
    leaves = np.array([r.leaf_count for r in rows], dtype=float)
    phi_evals = np.array([r.phi_eval_count for r in rows], dtype=float)
    times = np.array([r.time_s for r in rows], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # Panel A: outer_skin relative error vs depth (vs grid reference)
    for prop, color in zip(("area", "Ixx", "Iyy"), ("C0", "C1", "C2")):
        err = _series(rows, "outer_skin", "grid", prop)
        mask = np.isfinite(err) & (err > 0)
        if mask.any():
            ax_a.semilogy(depths[mask], err[mask], "o-", color=color, label=prop)
    ax_a.set_xlabel("max_depth")
    ax_a.set_ylabel("relative error vs grid")
    ax_a.set_title(f"outer_skin convergence ({tag}, twist={twist_deg:g} deg)")
    ax_a.grid(True, which="both", alpha=0.3)
    ax_a.legend(loc="best")

    # Panel B: error vs cost (log-log)
    for prop, color in zip(("area", "Ixx", "Iyy"), ("C0", "C1", "C2")):
        err = _series(rows, "outer_skin", "grid", prop)
        mask = np.isfinite(err) & (err > 0) & (leaves > 0)
        if mask.any():
            ax_b.loglog(leaves[mask], err[mask], "o-", color=color, label=f"{prop} vs leaves")
    ax_b.set_xlabel("leaf count")
    ax_b.set_ylabel("relative error vs grid")
    ax_b.set_title("accuracy vs cost (outer_skin)")
    ax_b.grid(True, which="both", alpha=0.3)
    ax_b.legend(loc="best", fontsize=8)

    # Panel C: cost vs depth
    ax_c.plot(depths, leaves, "o-", color="tab:blue", label="leaves")
    ax_c.set_xlabel("max_depth")
    ax_c.set_ylabel("leaf count", color="tab:blue")
    ax_c.tick_params(axis="y", labelcolor="tab:blue")
    ax_c2 = ax_c.twinx()
    ax_c2.plot(depths, phi_evals, "s--", color="tab:orange", label="phi evals")
    ax_c2.plot(depths, times, "^:", color="tab:green", label="time [s]")
    ax_c2.set_ylabel("phi evals / time [s]")
    ax_c.set_title("cost vs depth")
    ax_c.grid(True, alpha=0.3)
    lines1, labels1 = ax_c.get_legend_handles_labels()
    lines2, labels2 = ax_c2.get_legend_handles_labels()
    ax_c.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=8)

    # Panel D: thin feature (web_0) area error vs depth
    err_web = _series(rows, "shear_web_0", "grid", "area")
    mask_w = np.isfinite(err_web) & (err_web > 0)
    if mask_w.any():
        ax_d.semilogy(depths[mask_w], err_web[mask_w], "o-", color="tab:purple", label="vs grid")
    err_web_f = _series(rows, "shear_web_0", "finest", "area")
    mask_wf = np.isfinite(err_web_f) & (err_web_f > 0)
    if mask_wf.any():
        ax_d.semilogy(depths[mask_wf], err_web_f[mask_wf], "s--", color="tab:gray", label="vs finest qt")
    ax_d.set_xlabel("max_depth")
    ax_d.set_ylabel("relative area error")
    ax_d.set_title("web_0 (thin feature) area convergence")
    ax_d.grid(True, which="both", alpha=0.3)
    ax_d.legend(loc="best")

    g_skin = grid_ref.get("outer_skin", {})
    fig.suptitle(
        f"Quadtree convergence — {tag}  "
        f"(grid ref outer_skin area={g_skin.get('area', float('nan')):.6g} m²)",
        fontsize=11,
    )
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, out_png)
    plt.close(fig)


def rows_to_records(rows: list[DepthRow]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in rows:
        rec: dict[str, Any] = {
            "max_depth": row.max_depth,
            "leaf_count": row.leaf_count,
            "min_cell_width": row.min_cell_width,
            "phi_eval_count": row.phi_eval_count,
            "time_s": row.time_s,
        }
        for label in TRACK_LABELS:
            for k, v in row.properties.get(label, {}).items():
                rec[f"{label}.{k}"] = v
            for k, v in row.err_vs_grid.get(label, {}).items():
                rec[f"{label}.err_grid.{k}"] = v
            for k, v in row.err_vs_finest_qt.get(label, {}).items():
                rec[f"{label}.err_finest_qt.{k}"] = v
        records.append(rec)
    return records


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    fieldnames = list(records[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


def print_summary(tag: str, rows: list[DepthRow], grid_ref: dict[str, dict[str, float]]) -> None:
    print(f"\n=== {tag} ===")
    print(
        f"{'depth':>5} {'leaves':>8} {'min_w':>10} {'phi_eval':>10} {'time_s':>8}  "
        f"{'skin_A_err':>11} {'web_A_err':>11}"
    )
    for row in rows:
        skin_e = row.err_vs_grid.get("outer_skin", {}).get("area", float("nan"))
        web_e = row.err_vs_grid.get("shear_web_0", {}).get("area", float("nan"))
        print(
            f"{row.max_depth:5d} {row.leaf_count:8d} {row.min_cell_width:10.5f} "
            f"{row.phi_eval_count:10d} {row.time_s:8.2f}  "
            f"{skin_e:11.3e} {web_e:11.3e}"
        )
    g = grid_ref.get("outer_skin", {})
    print(
        f"Grid reference outer_skin: area={g.get('area', float('nan')):.6g}, "
        f"Ixx={g.get('Ixx', float('nan')):.6g}, Iyy={g.get('Iyy', float('nan')):.6g}"
    )


def run_case(
    airfoil_base: AirfoilSDF,
    web_x: tuple[float, ...],
    twist_deg: float,
    tag: str,
    *,
    chord: float,
    depths: list[int],
    max_cells: int,
    ref_nx: int,
    ref_ny: int,
    out_dir: Path,
) -> dict[str, Any]:
    twist_rad = float(np.deg2rad(twist_deg))
    section = make_section(airfoil_base, web_x, twist_rad)
    airfoil_b = airfoil_base.rotate(twist_rad) if abs(twist_rad) > 1e-12 else airfoil_base
    grid = SDFGrid.from_airfoil(airfoil_b, nx=int(ref_nx), ny=int(ref_ny))
    grid_ref = grid_reference_properties(section, grid)

    rows = sweep_depths(
        section,
        airfoil_b,
        chord=chord,
        twist_rad=twist_rad,
        depths=depths,
        max_cells=max_cells,
    )
    attach_errors(rows, grid_ref)

    records = rows_to_records(rows)
    payload = {
        "tag": tag,
        "twist_deg": twist_deg,
        "chord_m": chord,
        "depths": depths,
        "grid_reference": grid_ref,
        "rows": records,
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / f"convergence_{tag}.csv", records)
    write_json(out_dir / f"convergence_{tag}.json", payload)
    render_convergence_figure(
        rows,
        tag=tag,
        twist_deg=twist_deg,
        grid_ref=grid_ref,
        out_png=out_dir / f"convergence_{tag}.png",
    )
    print_summary(tag, rows, grid_ref)
    print(f"Wrote {out_dir / f'convergence_{tag}.png'}")
    print(f"Wrote {out_dir / f'convergence_{tag}.csv'}")
    print(f"Wrote {out_dir / f'convergence_{tag}.json'}")
    return payload


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--naca", default="2412")
    p.add_argument("--chord", type=float, default=1.0)
    p.add_argument("--min-depth", type=int, default=4)
    p.add_argument("--max-depth", type=int, default=9)
    p.add_argument("--max-cells", type=int, default=200000)
    p.add_argument("--ref-nx", type=int, default=1200)
    p.add_argument("--ref-ny", type=int, default=600)
    p.add_argument("--out-dir", type=Path, default=None)
    args = p.parse_args()

    if args.min_depth > args.max_depth:
        raise SystemExit("--min-depth must be <= --max-depth")

    out_dir = args.out_dir or (
        Path(__file__).resolve().parents[1] / "output" / "quadtree_convergence"
    )
    out_dir = Path(out_dir).resolve()
    chord = float(args.chord)
    depths = list(range(int(args.min_depth), int(args.max_depth) + 1))
    airfoil = AirfoilSDF.from_naca(str(args.naca).zfill(4), chord=chord)

    cases = [
        ((0.15 * chord, 0.50 * chord), 0.0, "2web_t0"),
        ((0.15 * chord, 0.50 * chord), 15.0, "2web_t15"),
    ]
    for web_x, twist_deg, tag in cases:
        run_case(
            airfoil,
            web_x,
            twist_deg,
            tag,
            chord=chord,
            depths=depths,
            max_cells=int(args.max_cells),
            ref_nx=int(args.ref_nx),
            ref_ny=int(args.ref_ny),
            out_dir=out_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
