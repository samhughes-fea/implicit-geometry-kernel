"""
visualisation/demos/quadtree_vs_shoelace_convergence_study.py
====================================================
Quadtree vs shoelace (exact polygon) section-property comparison.

Extends ``quadtree_convergence_study.py`` with a second integration method:
extract an explicit polygon boundary from the *same* implicit/SDF geometry
via marching squares, then compute area/centroid/second moments exactly via
the shoelace formula (Green's theorem,
``blade_geometry_model.shoelace_integration``).

Reports both methods' accuracy (vs the same high-resolution uniform-grid
reference used by the quadtree study) and wall time, and deliberately
separates two costs that get conflated if you only ask "is shoelace fast":

1. Shoelace's own per-call cost, once a polygon exists -- this is what was
   "instantaneous" under the retired *explicit* geometry kernel, because
   that kernel's geometry representation *was* a polygon already, at zero
   extra cost.
2. The SDF -> polygon extraction cost (marching squares) that the *current*
   implicit/SDF kernel would newly have to pay to get that polygon at all,
   since it represents geometry as signed-distance fields / CSG, not
   explicit vertex loops. This repo's only available extractor
   (``blade_geometry_model.root_blend_model.contouring.
   marching_squares_zero_contour``) is a pure-Python, unvectorized
   grid scan (no ``scikit-image`` in this environment) -- it is *not*
   representative of what a properly vectorized extractor would cost, and
   that caveat matters for interpreting the results below.

Run from the repository root:
    python src/blade_geometry_model/section_geometry/visualisation/demos/quadtree_vs_shoelace_convergence_study.py
"""

from __future__ import annotations

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
from shapely.geometry import Polygon as ShapelyPolygon

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

# Load the existing quadtree convergence study as a module so its
# section-building / grid-reference / depth-sweep logic is reused verbatim
# rather than duplicated (keeps the two studies apples-to-apples by
# construction).
import sys

_qcs_path = Path(__file__).resolve().parent / "quadtree_convergence_study.py"
_qcs_spec = importlib.util.spec_from_file_location("quadtree_convergence_study", _qcs_path)
qcs = importlib.util.module_from_spec(_qcs_spec)
assert _qcs_spec.loader is not None
sys.modules["quadtree_convergence_study"] = qcs  # dataclass() needs this in sys.modules
_qcs_spec.loader.exec_module(qcs)

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SDFGrid,
)
from blade_geometry_model.root_blend_model.contouring import marching_squares_zero_contour
from blade_geometry_model.shoelace_integration import section_properties_S

PROP_KEYS = qcs.PROP_KEYS
TRACK_LABELS = qcs.TRACK_LABELS
_rel_err = qcs._rel_err
_pick_props = qcs._pick_props

# Marching-squares extraction grid resolutions to sweep. Kept modest: the
# available extractor is an unvectorized pure-Python grid scan, and its cost
# grows roughly with nx*ny -- 400x200 already costs several seconds per
# component (see the SHOELACE_PER_CALL_REPS isolation below for contrast
# with shoelace's own near-instantaneous per-call cost).
MS_RESOLUTIONS = [(50, 25), (100, 50), (200, 100), (400, 200)]
SHOELACE_PER_CALL_REPS = 200


def _polygon_props_from_contour(contour: np.ndarray) -> dict[str, float]:
    poly = ShapelyPolygon(contour)
    A, (Cx, Cy), Ixx, Iyy, _Ixy = section_properties_S(poly)
    return {"area": A, "cx": Cx, "cy": Cy, "Ixx": Ixx, "Iyy": Iyy}


@dataclass
class ShoelaceRow:
    nx: int
    ny: int
    n_grid_cells: int
    extract_time_s: float
    shoelace_time_s: float
    total_time_s: float
    n_vertices: int
    properties: dict[str, dict[str, float]] = field(default_factory=dict)
    err_vs_grid: dict[str, dict[str, float]] = field(default_factory=dict)


def sweep_shoelace(
    section: Any,
    airfoil_b: Any,
    *,
    twist_rad: float,
    resolutions: list[tuple[int, int]],
) -> list[ShoelaceRow]:
    rows: list[ShoelaceRow] = []
    for nx, ny in resolutions:
        grid = SDFGrid.from_airfoil(airfoil_b, nx=nx, ny=ny)
        x1d = grid.X[0, :]
        y1d = grid.Y[:, 0]
        picked: dict[str, dict[str, float]] = {}
        extract_total = 0.0
        shoelace_total = 0.0
        n_vertices = 0
        for label in TRACK_LABELS:
            comp = None
            if label in section:
                comp = section[label]
            else:
                for lbl in section:
                    if qcs._label_key(lbl) == label:
                        comp = section[lbl]
                        break
            if comp is None:
                continue
            phi = grid.eval(comp)
            t0 = time.perf_counter()
            try:
                contour = marching_squares_zero_contour(x1d, y1d, phi)
            except ValueError:
                continue
            t1 = time.perf_counter()
            props = _polygon_props_from_contour(contour)
            t2 = time.perf_counter()
            extract_total += t1 - t0
            shoelace_total += t2 - t1
            n_vertices += len(contour)
            picked[label] = props
        rows.append(
            ShoelaceRow(
                nx=nx,
                ny=ny,
                n_grid_cells=nx * ny,
                extract_time_s=extract_total,
                shoelace_time_s=shoelace_total,
                total_time_s=extract_total + shoelace_total,
                n_vertices=n_vertices,
                properties=picked,
            )
        )
    return rows


def attach_shoelace_errors(rows: list[ShoelaceRow], grid_ref: dict[str, dict[str, float]]) -> None:
    for row in rows:
        for label in TRACK_LABELS:
            props = row.properties.get(label)
            if props is None:
                continue
            gref = grid_ref.get(label, {})
            row.err_vs_grid[label] = {
                k: _rel_err(props[k], gref[k]) for k in PROP_KEYS if k in props and k in gref
            }


def isolate_shoelace_per_call_cost(section: Any, airfoil_b: Any) -> dict[str, float]:
    """Time *only* ``section_properties_S`` on an already-extracted polygon.

    Isolates the cost that was "instantaneous" under the retired explicit
    kernel (a polygon evaluation, not an SDF->polygon extraction).
    """
    grid = SDFGrid.from_airfoil(airfoil_b, nx=400, ny=200)
    x1d = grid.X[0, :]
    y1d = grid.Y[:, 0]
    comp = section["outer_skin"] if "outer_skin" in section else None
    if comp is None:
        for lbl in section:
            if qcs._label_key(lbl) == "outer_skin":
                comp = section[lbl]
                break
    phi = grid.eval(comp)
    contour = marching_squares_zero_contour(x1d, y1d, phi)
    poly = ShapelyPolygon(contour)

    t0 = time.perf_counter()
    for _ in range(SHOELACE_PER_CALL_REPS):
        section_properties_S(poly)
    t1 = time.perf_counter()
    return {
        "n_vertices": len(contour),
        "n_reps": SHOELACE_PER_CALL_REPS,
        "total_s": t1 - t0,
        "per_call_s": (t1 - t0) / SHOELACE_PER_CALL_REPS,
    }


def render_comparison_figure(
    qt_rows: list[Any],
    sl_rows: list[ShoelaceRow],
    *,
    tag: str,
    twist_deg: float,
    out_png: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax_a, ax_b = axes

    qt_times = np.array([r.time_s for r in qt_rows], dtype=float)
    sl_times = np.array([r.total_time_s for r in sl_rows], dtype=float)
    sl_shoelace_only = np.array([r.shoelace_time_s for r in sl_rows], dtype=float)

    for ax, label, title in ((ax_a, "outer_skin", "outer_skin"), (ax_b, "shear_web_0", "web_0 (thin feature)")):
        qt_err = np.array(
            [r.err_vs_grid.get(label, {}).get("area", float("nan")) for r in qt_rows], dtype=float
        )
        sl_err = np.array(
            [r.err_vs_grid.get(label, {}).get("area", float("nan")) for r in sl_rows], dtype=float
        )
        mask_qt = np.isfinite(qt_err) & (qt_err > 0) & (qt_times > 0)
        mask_sl = np.isfinite(sl_err) & (sl_err > 0) & (sl_times > 0)
        if mask_qt.any():
            ax.loglog(qt_times[mask_qt], qt_err[mask_qt], "o-", color="tab:blue", label="quadtree (total)")
        if mask_sl.any():
            ax.loglog(
                sl_times[mask_sl], sl_err[mask_sl], "s-", color="tab:red",
                label="shoelace (extract+eval)",
            )
            ax.loglog(
                np.clip(sl_shoelace_only[mask_sl], 1e-9, None), sl_err[mask_sl], "^--",
                color="tab:orange", label="shoelace (eval only)",
            )
        ax.set_xlabel("wall time [s]")
        ax.set_ylabel("relative area error vs grid ref")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    fig.suptitle(f"Quadtree vs shoelace — {tag} (twist={twist_deg:g} deg)", fontsize=11)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, out_png)
    plt.close(fig)


def shoelace_rows_to_records(rows: list[ShoelaceRow]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in rows:
        rec: dict[str, Any] = {
            "nx": row.nx,
            "ny": row.ny,
            "n_grid_cells": row.n_grid_cells,
            "extract_time_s": row.extract_time_s,
            "shoelace_time_s": row.shoelace_time_s,
            "total_time_s": row.total_time_s,
            "n_vertices": row.n_vertices,
        }
        for label in TRACK_LABELS:
            for k, v in row.properties.get(label, {}).items():
                rec[f"{label}.{k}"] = v
            for k, v in row.err_vs_grid.get(label, {}).items():
                rec[f"{label}.err_grid.{k}"] = v
        records.append(rec)
    return records


def print_comparison_summary(
    tag: str,
    qt_rows: list[Any],
    sl_rows: list[ShoelaceRow],
    per_call: dict[str, float],
) -> None:
    print(f"\n=== {tag}: quadtree vs shoelace ===")
    print(f"{'depth':>5} {'qt_time_s':>10} {'qt_skin_A_err':>14}")
    for r in qt_rows:
        e = r.err_vs_grid.get("outer_skin", {}).get("area", float("nan"))
        print(f"{r.max_depth:5d} {r.time_s:10.3f} {e:14.3e}")
    print(f"\n{'nx':>5} {'ny':>5} {'extract_s':>10} {'shoelace_s':>11} {'total_s':>9} {'skin_A_err':>12}")
    for r in sl_rows:
        e = r.err_vs_grid.get("outer_skin", {}).get("area", float("nan"))
        print(f"{r.nx:5d} {r.ny:5d} {r.extract_time_s:10.4f} {r.shoelace_time_s:11.6f} {r.total_time_s:9.4f} {e:12.3e}")
    print(
        f"\nShoelace per-call cost once a polygon exists "
        f"({per_call['n_vertices']} vertices, {per_call['n_reps']} reps): "
        f"{per_call['per_call_s'] * 1e6:.2f} us/call"
    )
    if qt_rows:
        fastest_qt = min(qt_rows, key=lambda r: r.time_s)
        print(
            f"Fastest quadtree depth in this sweep: depth={fastest_qt.max_depth}, "
            f"time={fastest_qt.time_s:.3f}s -- "
            f"{fastest_qt.time_s / per_call['per_call_s']:.0f}x slower than one shoelace eval."
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
    section = qcs.make_section(airfoil_base, web_x, twist_rad)
    airfoil_b = airfoil_base.rotate(twist_rad) if abs(twist_rad) > 1e-12 else airfoil_base
    grid = SDFGrid.from_airfoil(airfoil_b, nx=int(ref_nx), ny=int(ref_ny))
    grid_ref = qcs.grid_reference_properties(section, grid)

    qt_rows = qcs.sweep_depths(
        section, airfoil_b, chord=chord, twist_rad=twist_rad, depths=depths, max_cells=max_cells,
    )
    qcs.attach_errors(qt_rows, grid_ref)

    sl_rows = sweep_shoelace(section, airfoil_b, twist_rad=twist_rad, resolutions=MS_RESOLUTIONS)
    attach_shoelace_errors(sl_rows, grid_ref)

    per_call = isolate_shoelace_per_call_cost(section, airfoil_b)

    out_dir.mkdir(parents=True, exist_ok=True)
    sl_records = shoelace_rows_to_records(sl_rows)
    with (out_dir / f"shoelace_{tag}.csv").open("w", newline="", encoding="utf-8") as fh:
        if sl_records:
            fieldnames: list[str] = []
            for rec in sl_records:
                for k in rec:
                    if k not in fieldnames:
                        fieldnames.append(k)
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(sl_records)

    payload = {
        "tag": tag,
        "twist_deg": twist_deg,
        "chord_m": chord,
        "grid_reference": grid_ref,
        "quadtree_rows": qcs.rows_to_records(qt_rows),
        "shoelace_rows": sl_records,
        "shoelace_per_call_cost": per_call,
    }
    with (out_dir / f"comparison_{tag}.json").open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)

    render_comparison_figure(
        qt_rows, sl_rows, tag=tag, twist_deg=twist_deg, out_png=out_dir / f"comparison_{tag}.png",
    )
    print_comparison_summary(tag, qt_rows, sl_rows, per_call)
    print(f"Wrote {out_dir / f'comparison_{tag}.png'}")
    print(f"Wrote {out_dir / f'comparison_{tag}.csv'}")
    print(f"Wrote {out_dir / f'comparison_{tag}.json'}")
    return payload


def main() -> int:
    chord = 1.0
    depths = list(range(4, 8))  # narrower than the base study (4..9) to bound total runtime
    max_cells = 200_000
    ref_nx, ref_ny = 1200, 600
    out_dir = (Path(__file__).resolve().parents[1] / "output" / "quadtree_vs_shoelace_convergence").resolve()

    airfoil = AirfoilSDF.from_naca("2412", chord=chord)
    cases = [
        ((0.15 * chord, 0.50 * chord), 0.0, "2web_t0"),
        ((0.15 * chord, 0.50 * chord), 15.0, "2web_t15"),
    ]
    for web_x, twist_deg, tag in cases:
        run_case(
            airfoil, web_x, twist_deg, tag,
            chord=chord, depths=depths, max_cells=max_cells,
            ref_nx=ref_nx, ref_ny=ref_ny, out_dir=out_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
