"""
visualisation/demos/robustness_monte_carlo_demo.py
=======================================
Randomized robustness sweep for implicit section geometry.

Writes under ``visualisation/output/robustness_monte_carlo/``:

* ``robustness_monte_carlo.png`` - every successful section overlaid, and the
  median-normalised spread of area, Iyy and Ixx.
* ``robustness_monte_carlo.json`` - per-sample metadata and aggregate statistics.
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

from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    RED,
    TEAL,
    strip_axes_chrome,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot import (
    plot_section_hybrid,
)

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    build_section_view,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_structure_utilities.geometry_utilities.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)

_DPI = DEFAULT_DPI


def _safe_float(v: float) -> float | None:
    return float(v) if np.isfinite(v) else None


def _totals(props: dict[str, dict[str, float]]) -> dict[str, float]:
    area = float(np.nansum([float(row.get("area", np.nan)) for row in props.values()]))
    ixx = float(np.nansum([float(row.get("Ixx", np.nan)) for row in props.values()]))
    iyy = float(np.nansum([float(row.get("Iyy", np.nan)) for row in props.values()]))
    return {"area": area, "Ixx": ixx, "Iyy": iyy}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--key", default="2D-F")
    p.add_argument("--naca", default="2412")
    p.add_argument("--chord", type=float, default=1.0)
    p.add_argument("--twist-deg", type=float, default=15.0)
    p.add_argument("--n-samples", type=int, default=120)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--presentation", action="store_true", help="Deck styling.")
    args = p.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = args.out_dir or (Path(__file__).resolve().parents[1] / "output" / "robustness_monte_carlo")
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(int(args.seed))
    twist_nom_deg = float(args.twist_deg)
    layout = resolve_system_type(str(args.key))
    airfoil = AirfoilSDF.from_naca(str(args.naca), chord=float(args.chord), n_points=170)
    spec = quadtree_spec_for_airfoil(airfoil, chord=float(args.chord), max_depth=7, max_cells=140000)

    samples: list[dict[str, object]] = []
    ok_area: list[float] = []
    ok_ixx: list[float] = []
    ok_iyy: list[float] = []
    ghost_sections: list = []

    for i in range(int(args.n_samples)):
        skin = float(np.clip(rng.normal(0.0035, 0.0008), 0.0015, 0.0065))
        cap = float(np.clip(rng.normal(0.0120, 0.0025), 0.0060, 0.0200))
        web = float(np.clip(rng.normal(0.0040, 0.0010), 0.0015, 0.0080))
        twist_deg = float(np.clip(rng.normal(float(args.twist_deg), 4.0), -2.0, 30.0))
        twist_rad = float(np.deg2rad(twist_deg))

        row: dict[str, object] = {
            "i": i,
            "skin_thickness_m": skin,
            "spar_cap_thickness_m": cap,
            "shear_web_thickness_m": web,
            "twist_deg": twist_deg,
            "twist_dev_deg": twist_deg - twist_nom_deg,
            "ok": False,
            "error": None,
            "area": None,
            "Ixx": None,
            "Iyy": None,
        }
        try:
            section = build_section_view(
                airfoil,
                layout,
                twist_rad=twist_rad,
                skin_thickness=skin,
                spar_cap_thickness=cap,
                shear_web_thickness=web,
            )
            props, _meta = compute_section_properties(
                section,
                backend="hybrid",
                quadtree_spec=spec,
                chord=float(args.chord),
                twist_rad=twist_rad,
                airfoil_for_spec=airfoil,
            )
            t = _totals(props)
            good = all(np.isfinite([t["area"], t["Ixx"], t["Iyy"]])) and t["area"] > 0.0
            row["ok"] = bool(good)
            row["area"] = _safe_float(t["area"])
            row["Ixx"] = _safe_float(t["Ixx"])
            row["Iyy"] = _safe_float(t["Iyy"])
            if good:
                ok_area.append(float(t["area"]))
                ok_ixx.append(float(t["Ixx"]))
                ok_iyy.append(float(t["Iyy"]))
                ghost_sections.append(section)
        except Exception as exc:
            row["error"] = str(exc)

        samples.append(row)

    ok_mask = np.array([bool(r["ok"]) for r in samples], dtype=bool)
    area_vals = np.array([np.nan if r["area"] is None else float(r["area"]) for r in samples], dtype=float)
    fail_count = int(np.sum(~ok_mask))
    ok_count = int(np.sum(ok_mask))
    success_rate = float(ok_count / max(len(samples), 1))

    # Left: every successful section drawn on top of the others at low opacity,
    # so the envelope of the whole sweep is visible at once. Right: the spread of
    # each response about its own median (area is m², inertias m⁴, so raw values
    # cannot share an axis).
    fig = plt.figure(figsize=(14.0, 5.2))
    grid_spec = fig.add_gridspec(1, 2, width_ratios=(1.9, 1.0), wspace=0.12)
    ax_ghost = fig.add_subplot(grid_spec[0, 0])
    ax_spread = fig.add_subplot(grid_spec[0, 1])

    ghost_alpha = float(np.clip(6.0 / max(ok_count, 1), 0.01, 0.2))
    for section in ghost_sections:
        plot_section_hybrid(
            section,
            residual_mesh=None,
            ax=ax_ghost,
            alpha=ghost_alpha,
            core_alpha=0.0,
            omit_cores=True,
            title="",
            show_legend=False,
        )
    ax_ghost.set_aspect("equal", adjustable="datalim")
    strip_axes_chrome(ax_ghost)
    ax_ghost.set_title(
        f"{ok_count} randomised {args.key} sections, overlaid",
        loc="left",
    )
    ax_ghost.text(
        0.99, 0.04,
        f"{ok_count} / {len(samples)} built · {fail_count} failures",
        transform=ax_ghost.transAxes, ha="right", va="bottom",
        fontsize=16, fontweight="bold", color=TEAL if fail_count == 0 else RED,
    )

    responses = {
        "area": np.array(ok_area, dtype=float),
        "Iyy": np.array(ok_iyy, dtype=float),
        "Ixx": np.array(ok_ixx, dtype=float),
    }
    colours = {"area": MUTED, "Iyy": MUTED, "Ixx": RED}
    jitter = np.random.default_rng(0)
    for position, (name, values) in enumerate(responses.items()):
        if values.size == 0:
            continue
        normalised = values / np.median(values)
        ax_spread.scatter(
            position + jitter.uniform(-0.18, 0.18, normalised.size), normalised,
            s=9, alpha=0.35, color=colours[name], linewidths=0,
        )
        low, high = np.percentile(normalised, [5, 95])
        ax_spread.plot([position, position], [low, high], color=INK, lw=2.4)
        ax_spread.text(
            position + 0.28, high, f"×{high / low:.2f}", va="center", fontsize=12,
            color=colours[name], fontweight="bold",
        )
    ax_spread.axhline(1.0, color=MUTED, lw=0.9, ls="--", zorder=0)
    ax_spread.set_xticks(range(len(responses)))
    ax_spread.set_xticklabels(list(responses))
    ax_spread.set_xlim(-0.6, len(responses) - 0.2)
    ax_spread.set_ylabel("response / median")
    ax_spread.set_title("Spread (5–95 %): only Ixx really moves", loc="left")
    ax_spread.grid(True, axis="y", alpha=0.4)

    out_png = out_dir / "robustness_monte_carlo.png"
    save_figure(fig, out_png, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {out_png}")

    payload = {
        "key": str(args.key),
        "naca": str(args.naca),
        "chord_m": float(args.chord),
        "seed": int(args.seed),
        "n_samples": int(args.n_samples),
        "ok_count": ok_count,
        "fail_count": fail_count,
        "success_rate": success_rate,
        "area_mean": _safe_float(float(np.nanmean(area_vals)) if np.any(np.isfinite(area_vals)) else np.nan),
        "area_std": _safe_float(float(np.nanstd(area_vals)) if np.any(np.isfinite(area_vals)) else np.nan),
        "samples": samples,
        "png": out_png.name,
    }
    out_json = out_dir / "robustness_monte_carlo.json"
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

