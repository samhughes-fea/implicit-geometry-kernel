"""
visualisation/demos/sdf_metric_quality_demo.py
===================================
Computational-geometry demo: metric SDF quality and CSG-seam ridges.

True SDFs satisfy |∇φ| ≈ 1 almost everywhere (Eikonal equation). Boolean CSG
and offset/slab recipes produce *signed-distance-like* fields with ridges
where |∇φ| collapses. This demo makes that visible on:

* airfoil polyline SDF (near-metric)
* outer-skin half-annulus slab (SS) after Intersect/camber clip
* multi-cell CSG union of wall solids

Writes:

* ``metric_quality_fields.png`` — φ and |∇φ| panels for the three fields
* ``eikonal_residual_hist.png`` — residual histograms + ridge scatter
* ``sdf_metric_quality.json`` — residual percentiles

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/sdf_metric_quality_demo.py
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

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    demo_output_dir,
    save_figure,
)

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
)
from blade_geometry_model.section_geometry.visualisation.plot_sdf import (
    plot_sdf_field,
)
from blade_geometry_model.section_geometry.sdf.csg import union_all

_DPI = DEFAULT_DPI
_BAND_M = 0.02  # exclude a thin band around φ=0 from residual stats
_RIDGE_GRAD = 0.55


def _build_section(af: AirfoilSDF) -> MultiCellSection:
    return MultiCellSection.twin_web(
        af,
        x_fore=0.20,
        x_aft=0.50,
        skin_thickness=0.005,
        web_thickness=0.006,
        cap_height=(0.014, 0.012),
        web_alignment="chord_normal",
        twist_angle=0.0,
        structural_family="D",
        core_enabled=False,
    )


def _residual_stats(gm: np.ndarray, phi: np.ndarray, *, band_m: float) -> dict[str, float]:
    mask = np.abs(phi) >= float(band_m)
    resid = np.abs(gm[mask] - 1.0)
    if resid.size == 0:
        return {
            "n": 0,
            "p50": float("nan"),
            "p90": float("nan"),
            "p99": float("nan"),
            "mean": float("nan"),
            "max": float("nan"),
        }
    return {
        "n": int(resid.size),
        "p50": float(np.percentile(resid, 50)),
        "p90": float(np.percentile(resid, 90)),
        "p99": float(np.percentile(resid, 99)),
        "mean": float(np.mean(resid)),
        "max": float(np.max(resid)),
    }


def _ridge_xy(grid: SDFGrid, gm: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Sample points where |∇φ| collapses near material (CSG / offset ridges)."""
    mask = (gm < _RIDGE_GRAD) & (phi < 0.02) & (phi > -0.08)
    ys, xs = np.nonzero(mask)
    if xs.size == 0:
        return np.zeros((0, 2), dtype=float)
    # Subsample for scatter clarity.
    step = max(1, xs.size // 2500)
    return np.column_stack(
        [grid.X[ys[::step], xs[::step]], grid.Y[ys[::step], xs[::step]]]
    )


def _save_fig(fig, path: Path) -> None:
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "sdf_metric_quality"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    args = parser.parse_args(argv)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", chord=1.0, n_points=160)
    section = _build_section(af)
    grid = SDFGrid.from_airfoil(af, padding=0.08, nx=500, ny=220)

    skin_ss = section["skin_suction_surface"]
    wall_solids = [
        section[lab]
        for lab in section.labels
        if str(getattr(lab, "value", lab)).startswith(
            ("outer_skin", "web_", "spar_cap_")
        )
    ]
    union_fn = union_all(wall_solids)

    # Carry the stats key alongside each field. Deriving it from the title by
    # substring match silently mislabelled the middle field and left its entry
    # missing, which crashed the residual-histogram figure below.
    fields = [
        ("airfoil SDF", af, 0.08, "airfoil"),
        ("skin_suction_surface (slab ∩ camber)", skin_ss, 0.04, "skin_suction_surface"),
        ("wall-solid CSG union", union_fn, 0.04, "wall_solid_union"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 7.5))
    stats: dict[str, dict] = {}
    ridge_store: dict[str, np.ndarray] = {}

    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
        rotated_naca_grid_to_section_plot_axes,
    )

    # The top row goes through plot_sdf_field, which maps the grid into the B plot
    # frame. The bottom row draws the same grid, so it has to be mapped too —
    # otherwise apply_structural_section_axes labels raw chord coordinates as Y_B
    # and the two rows disagree about where the section is.
    grid_h, grid_v = rotated_naca_grid_to_section_plot_axes(
        grid.X, grid.Y, twist_rad=0.0, chord_m=float(af.chord)
    )

    for col, (title, sdf_fn, vmax, key) in enumerate(fields):
        phi = grid.eval(sdf_fn)
        gm = grid.grad_magnitude(phi)
        plot_sdf_field(
            phi,
            grid,
            ax=axes[0, col],
            title=f"φ — {title}",
            vmax=vmax,
            colorbar=(col == 2),
        )
        im = axes[1, col].pcolormesh(
            grid_h,
            grid_v,
            gm,
            cmap="viridis",
            vmin=0.0,
            vmax=1.6,
            shading="auto",
        )
        axes[1, col].contour(
            grid_h, grid_v, phi, levels=[0.0], colors="w", linewidths=1.0
        )
        ridges = _ridge_xy(grid, gm, phi)
        ridge_store[title] = ridges
        if ridges.size:
            ridge_h, ridge_v = rotated_naca_grid_to_section_plot_axes(
                ridges[:, 0], ridges[:, 1], twist_rad=0.0, chord_m=float(af.chord)
            )
            axes[1, col].scatter(
                ridge_h,
                ridge_v,
                s=2,
                c="#FF5252",
                alpha=0.55,
                zorder=5,
                label=f"|∇φ|<{_RIDGE_GRAD:g} ridges",
            )
            axes[1, col].legend(fontsize=7, loc="upper right")
        axes[1, col].set_title(f"|∇φ| — {title}")
        apply_structural_section_axes(axes[1, col])
        if col == 2:
            fig.colorbar(im, ax=axes[1, col], label="|∇φ|")
        stats[key] = {
            "title": title,
            "eikonal_residual_excl_band_m": _BAND_M,
            **_residual_stats(gm, phi, band_m=_BAND_M),
            "n_ridge_samples": int(ridges.shape[0]),
        }

    fig.suptitle(
        "Metric SDF quality — Eikonal |∇φ|≈1; red points mark CSG/offset ridges",
        fontsize=12,
    )
    _save_fig(fig, out_dir / "metric_quality_fields.png")

    fig2, axes2 = plt.subplots(1, 3, figsize=(14.0, 4.0))
    for ax, (title, sdf_fn, _vmax, key) in zip(axes2, fields):
        phi = grid.eval(sdf_fn)
        gm = grid.grad_magnitude(phi)
        mask = np.abs(phi) >= _BAND_M
        resid = np.abs(gm[mask] - 1.0)
        ax.hist(resid, bins=60, range=(0.0, 1.5), color="#455A64", alpha=0.9)
        st = stats[key]
        ax.axvline(st["p50"], color="#FFC107", lw=1.5, label=f"p50={st['p50']:.3f}")
        ax.axvline(st["p90"], color="#FF7043", lw=1.5, label=f"p90={st['p90']:.3f}")
        ax.axvline(st["p99"], color="#E53935", lw=1.5, label=f"p99={st['p99']:.3f}")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(r"$||\nabla\phi|-1|$")
        ax.set_ylabel("count")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    fig2.suptitle(
        f"Eikonal residual histograms (excluding |φ| < {_BAND_M:g} m band)",
        fontsize=11,
    )
    _save_fig(fig2, out_dir / "eikonal_residual_hist.png")

    payload = {
        "airfoil": "NACA 2412",
        "chord_m": float(af.chord),
        "grid": {"nx": grid.nx, "ny": grid.ny, "dx": grid.dx, "dy": grid.dy},
        "ridge_grad_threshold": _RIDGE_GRAD,
        "fields": stats,
        "notes": (
            "Airfoil polyline SDF stays near-metric; slab/CSG unions elevate "
            "Eikonal residual and create |∇φ| ridges at Boolean seams."
        ),
    }
    json_path = out_dir / "sdf_metric_quality.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    for k, st in stats.items():
        print(
            f"  {k}: residual p50={st['p50']:.4f} p90={st['p90']:.4f} "
            f"p99={st['p99']:.4f} ridges={st['n_ridge_samples']}"
        )
    return out_dir


if __name__ == "__main__":
    main()
