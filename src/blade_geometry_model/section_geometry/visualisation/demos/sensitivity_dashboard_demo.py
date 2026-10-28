"""
visualisation/demos/sensitivity_dashboard_demo.py
=====================================
Finite-difference sensitivity dashboard for implicit section geometry.

Writes under ``visualisation/output/sensitivity_dashboard/``:

* ``sensitivity_dashboard.png`` - derivative bars, local sweep responses, and base geometry.
* ``sensitivity_dashboard.json`` - derivatives and sweep data.
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
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SDFGrid,
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.visualisation.plot import (
    plot_section,
)
from blade_geometry_model.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    NAVY,
    RED,
    TEAL,
)

_DPI = DEFAULT_DPI


def _totals(props: dict[str, dict[str, float]]) -> dict[str, float]:
    return {
        "area": float(np.nansum([float(v.get("area", np.nan)) for v in props.values()])),
        "Ixx": float(np.nansum([float(v.get("Ixx", np.nan)) for v in props.values()])),
        "Iyy": float(np.nansum([float(v.get("Iyy", np.nan)) for v in props.values()])),
    }


def _evaluate(layout, airfoil, spec, chord_m: float, pars: dict[str, float]) -> dict[str, float]:
    section = build_section_view(
        airfoil,
        layout,
        twist_rad=float(np.deg2rad(pars["twist_deg"])),
        skin_thickness=float(pars["skin_thickness"]),
        spar_cap_thickness=float(pars["spar_cap_thickness"]),
        shear_web_thickness=float(pars["shear_web_thickness"]),
    )
    props, _meta = compute_section_properties(
        section,
        backend="hybrid",
        quadtree_spec=spec,
        chord=float(chord_m),
        twist_rad=float(np.deg2rad(pars["twist_deg"])),
        airfoil_for_spec=airfoil,
    )
    return _totals(props)


#: Slide colours, one per design variable, in the deck's palette.
_PARAMETER_COLOURS = {
    "skin_thickness": NAVY,
    "spar_cap_thickness": RED,
    "shear_web_thickness": TEAL,
    "twist_deg": "#9E86C8",
}


def _plot_sensitivity_slide(
    base: dict[str, float],
    sweep: dict[str, list[dict[str, float]]],
    deriv: dict[str, dict[str, float]],
    y0: dict[str, float],
    path: Path,
    *,
    presentation: bool,
) -> None:
    """The slide cut of the dashboard: normalised sweeps and elasticities.

    Both panels are deliberately dimensionless. The dashboard plots every
    parameter against one shared axis in its own units, which puts three
    millimetre-scale gauges on top of each other at x≈0 and a 15-degree twist
    fifteen units away; and it compares raw derivatives whose units are m²/m
    against m²/deg. Neither is readable as a comparison. Normalising the sweep
    by each parameter's own base value, and the derivative into an elasticity,
    puts all four on one scale where they can actually be ranked.
    """
    figure, (axis_sweep, axis_bars) = plt.subplots(1, 2, figsize=(14.5, 5.4))

    for name, rows in sweep.items():
        base_value = float(base[name])
        if base_value == 0.0:
            continue
        parameter_change = [100.0 * (float(r["x"]) / base_value - 1.0) for r in rows]
        area_change = [100.0 * (float(r["area"]) / float(y0["area"]) - 1.0) for r in rows]
        axis_sweep.plot(
            parameter_change,
            area_change,
            marker="o",
            markersize=5.5,
            linewidth=2.0,
            color=_PARAMETER_COLOURS.get(name, INK),
            label=name.replace("_", " "),
        )

    axis_sweep.axhline(0.0, color=MUTED, linewidth=0.8, linestyle=":")
    axis_sweep.axvline(0.0, color=MUTED, linewidth=0.8, linestyle=":")
    axis_sweep.set_xlabel("change in design variable [% of base]")
    axis_sweep.set_ylabel("change in section area [%]")
    axis_sweep.set_title("Local response: straight and monotone in every variable")
    axis_sweep.grid(True, alpha=0.25)
    axis_sweep.legend(fontsize=9)

    # Elasticity: dimensionless, so gauges in metres and twist in degrees can be
    # ranked against each other. This is the number the optimiser experiences.
    names = list(sweep.keys())
    positions = np.arange(len(names))
    responses = ("area", "Ixx", "Iyy")
    bar_width = 0.26
    for offset_index, response in enumerate(responses):
        elasticity = [
            deriv[response][name] * float(base[name]) / float(y0[response])
            if float(y0[response]) != 0.0
            else 0.0
            for name in names
        ]
        axis_bars.bar(
            positions + (offset_index - 1) * bar_width,
            elasticity,
            width=bar_width,
            label=f"d ln {response} / d ln x",
            color=(NAVY, RED, TEAL)[offset_index],
        )

    axis_bars.axhline(0.0, color=INK, linewidth=0.9)
    axis_bars.set_xticks(positions)
    axis_bars.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=9)
    axis_bars.set_ylabel("elasticity (dimensionless)")
    axis_bars.set_title("What the optimiser sees: response elasticities")
    axis_bars.grid(True, axis="y", alpha=0.25)
    axis_bars.legend(fontsize=9)

    figure.suptitle(
        "Section response is smooth in the design variables — "
        "nothing discontinuous for a gradient step to trip on",
        fontsize=13,
    )
    figure.tight_layout()
    save_figure(figure, path, presentation=presentation)
    plt.close(figure)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--key", default="2D-F")
    p.add_argument("--naca", default="2412")
    p.add_argument("--chord", type=float, default=1.0)
    p.add_argument("--twist-deg", type=float, default=15.0)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )
    args = p.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = args.out_dir or (Path(__file__).resolve().parents[1] / "output" / "sensitivity_dashboard")
    out_dir.mkdir(parents=True, exist_ok=True)

    layout = resolve_system_type(str(args.key))
    airfoil = AirfoilSDF.from_naca(str(args.naca), chord=float(args.chord), n_points=170)
    spec = quadtree_spec_for_airfoil(airfoil, chord=float(args.chord), max_depth=7, max_cells=140000)

    base = {
        "skin_thickness": 0.0035,
        "spar_cap_thickness": 0.0120,
        "shear_web_thickness": 0.0040,
        "twist_deg": float(args.twist_deg),
    }
    deltas = {
        "skin_thickness": 2.5e-4,
        "spar_cap_thickness": 8.0e-4,
        "shear_web_thickness": 3.5e-4,
        "twist_deg": 0.5,
    }

    y0 = _evaluate(layout, airfoil, spec, float(args.chord), base)

    deriv: dict[str, dict[str, float]] = {"area": {}, "Ixx": {}, "Iyy": {}}
    for name, h in deltas.items():
        xp = dict(base)
        xm = dict(base)
        xp[name] = float(xp[name] + h)
        xm[name] = float(xm[name] - h)
        yp = _evaluate(layout, airfoil, spec, float(args.chord), xp)
        ym = _evaluate(layout, airfoil, spec, float(args.chord), xm)
        for out_name in deriv:
            deriv[out_name][name] = float((yp[out_name] - ym[out_name]) / (2.0 * h))

    sweep: dict[str, list[dict[str, float]]] = {}
    for name, h in deltas.items():
        rows: list[dict[str, float]] = []
        for t in [-2.0, -1.0, 0.0, 1.0, 2.0]:
            x = dict(base)
            x[name] = float(x[name] + t * h)
            y = _evaluate(layout, airfoil, spec, float(args.chord), x)
            rows.append({"x": x[name], "area": y["area"], "Ixx": y["Ixx"], "Iyy": y["Iyy"]})
        sweep[name] = rows

    section_base = build_section_view(
        airfoil,
        layout,
        twist_rad=float(np.deg2rad(base["twist_deg"])),
        skin_thickness=base["skin_thickness"],
        spar_cap_thickness=base["spar_cap_thickness"],
        shear_web_thickness=base["shear_web_thickness"],
    )
    airfoil_b = airfoil.rotate(float(np.deg2rad(base["twist_deg"])))
    grid = SDFGrid.from_airfoil(airfoil_b, nx=520, ny=220)

    fig, axs = plt.subplots(2, 2, figsize=(16, 10.5))
    ax0, ax1, ax2, ax3 = axs.ravel()

    names = list(deltas.keys())
    xloc = np.arange(len(names))
    ax0.bar(xloc - 0.25, [deriv["area"][k] for k in names], width=0.25, label="dA/dx")
    ax0.bar(xloc, [deriv["Ixx"][k] for k in names], width=0.25, label="dIxx/dx")
    ax0.bar(xloc + 0.25, [deriv["Iyy"][k] for k in names], width=0.25, label="dIyy/dx")
    ax0.set_xticks(xloc)
    ax0.set_xticklabels(names, rotation=20, ha="right")
    ax0.set_title("Finite-difference sensitivities")
    ax0.grid(True, axis="y", alpha=0.25)
    ax0.legend(loc="best")

    for name in names:
        rows = sweep[name]
        xs = [r["x"] for r in rows]
        ys = [r["area"] for r in rows]
        ax1.plot(xs, ys, marker="o", lw=1.2, label=name)
    ax1.set_title("Local sweep: area response")
    ax1.set_xlabel("parameter value")
    ax1.set_ylabel("area [m²]")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best", fontsize=8)

    for name in names:
        rows = sweep[name]
        xs = [r["x"] for r in rows]
        yx = [r["Ixx"] for r in rows]
        yy = [r["Iyy"] for r in rows]
        ax2.plot(xs, yx, marker="o", lw=1.0, alpha=0.9, label=f"{name}: Ixx")
        ax2.plot(xs, yy, marker="s", lw=1.0, alpha=0.9, linestyle="--", label=f"{name}: Iyy")
    ax2.set_title("Local sweep: inertia responses")
    ax2.set_xlabel("parameter value")
    ax2.set_ylabel("inertia [m⁴]")
    ax2.grid(True, alpha=0.25)

    plot_section(
        section_base,
        grid,
        ax=ax3,
        show_airfoil=False,
        show_legend=False,
        auto_limits=True,
        alpha=0.52,
        title=f"Base geometry - {args.key} ({base['twist_deg']:.1f} deg)",
    )

    fig.suptitle("Sensitivity dashboard (implicit CSG section geometry)", fontsize=14)
    fig.tight_layout()
    out_png = out_dir / "sensitivity_dashboard.png"
    save_figure(fig, out_png, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {out_png}")

    _plot_sensitivity_slide(
        base,
        sweep,
        deriv,
        y0,
        out_dir / "sensitivity_slide.png",
        presentation=bool(args.presentation),
    )

    payload = {
        "key": str(args.key),
        "naca": str(args.naca),
        "chord_m": float(args.chord),
        "base": base,
        "base_outputs": y0,
        "deltas": deltas,
        "derivatives": deriv,
        "sweeps": sweep,
        "png": out_png.name,
    }
    out_json = out_dir / "sensitivity_dashboard.json"
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
