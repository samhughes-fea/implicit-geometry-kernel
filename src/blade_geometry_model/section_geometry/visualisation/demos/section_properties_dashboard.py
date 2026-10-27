"""
visualisation/demos/section_properties_dashboard.py
========================================
Section-property analysis dashboard for the hybrid geometry kernel.

For a twin-web (2D-CN) NACA 2412 section, writes:

* ``props_dashboard.png`` — hybrid geometry + area bars + I contribution stacks
  + S→B second-moment transform under twist
* ``props_dashboard.json`` — per-label properties and totals

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/section_properties_dashboard.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

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

from blade_geometry_model.coordinate_frames import (
    CoordinateFrameTransforms,
)
from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.visualisation.plot import (
    plot_section_hybrid,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)

_DPI = DEFAULT_DPI
_TWIST_DEG = 15.0
_LABEL_ORDER = (
    "skin_suction_surface",
    "skin_pressure_surface",
    "spar_cap_suction_surface",
    "spar_cap_pressure_surface",
    "shear_web_0",
    "shear_web_1",
    "sandwich_core_0",
    "sandwich_core_1",
    "sandwich_core_2",
)
_COLOURS = {
    "skin_suction_surface": "#1976D2",
    "skin_pressure_surface": "#64B5F6",
    "spar_cap_suction_surface": "#E53935",
    "spar_cap_pressure_surface": "#EF9A9A",
    "shear_web_0": "#43A047",
    "shear_web_1": "#81C784",
    "sandwich_core_0": "#FFB300",
    "sandwich_core_1": "#FFD54F",
    "sandwich_core_2": "#FFE082",
}


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


def _ordered_labels(props: dict[str, dict[str, float]]) -> list[str]:
    ordered = [lab for lab in _LABEL_ORDER if lab in props]
    for lab in props:
        if lab not in ordered:
            ordered.append(lab)
    return ordered


def _sum_props(props: dict[str, dict[str, float]], labels: list[str]) -> dict[str, float]:
    """Parallel-axis assemble selected labels via the section_properties worker."""
    from blade_geometry_model.section_properties import (
        aggregate_solid_label_properties,
    )

    filtered = {lab: props[lab] for lab in labels if lab in props}
    sp = aggregate_solid_label_properties(
        filtered, jt="skipped", exclude_label_prefix=None
    )
    if sp.A <= 0.0:
        return {
            "area": 0.0,
            "cx": float("nan"),
            "cy": float("nan"),
            "Ixx": float("nan"),
            "Iyy": float("nan"),
            "Ixy": float("nan"),
        }
    return {
        "area": float(sp.A),
        "cx": float(sp.Cx),
        "cy": float(sp.Cy),
        "Ixx": float(sp.Ixx),
        "Iyy": float(sp.Iyy),
        "Ixy": float(sp.Ixy),
    }


def _save_fig(fig, path: Path) -> None:
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "section_properties_dashboard"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    args = parser.parse_args(argv)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", n_points=140, chord=1.0)
    twist = float(np.deg2rad(_TWIST_DEG))
    sec0 = _build(af, twist_rad=0.0)
    sec_t = _build(af, twist_rad=twist)
    spec = quadtree_spec_for_airfoil(af, chord=1.0, max_depth=7, max_cells=120000)

    props0, meta0, mesh0 = compute_section_properties(
        sec0,
        backend="hybrid",
        quadtree_spec=spec,
        airfoil_for_spec=af,
        twist_rad=0.0,
        with_mesh=True,
    )
    props_t, _meta_t, _mesh_t = compute_section_properties(
        sec_t,
        backend="hybrid",
        quadtree_spec=spec,
        airfoil_for_spec=af,
        twist_rad=twist,
        with_mesh=True,
    )
    labels = _ordered_labels(props0)
    totals0 = _sum_props(props0, labels)
    totals_t = _sum_props(props_t, labels)

    # S (twist=0 hybrid) → B via frame transforms; compare to twisted-section totals.
    xform = CoordinateFrameTransforms()
    Ixx_B, Iyy_B, Ixy_B = xform.transform_second_moments_S_to_B(
        totals0["Ixx"], totals0["Iyy"], totals0["Ixy"], _TWIST_DEG
    )

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5))

    # A — hybrid geometry
    plot_section_hybrid(
        sec0,
        residual_mesh=mesh0,
        ax=axes[0, 0],
        title="Hybrid geometry (constructive + residual) — twist 0°",
        show_legend=True,
    )

    # B — area bars
    areas = [float(props0[lab]["area"]) for lab in labels]
    colours = [_COLOURS.get(lab, "#78909C") for lab in labels]
    axes[0, 1].barh(labels[::-1], areas[::-1], color=colours[::-1], edgecolor="k", lw=0.3)
    axes[0, 1].set_xlabel(r"area [m$^2$]")
    axes[0, 1].set_title(f"Per-label area  (ΣA = {totals0['area']:.5f} m²)")
    axes[0, 1].grid(True, axis="x", alpha=0.35)

    # C — I contribution stacks (centroidal about section centroid via parallel axis)
    # Use fraction of |Ixx|, |Iyy| after shifting each label to section centroid.
    cx_t, cy_t = totals0["cx"], totals0["cy"]
    ixx_parts = []
    iyy_parts = []
    for lab in labels:
        p = props0[lab]
        a = float(p["area"])
        if a <= 0.0:
            ixx_parts.append(0.0)
            iyy_parts.append(0.0)
            continue
        dx = float(p["cx"]) - cx_t
        dy = float(p["cy"]) - cy_t
        ixx_parts.append(float(p["Ixx"]) + a * dy * dy)
        iyy_parts.append(float(p["Iyy"]) + a * dx * dx)
    x = np.arange(2)
    bottom = np.zeros(2)
    for lab, ix, iy, col in zip(labels, ixx_parts, iyy_parts, colours):
        vals = np.array([ix, iy])
        axes[1, 0].bar(x, vals, bottom=bottom, color=col, edgecolor="k", lw=0.3, label=lab)
        bottom = bottom + vals
    axes[1, 0].set_xticks(x, [r"$I_{xx}$ (flap)", r"$I_{yy}$ (edge)"])
    axes[1, 0].set_ylabel(r"second moment [m$^4$]")
    axes[1, 0].set_title("Label contributions about section centroid (twist 0°)")
    axes[1, 0].legend(fontsize=6, loc="upper right", ncol=2)
    axes[1, 0].grid(True, axis="y", alpha=0.35)

    # D — S vs B moments
    cats = ["Ixx", "Iyy", "Ixy"]
    s_vals = [totals0["Ixx"], totals0["Iyy"], totals0["Ixy"]]
    b_xform = [float(Ixx_B), float(Iyy_B), float(Ixy_B)]
    b_geom = [totals_t["Ixx"], totals_t["Iyy"], totals_t["Ixy"]]
    xpos = np.arange(len(cats))
    w = 0.25
    axes[1, 1].bar(xpos - w, s_vals, w, label="S (twist 0° hybrid)", color="#1565C0")
    axes[1, 1].bar(
        xpos,
        b_xform,
        w,
        label=f"B via S→B transform ({_TWIST_DEG:g}°)",
        color="#2E7D32",
    )
    axes[1, 1].bar(
        xpos + w,
        b_geom,
        w,
        label=f"B geometry (section twist {_TWIST_DEG:g}°)",
        color="#EF6C00",
    )
    axes[1, 1].set_xticks(xpos, cats)
    axes[1, 1].set_ylabel(r"[m$^4$]")
    axes[1, 1].set_title("Second moments: S vs B (transform vs twisted outline)")
    axes[1, 1].legend(fontsize=7, loc="best")
    axes[1, 1].grid(True, axis="y", alpha=0.35)

    fig.suptitle(
        "Section properties dashboard — NACA 2412 twin-web (hybrid backend)",
        fontsize=13,
    )
    _save_fig(fig, out_dir / "props_dashboard.png")

    per_label: dict[str, Any] = {}
    for lab in labels:
        p = props0[lab]
        per_label[lab] = {
            "area": float(p["area"]),
            "cx": float(p["cx"]),
            "cy": float(p["cy"]),
            "Ixx": float(p["Ixx"]),
            "Iyy": float(p["Iyy"]),
            "Ixy": float(p["Ixy"]),
        }

    payload = {
        "airfoil": "NACA 2412",
        "layout": "twin_web family D",
        "backend": "hybrid",
        "twist_deg_for_B_compare": _TWIST_DEG,
        "per_label_twist_0": per_label,
        "totals_twist_0_S": totals0,
        "totals_twist_15_B_geometry": totals_t,
        "totals_B_from_S_transform": {
            "Ixx": float(Ixx_B),
            "Iyy": float(Iyy_B),
            "Ixy": float(Ixy_B),
        },
        "notes": (
            "Contribution stacks use parallel-axis shift to the section centroid. "
            "S→B moments use CoordinateFrameTransforms with +Y→TE frame sense; "
            "twisted-outline totals are from a section built with twist_angle."
        ),
    }
    json_path = out_dir / "props_dashboard.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(
        f"  sum(A)={totals0['area']:.6f}  Ixx_S={totals0['Ixx']:.3e}  "
        f"Iyy_S={totals0['Iyy']:.3e}  -> Ixx_B(xform)={float(Ixx_B):.3e}"
    )
    return out_dir


if __name__ == "__main__":
    main()
