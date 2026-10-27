"""
visualisation/demos/thin_feature_robustness_demo.py
========================================
Computational-geometry demo: knife-edge web capsules under twist.

Shear webs are capsule ∩ cavity solids. Adaptive quadtree refinement must
resolve them where a coarse uniform grid would miss area. This demo shows:

* hybrid geometry + QT leaf overlay zoomed on ``web_0`` at 0° and 15° twist
* ``web_0`` area vs quadtree depth (and constructive shoelace reference)
* area Δ (QT − constructive) as a stability metric

Writes:

* ``thin_web_zoom_twist.png`` — 0° / 15° hybrid + QT zooms
* ``web_area_vs_depth.png`` — depth sweep
* ``thin_feature_robustness.json`` — areas and deltas

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/thin_feature_robustness_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

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
    apply_presentation_style,
    save_figure,
)
from blade_geometry_model.section_geometry.visualisation.palette import (
    strip_axes_chrome,
)

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
)
from blade_geometry_model.section_geometry.props_emit.area_mesh import (
    section_properties_via_quadtree,
)
from blade_geometry_model.section_geometry.constructive import (
    emit_wall_solid_polygons,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.visualisation.plot import (
    plot_section_hybrid,
)
from blade_geometry_model.section_geometry.visualisation.plot_section import (
    _naca_xy_to_plot_hv,
    apply_structural_section_axes,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.area_2d.quadrature import label_boundary_segments

_DPI = DEFAULT_DPI
# Floor is MIN_REALISTIC_WEB_WALL_THICKNESS_M (4 mm); request that minimum.
_WEB_T = 0.004
_DEPTHS = (4, 5, 6, 7, 8)


def _build_section(af: AirfoilSDF, *, twist_rad: float) -> MultiCellSection:
    return MultiCellSection.twin_web(
        af,
        x_fore=0.22,
        x_aft=0.48,
        skin_thickness=0.004,
        web_thickness=_WEB_T,
        cap_height=(0.012, 0.010),
        web_alignment="chord_normal",
        twist_angle=float(twist_rad),
        structural_family="D",
        core_enabled=False,
    )


def _constructive_web_area(section: MultiCellSection) -> float:
    polys = emit_wall_solid_polygons(section)
    geom = polys.get("shear_web_0")
    if geom is None or geom.is_empty:
        return 0.0
    return float(geom.area)


def _web_centroid(section: MultiCellSection) -> np.ndarray:
    polys = emit_wall_solid_polygons(section)
    geom = polys["shear_web_0"]
    c = geom.centroid
    return np.array([float(c.x), float(c.y)], dtype=float)


def _save_fig(fig, path: Path, *, presentation: bool = False) -> None:
    fig.tight_layout()
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def _plot_zoom_panel(
    ax,
    section,
    mesh,
    centre_naca: np.ndarray,
    half: float,
    title: str,
    *,
    twist_rad: float,
    chord_m: float,
) -> None:
    """Hybrid fill + QT overlay in canonical B (pitch at origin, LE on RHS)."""
    twist = float(twist_rad)
    chord = float(chord_m)
    plot_section_hybrid(
        section,
        residual_mesh=mesh,
        ax=ax,
        alpha=0.5,
        title=title,
        show_legend=False,
        structural_axes=True,
        twist_rad=twist,
        chord_m=chord,
    )
    centre_hv = _naca_xy_to_plot_hv(
        np.asarray(centre_naca, dtype=float).reshape(1, 2),
        twist_rad=twist,
        chord_m=chord,
    )[0]
    for cell in mesh.leaf_cells():
        corners_naca = np.array(
            [
                [cell.x_min, cell.y_min],
                [cell.x_max, cell.y_min],
                [cell.x_max, cell.y_max],
                [cell.x_min, cell.y_max],
            ],
            dtype=float,
        )
        corners_hv = _naca_xy_to_plot_hv(
            corners_naca, twist_rad=twist, chord_m=chord
        )
        # AABB filter in plot space around the web centroid.
        if (
            float(corners_hv[:, 0].max()) < centre_hv[0] - half
            or float(corners_hv[:, 0].min()) > centre_hv[0] + half
            or float(corners_hv[:, 1].max()) < centre_hv[1] - half
            or float(corners_hv[:, 1].min()) > centre_hv[1] + half
        ):
            continue
        ax.add_patch(
            Polygon(
                corners_hv,
                closed=True,
                fill=False,
                edgecolor="0.35",
                linewidth=0.35,
                alpha=0.7,
                zorder=6,
            )
        )
    segs = label_boundary_segments(mesh, "shear_web_0")
    if len(segs):
        segs_hv = [
            _naca_xy_to_plot_hv(np.asarray(s, dtype=float), twist_rad=twist, chord_m=chord)
            for s in segs
        ]
        ax.add_collection(
            LineCollection(segs_hv, colors=["#1B5E20"], linewidths=1.6, zorder=8)
        )
    x0, x1 = float(centre_hv[0] - half), float(centre_hv[0] + half)
    y0, y1 = float(centre_hv[1] - half), float(centre_hv[1] + half)
    apply_structural_section_axes(ax)
    if ax.xaxis_inverted():
        ax.set_xlim(x1, x0)
    else:
        ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal", adjustable="box")
    strip_axes_chrome(ax, keep_frame=True)


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "thin_feature_robustness"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    parser.add_argument("--presentation", action="store_true", help="Deck styling.")
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", chord=1.0, n_points=140)
    cases = [("twist_0deg", 0.0), ("twist_15deg", float(np.deg2rad(15.0)))]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 5.2))
    zoom_meta: dict[str, dict] = {}
    for ax, (tag, twist) in zip(axes, cases):
        section = _build_section(af, twist_rad=twist)
        # QT domain must cover rotated-NACA solids (same as convergence demos).
        af_for_spec = af.rotate(float(twist)) if abs(float(twist)) > 1e-12 else af
        spec = quadtree_spec_for_airfoil(
            af_for_spec, chord=1.0, max_depth=7, max_cells=120000
        )
        props, _meta, mesh = compute_section_properties(
            section,
            backend="quadtree",
            quadtree_spec=spec,
            airfoil_for_spec=af_for_spec,
            twist_rad=float(twist),
            with_mesh=True,
        )
        assert mesh is not None
        centre = _web_centroid(section)
        a_qt = float(props["shear_web_0"]["area"])
        a_c = _constructive_web_area(section)
        _plot_zoom_panel(
            ax,
            section,
            mesh,
            centre,
            half=0.06,
            title=(
                f"{np.rad2deg(twist):.0f}° twist · {_WEB_T * 1e3:g} mm web\n"
                f"quadtree area {abs(a_qt - a_c) / a_c:.2%} from exact"
            ),
            twist_rad=float(twist),
            chord_m=1.0,
        )
        centre_hv = _naca_xy_to_plot_hv(
            centre.reshape(1, 2), twist_rad=float(twist), chord_m=1.0
        )[0]
        zoom_meta[tag] = {
            "twist_rad": float(twist),
            "web_0_area_quadtree_m2": a_qt,
            "web_0_area_constructive_m2": a_c,
            "area_delta_qt_minus_constr_m2": a_qt - a_c,
            "n_leaves": int(len(mesh.leaf_cells())),
            "centre_rotated_naca": [float(centre[0]), float(centre[1])],
            "centre_plot_YB_XB": [float(centre_hv[0]), float(centre_hv[1])],
        }
    _save_fig(fig, out_dir / "thin_web_zoom_twist.png", presentation=presentation)

    # Depth sweep at 0° twist — clean QT ↔ constructive convergence on the
    # same wall solid. Zoom panels above still show the web under 15° twist.
    section = _build_section(af, twist_rad=0.0)
    a_c = _constructive_web_area(section)
    depth_rows: list[dict] = []
    for depth in _DEPTHS:
        spec = quadtree_spec_for_airfoil(
            af, chord=1.0, max_depth=int(depth), max_cells=200000
        )
        props, mesh, qt_stats = section_properties_via_quadtree(
            section, spec, twist_rad=0.0
        )
        a_qt = float(props["shear_web_0"]["area"])
        n_leaves = int(qt_stats.get("leaf_count", len(mesh.leaf_cells())))
        depth_rows.append(
            {
                "max_depth": int(depth),
                "web_0_area_quadtree_m2": a_qt,
                "web_0_area_constructive_m2": a_c,
                "area_delta_m2": a_qt - a_c,
                "rel_err_vs_constructive": abs(a_qt - a_c) / max(a_c, 1e-30),
                "n_leaves": n_leaves,
            }
        )

    fig2, ax2 = plt.subplots(figsize=(7.5, 4.5))
    depths = [r["max_depth"] for r in depth_rows]
    areas = [r["web_0_area_quadtree_m2"] for r in depth_rows]
    ax2.plot(depths, areas, "o-", color="#2E7D32", lw=2, label="quadtree web_0 area")
    ax2.axhline(a_c, color="#1565C0", ls="--", lw=1.5, label="constructive shoelace")
    ax2.set_xlabel("quadtree max_depth")
    ax2.set_ylabel(r"web_0 area [m$^2$]")
    ax2.set_title("Thin web area vs adaptive depth (twist = 0 deg)")
    ax2.grid(True, alpha=0.35)
    ax2.legend(fontsize=9)
    _save_fig(fig2, out_dir / "web_area_vs_depth.png")

    payload = {
        "airfoil": "NACA 2412",
        "chord_m": 1.0,
        "web_thickness_m": _WEB_T,
        "zoom_cases": zoom_meta,
        "depth_sweep_twist_0deg": depth_rows,
        "notes": (
            "Minimum-thickness webs remain finite under twist (zoom panels); "
            "at zero twist, QT area approaches constructive shoelace with depth."
        ),
    }
    json_path = out_dir / "thin_feature_robustness.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    for r in depth_rows:
        print(
            f"  depth={r['max_depth']}: A_qt={r['web_0_area_quadtree_m2']:.6e} "
            f"rel_err={r['rel_err_vs_constructive']:.3e}"
        )
    return out_dir


if __name__ == "__main__":
    main()
