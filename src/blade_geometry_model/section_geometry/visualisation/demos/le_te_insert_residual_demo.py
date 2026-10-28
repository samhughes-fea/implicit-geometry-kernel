"""
visualisation/demos/le_te_insert_residual_demo.py
======================================
Demonstrate LE/TE inserts as first-class residual-φ consumers.

Builds a twin-web NACA 2412 section with ``leading_edge_insert_x_m`` / ``trailing_edge_insert_x_m``,
then plots constructive wall solids + residual inserts (hybrid / quadtree).

Writes:

* ``inserts_hybrid_quadtree.png`` — full hybrid geometry + union quadtree
* ``inserts_le_te_zooms.png`` — zoomed LE insert and TE insert panels
* ``le_trailing_edge_insert_residual.json`` — label classification + optional areas

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/le_trailing_edge_insert_residual_demo.py

Outputs default to ``visualisation/output/le_trailing_edge_insert_residual/``.
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
    MultiCellSection,
    SDFGrid,
    SectionPropertiesReport,
)
from blade_geometry_model.section_geometry.constructive import (
    classify_property_labels,
    is_residual_label,
)
from blade_geometry_model.section_geometry.visualisation.plot import (
    plot_section_hybrid,
    plot_section_with_quadtree,
)
from blade_geometry_model.section_geometry.sdf import offset

_DPI = DEFAULT_DPI


_SKIN_T = 0.003


def _build_section(af: AirfoilSDF) -> MultiCellSection:
    return MultiCellSection.twin_web(
        af,
        x_fore=0.20,
        x_aft=0.50,
        skin_thickness=_SKIN_T,
        web_thickness=0.004,
        cap_height=(0.014, 0.012),
        web_alignment="chord_normal",
        twist_angle=0.0,
        trailing_edge_insert_x_m=0.75,
        leading_edge_insert_x_m=0.10,
        core_enabled=True,
    )


def _overlay_skin_interfaces(ax, af: AirfoilSDF, skin_t: float, grid: SDFGrid) -> None:
    """Draw inner-skin (cavity) and mid-surface loci for flush checks."""
    phi_inner = grid.eval(offset(af, -skin_t))
    phi_mid = grid.eval(offset(af, -0.5 * skin_t))
    ax.contour(
        grid.X,
        grid.Y,
        phi_inner,
        levels=[0.0],
        colors=["#111111"],
        linewidths=1.4,
        linestyles=["--"],
        zorder=25,
    )
    ax.contour(
        grid.X,
        grid.Y,
        phi_mid,
        levels=[0.0],
        colors=["#9E9E9E"],
        linewidths=0.9,
        linestyles=[":"],
        zorder=24,
    )
    # Legend proxies (contour artists are not labelled cleanly).
    ax.plot([], [], color="#111111", ls="--", lw=1.4, label="inner skin (cavity)")
    ax.plot([], [], color="#9E9E9E", ls=":", lw=0.9, label="skin mid-surface")


def _zoom(ax, x0, x1, y0, y1) -> None:
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)


def _save_fig(fig, path: Path, *, presentation: bool = False) -> None:
    fig.tight_layout()
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "le_trailing_edge_insert_residual"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=default_out,
        help=f"Output directory (default: {default_out})",
    )
    parser.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", chord=1.0)
    section = _build_section(af)
    report = SectionPropertiesReport(
        section,
        backend="hybrid",
        chord=float(af.chord),
        twist_rad=0.0,
        airfoil_for_spec=af,
    )
    mesh = report.mesh
    if mesh is None:
        raise RuntimeError("Hybrid props report did not expose a residual quadtree mesh.")

    labels = list(section.labels)
    constructive, residual = classify_property_labels(labels)

    # --- full hybrid + quadtree ---
    # Right panel draws from residual mesh phi only (leading_edge_insert / trailing_edge_insert).
    fig, _axes = plot_section_with_quadtree(
        mesh,
        labels=residual or list(mesh.phi_by_label.keys()),
        section=section,
        title="Hybrid geometry (constructive wall solids + residual LE/TE inserts)",
        mesh_title="Residual quadtree phi path for LE/TE inserts",
        figsize=(16.0, 5.5),
        alpha=0.5,
    )
    fig.suptitle(
        "LE/TE inserts: residual-phi consumers, flush with inner-skin / core interface",
        fontsize=11,
        y=1.02,
    )
    _save_fig(fig, out_dir / "inserts_hybrid_quadtree.png", presentation=presentation)

    # --- LE / TE zooms ---
    le = np.asarray(af.leading_edge, dtype=float)
    te = np.asarray(af.trailing_edge, dtype=float)
    grid_zoom = SDFGrid.from_airfoil(af, padding=0.05, nx=520, ny=240)
    fig2, axes2 = plt.subplots(1, 2, figsize=(12.5, 5.2))
    for ax, title in (
        (axes2[0], "LE insert zoom (flush with inner skin)"),
        (axes2[1], "TE insert zoom (flush with inner skin)"),
    ):
        plot_section_hybrid(
            section,
            residual_mesh=mesh,
            ax=ax,
            alpha=0.55,
            title=title,
            show_legend=True,
        )
        _overlay_skin_interfaces(ax, af, _SKIN_T, grid_zoom)
    _zoom(axes2[0], le[0] - 0.02, le[0] + 0.14, le[1] - 0.07, le[1] + 0.07)
    _zoom(axes2[1], te[0] - 0.28, te[0] + 0.03, te[1] - 0.06, te[1] + 0.06)
    axes2[0].plot(le[0], le[1], "k*", ms=10)
    axes2[1].plot(te[0], te[1], "k*", ms=10)
    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
    )

    for ax in axes2:
        apply_structural_section_axes(ax)
        ax.legend(fontsize=7, loc="best")
    fig2.suptitle(
        "Inserts use Offset(airfoil, -skin_thickness) — cavity interface, not skin mid-surface",
        fontsize=11,
    )
    _save_fig(fig2, out_dir / "inserts_le_te_zooms.png", presentation=presentation)

    areas: dict[str, float] = {}
    for lab, props in report.to_dict().items():
        if isinstance(props, dict) and "area" in props:
            areas[str(lab)] = float(props["area"])

    residual_areas = {k: v for k, v in areas.items() if is_residual_label(k)}
    constructive_areas = {k: v for k, v in areas.items() if k not in residual_areas}

    payload = {
        "airfoil": "NACA 2412",
        "chord_m": float(af.chord),
        "leading_edge_insert_x_m": 0.10,
        "trailing_edge_insert_x_m": 0.75,
        "labels": [str(x) for x in labels],
        "constructive_labels": constructive,
        "residual_labels": residual,
        "residual_maturity": {
            "leading_edge_insert": "residual phi; geometry less mature than core",
            "trailing_edge_insert": "residual phi; geometry less mature than core",
        },
        "areas_m2": areas,
        "residual_area_sum_m2": float(sum(residual_areas.values())) if residual_areas else None,
        "constructive_area_sum_m2": (
            float(sum(constructive_areas.values())) if constructive_areas else None
        ),
        "skin_thickness_m": _SKIN_T,
        "notes": (
            "Dual-rep model: wall solids/cores use constructive polygons; "
            "leading_edge_insert/trailing_edge_insert stay on the shared adaptive quadtree over phi. "
            "Insert cavity envelope is Offset(airfoil, -skin_thickness), flush with the "
            "inner-skin / core interface (not the skin mid-surface at -t/2)."
        ),
    }
    json_path = out_dir / "le_trailing_edge_insert_residual.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(
        f"  labels={len(labels)}; residual={residual}; "
        f"residual area sum={payload['residual_area_sum_m2']}"
    )
    return out_dir


if __name__ == "__main__":
    main()
