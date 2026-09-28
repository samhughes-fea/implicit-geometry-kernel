"""
visualisation/demos/le_te_exactness_demo.py
================================
Demonstrate LE/TE geometric exactness of the implicit / constructive kernel.

Shows:

1. TE tip — dense-grid ``plot_section`` vs constructive ``plot_section_hybrid``
   (hybrid polygons are not marching-squares limited).
2. Skin mid-surface ending short of the aero TE (zero-area spur avoidance).
3. LE nose — polyhedral SDF fidelity (|φ| ≈ 0 on airfoil vertices).
4. Optional 15° twist case with B-frame LE/TE markers.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/le_te_exactness_demo.py

Outputs default to ``visualisation/output/le_te_exactness/``.
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
    NAVY,
    RED,
    TEAL,
)

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
    build_section_geometry,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot import (
    plot_section,
    plot_section_hybrid,
    plot_sdf_field,
)
from blade_structure_utilities.geometry_utilities.coordinate_frames import (
    naca_to_twisted_pitch_axis_blade_frame,
    section_plot_axes,
)

_DPI = DEFAULT_DPI


def _build_section(af: AirfoilSDF, *, twist_rad: float = 0.0) -> MultiCellSection:
    return MultiCellSection.twin_web(
        af,
        x_fore=0.20,
        x_aft=0.50,
        skin_thickness=0.005,
        web_thickness=0.006,
        cap_height=(0.014, 0.012),
        web_alignment="chord_normal",
        twist_angle=float(twist_rad),
        structural_family="D",
        core_enabled=False,
    )


def _zoom(ax, x0, x1, y0, y1) -> None:
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)


def _phi_on_vertices(af: AirfoilSDF, *, x_max: float | None = None, x_min: float | None = None) -> dict:
    verts = np.asarray(af.vertices, dtype=float)
    mask = np.ones(len(verts), dtype=bool)
    if x_max is not None:
        mask &= verts[:, 0] <= float(x_max)
    if x_min is not None:
        mask &= verts[:, 0] >= float(x_min)
    sel = verts[mask]
    if sel.size == 0:
        return {"n_vertices": 0, "max_abs_phi_m": None, "mean_abs_phi_m": None}
    phi = af(sel[:, 0], sel[:, 1])
    abs_phi = np.abs(np.asarray(phi, dtype=float))
    return {
        "n_vertices": int(sel.shape[0]),
        "max_abs_phi_m": float(np.max(abs_phi)),
        "mean_abs_phi_m": float(np.mean(abs_phi)),
    }


def _skin_midline_te_gap_from_mid(mid: np.ndarray, af: AirfoilSDF) -> dict:
    te = np.asarray(af.trailing_edge, dtype=float)
    mid = np.asarray(mid, dtype=float)
    # Midline is closed; points nearest TE are the open TE gap ends
    d = np.linalg.norm(mid - te.reshape(1, 2), axis=1)
    i_near = int(np.argmin(d))
    gap = float(d[i_near])
    # Also distance from max-x midline points (TE-side)
    te_side = mid[mid[:, 0] >= float(np.percentile(mid[:, 0], 95))]
    gap_te_side = float(np.min(np.linalg.norm(te_side - te.reshape(1, 2), axis=1)))
    return {
        "aero_te_s": [float(te[0]), float(te[1])],
        "nearest_midline_point_s": [float(mid[i_near, 0]), float(mid[i_near, 1])],
        "gap_nearest_m": gap,
        "gap_te_side_m": gap_te_side,
        "midline_max_x_m": float(np.max(mid[:, 0])),
        "aero_te_x_m": float(te[0]),
    }


def _storage_to_plot(pts: np.ndarray, chord: float) -> np.ndarray:
    """Storage (x_chord, y_thickness) -> the untwisted plot frame used by plot_section.

    ``plot_section`` draws in the canonical B frame (horizontal = Y_B, measured from
    the pitch axis), so a zoom window expressed in storage coordinates lands nowhere
    near the section.
    """
    pts = np.asarray(pts, dtype=float).reshape(-1, 2)
    b = naca_to_twisted_pitch_axis_blade_frame(pts, 0.0, float(chord))
    h, v = section_plot_axes(b)
    return np.column_stack([np.asarray(h, dtype=float), np.asarray(v, dtype=float)])


def _sampled_nose_deviation(af: AirfoilSDF, nx: int, *, window_x: float = 0.05, padding: float = 0.02):
    """Contour the airfoil field on an nx-wide square-celled grid over the nose.

    Returns the worst deviation of that sampled boundary from the true one. Because
    ``af`` is itself the exact signed distance, ``|af(p)|`` at a sampled contour
    point *is* that point's distance to the real boundary — no nearest-point search
    and no discretisation of the reference.
    """
    verts = np.asarray(af.vertices, dtype=float)
    x_min, x_max = verts[:, 0].min() - padding, verts[:, 0].max() + padding
    y_min, y_max = verts[:, 1].min() - padding, verts[:, 1].max() + padding
    cell = (x_max - x_min) / float(nx - 1)
    ny = max(4, int(round((y_max - y_min) / cell)) + 1)

    grid = SDFGrid.from_bbox(x_min, x_max, y_min, y_max, nx=int(nx), ny=int(ny))
    segments = grid.zero_contour(grid.eval(af))

    deviation = 0.0
    for seg in segments:
        in_window = seg[:, 0] <= float(window_x)
        if in_window.any():
            phi_here = np.abs(np.asarray(af(seg[in_window, 0], seg[in_window, 1]), dtype=float))
            deviation = max(deviation, float(np.max(phi_here)))
    return deviation, float(cell), segments


def _save_fig(fig, path: Path, *, presentation: bool = False) -> None:
    fig.tight_layout()
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "le_te_exactness"
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
        help="Deck styling for the slide-facing nose figure.",
    )
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", chord=1.0, n_points=240)
    section = _build_section(af, twist_rad=0.0)
    # Intentionally moderate grid — TE tip should look mesh-limited vs hybrid
    grid = SDFGrid.from_airfoil(af, padding=0.04, nx=280, ny=120)

    # --- 1. TE tip: grid vs hybrid ---
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    plot_section(
        section,
        grid,
        ax=axes[0],
        alpha=0.55,
        title="Dense-grid fill (mesh-limited TE tip)",
        show_legend=False,
        omit_cores=True,
        auto_limits=False,
    )
    plot_section_hybrid(
        section,
        residual_mesh=None,
        ax=axes[1],
        alpha=0.55,
        title="Constructive hybrid polygons\n(sharp TE — not marching-squares limited)",
        show_legend=False,
    )
    te = np.asarray(af.trailing_edge, dtype=float)
    from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
    )

    # plot_section draws in the plot frame, so the zoom and the marker have to be
    # taken there too — in storage coordinates this window misses the section.
    te_plot = _storage_to_plot(te, af.chord)[0]
    for ax in axes:
        _zoom(ax, te_plot[0] - 0.02, te_plot[0] + 0.08, te_plot[1] - 0.04, te_plot[1] + 0.04)
        ax.plot(te_plot[0], te_plot[1], "k*", ms=10, label="aero TE")
        ax.legend(fontsize=8, loc="upper left")
        apply_structural_section_axes(ax)
    fig.suptitle("TE tip exactness — NACA 2412 twin-web (family D)", fontsize=12)
    _save_fig(fig, out_dir / "te_tip_grid_vs_hybrid.png")

    # --- 2. Skin mid-surface vs aero tip ---
    # Full closed mid-ring for TE-gap demo (union of SS/PS open half-rings).
    from blade_structure_utilities.geometry_utilities.section_geometry.recipes.sdf_ir.outer_skin import (
        OuterSkin,
    )

    skin_full = OuterSkin(af, 0.005)
    mid = np.asarray(skin_full.midline_polyline(), dtype=float)
    te_gap = _skin_midline_te_gap_from_mid(mid, af)
    fig2, ax2 = plt.subplots(figsize=(8.5, 5.5))
    verts = np.asarray(af.vertices, dtype=float)
    ax2.plot(
        np.r_[verts[:, 0], verts[0, 0]],
        np.r_[verts[:, 1], verts[0, 1]],
        "k-",
        lw=1.0,
        label="aero contour",
    )
    ax2.plot(mid[:, 0], mid[:, 1], color="#1565C0", lw=1.8, label="skin mid-surface")
    ax2.plot(te[0], te[1], "r*", ms=12, label="aero TE")
    near = np.asarray(te_gap["nearest_midline_point_s"], dtype=float)
    ax2.plot(near[0], near[1], "o", color="#FF6F00", ms=7, label="nearest mid-surface pt")
    ax2.annotate(
        f"gap ≈ {te_gap['gap_te_side_m']:.4f} m\n(no zero-area TE spur)",
        xy=near,
        xytext=(te[0] - 0.06, te[1] + 0.025),
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color="0.3"),
    )
    _zoom(ax2, te[0] - 0.10, te[0] + 0.03, te[1] - 0.045, te[1] + 0.045)
    ax2.set_title("OuterSkin mid-surface ends short of aero TE by construction")
    ax2.legend(fontsize=8, loc="upper left")
    apply_structural_section_axes(
        ax2,
        xlabel=r"$Y_S$ (edgewise) [m]",
        ylabel=r"$X_S$ (flapwise) [m]",
    )
    _save_fig(fig2, out_dir / "skin_midline_vs_aero_te.png")

    # --- 3. LE nose SDF fidelity ---
    le = np.asarray(af.leading_edge, dtype=float)
    grid_le = SDFGrid.from_bbox(
        le[0] - 0.04,
        le[0] + 0.08,
        le[1] - 0.05,
        le[1] + 0.05,
        nx=320,
        ny=320,
    )
    phi_le = grid_le.eval(af)
    fig3, axes3 = plt.subplots(1, 2, figsize=(12.0, 5.2))
    plot_sdf_field(
        phi_le,
        grid_le,
        ax=axes3[0],
        title="Airfoil φ near LE (exact polyhedral SDF)",
        vmax=0.02,
    )
    axes3[0].plot(le[0], le[1], "k*", ms=10, label="LE point")
    # Overlay contour vertices in LE window
    le_mask = verts[:, 0] <= 0.12
    axes3[0].scatter(
        verts[le_mask, 0],
        verts[le_mask, 1],
        s=8,
        c="k",
        alpha=0.7,
        label="vertices",
        zorder=5,
    )
    axes3[0].legend(fontsize=8, loc="best")

    # Sample |φ| along LE vertices as a stem plot
    le_verts = verts[le_mask]
    phi_v = np.abs(af(le_verts[:, 0], le_verts[:, 1]))
    axes3[1].semilogy(le_verts[:, 0], np.maximum(phi_v, 1e-16), "o-", ms=4, color="#C62828")
    axes3[1].axhline(1e-12, color="0.5", ls="--", lw=0.9, label="1e-12 m")
    axes3[1].set_xlabel(r"$Y_B$ (edgewise) [m]")
    axes3[1].set_ylabel(r"$|\varphi|$ on vertices [m]")
    axes3[1].set_title("Polyhedral SDF: |φ| ≈ 0 on airfoil vertices (LE window)")
    axes3[1].grid(True, which="both", alpha=0.3)
    axes3[1].legend(fontsize=8)
    fig3.suptitle("LE nose exactness — AirfoilSDF polyhedral kernel", fontsize=12)
    _save_fig(fig3, out_dir / "le_nose_sdf_fidelity.png")

    # --- 3b. What sampling costs you at the nose ---
    # The flat 1e-16 result above proves the absence of a problem the reader has
    # never seen. Draw the problem: contour the same field off a grid and watch
    # the nose facet away.
    resolutions = (48, 96, 192, 384)
    nose_rows: list[dict] = []
    fig3b, axes3b = plt.subplots(
        1, 2, figsize=(14.0, 5.6), gridspec_kw={"width_ratios": (1.0, 1.25), "wspace": 0.28}
    )

    ax_nose, ax_err = axes3b
    # Coarse grids get the strongest colour: they are the ones that visibly fail.
    colours = [RED, "#e46a86", "#eea3b4", "#f5ccd6"]

    for colour, nx in zip(colours, resolutions):
        deviation, cell, segments = _sampled_nose_deviation(af, nx)
        if not segments:
            continue
        nose_rows.append({"nx": int(nx), "cell_m": cell, "nose_deviation_m": deviation})
        for k, seg in enumerate(segments):
            ax_nose.plot(
                seg[:, 0] * 1e3,
                seg[:, 1] * 1e3,
                color=colour,
                lw=2.0,
                label=f"sampled, {cell * 1e3:.1f} mm cell" if k == 0 else None,
                zorder=3,
            )

    ax_nose.plot(
        np.r_[verts[:, 0], verts[0, 0]] * 1e3,
        np.r_[verts[:, 1], verts[0, 1]] * 1e3,
        color=NAVY,
        lw=2.8,
        label="exact polyhedral SDF",
        zorder=5,
    )
    ax_nose.set_xlim((le[0] - 0.008) * 1e3, (le[0] + 0.05) * 1e3)
    ax_nose.set_ylim((le[1] - 0.030) * 1e3, (le[1] + 0.030) * 1e3)
    ax_nose.set_aspect("equal", adjustable="box")
    ax_nose.set_xlabel("chordwise [mm]")
    ax_nose.set_ylabel("thickness [mm]")
    ax_nose.set_title("Sample the field and the nose facets away", loc="left")
    ax_nose.legend(loc="center right", fontsize=11)

    if nose_rows:
        cells_mm = np.array([r["cell_m"] for r in nose_rows], dtype=float) * 1e3
        errors_mm = np.array([max(r["nose_deviation_m"], 1e-18) for r in nose_rows], dtype=float) * 1e3
        # The exact path's error is what |phi| already measured on the vertex loop.
        exact_floor_mm = max(float(np.max(phi_v)), 1e-17) * 1e3

        ax_err.loglog(cells_mm, errors_mm, "o-", color=RED, ms=9, lw=2.4, zorder=3)
        for cell_mm, error_mm in zip(cells_mm, errors_mm):
            ax_err.annotate(
                f"{error_mm:.2g} mm",
                xy=(cell_mm, error_mm),
                xytext=(0, 12),
                textcoords="offset points",
                ha="center",
                fontsize=11,
                color=RED,
                fontweight="bold",
            )
        exact_is_zero = float(np.max(phi_v)) == 0.0
        ax_err.axhline(exact_floor_mm, color=TEAL, ls="--", lw=2.2, zorder=2)
        ax_err.text(
            cells_mm.max() * 1.4,
            exact_floor_mm * 8,
            (
                "exact polyhedral SDF: |φ| = 0 on the vertex loop\n(drawn at 1e-14 mm only so a log axis can show it)"
                if exact_is_zero
                else f"exact polyhedral SDF: max |φ| = {exact_floor_mm:.1e} mm on the vertex loop"
            ),
            ha="right",
            va="bottom",
            fontsize=11.5,
            color=TEAL,
            fontweight="bold",
        )

        ax_err.set_xticks(cells_mm)
        ax_err.set_xticklabels([f"{value:.1f}" for value in cells_mm])
        ax_err.minorticks_off()
        ax_err.set_xlim(cells_mm.min() * 0.35, cells_mm.max() * 1.5)
        ax_err.set_ylim(exact_floor_mm * 0.05, errors_mm.max() * 20)
        ax_err.set_xlabel("grid cell size [mm]")
        ax_err.set_ylabel("worst deviation from the true nose [mm]")
        ax_err.set_title("The nose is only ever as sharp as the cell", loc="left")
        ax_err.grid(True, which="major", axis="y", alpha=0.5)

    _save_fig(fig3b, out_dir / "le_nose_grid_vs_exact.png", presentation=presentation)

    # --- 4. Twist case: B-frame LE/TE ---
    twist = float(np.radians(15.0))
    geom = build_section_geometry(
        resolve_system_type("2D-F"),
        af,
        skin_thickness=0.005,
        spar_cap_thickness=0.014,
        shear_web_thickness=0.006,
        twist_rad=twist,
        outputs=("midlines",),
    )
    section_tw = geom.section
    grid_tw = SDFGrid.from_airfoil(af, padding=0.06, nx=400, ny=180)
    # Evaluate section in chord frame; rotate airfoil outline for overlay
    fig4, ax4 = plt.subplots(figsize=(11.0, 5.0))
    plot_section(
        section_tw,
        grid_tw,
        ax=ax4,
        alpha=0.45,
        title=f"Twin-web @ {np.degrees(twist):.0f}° twist — B-frame LE/TE markers",
        show_legend=True,
        omit_cores=True,
    )
    # Overlay LE/TE in canonical B so markers match plot_section (Y_B horizontal).
    le_b = naca_to_twisted_pitch_axis_blade_frame(
        le.reshape(1, 2), twist, float(af.chord)
    )[0]
    te_b = naca_to_twisted_pitch_axis_blade_frame(
        te.reshape(1, 2), twist, float(af.chord)
    )[0]
    h_le, v_le = section_plot_axes(le_b.reshape(1, 2))
    h_te, v_te = section_plot_axes(te_b.reshape(1, 2))
    ax4.plot(h_le, v_le, "c*", ms=12, label="LE (B)")
    ax4.plot(h_te, v_te, "m*", ms=12, label="TE (B)")
    ax4.annotate(
        f"LE_B=({le_b[0]:.3f},{le_b[1]:.3f})\nTE_B=({te_b[0]:.3f},{te_b[1]:.3f})",
        xy=(0.02, 0.98),
        xycoords="axes fraction",
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.6"),
    )
    ax4.legend(fontsize=7, loc="lower right")
    _save_fig(fig4, out_dir / "twist15_le_te_markers.png")

    phi_le_stats = _phi_on_vertices(af, x_max=0.12)
    phi_te_stats = _phi_on_vertices(af, x_min=0.90)
    payload = {
        "airfoil": "NACA 2412",
        "chord_m": float(af.chord),
        "skin_thickness_m": 0.005,
        "te_tip_gap": te_gap,
        "phi_on_vertices_le_window": phi_le_stats,
        "phi_on_vertices_te_window": phi_te_stats,
        "twist_rad": twist,
        "leading_edge_twisted_pitch_axis_blade_frame": [float(le_b[0]), float(le_b[1])],
        "trailing_edge_twisted_pitch_axis_blade_frame": [float(te_b[0]), float(te_b[1])],
        "notes": {
            "hybrid_te": "plot_section_hybrid uses constructive wall-solid polygons; TE tip is not mesh-limited.",
            "skin_midline": (
                "OuterSkin.midline_polyline ends short of the aero tip where local "
                "thickness < t — avoids a zero-area TE spur."
            ),
            "airfoil_sdf": "AirfoilSDF builds an exact polyhedral SDF via sdf_polygon.",
        },
    }
    json_path = out_dir / "le_te_exactness.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(
        f"  TE mid-surface gap={te_gap['gap_te_side_m']:.4e} m; "
        f"max |phi| LE verts={phi_le_stats['max_abs_phi_m']:.3e} m"
    )
    return out_dir


if __name__ == "__main__":
    main()
