"""Labelled anatomy figure for layout 2D-F (teaching companion to the gallery).

Writes two figures to visualisation/output/anatomy_2d_f/:

* ``2D-F_labelled.png`` — the section with the standard component legend, so docs
  can name suction/pressure skin, webs, and spar caps.
* ``2D-F_anatomy_annotated.png`` — the same section with the names written *on*
  the regions via leader lines, plus the datums (OML, inner mould line, chord,
  pitch axis, LE/TE) and the skin gauge. For an audience meeting the object for
  the first time, a legend is an extra decoding step; leader lines are not.

Gallery tiles from structural_family_grids_demo.py stay legend-less.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/anatomy_2d_f_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_geometry import (
    SectionGeometryRunner,
)
from blade_structure_utilities.geometry_utilities.section_geometry.api.factory import (
    build_section_view,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.airfoil import AirfoilSDF
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    strip_axes_chrome,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot import (
    plot_section_hybrid,
)
from blade_structure_utilities.geometry_utilities.section_geometry.constructive import (
    emit_wall_solid_polygons,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
    apply_structural_section_axes,
)
from blade_structure_utilities.geometry_utilities.coordinate_frames import (
    naca_to_twisted_pitch_axis_blade_frame,
    section_plot_axes,
)

_FIGSIZE = (12.0, 5.0)
_FIGSIZE_ANNOTATED = (13.5, 6.6)
_LAYOUT_KEY = "2D-F"
_SKIN_THICKNESS_M = 0.006
_SPAR_CAP_THICKNESS_M = 0.014
_SHEAR_WEB_THICKNESS_M = 0.006

# Which side of the section each wall's label sits on. Positions are worked out
# from the geometry at draw time, not hard-coded.
_WALL_SIDE = {
    "skin_suction_surface": ("Suction-side outer skin", "up"),
    "spar_cap_suction_surface": ("Spar cap — suction", "up"),
    "shear_web_0": ("Shear web 0 (fore)", "up"),
    "shear_web_1": ("Shear web 1 (aft)", "up"),
    "skin_pressure_surface": ("Pressure-side outer skin", "down"),
    "spar_cap_pressure_surface": ("Spar cap — pressure", "down"),
}

# The bays between the walls. emit_wall_solid_polygons returns wall solids only, and
# the sandwich_core_* regions overlap rather than partitioning the cavity, so the
# cells are located from the resolved web stations instead.
_CELL_NAMES = ("Leading-edge cell", "Box cell (between webs)", "Trailing-edge cell")

_REGION_COLOUR = {
    "skin_suction_surface": "#1565C0",
    "skin_pressure_surface": "#64B5F6",
    "spar_cap_suction_surface": "#E53935",
    "spar_cap_pressure_surface": "#FB8C00",
    "shear_web_0": "#2E7D32",
    "shear_web_1": "#26A69A",
    "sandwich_core_0": "#CE93D8",
    "sandwich_core_1": "#8E24AA",
    "sandwich_core_2": "#BA68C8",
}


def _to_plot(points, chord: float, twist_rad: float = 0.0) -> np.ndarray:
    """Rotated-NACA storage coordinates -> the canonical B plot frame.

    ``apply_structural_section_axes`` labels the axes Y_B / X_B but does not
    transpose data, so everything drawn here is converted first.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 2)
    blade = naca_to_twisted_pitch_axis_blade_frame(pts, float(twist_rad), float(chord))
    h, v = section_plot_axes(blade)
    return np.column_stack([np.asarray(h, dtype=float), np.asarray(v, dtype=float)])


def _polygon_parts(geom):
    """Every simple polygon in a possibly-multi geometry, largest first."""
    parts = list(getattr(geom, "geoms", [geom]))
    return sorted(parts, key=lambda g: float(g.area), reverse=True)


def _place_labels(ax, items, *, text_row: float, va: str) -> None:
    """Lay leader-line labels out across the span in anchor order.

    Sorting the labels by their anchor and distributing the text evenly in the
    same order is what keeps the leaders from crossing — hand-tuned offsets stop
    working the moment the geometry changes.
    """
    if not items:
        return
    ordered = sorted(items, key=lambda it: -float(it[1][0]))
    span_hi, span_lo = 0.62, -0.30
    n = len(ordered)
    text_positions = (
        [0.5 * (span_hi + span_lo)]
        if n == 1
        else np.linspace(span_hi, span_lo, n)
    )
    for (name, anchor, colour), text_x in zip(ordered, text_positions):
        ax.annotate(
            name,
            xy=(float(anchor[0]), float(anchor[1])),
            xytext=(float(text_x), text_row),
            fontsize=10,
            ha="center",
            va=va,
            zorder=8,
            arrowprops=dict(arrowstyle="-", color="#37474F", lw=0.9, shrinkA=0.0, shrinkB=3.0),
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=colour, lw=1.4, alpha=0.96),
        )


def _draw_annotated(section, airfoil, web_x_m, *, twist_rad: float, out_dir: Path,
                    presentation: bool) -> Path:
    """The section with its parts named on the drawing, plus the datums."""
    chord = float(airfoil.chord)
    polygons = emit_wall_solid_polygons(section)

    fig, ax = plt.subplots(figsize=_FIGSIZE_ANNOTATED)

    up_items: list[tuple[str, np.ndarray, str]] = []
    down_items: list[tuple[str, np.ndarray, str]] = []
    for label, geom in polygons.items():
        if label not in _WALL_SIDE:
            continue
        name, side = _WALL_SIDE[label]
        colour = _REGION_COLOUR.get(label, "#90A4AE")
        parts = _polygon_parts(geom)
        for part in parts:
            ring = _to_plot(np.asarray(part.exterior.coords, dtype=float), chord, twist_rad)
            ax.fill(ring[:, 0], ring[:, 1], color=colour, alpha=0.85, lw=0.0, zorder=2)
        # representative_point is guaranteed inside, unlike a centroid on a
        # crescent-shaped skin band.
        rep = parts[0].representative_point()
        anchor = _to_plot([[rep.x, rep.y]], chord, twist_rad)[0]
        (up_items if side == "up" else down_items).append((name, anchor, colour))

    outer = _to_plot(np.asarray(airfoil.vertices, dtype=float), chord, twist_rad)
    ax.plot(
        np.r_[outer[:, 0], outer[0, 0]],
        np.r_[outer[:, 1], outer[0, 1]],
        color="#0D1B2A",
        lw=2.0,
        zorder=4,
    )

    # Bays, located between the resolved web stations rather than from the
    # overlapping sandwich_core_* regions.
    stations = [0.0, *[float(x) for x in web_x_m], chord]
    for name, x_lo, x_hi in zip(_CELL_NAMES, stations[:-1], stations[1:]):
        x_mid = 0.5 * (x_lo + x_hi)
        camber_y = float(np.interp(x_mid, [0.0, chord], [airfoil.leading_edge[1], airfoil.trailing_edge[1]]))
        anchor = _to_plot([[x_mid, camber_y]], chord, twist_rad)[0]
        down_items.append((name, anchor, "#8E24AA"))

    le = _to_plot([airfoil.leading_edge], chord, twist_rad)[0]
    te = _to_plot([airfoil.trailing_edge], chord, twist_rad)[0]
    ax.plot([le[0], te[0]], [le[1], te[1]], ls=":", lw=1.2, color="#546E7A", zorder=3)
    for point in (le, te):
        ax.plot(point[0], point[1], marker="*", ms=16, color="#0D1B2A", zorder=6)
    ax.annotate("leading\nedge", xy=(le[0], le[1]), xytext=(le[0] - 0.055, le[1] + 0.02),
                fontsize=9, ha="center", va="center", zorder=7)
    ax.annotate("trailing\nedge", xy=(te[0], te[1]), xytext=(te[0] + 0.055, te[1] + 0.02),
                fontsize=9, ha="center", va="center", zorder=7)

    # The pitch axis is the origin of the blade frame, so it is exactly Y_B = 0.
    ax.axvline(0.0, ls="--", lw=1.2, color="#B0BEC5", zorder=1)

    # Explicit rows, not values read back from the axes: with equal aspect the
    # y-limits are only settled at draw time, so reading them here places the
    # labels on top of the section.
    row_up, row_down = 0.150, -0.150

    _place_labels(ax, up_items, text_row=row_up, va="bottom")
    _place_labels(ax, down_items, text_row=row_down, va="top")

    # Anchor on the crest of the outer loop; the loop is closed, so interpolating
    # against its x would pick a point on whichever surface came first.
    crest = outer[int(np.argmax(outer[:, 1]))]
    ax.annotate(
        f"exterior / interior skin-shell surface\nskin gauge {_SKIN_THICKNESS_M * 1e3:.0f} mm — one scalar in the tree",
        xy=(float(crest[0]), float(crest[1])),
        # Kept clear of the upper label row in x, but pulled down in y: the long
        # leader read as pointing at empty space rather than at the two lines.
        xytext=(float(crest[0]) + 0.16, row_up + 0.048),
        fontsize=9,
        color="#37474F",
        ha="center",
        zorder=8,
        arrowprops=dict(arrowstyle="->", color="#37474F", lw=1.0),
    )
    ax.annotate("pitch axis", xy=(0.0, row_down - 0.030),
                fontsize=9, color="#607D8B", ha="center", va="top", zorder=7)

    ax.set_title(
        "One station of a tidal blade: nine named regions, every wall gauge a design variable\n"
        f"layout {_LAYOUT_KEY} — NACA 2412, chord {chord:g} m",
        fontsize=12,
    )
    apply_structural_section_axes(ax)
    ax.grid(False)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(outer[:, 0].max() + 0.09, outer[:, 0].min() - 0.09)
    ax.set_ylim(row_down - 0.075, row_up + 0.135)
    if presentation:
        # On a slide the tick numbers carry nothing — every dimension a viewer
        # needs is already called out by a leader — and dropping them plus the
        # frame hands their space back to the section. The docs build keeps them.
        strip_axes_chrome(ax)

    fig.tight_layout()
    out_png = out_dir / "2D-F_anatomy_annotated.png"
    save_figure(fig, out_png, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {out_png}")
    return out_png


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="anatomy_2d_f")
    parser.add_argument("--twist-deg", type=float, default=0.0)
    parser.add_argument("--naca", default="2412")
    parser.add_argument("--chord", type=float, default=1.0)
    args = parser.parse_args(argv)

    apply_presentation_style(bool(args.presentation))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    airfoil = AirfoilSDF.from_naca(args.naca, n_points=140, chord=float(args.chord))
    layout = resolve_system_type(_LAYOUT_KEY)
    twist_rad = float(np.deg2rad(float(args.twist_deg)))
    section = build_section_view(
        airfoil,
        layout,
        twist_angle_rad=twist_rad,
        skin_thickness=_SKIN_THICKNESS_M,
        spar_cap_thickness=_SPAR_CAP_THICKNESS_M,
        shear_web_thickness=_SHEAR_WEB_THICKNESS_M,
    )

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    plot_section_hybrid(
        section,
        residual_mesh=None,
        ax=ax,
        alpha=0.5,
        title=(
            f"{_LAYOUT_KEY} anatomy - suction/pressure skin, spar caps, shear webs\n"
            f"(NACA {args.naca}, chord={args.chord:g} m, twist={args.twist_deg:g} deg)"
        ),
        show_legend=True,
    )
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()

    out_png = out_dir / "2D-F_labelled.png"
    save_figure(fig, out_png, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {out_png}")

    resolution = (
        SectionGeometryRunner(
            layout=layout,
            airfoil=airfoil,
            skin_thickness=_SKIN_THICKNESS_M,
            spar_cap_thickness=_SPAR_CAP_THICKNESS_M,
            shear_web_thickness=_SHEAR_WEB_THICKNESS_M,
            twist_rad=twist_rad,
        )
        .resolve()
        .get_layout_resolution()
    )

    _draw_annotated(
        section,
        airfoil,
        resolution.web_x_positions_m,
        twist_rad=twist_rad,
        out_dir=out_dir,
        presentation=bool(args.presentation),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
