"""
visualisation/demos/naca2412_section.py
============================
End-to-end demonstration of the section_geometry package.

Demonstrates:
  1.  NACA 2412 airfoil SDF
  2.  Classic twin-web box spar (corrected curved spar caps)
  3.  3-web torsion box at 15° twist — separate all-flapwise (F) and
      all-chord-normal (CN) layouts (not mixed per-web alignment)
  4.  Worker path via ``build_section_geometry`` (props + constructive midlines)
  5.  Section properties report
  6.  Per-component SDF fields

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/naca2412_section.py

Outputs under ``visualisation/output/naca2412_section/``.

For thickness-only geometry (no properties sweep), see
``geometry_from_design_vector.py``.
"""

from __future__ import annotations

from pathlib import Path
import importlib.util

# When run as a script, load repo root onto sys.path so ``import blade_structure_utilities`` works.
_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
if not _erp.is_file():
    raise RuntimeError(f"Missing helper {_erp.name}; run from the repository checkout.")
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)

import numpy as np
import matplotlib.pyplot as plt

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
    SectionPropertiesReport,
    build_section_geometry,
    export_section_json,
    plot_section,
    resolve_system_type,
)
from blade_geometry_model.coordinate_frames import (
    naca_to_twisted_pitch_axis_blade_frame,
    section_plot_axes,
)
from blade_geometry_model.section_geometry.visualisation.plot_section import (
    apply_structural_section_axes,
    blade_plot_sampling_grid,
)

import argparse

_parser = argparse.ArgumentParser(description=__doc__)
add_output_and_presentation_args(_parser, default_stem="naca2412_section")
_args = _parser.parse_args()
apply_presentation_style(bool(_args.presentation))
_PRESENTATION = bool(_args.presentation)
OUTPUT = Path(_args.out_dir)
OUTPUT.mkdir(parents=True, exist_ok=True)


def _build_torsion_box(af: AirfoilSDF, *, twist_rad: float, web_alignment: str) -> MultiCellSection:
    """Uniform web alignment for the 3-web torsion-box demo (F or CN)."""
    return MultiCellSection(
        airfoil_sdf=af,
        web_x_positions=[0.15, 0.35, 0.55],
        web_thickness=[0.005, 0.004, 0.005],
        web_alignment=str(web_alignment),
        cap_height=(0.014, 0.012),
        skin_thickness=0.003,
        twist_angle=float(twist_rad),
        trailing_edge_insert_x_m=0.80,
        leading_edge_insert_x_m=0.08,
        core_enabled=True,
    )


def _save_torsion_section_plot(
    section: MultiCellSection,
    grid: SDFGrid,
    *,
    af: AirfoilSDF,
    xc: np.ndarray,
    yc: np.ndarray,
    twist_rad: float,
    title: str,
    out_name: str,
) -> None:
    fig, ax = plot_section(section, grid, title=title)
    camber_b = naca_to_twisted_pitch_axis_blade_frame(
        np.column_stack([xc, yc]), float(twist_rad), float(af.chord)
    )
    h_c, v_c = section_plot_axes(camber_b)
    ax.plot(h_c, v_c, "k--", lw=1.2, label="camber line (analytical)")
    ax.plot(0.0, 0.0, "k+", ms=10, mew=1.5, label="pitch axis (origin)")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    path = OUTPUT / out_name
    save_figure(fig, path, presentation=_PRESENTATION)
    plt.close(fig)
    print(f"  Saved: {path}")


# =============================================================================
# 1. Airfoil
# =============================================================================

print("=" * 60)
print("  NACA 2412 Blade Section Geometry Demo")
print("=" * 60)

af = AirfoilSDF.from_naca("2412", chord=1.0)
print(f"\nAirfoil: {af}")

xc, yc = af.camber_line(n_points=100)
xc_t, t_dist = af.thickness_distribution(n_points=100)
print(f"  Max thickness : {t_dist.max():.4f}c at x/c = {xc_t[np.argmax(t_dist)]:.3f}")
print(f"  Max camber    : {yc.max():.4f}c  at x/c = {xc[np.argmax(yc)]:.3f}")


# =============================================================================
# 2. Twin-web box spar (corrected curved spar caps)
# =============================================================================

print("\n--- Twin-web box spar (0 twist) -------------------------")

twin = MultiCellSection.twin_web(
    af,
    x_fore         = 0.20,
    x_aft          = 0.50,
    skin_thickness = 0.003,
    web_thickness  = 0.004,
    cap_height     = (0.014, 0.012),    # (upper, lower)
    web_alignment  = "chord_normal",
    twist_angle    = 0.0,
    trailing_edge_insert_x_m    = 0.75,
    leading_edge_insert_x_m    = 0.10,
    core_enabled   = True,
)
print(f"  Components: {twin.labels}")
print(f"  n_webs={twin.n_webs}, n_cells={twin.n_cells}")

grid_twin = SDFGrid.from_airfoil(af, padding=0.08, nx=600, ny=240)

fig, ax = plot_section(twin, grid_twin,
                       title="NACA 2412 — twin-web box spar (corrected curved caps)")
camber_b = naca_to_twisted_pitch_axis_blade_frame(
    np.column_stack([xc, yc]), 0.0, float(af.chord)
)
h_c, v_c = section_plot_axes(camber_b)
ax.plot(h_c, v_c, "k--", lw=1.2, label="camber line")
ax.plot(0.0, 0.0, "k+", ms=10, mew=1.5, label="pitch axis (origin)")
ax.legend(fontsize=7, loc="upper right")
fig.tight_layout()
save_figure(fig, OUTPUT / "twin_web_section.png", presentation=_PRESENTATION)
plt.close(fig)
print(f"  Saved: {OUTPUT / 'twin_web_section.png'}")


# =============================================================================
# 3. Three-web torsion box @ 15° — all-F then all-CN (uniform alignment)
# =============================================================================

print("\n--- 3-web torsion box (15 deg twist, F then CN) ----------")

TWIST = np.radians(15)
grid_tors = SDFGrid.from_airfoil(af, padding=0.08, nx=700, ny=280)

torsion_f = _build_torsion_box(af, twist_rad=float(TWIST), web_alignment="flapwise")
torsion_cn = _build_torsion_box(af, twist_rad=float(TWIST), web_alignment="chord_normal")
# Downstream props / SDF panels use CN (webs tilt with the section).
torsion = torsion_cn

for tag, section, align_note, out_name in (
    (
        "F",
        torsion_f,
        "all webs flapwise (B-vertical)",
        "torsion_box_section_flapwise.png",
    ),
    (
        "CN",
        torsion_cn,
        "all webs chord-normal (tilt with section)",
        "torsion_box_section_chord_normal.png",
    ),
):
    print(f"  [{tag}] Components : {section.labels}")
    print(f"  [{tag}] n_webs={section.n_webs}, n_cells={section.n_cells}")
    _save_torsion_section_plot(
        section,
        grid_tors,
        af=af,
        xc=xc,
        yc=yc,
        twist_rad=float(TWIST),
        title=(
            f"NACA 2412 — 3-web torsion box @ {np.degrees(TWIST):.0f}° twist\n"
            f"({align_note})"
        ),
        out_name=out_name,
    )

# Legacy alias: same as CN (replaces the old mixed-alignment figure).
_save_torsion_section_plot(
    torsion_cn,
    grid_tors,
    af=af,
    xc=xc,
    yc=yc,
    twist_rad=float(TWIST),
    title=(
        f"NACA 2412 — 3-web torsion box @ {np.degrees(TWIST):.0f}° twist\n"
        "(all webs chord-normal)"
    ),
    out_name="torsion_box_section.png",
)
print(f"  Twist      : {np.degrees(TWIST):.1f}°")


# =============================================================================
# 4. Worker path: build_section_geometry
# =============================================================================

print("\n--- build_section_geometry (2D-F layout) -----------------")

geom = build_section_geometry(
    resolve_system_type("2D-F"),
    af,
    skin_thickness=0.003,
    spar_cap_thickness=0.012,
    shear_web_thickness=0.004,
    twist_rad=float(TWIST),
    outputs=("properties", "midlines"),
)
n_strips = len(geom.midline_strips or ())
n_labels = len(geom.properties_map or {})
print(f"  Midline strips : {n_strips}")
print(f"  Property labels: {n_labels}")
print(f"  Backend        : {geom.properties_meta.backend if geom.properties_meta else None}")


# =============================================================================
# 5. Section properties
# =============================================================================

print("\n--- Section properties -----------------------------------")

report_twin   = SectionPropertiesReport(twin,    grid_twin)
report_torsion = SectionPropertiesReport(torsion, grid_tors)

print("\nTwin-web box spar:")
print(report_twin.summary())

print("\n3-web torsion box (15° twist, chord-normal webs):")
print(report_torsion.summary())


# =============================================================================
# 6. Per-component SDF grid (CN torsion box — 3 webs, 4 full-span core bays)
# =============================================================================

labels_to_plot = [
    "skin_suction_surface",
    "skin_pressure_surface",
    "spar_cap_suction_surface",
    "spar_cap_pressure_surface",
    "shear_web_0",
    "shear_web_1",
    "shear_web_2",
    "sandwich_core_0",
    "sandwich_core_1",
    "sandwich_core_2",
]
fig, axes = plt.subplots(2, 5, figsize=(20, 7))
Y_plot, X_plot, QX, QY = blade_plot_sampling_grid(
    af,
    twist_rad=float(TWIST),
    chord_m=float(af.chord),
    nx=int(grid_tors.X.shape[1]),
    ny=int(grid_tors.X.shape[0]),
)
for ax, lbl in zip(axes.ravel(), labels_to_plot):
    phi = np.asarray(torsion[lbl](QX, QY), dtype=float)
    im = ax.pcolormesh(
        Y_plot,
        X_plot,
        phi,
        cmap="RdBu_r",
        vmin=-0.05,
        vmax=0.05,
        shading="auto",
    )
    ax.contour(Y_plot, X_plot, phi, levels=[0.0], colors="k", linewidths=1.0)
    ax.set_title(lbl, fontsize=9)
    # Keep data limits (airfoil B bbox); equal-aspect with datalim would pad
    # huge empty flapwise margins and look like a clipped thin band.
    apply_structural_section_axes(ax)
    ax.set_aspect("equal", adjustable="box")
    fig.colorbar(im, ax=ax, label="φ")
fig.suptitle(
    "Per-component SDF fields — NACA 2412 torsion box @ 15° twist (all chord-normal)",
    fontsize=11,
)
fig.tight_layout()
save_figure(fig, OUTPUT / "component_sdfs_torsion.png", presentation=_PRESENTATION)
plt.close(fig)
print(f"  Saved: {OUTPUT / 'component_sdfs_torsion.png'}")


# =============================================================================
# 7. Export
# =============================================================================

print("\n--- Exporting data ---------------------------------------")

report_torsion.to_json(str(OUTPUT / "section_properties_torsion.json"))
report_torsion.to_csv(str(OUTPUT / "section_properties_torsion.csv"))
export_section_json(
    torsion,
    grid_tors,
    {},
    str(OUTPUT / "section_full_torsion.json"),
)

print("  Exported: section_properties_torsion.json/csv")
print("  Exported: section_full_torsion.json")
print("\nDone.")
