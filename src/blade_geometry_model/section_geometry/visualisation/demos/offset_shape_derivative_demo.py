"""
visualisation/demos/offset_shape_derivative_demo.py
===================================================
The derivative the field gives for free: d/dt of an offset is a boundary integral.

Thickness is ``Offset(aerofoil, -t)`` (slide "Thickness is a field operation").
For an exact signed-distance field, the inner face ``{phi = -t}`` moves inward
at unit normal speed as ``t`` grows, so by the level-set (Hadamard) shape
derivative the skin region ``{-t < phi < 0}`` changes by

    dA/dt   = ∮_{phi=-t} ds          (the length of the inner face)
    dIxx/dt = ∮_{phi=-t} y² ds       (its second moment about y = 0)

No re-solve, no step size: one contour of the field already evaluated.

The check is against the kernel's own property path — central finite
differences of the skin-only ``0A`` layout through
``build_section_view`` → ``compute_section_properties(backend="hybrid")``, the
path slide "What the sizing search sees" differentiates numerically.

The identity needs ``|grad phi| = 1`` on the contour. That holds for the
aerofoil field, which is why the package offsets the aerofoil and never an
assembled union (slide "…and exactly where that stops being true").

Sweep starts at 3 mm: ``layout/thickness_limits.MIN_REALISTIC_SKIN_WALL_THICKNESS_M``
(2.5 mm) floors the requested skin, so below it the finite difference is zero
by construction, not by geometry.

Writes under ``visualisation/output/offset_shape_derivative/``:

* ``offset_shape_derivative.png``
* ``offset_shape_derivative.json``

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/offset_shape_derivative_demo.py --presentation
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

from blade_geometry_model.section_geometry import (
    AirfoilSDF,
    SDFGrid,
    build_section_view,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.layout.thickness_limits import (
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M,
)
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    quadtree_spec_for_airfoil,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import (
    compute_section_properties,
)
from blade_geometry_model.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_geometry_model.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    NAVY,
    RED,
    TEAL,
    strip_axes_chrome,
)

_FD_STEP_M = 2.0e-5
_TOLERANCE = 1.0e-2


def _skin_totals(airfoil, layout, spec, thickness: float) -> tuple[float, float]:
    """Skin area and second moment about y = 0 from the kernel's property path."""
    section = build_section_view(
        airfoil, layout, skin_thickness=thickness, spar_cap_thickness=0.0, shear_web_thickness=0.0
    )
    props, _ = compute_section_properties(section, backend="hybrid", quadtree_spec=spec, airfoil_for_spec=airfoil)
    labels = [label for label in props if str(label).startswith("skin")]
    if not labels:
        raise RuntimeError("0A section produced no skin regions")
    area = sum(float(props[label]["area"]) for label in labels)
    # Centroidal Ixx back to the fixed axis y = 0, so it can be differentiated.
    moment = sum(float(props[label]["Ixx"]) + float(props[label]["area"]) * float(props[label]["cy"]) ** 2
                 for label in labels)
    return area, moment


def _boundary_integrals(grid: SDFGrid, phi: np.ndarray, thickness: float) -> tuple[float, float, list]:
    """∮ ds and ∮ y² ds on the contour phi = -thickness (exact per straight segment)."""
    loops = grid.zero_contour(phi + thickness)
    length = 0.0
    moment = 0.0
    for loop in loops:
        points = np.asarray(loop, dtype=float)
        segments = np.diff(points, axis=0)
        lengths = np.linalg.norm(segments, axis=1)
        mid_y = 0.5 * (points[1:, 1] + points[:-1, 1])
        length += float(lengths.sum())
        moment += float((lengths * (mid_y**2 + segments[:, 1] ** 2 / 12.0)).sum())
    return length, moment, loops


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="offset_shape_derivative")
    parser.add_argument("--t-min", type=float, default=0.003)
    parser.add_argument("--t-max", type=float, default=0.030)
    parser.add_argument("--samples", type=int, default=14)
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))
    if float(args.t_min) - _FD_STEP_M <= MIN_REALISTIC_SKIN_WALL_THICKNESS_M:
        raise SystemExit(
            f"--t-min must exceed the {MIN_REALISTIC_SKIN_WALL_THICKNESS_M * 1e3:g} mm skin floor plus the FD step"
        )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    airfoil = AirfoilSDF.from_naca("2412", chord=1.0, n_points=240)
    layout = resolve_system_type("0A")
    spec = quadtree_spec_for_airfoil(airfoil, chord=1.0, max_depth=7, max_cells=120000)
    grid = SDFGrid.from_airfoil(airfoil, nx=3000, ny=1300, padding=0.02)
    phi = np.asarray(grid.eval(airfoil), dtype=float)

    thicknesses = np.linspace(float(args.t_min), float(args.t_max), int(args.samples))
    rows = []
    for thickness in thicknesses:
        area_plus, moment_plus = _skin_totals(airfoil, layout, spec, thickness + _FD_STEP_M)
        area_minus, moment_minus = _skin_totals(airfoil, layout, spec, thickness - _FD_STEP_M)
        length, moment, _ = _boundary_integrals(grid, phi, thickness)
        fd_area = (area_plus - area_minus) / (2.0 * _FD_STEP_M)
        fd_moment = (moment_plus - moment_minus) / (2.0 * _FD_STEP_M)
        rows.append({
            "t_m": float(thickness),
            "dA_dt_analytic": length,
            "dA_dt_fd": fd_area,
            "dI_dt_analytic": moment,
            "dI_dt_fd": fd_moment,
            "rel_err_A": abs(fd_area - length) / length,
            "rel_err_I": abs(fd_moment - moment) / moment,
        })
    worst = max(max(row["rel_err_A"], row["rel_err_I"]) for row in rows)
    if worst > _TOLERANCE:
        raise RuntimeError(f"analytic and finite-difference derivatives disagree by {worst:.2e} (> {_TOLERANCE})")

    fig = plt.figure(figsize=(14.0, 6.2))
    grid_spec = fig.add_gridspec(2, 2, height_ratios=(0.62, 1.0), hspace=0.42, wspace=0.2)
    ax_section = fig.add_subplot(grid_spec[0, :])
    ax_value = fig.add_subplot(grid_spec[1, 0])
    ax_error = fig.add_subplot(grid_spec[1, 1])

    # Left: the moving boundary. Inner faces at several gauges, and the unit
    # normal speed drawn on one of them.
    outline = np.asarray(airfoil.vertices, dtype=float)
    ax_section.fill(outline[:, 0], outline[:, 1], color="#eef1f5", zorder=0)
    ax_section.plot(np.r_[outline[:, 0], outline[0, 0]], np.r_[outline[:, 1], outline[0, 1]], color=INK, lw=1.8)
    shown = (0.005, 0.015, 0.025)
    shades = (TEAL, NAVY, MUTED)
    for thickness, shade in zip(shown, shades):
        _, _, loops = _boundary_integrals(grid, phi, thickness)
        for index, loop in enumerate(loops):
            loop = np.asarray(loop, dtype=float)
            ax_section.plot(loop[:, 0], loop[:, 1], color=shade, lw=1.6,
                            label=f"inner face, t = {thickness * 1e3:g} mm" if index == 0 else None)
    _, _, loops = _boundary_integrals(grid, phi, shown[0])
    face = np.asarray(max(loops, key=len), dtype=float)
    gradient_x, gradient_y = np.gradient(phi, grid.Y[:, 0], grid.X[0, :])[1], np.gradient(phi, grid.Y[:, 0], grid.X[0, :])[0]
    # The trailing-edge wedge is where the inner face meets the medial axis and
    # the normal flips; keep the arrows to where the face is smooth.
    face = face[face[:, 0] < 0.85]
    picks = face[:: max(1, len(face) // 30)]
    columns = np.clip(np.searchsorted(grid.X[0, :], picks[:, 0]), 0, grid.X.shape[1] - 1)
    rows_index = np.clip(np.searchsorted(grid.Y[:, 0], picks[:, 1]), 0, grid.X.shape[0] - 1)
    normal_x = -gradient_x[rows_index, columns]
    normal_y = -gradient_y[rows_index, columns]
    ax_section.quiver(picks[:, 0], picks[:, 1], normal_x, normal_y, color=RED, scale=38, width=0.0028,
                      headwidth=4, zorder=5)
    ax_section.set_aspect("equal")
    ax_section.set_xlim(-0.02, 1.02)
    strip_axes_chrome(ax_section)
    ax_section.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=10.5)
    ax_section.set_title("Offset moves the inner face at unit speed (red): dA/dt is its length", loc="left",
                         fontsize=13, color=INK)

    t_mm = np.array([row["t_m"] for row in rows]) * 1e3
    ax_value.plot(t_mm, [row["dA_dt_analytic"] for row in rows], color=TEAL, lw=2.6, label="∮ ds on one contour")
    ax_value.plot(t_mm, [row["dA_dt_fd"] for row in rows], "o", color=INK, ms=6, mfc="none", mew=1.6,
                  label="finite difference, 2 kernel solves each")
    ax_value.set_ylabel("dA/dt  [m²/m]")
    ax_value.set_xlabel("skin gauge t [mm]")
    ax_value.set_title("Skin area sensitivity to gauge", loc="left", fontsize=13, color=INK)
    ax_value.legend(frameon=False, fontsize=10.5, loc="upper right")
    ax_value.grid(True, alpha=0.4)

    ax_error.semilogy(t_mm, [row["rel_err_A"] for row in rows], "o-", color=TEAL, lw=2.0, ms=5, label="dA/dt")
    ax_error.semilogy(t_mm, [row["rel_err_I"] for row in rows], "s-", color=NAVY, lw=2.0, ms=5,
                      label="dIxx/dt  (∮ y² ds)")
    ax_error.set_xlabel("skin gauge t [mm]")
    ax_error.set_title("Analytic vs finite difference", loc="left", fontsize=13, color=INK)
    ax_error.set_ylabel("relative difference")
    ax_error.set_ylim(1e-7, 1e-2)
    ax_error.legend(frameon=False, fontsize=10.5, loc="lower right", ncol=2)
    ax_error.grid(True, which="major", alpha=0.4)

    png_path = out_dir / "offset_shape_derivative.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    payload = {
        "airfoil": "NACA 2412, chord 1 m",
        "layout": "0A (skin only)",
        "fd_step_m": _FD_STEP_M,
        "skin_floor_m": MIN_REALISTIC_SKIN_WALL_THICKNESS_M,
        "max_rel_err_A": max(row["rel_err_A"] for row in rows),
        "max_rel_err_I": max(row["rel_err_I"] for row in rows),
        "rows": rows,
    }
    json_path = out_dir / "offset_shape_derivative.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(f"max rel err A {payload['max_rel_err_A']:.2e}  I {payload['max_rel_err_I']:.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
