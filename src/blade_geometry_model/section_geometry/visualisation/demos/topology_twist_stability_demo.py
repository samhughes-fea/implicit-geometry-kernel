"""
visualisation/demos/topology_twist_stability_demo.py
=========================================
Computational-geometry demo: shell topology stability under twist.

``build_shell_strip_topology`` snaps web feet onto skin and inserts cap–web
T-junction vertices. This demo sweeps section twist and shows that junction
gaps stay near geometry tolerance — topology ops are stable under rigid motion
of the IR section.

Writes:

* ``junction_before_after.png`` — the worst raw web foot before and after topology
  preparation, beside the gap across the twist sweep (the slide figure)
* ``gap_vs_twist.png`` — max web–skin gap (raw vs topology) vs twist
* ``junction_zoom_worst_twist.png`` — junction zooms at the worst twist
* ``topology_twist_stability.json`` — per-twist metrics

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/topology_twist_stability_demo.py
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

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
)
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.ring_geometry import (
    point_to_open_polyline_nearest,
)
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_structure_utilities.geometry_utilities.section_geometry.topology.shell_strip_topology import (
    build_shell_strip_topology,
)
from blade_structure_utilities.geometry_utilities.coordinate_frames import (
    section_plot_axes,
)

_DPI = DEFAULT_DPI
_TWIST_DEG = (0.0, 5.0, 10.0, 15.0, 30.0)


def _build_section(af: AirfoilSDF, *, twist_rad: float) -> MultiCellSection:
    return MultiCellSection.torsion_box(
        af,
        x_fore=0.15,
        x_mid=0.35,
        x_aft=0.55,
        skin_thickness=0.005,
        web_thickness=0.006,
        cap_height=(0.014, 0.012),
        web_alignment="chord_normal",
        twist_angle=float(twist_rad),
        structural_family="D",
        core_enabled=False,
    )


def _nearest_point_distance(pt: np.ndarray, poly: np.ndarray) -> float:
    return float(np.min(np.linalg.norm(poly - pt.reshape(1, 2), axis=1)))


def _web_skin_gaps(strips) -> list[dict]:
    skins = [np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float) for s in strips if s.kind == "skin"]
    if not skins:
        return []
    skin = np.vstack(skins)
    gaps: list[dict] = []
    for s in strips:
        if s.kind != "web":
            continue
        pts = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if pts.shape[0] < 2:
            continue
        for end_name, row in (("top", 0), ("bot", -1)):
            pt = pts[row]
            gaps.append(
                {
                    "web_label": str(s.label),
                    "end": end_name,
                    "gap_to_skin_vertex_m": _nearest_point_distance(pt, skin),
                }
            )
    return gaps


def _cap_web_shared_count(strips, *, tol_m: float = 2.5e-3) -> int:
    webs = [
        (s.label, np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float))
        for s in strips
        if s.kind == "web" and np.asarray(s.midline_twisted_pitch_axis_blade_frame).shape[0] >= 2
    ]
    n = 0
    for s in strips:
        if s.kind != "cap":
            continue
        pts = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for pt in pts:
            best = 1e300
            for _wlab, wpts in webs:
                gap, _q = point_to_open_polyline_nearest(pt, wpts)
                best = min(best, float(gap))
            if best <= tol_m:
                n += 1
    return n


def _overlay_strips(ax, strips, *, lw: float = 1.8) -> None:
    kind_colour = {"skin": "#2196F3", "web": "#4CAF50", "cap": "#F44336"}
    for strip in strips:
        pts = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if pts.shape[0] < 2:
            continue
        h, v = section_plot_axes(pts)
        ax.plot(
            h,
            v,
            color=kind_colour.get(str(strip.kind), "#607D8B"),
            lw=lw,
            alpha=0.95,
        )


def _save_fig(fig, path: Path) -> None:
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)
    print(f"Wrote {path}")


def _worst_raw_foot(strips) -> tuple[str, int, np.ndarray]:
    """The web end furthest from any skin vertex in the raw midlines."""
    skin = np.vstack(
        [np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float) for s in strips if s.kind == "skin"]
    )
    worst: tuple[float, str, int, np.ndarray] | None = None
    for strip in strips:
        if strip.kind != "web":
            continue
        points = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        for row in (0, -1):
            gap = _nearest_point_distance(points[row], skin)
            if worst is None or gap > worst[0]:
                worst = (gap, str(strip.label), row, points[row].copy())
    if worst is None:
        raise RuntimeError("No web strips to inspect.")
    return worst[1], worst[2], worst[3]


def _draw_junction_zoom(ax, strips, centre: np.ndarray, *, half: float, show_nodes: bool) -> None:
    from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
        region_color,
        strip_axes_chrome,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
    )

    kind_colour = {
        "skin": region_color("skin_suction_surface"),
        "web": region_color("shear_web"),
        "cap": region_color("spar_cap_suction_surface"),
    }
    for strip in strips:
        points = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if points.shape[0] < 2:
            continue
        horizontal, vertical = section_plot_axes(points)
        colour = kind_colour.get(str(strip.kind), "#607D8B")
        ax.plot(horizontal, vertical, color=colour, lw=3.0, zorder=2, solid_capstyle="round")
        if show_nodes:
            ax.plot(horizontal, vertical, "o", ms=5, color=colour, mec="white", mew=0.8, zorder=3)
    centre_h, centre_v = section_plot_axes(np.asarray(centre, dtype=float).reshape(1, 2))
    apply_structural_section_axes(ax)
    ax.set_xlim(float(centre_h[0]) - half, float(centre_h[0]) + half)
    ax.set_ylim(float(centre_v[0]) - half, float(centre_v[0]) + half)
    ax.set_aspect("equal", adjustable="box")
    strip_axes_chrome(ax, keep_frame=True)


def _draw_junction_story(rows, raw_by_twist, topo_by_twist, out_path: Path, *, presentation: bool) -> None:
    """Before / after at the worst web foot, then the gap across the twist sweep."""
    from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
        INK,
        MUTED,
        RED,
        TEAL,
    )

    story_twist = 15.0 if 15.0 in raw_by_twist else rows[0]["twist_deg"]
    raw = raw_by_twist[story_twist]
    topo = topo_by_twist[story_twist]
    web_label, row_index, foot = _worst_raw_foot(raw)
    skin = np.vstack(
        [np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float) for s in raw if s.kind == "skin"]
    )
    nearest_skin = skin[np.argmin(np.linalg.norm(skin - foot.reshape(1, 2), axis=1))]
    raw_gap = float(np.linalg.norm(nearest_skin - foot))
    prepared_foot = next(
        np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)[row_index]
        for s in topo
        if s.kind == "web" and str(s.label) == web_label
    )
    half = max(4.0 * raw_gap, 0.012)

    fig = plt.figure(figsize=(14.0, 5.0))
    grid_spec = fig.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 1.25), wspace=0.3)
    ax_raw = fig.add_subplot(grid_spec[0, 0])
    ax_prepared = fig.add_subplot(grid_spec[0, 1])
    ax_sweep = fig.add_subplot(grid_spec[0, 2])

    _draw_junction_zoom(ax_raw, raw, foot, half=half, show_nodes=True)
    (foot_h,), (foot_v,) = section_plot_axes(foot.reshape(1, 2))
    (skin_h,), (skin_v,) = section_plot_axes(nearest_skin.reshape(1, 2))
    ax_raw.annotate(
        "", xy=(skin_h, skin_v), xytext=(foot_h, foot_v),
        arrowprops={"arrowstyle": "<->", "color": RED, "lw": 2.0}, zorder=5,
    )
    ax_raw.scatter([foot_h], [foot_v], s=260, facecolors="none", edgecolors=RED, linewidths=2.4, zorder=6)
    ax_raw.text(
        0.04, 0.96, f"web foot lands {raw_gap * 1e3:.1f} mm from\nthe nearest skin node: a hinge",
        transform=ax_raw.transAxes, fontsize=12.5, color=RED, fontweight="bold", va="top",
    )
    ax_raw.set_title("Raw midlines", loc="left", fontsize=14, color=INK)

    _draw_junction_zoom(ax_prepared, topo, foot, half=half, show_nodes=True)
    (shared_h,), (shared_v,) = section_plot_axes(prepared_foot.reshape(1, 2))
    ax_prepared.scatter(
        [shared_h], [shared_v], s=260, facecolors="none", edgecolors=TEAL, linewidths=2.4, zorder=6
    )
    ax_prepared.text(
        0.04, 0.96, "skin split at the foot:\none shared node",
        transform=ax_prepared.transAxes, fontsize=12.5, color=TEAL, fontweight="bold", va="top",
    )
    ax_prepared.set_title("After topology preparation", loc="left", fontsize=14, color=INK)

    degrees = np.array([r["twist_deg"] for r in rows], dtype=float)
    raw_gaps = np.array([max(r["max_web_skin_gap_raw_m"], 1e-16) for r in rows]) * 1e3
    prepared_gaps = np.array([max(r["max_web_skin_gap_topology_m"], 1e-16) for r in rows]) * 1e3
    ax_sweep.semilogy(degrees, raw_gaps, "o-", color=RED, lw=2.4, ms=8)
    ax_sweep.semilogy(degrees, prepared_gaps, "o-", color=TEAL, lw=2.4, ms=8)
    ax_sweep.text(degrees[-1], raw_gaps[-1] * 4, f"raw  {raw_gaps.max():.1f} mm", ha="right",
                  va="bottom", fontsize=12, color=RED, fontweight="bold")
    ax_sweep.text(degrees[-1], prepared_gaps[-1] * 4, f"prepared  {prepared_gaps.max():.0e} mm",
                  ha="right", va="bottom", fontsize=12, color=TEAL, fontweight="bold")
    orders = int(np.floor(np.log10(raw_gaps.min() / prepared_gaps.max())))
    middle = 0.25 * (degrees[0] + degrees[-1])
    ax_sweep.annotate(
        "", xy=(middle, prepared_gaps.max() * 1.5), xytext=(middle, raw_gaps.min() / 1.5),
        arrowprops={"arrowstyle": "<->", "color": INK, "lw": 1.4},
    )
    ax_sweep.text(middle + 1.0, np.sqrt(raw_gaps.min() * prepared_gaps.max()), f"{orders} orders",
                  fontsize=13, color=INK, fontweight="bold", va="center")
    ax_sweep.set_xlabel("section twist [deg]")
    ax_sweep.set_ylabel("worst web-foot gap [mm]")
    ax_sweep.set_ylim(prepared_gaps.min() / 100, raw_gaps.max() * 100)
    ax_sweep.set_title("Every foot, every twist", loc="left", fontsize=14, color=INK)
    ax_sweep.grid(True, which="major", axis="y", alpha=0.4)
    ax_sweep.tick_params(axis="y", colors=MUTED)

    save_figure(fig, out_path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "topology_twist_stability"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    parser.add_argument("--presentation", action="store_true", help="Deck styling for the junction story.")
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", chord=1.0)
    rows: list[dict] = []
    topo_by_twist: dict[float, object] = {}
    raw_by_twist: dict[float, object] = {}

    for deg in _TWIST_DEG:
        twist = float(np.deg2rad(deg))
        section = _build_section(af, twist_rad=twist)
        raw = build_shell_midline_strips(section, twist_rad=twist)
        topo = build_shell_strip_topology(raw)
        topo_by_twist[deg] = topo
        raw_by_twist[deg] = raw
        gaps_raw = _web_skin_gaps(raw)
        gaps_topo = _web_skin_gaps(topo)
        row = {
            "twist_deg": float(deg),
            "twist_rad": twist,
            "max_web_skin_gap_raw_m": max(
                (g["gap_to_skin_vertex_m"] for g in gaps_raw), default=float("nan")
            ),
            "max_web_skin_gap_topology_m": max(
                (g["gap_to_skin_vertex_m"] for g in gaps_topo), default=float("nan")
            ),
            "mean_web_skin_gap_topology_m": float(
                np.mean([g["gap_to_skin_vertex_m"] for g in gaps_topo])
            )
            if gaps_topo
            else float("nan"),
            "n_cap_web_shared_vertices": _cap_web_shared_count(topo),
            "n_skin_strips": sum(1 for s in topo if s.kind == "skin"),
            "n_web_strips": sum(1 for s in topo if s.kind == "web"),
        }
        rows.append(row)

    _draw_junction_story(
        rows, raw_by_twist, topo_by_twist, out_dir / "junction_before_after.png",
        presentation=presentation,
    )

    # Gap vs twist
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    degs = [r["twist_deg"] for r in rows]
    raw_g = [max(r["max_web_skin_gap_raw_m"], 1e-16) for r in rows]
    topo_g = [max(r["max_web_skin_gap_topology_m"], 1e-16) for r in rows]
    ax.semilogy(degs, raw_g, "o--", color="#EF6C00", lw=1.8, label="max gap raw")
    ax.semilogy(
        degs, topo_g, "s-", color="#2E7D32", lw=2.0, label="max gap after topology"
    )
    ax.set_xlabel("section twist [deg]")
    ax.set_ylabel("max web–skin foot gap [m]")
    ax.set_title("Topology snap stability under twist (NACA 2412 3-web box)")
    ax.grid(True, which="both", alpha=0.35)
    ax.legend(fontsize=9)
    _save_fig(fig, out_dir / "gap_vs_twist.png")

    # Worst topology gap (should still be tiny); zoom that case.
    worst = max(rows, key=lambda r: r["max_web_skin_gap_topology_m"])
    topo = topo_by_twist[worst["twist_deg"]]
    gaps_topo = _web_skin_gaps(topo)
    # Pick a few feet for zoom panels.
    centres = []
    for g in gaps_topo:
        # Reconstruct point from strips for zoom centre.
        pass
    # Use web endpoints from topo strips.
    for s in topo:
        if s.kind != "web":
            continue
        pts = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if pts.shape[0] < 2:
            continue
        centres.append((f"{s.label} top @ {worst['twist_deg']:g}°", pts[0].copy()))
        centres.append((f"{s.label} bot @ {worst['twist_deg']:g}°", pts[-1].copy()))
    centres = centres[:6]

    n = len(centres)
    ncols = 3 if n > 3 else max(n, 1)
    nrows = int(np.ceil(n / ncols)) if n else 1
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(13.0, 3.8 * nrows))
    axes_flat = np.atleast_1d(axes2).ravel()
    half = 0.035
    for ax, (title, centre) in zip(axes_flat, centres):
        _overlay_strips(ax, topo, lw=2.0)
        # centre is canonical B (X_B, Y_B); plot horizontal = Y_B.
        h_c, v_c = section_plot_axes(np.asarray(centre, dtype=float).reshape(1, 2))
        ax.scatter(
            [float(h_c[0])],
            [float(v_c[0])],
            s=40,
            c="#FFC107",
            edgecolors="k",
            linewidths=0.4,
            zorder=10,
        )
        ax.set_xlim(float(h_c[0]) - half, float(h_c[0]) + half)
        ax.set_ylim(float(v_c[0]) - half, float(v_c[0]) + half)
        ax.grid(True, alpha=0.3)
        ax.set_title(title, fontsize=9)
        from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
            apply_structural_section_axes,
        )

        apply_structural_section_axes(ax)
    for ax in axes_flat[n:]:
        ax.axis("off")
    fig2.suptitle(
        f"Junction zooms at worst twist "
        f"(topo max gap={worst['max_web_skin_gap_topology_m']:.3e} m)",
        fontsize=11,
    )
    _save_fig(fig2, out_dir / "junction_zoom_worst_twist.png")

    payload = {
        "airfoil": "NACA 2412",
        "layout": "torsion_box family D (x=0.15/0.35/0.55)",
        "twist_deg_sweep": list(_TWIST_DEG),
        "rows": rows,
        "worst_twist_deg": float(worst["twist_deg"]),
        "max_topology_gap_over_sweep_m": max(
            r["max_web_skin_gap_topology_m"] for r in rows
        ),
        "notes": (
            "Web–skin snap gaps after topology remain near machine/geometry "
            "tolerance across the twist sweep."
        ),
    }
    json_path = out_dir / "topology_twist_stability.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    for r in rows:
        print(
            f"  twist={r['twist_deg']:5.1f}°  raw={r['max_web_skin_gap_raw_m']:.3e}  "
            f"topo={r['max_web_skin_gap_topology_m']:.3e}  "
            f"cap-web nodes={r['n_cap_web_shared_vertices']}"
        )
    return out_dir


if __name__ == "__main__":
    main()
