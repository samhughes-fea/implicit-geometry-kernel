"""
visualisation/demos/shell_junction_topology_demo.py
=======================================
Demonstrate constructive shell midline junction topology.

Builds a three-web NACA 2412 box section, exports B-frame midlines, then
applies ``prepare_shell_midline_junctions`` (skin split at web/cap feet +
cap-web T-junction vertices). Writes:

* ``topology_before_after.png`` — raw midlines vs topology-prepared strips
* ``topology_junction_zooms.png`` — web–skin feet and cap–web T-junction zooms
* ``shell_junction_topology.json`` — strip counts + endpoint gap metrics

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/shell_junction_topology_demo.py

Outputs default to ``visualisation/output/shell_junction_topology/``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    DEFAULT_DPI,
    add_output_and_presentation_args,
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot import (
    plot_shell_midline_strip_topology,
)
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.ring_geometry import (
    point_to_open_polyline_nearest,
)
from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import (
    prepare_shell_midline_junctions,
)
from blade_structure_utilities.geometry_utilities.coordinate_frames import (
    section_plot_axes,
)

_DPI = DEFAULT_DPI
_FIGSIZE_FULL = (14.0, 5.5)
_FIGSIZE_ZOOM = (14.0, 10.0)

_KIND_COLOUR = {"skin": "#2196F3", "web": "#4CAF50", "cap": "#F44336"}
_JUNCTION_COLOUR = "#FFC107"


def _build_section(af: AirfoilSDF) -> MultiCellSection:
    # Three-web family D: outermost webs define the box; the mid web crosses
    # the caps and creates interior cap-web T-junction nodes.
    return MultiCellSection.torsion_box(
        af,
        x_fore=0.15,
        x_mid=0.35,
        x_aft=0.55,
        skin_thickness=0.005,
        web_thickness=0.006,
        cap_height=(0.014, 0.012),
        web_alignment="chord_normal",
        twist_angle=0.0,
        structural_family="D",
        core_enabled=False,
    )


def _strip_counts(strips) -> dict[str, int]:
    counts: dict[str, int] = {}
    for s in strips:
        counts[s.kind] = counts.get(s.kind, 0) + 1
    return counts


def _nearest_point_distance(pt: np.ndarray, poly: np.ndarray) -> float:
    d = np.linalg.norm(poly - pt.reshape(1, 2), axis=1)
    return float(np.min(d))


def _web_skin_gaps(strips) -> list[dict]:
    """Distance from each web endpoint to nearest skin vertex."""
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
                    "web_label": s.label,
                    "end": end_name,
                    "point_twisted_pitch_axis_blade_frame": [float(pt[0]), float(pt[1])],
                    "gap_to_skin_vertex_m": _nearest_point_distance(pt, skin),
                }
            )
    return gaps


def _cap_web_shared_vertices(strips, *, tol_m: float = 5e-4) -> list[dict]:
    """Cap vertices near a web polyline (interior T-nodes or bay-end joints)."""
    webs = [
        (s.label, np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float))
        for s in strips
        if s.kind == "web" and np.asarray(s.midline_twisted_pitch_axis_blade_frame).shape[0] >= 2
    ]
    hits: list[dict] = []
    for s in strips:
        if s.kind != "cap":
            continue
        pts = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        n = int(pts.shape[0])
        for i, pt in enumerate(pts):
            best_gap = float("inf")
            best_web = None
            best_q = None
            for wlab, wpts in webs:
                gap, q = point_to_open_polyline_nearest(pt, wpts)
                if gap < best_gap:
                    best_gap = float(gap)
                    best_web = str(wlab)
                    best_q = np.asarray(q, dtype=float).reshape(2)
            if best_web is None or best_gap > tol_m:
                continue
            role = "endpoint" if i in (0, n - 1) else "interior"
            hits.append(
                {
                    "cap_label": str(getattr(s.label, "value", s.label)),
                    "web_label": best_web,
                    "vertex_index": int(i),
                    "role": role,
                    "point_twisted_pitch_axis_blade_frame": [float(pt[0]), float(pt[1])],
                    "web_point_twisted_pitch_axis_blade_frame": [float(best_q[0]), float(best_q[1])],
                    "gap_to_web_m": best_gap,
                }
            )
    return hits


def _zoom_axes(ax, centre: np.ndarray, half_span: float) -> None:
    from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
    )

    # centre is canonical B (X_B, Y_B); plot horizontal = Y_B.
    h, v = section_plot_axes(np.asarray(centre, dtype=float).reshape(1, 2))
    cx, cy = float(h[0]), float(v[0])
    ax.set_xlim(cx - half_span, cx + half_span)
    ax.set_ylim(cy - half_span, cy + half_span)
    ax.grid(True, alpha=0.3)
    apply_structural_section_axes(ax)


def _pretty_label(label: object) -> str:
    text = str(getattr(label, "value", label))
    if text.startswith("ComponentLabel."):
        text = text.split(".", 1)[1].lower()
    # Collapse split skin segments to a single legend entry
    if ":seg" in text:
        text = text.split(":seg", 1)[0]
    return text


def _overlay_strips(ax, strips, *, lw: float = 1.6, label_strips: bool = True) -> None:
    seen: set[str] = set()
    for strip in strips:
        pts = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=float)
        if pts.shape[0] < 2:
            continue
        colour = _KIND_COLOUR.get(str(strip.kind), "#607D8B")
        pretty = _pretty_label(strip.label)
        legend = pretty if label_strips and pretty not in seen else None
        seen.add(pretty)
        _h, _v = section_plot_axes(pts)  # canonical B -> horizontal = Y_B
        ax.plot(_h, _v, color=colour, lw=lw, label=legend, alpha=0.95)
        if strip.closed:
            closed = np.vstack([pts[-1], pts[0]])
            ch, cv = section_plot_axes(closed)
            ax.plot(
                ch,
                cv,
                color=colour,
                lw=lw,
                alpha=0.95,
            )


def _mark_points(ax, points: list[np.ndarray], *, colour: str, label: str | None) -> None:
    if not points:
        return
    xy = np.asarray(points, dtype=float).reshape(-1, 2)
    h, v = section_plot_axes(xy)
    ax.scatter(
        h,
        v,
        s=36,
        c=colour,
        edgecolors="k",
        linewidths=0.4,
        zorder=10,
        label=label,
    )


def _draw_shared_key(fig, spare_axes) -> None:
    """Draw one key for every zoom panel: three strip kinds plus the junction marker.

    Hosted in the first unused grid cell when the panel count leaves one, so
    ``tight_layout`` cannot push it over the data; otherwise it sits under the
    figure title.
    """
    handles = [
        Line2D([], [], color=colour, lw=2.0, label=kind)
        for kind, colour in _KIND_COLOUR.items()
    ]
    handles.append(
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=7,
            markerfacecolor=_JUNCTION_COLOUR,
            markeredgecolor="k",
            markeredgewidth=0.4,
            label="junction (shared vertex)",
        )
    )

    spare = list(spare_axes)
    if spare:
        spare[0].legend(handles=handles, loc="center", fontsize=11, frameon=False)
    else:
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), fontsize=10)


def _save_fig(fig, path: Path, *, presentation: bool = False) -> None:
    fig.tight_layout()
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="shell_junction_topology")
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))
    out_dir: Path = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    af = AirfoilSDF.from_naca("2412", chord=1.0)
    section = _build_section(af)
    raw = build_shell_midline_strips(section, twist_rad=0.0)
    topo = prepare_shell_midline_junctions(raw)

    # --- full before / after ---
    fig, axes = plt.subplots(1, 2, figsize=_FIGSIZE_FULL)
    plot_shell_midline_strip_topology(
        raw,
        ax=axes[0],
        title="Raw constructive midlines (pre-topology)",
    )
    plot_shell_midline_strip_topology(
        topo,
        ax=axes[1],
        title="After prepare_shell_midline_junctions\n(skin split + cap–web T nodes)",
    )
    _save_fig(fig, out_dir / "topology_before_after.png", presentation=bool(args.presentation))

    # --- junction metrics ---
    gaps_raw = _web_skin_gaps(raw)
    gaps_topo = _web_skin_gaps(topo)
    # Interior mid-web T-nodes are sub-mm; bay-end cap/web feet are ~1-2 mm.
    t_junctions = _cap_web_shared_vertices(topo, tol_m=2.5e-3)

    fore_feet = [g for g in gaps_topo if g["web_label"] == "shear_web_0"]
    aft_feet = [g for g in gaps_topo if g["web_label"] == "shear_web_2"]
    mid_feet = [g for g in gaps_topo if g["web_label"] == "shear_web_1"]
    def _cap_side(label: str) -> str | None:
        key = str(label).lower()
        if "upper" in key:
            return "upper"
        if "lower" in key:
            return "lower"
        return None

    interior_t = [h for h in t_junctions if h.get("role") == "interior"]
    upper_t = [
        np.asarray(h["point_twisted_pitch_axis_blade_frame"], dtype=float)
        for h in interior_t
        if _cap_side(h["cap_label"]) == "upper"
    ]
    lower_t = [
        np.asarray(h["point_twisted_pitch_axis_blade_frame"], dtype=float)
        for h in interior_t
        if _cap_side(h["cap_label"]) == "lower"
    ]

    centres: list[tuple[str, np.ndarray]] = []
    if fore_feet:
        centres.append(("Fore web-skin (shear_web_0 top)", np.asarray(fore_feet[0]["point_twisted_pitch_axis_blade_frame"])))
        centres.append(("Fore web-skin (shear_web_0 bot)", np.asarray(fore_feet[-1]["point_twisted_pitch_axis_blade_frame"])))
    if aft_feet:
        centres.append(("Aft web-skin (web_2 top)", np.asarray(aft_feet[0]["point_twisted_pitch_axis_blade_frame"])))
    if upper_t:
        centres.append(("Cap-web T (upper mid web)", upper_t[0]))
    if lower_t:
        centres.append(("Cap-web T (lower mid web)", lower_t[0]))
    elif mid_feet:
        centres.append(("Mid web-skin (web_1 top)", np.asarray(mid_feet[0]["point_twisted_pitch_axis_blade_frame"])))
    centres = centres[:6]

    n = len(centres)
    ncols = 3 if n > 3 else max(n, 1)
    nrows = int(np.ceil(n / ncols)) if n else 1
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=_FIGSIZE_ZOOM)
    axes_flat = np.atleast_1d(axes2).ravel()
    half = 0.035
    for ax, (title, centre) in zip(axes_flat, centres):
        # One entry per strip label would be a 13-row legend sitting on top of the
        # junction each panel exists to show; the shared key below carries it instead.
        _overlay_strips(ax, topo, lw=2.0, label_strips=False)
        _mark_points(ax, [centre], colour=_JUNCTION_COLOUR, label=None)
        _zoom_axes(ax, centre, half)
        ax.set_title(title, fontsize=10)
    for ax in axes_flat[n:]:
        ax.axis("off")
    _draw_shared_key(fig2, axes_flat[n:])
    fig2.suptitle(
        "Shell junction zooms - shared vertices after topology prep (NACA 2412 3-web box)",
        fontsize=12,
    )
    _save_fig(fig2, out_dir / "topology_junction_zooms.png", presentation=bool(args.presentation))

    payload = {
        "airfoil": "NACA 2412",
        "chord_m": float(af.chord),
        "layout": "torsion_box family D (x=0.15/0.35/0.55)",
        "strip_counts_raw": _strip_counts(raw),
        "strip_counts_topology": _strip_counts(topo),
        "web_skin_gaps_raw_m": gaps_raw,
        "web_skin_gaps_topology_m": gaps_topo,
        "max_web_skin_gap_raw_m": max((g["gap_to_skin_vertex_m"] for g in gaps_raw), default=None),
        "max_web_skin_gap_topology_m": max(
            (g["gap_to_skin_vertex_m"] for g in gaps_topo), default=None
        ),
        "cap_web_shared_vertices": t_junctions,
        "n_cap_web_shared_vertices": len(t_junctions),
        "notes": (
            "Topology snaps web feet onto skin vertices and inserts shared "
            "cap–web T-junction nodes for mesh continuity."
        ),
    }
    json_path = out_dir / "shell_junction_topology.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(
        f"  max web-skin gap raw={payload['max_web_skin_gap_raw_m']:.3e} m -> "
        f"topo={payload['max_web_skin_gap_topology_m']:.3e} m; "
        f"cap-web shared nodes={payload['n_cap_web_shared_vertices']}"
    )
    return out_dir


if __name__ == "__main__":
    main()
