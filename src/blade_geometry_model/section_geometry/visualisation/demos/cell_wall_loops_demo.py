"""
visualisation/demos/cell_wall_loops_demo.py
=====================================
Plot constructive cell-wall loops on NACA midlines (0…5 webs).

Builds ``MultiCellSection``, exports B-frame midline strips, then
``build_cell_wall_loops``. Writes **one PNG per case** (no subplot grids):

* ``cell_wall_loops_{name}.png`` — filled cell loops over midlines (family A)
* ``cell_wall_loops_on_section_{name}.png`` — structural ``plot_section``
  with wall midlines + cell-wall loops as thin dashed centreline rings
* ``cell_wall_loops.json`` — serializable strips + loops for all cases

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/cell_wall_loops_demo.py

Outputs default to ``visualisation/output/cell_wall_loops/``.
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
    apply_presentation_style,
    demo_output_dir,
    save_figure,
)

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    MultiCellSection,
    SDFGrid,
)
# ``interface/`` is a re-export shim and carries none of these any more: the
# midline builder moved to ``midlines/`` and the loop builder to
# ``cell_wall_loops/``, alongside the rest of their own stage.
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.shell_midline_export import (
    build_shell_midline_strips,
)
from blade_structure_utilities.geometry_utilities.section_geometry.cell_wall_loops.cell_wall_loop_builders import (
    build_cell_wall_loops,
)
from blade_structure_utilities.geometry_utilities.section_geometry import plot_section
from blade_structure_utilities.geometry_utilities.coordinate_frames import (
    section_plot_axes,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.demos.cell_wall_loops_case_payload import (
    cell_wall_loops_case_payload,
)

_CASE_SPECS: tuple[tuple[str, tuple[float, ...]], ...] = (
    ("n_webs_0", ()),
    ("n_webs_1", (0.35,)),
    ("n_webs_2", (0.25, 0.55)),
    ("n_webs_3", (0.20, 0.40, 0.65)),
    ("n_webs_4", (0.15, 0.35, 0.55, 0.75)),
    ("n_webs_5", (0.12, 0.28, 0.44, 0.60, 0.78)),
)

# High-contrast loop colours (fill + edge); enough for n_webs=5 â†’ 6 cells
_LOOP_STYLES: tuple[tuple[str, str], ...] = (
    ("#e41a1c", "#fbb4ae"),  # red
    ("#377eb8", "#b3cde3"),  # blue
    ("#4daf4a", "#ccebc5"),  # green
    ("#ff7f00", "#fed9a6"),  # orange
    ("#984ea3", "#decbe4"),  # purple
    ("#a65628", "#e5d8bd"),  # brown
)

_FIGSIZE = (14.5, 5.8)
_DPI = DEFAULT_DPI


def _build_section(
    af: AirfoilSDF,
    web_fracs: tuple[float, ...],
    *,
    structural_family: str = "A",
) -> MultiCellSection:
    chord = float(af.chord)
    xs = [f * chord for f in web_fracs]
    # Family D's box spar spans the outermost webs, so it needs two. Follow the
    # resolver's own step-down: one web hands D to C (a cap per web), none to A.
    family = structural_family
    if not xs and family not in ("A", "B"):
        family = "A"
    elif len(xs) == 1 and family == "D":
        family = "C"
    # Exaggerated member widths for demo readability (not design values).
    return MultiCellSection(
        af,
        web_x_positions=xs,
        web_thickness=0.022,
        skin_thickness=0.016,
        web_alignment="chord_normal",
        twist_angle=0.0,
        structural_family=family,
        cap_height=0.040,
        core_enabled=False,
    )


def _axis_limits(strips, pad_frac: float = 0.06) -> tuple[float, float, float, float]:
    """Tight padded box in section-plot coords (h=Y_B, v=X_B)."""
    pts = [np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float) for s in strips if s.kind in ("skin", "web")]
    all_xy = np.vstack(pts)
    h, v = section_plot_axes(all_xy)
    xmin, xmax = float(h.min()), float(h.max())
    ymin, ymax = float(v.min()), float(v.max())
    dx = max(xmax - xmin, 1e-9)
    dy = max(ymax - ymin, 1e-9)
    return (
        xmin - pad_frac * dx,
        xmax + pad_frac * dx,
        ymin - pad_frac * dy,
        ymax + pad_frac * dy,
    )


def _plot_midline_loops(ax, strips, loops, title: str) -> None:
    legend_handles: list[Line2D] = []

    for s in strips:
        if s.kind == "skin":
            xy = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
            h, v = section_plot_axes(xy)
            (ln,) = ax.plot(
                h,
                v,
                color="#444444",
                lw=1.4,
                ls=":",
                zorder=1,
                label="skin midline",
            )
            legend_handles.append(ln)
        elif s.kind == "web":
            xy = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
            h, v = section_plot_axes(xy)
            (ln,) = ax.plot(
                h,
                v,
                color="#111111",
                lw=1.8,
                ls="--",
                zorder=5,
                label="web midline",
            )
            if not any(handle.get_label() == "web midline" for handle in legend_handles):
                legend_handles.append(ln)

    for i, lp in enumerate(loops):
        xy = np.asarray(lp.polyline_twisted_pitch_axis_blade_frame, dtype=float)
        edge, fill = _LOOP_STYLES[i % len(_LOOP_STYLES)]
        h, v = section_plot_axes(xy)  # polyline is canonical B
        ax.fill(
            h,
            v,
            facecolor=fill,
            edgecolor=edge,
            linewidth=2.2,
            alpha=0.55,
            zorder=3,
        )
        ax.plot(h, v, color=edge, lw=2.2, solid_capstyle="round", zorder=4)
        legend_handles.append(
            Line2D([0], [0], color=edge, lw=2.2, label=f"cell_{lp.id}  A={lp.area_m2:.4g} mÂ²")
        )

    from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (
        apply_structural_section_axes,
    )

    xmin, xmax, ymin, ymax = _axis_limits(strips)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    apply_structural_section_axes(ax)
    ax.set_title(title, fontsize=13, pad=8)
    ax.tick_params(labelsize=11)
    ax.set_facecolor("#f7f7f7")
    ax.grid(True, which="major", color="0.82", lw=0.9, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("0.3")
        spine.set_linewidth(1.1)

    seen: dict[str, Line2D] = {}
    for h in legend_handles:
        lab = h.get_label()
        if lab and lab not in seen and not lab.startswith("_"):
            seen[lab] = h
    ax.legend(
        seen.values(),
        seen.keys(),
        fontsize=9,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=True,
        fancybox=False,
        edgecolor="0.55",
        framealpha=0.97,
    )


def _draw_component_midlines(ax, strips) -> list[Line2D]:
    """Faint solid midlines so dashed cell loops can be checked against centres."""
    handles: list[Line2D] = []
    skin_drawn = False
    web_drawn = False
    for s in strips:
        xy = np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float)
        h, v = section_plot_axes(xy)
        if s.kind == "skin":
            (ln,) = ax.plot(
                h,
                v,
                color="#1a1a1a",
                lw=0.9,
                ls="-",
                alpha=0.55,
                zorder=15,
                label="skin midline",
            )
            if not skin_drawn:
                handles.append(ln)
                skin_drawn = True
        elif s.kind == "web":
            (ln,) = ax.plot(
                h,
                v,
                color="#111111",
                lw=1.0,
                ls="-",
                alpha=0.65,
                zorder=16,
                label="web midline",
            )
            if not web_drawn:
                handles.append(ln)
                web_drawn = True
    return handles


def _overlay_loops_dashed(ax, loops, *, lw: float = 1.35) -> list[Line2D]:
    """Draw each cell wall loop as a thin dashed centreline over a section plot."""
    handles: list[Line2D] = []
    for i, lp in enumerate(loops):
        xy = np.asarray(lp.polyline_twisted_pitch_axis_blade_frame, dtype=float)
        edge, _fill = _LOOP_STYLES[i % len(_LOOP_STYLES)]
        h, v = section_plot_axes(xy)
        (ln,) = ax.plot(
            h,
            v,
            color=edge,
            lw=lw,
            ls=(0, (2.8, 1.6)),
            solid_capstyle="round",
            zorder=20,
            label=f"cell_{lp.id}  A={lp.area_m2:.4g} mÂ²",
        )
        handles.append(ln)
    return handles


def _plot_section_with_loops(ax, section, grid, strips, loops, title: str) -> None:
    plot_section(
        section,
        grid,
        ax=ax,
        alpha=0.42,
        show_airfoil=False,
        title=title,
        show_legend=True,
        omit_cores=True,
        auto_limits=True,
        limit_pad_frac=0.06,
        linewidths=0.85,
    )
    mid_handles = _draw_component_midlines(ax, strips)
    loop_handles = _overlay_loops_dashed(ax, loops)
    # Append cell-loop legend below the component legend
    ax.add_artist(ax.get_legend())
    ax.legend(
        mid_handles + loop_handles,
        [h.get_label() for h in mid_handles + loop_handles],
        fontsize=9,
        loc="center left",
        bbox_to_anchor=(1.01, 0.35),
        frameon=True,
        fancybox=False,
        edgecolor="0.55",
        framealpha=0.97,
        title="centrelines / cell loops",
        title_fontsize=9,
    )


def _save_fig(fig, path: Path, *, presentation: bool = False) -> None:
    fig.tight_layout()
    save_figure(fig, path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {path}")


def main(argv: list[str] | None = None) -> Path:
    default_out = Path(__file__).resolve().parents[1] / "output" / "cell_wall_loops"
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

    # Remove legacy combined subplot PNGs if present
    for legacy in (
        "cell_wall_loops.png",
        "cell_wall_loops_on_section.png",
    ):
        legacy_path = out_dir / legacy
        if legacy_path.is_file():
            legacy_path.unlink()

    af = AirfoilSDF.from_naca("2412", chord=1.0)
    # Finer grid for single-panel section renders
    grid = SDFGrid.from_airfoil(af, padding=0.05, nx=720, ny=300)
    cases: list[dict] = []

    for name, fracs in _CASE_SPECS:
        n_webs = len(fracs)

        # Family A â€” midline / filled-loop detail
        section_a = _build_section(af, fracs, structural_family="A")
        strips_a = build_shell_midline_strips(section_a, twist_rad=0.0)
        loops_a = build_cell_wall_loops(strips_a)
        cases.append(cell_wall_loops_case_payload(strips_a, name=name, loops=loops_a))
        if len(loops_a) != n_webs + 1:
            raise RuntimeError(f"{name}: expected {n_webs + 1} loops, got {len(loops_a)}")

        fig, ax = plt.subplots(figsize=_FIGSIZE)
        _plot_midline_loops(
            ax,
            strips_a,
            loops_a,
            title=(
                f"Void-cell wall loops â€” {name}  "
                f"(NACA 2412, family A, {n_webs} web(s) â†’ {len(loops_a)} loop(s))"
            ),
        )
        _save_fig(fig, out_dir / f"cell_wall_loops_{name}.png", presentation=presentation)

        # Family D (A if 0 webs) â€” structural section + dashed centreline loops
        section_d = _build_section(af, fracs, structural_family="D")
        strips_d = build_shell_midline_strips(section_d, twist_rad=0.0)
        loops_d = build_cell_wall_loops(strips_d)
        fam = "A" if n_webs == 0 else ("C" if n_webs == 1 else "D")

        fig2, ax2 = plt.subplots(figsize=_FIGSIZE)
        _plot_section_with_loops(
            ax2,
            section_d,
            grid,
            strips_d,
            loops_d,
            title=(
                f"Section geometry + cell-wall loops — {name}  "
                f"(NACA 2412, family {fam}, {n_webs} web(s) → {len(loops_d)} loop(s))"
            ),
        )
        _save_fig(fig2, out_dir / f"cell_wall_loops_on_section_{name}.png", presentation=presentation)
        print(f"  {name}: n_webs={n_webs}, n_loops={len(loops_d)}, family={fam}")

    json_path = out_dir / "cell_wall_loops.json"
    json_path.write_text(json.dumps({"cases": cases}, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return out_dir


if __name__ == "__main__":
    main()

