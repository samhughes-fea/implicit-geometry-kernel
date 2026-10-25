"""
visualisation/demos/spanwise_yz_plan_flyover.py
===============================================
Fly outboard along the blade and watch where the shear webs go.

A section is a cut; a blade is a taper. This animates the consequence, comparing
two ways of placing webs down the span:

* **station-relative** — each web holds a constant fraction of the *local* chord,
  so the web pattern fans with the planform. This is what the package does
  (``api/factory`` multiplies ``layout.web_chord_fracs`` by the local chord).
* **root-anchored** — each web keeps the absolute ``Y_B`` offset it had at the
  root. As the chord shrinks the webs march toward the edges in ``x/c`` until the
  local cavity can no longer anchor them, and they are dropped.

The point of the comparison is that a web which no longer fits is *dropped, not
slid*: sliding it would silently move a load path.

Panels, top to bottom:

1. the section at the current station, both models side by side;
2. the plan view ``(Z_R, Y_B)`` — looking down on the blade, with the web tracks
   running out along the span inside the planform;
3. web station in ``x/c`` against radius, dashed for the station-relative target
   and solid for where root-anchoring actually puts each web.

Writes ``spanwise_yz_plan_flyover.gif`` plus its frames and a JSON sidecar to
``visualisation/output/spanwise_yz_plan_flyover/``.

Requires ImageMagick (``magick``) on PATH for GIF encoding; ``--no-gif`` renders
frames and the static contact-sheet frame only.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/spanwise_yz_plan_flyover.py

Provenance note: this script is a reconstruction. An earlier misspelled GIF
(``spawnise_YZ_plan_flyover.gif``) sat orphaned with no generator in the repo.
This rebuilds the figure from the shared spanwise station loader and the same
``resolve_web_anchors`` / ``inner_y_at_x`` calls that ``spanwise_web_anchors.py``
probes with. The deck tracks ``spanwise_yz_plan_flyover.gif`` only.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_REPO = _ensure_mod.ensure_repo_path()

from spanwise_web_anchors_shared import build_stations as _build_spanwise_stations  # noqa: E402
from blade_structure_precompute.module.section_selection.processing.adapters.scale_adapter import (  # noqa: E402
    airfoil_sdf_from_station,
)
from blade_structure_utilities.geometry_utilities.section_geometry import (  # noqa: E402
    MultiCellSection,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sections.web_layout import (  # noqa: E402
    inner_y_at_x,
    resolve_web_anchors,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (  # noqa: E402
    apply_presentation_style,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot import (  # noqa: E402
    plot_section_hybrid,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.plot_section import (  # noqa: E402
    apply_structural_section_axes,
)

import feature_flythrough_render as _render  # noqa: E402
from feature_flythrough_specs import DEFAULT_PROFILE, PROFILES  # noqa: E402

logging.disable(logging.INFO)

ROTOR = "new_3.2m_rotor"
MATERIALS = {"skin": 17, "spar": 17, "web": 17}
LAYOUT_KEY = "5A-F"

# Reference gauges. The point of this figure is where the webs sit, not the sizing.
T_SKIN = 0.001
T_WEB = 0.0015

# MultiCellSection's own default, and the frame origin: Y_B is measured from the
# pitch axis, so a web at this chord fraction sits at Y_B = 0 and is the one web
# whose x/c does not move under root anchoring.
PIVOT_FRACTION = 1.0 / 3.0

# LE-most to TE-most. Picked to stay distinguishable in a GIF's 256-colour palette.
WEB_COLOURS = ("#43A047", "#00897B", "#8BC34A", "#2E7D32", "#1B5E20")

PLANFORM_FILL = "#90CAF9"
PLANFORM_EDGE = "#42A5F5"
ROOT_BAND = "#ECEFF1"
STATION_MARKER = "#E64A19"
DROP_MARKER = "#2E7D32"

# The loader reports each station as zL in [0, 1] along the blade, not as a radius,
# so the radial axis needs the rotor's root cut-out. 0.125 of tip radius reproduces
# the 3.2 m rotor's 0.2 m root (R = 1.6 m, blade length 1.4 m), which is what the
# Scale stage reports for these benchmark rotors.
ROOT_RADIUS_FRACTION = 0.125

OUT_DIR = Path(__file__).resolve().parents[1] / "output" / "spanwise_yz_plan_flyover"


def _web_colour(index: int) -> str:
    return WEB_COLOURS[index % len(WEB_COLOURS)]


def build_span(
    *,
    rotor: str,
    layout_key: str,
    tip_radius_m: float,
    skin_thickness_m: float = T_SKIN,
) -> dict[str, Any]:
    """Resolve both web-placement models at every station along the span.

    Returns the per-station tracks the animation draws, plus the drop events that
    mark where root anchoring loses a web.
    """
    stations = _build_spanwise_stations(
        repo_root=_REPO,
        rotor=rotor,
        materials=MATERIALS,
        structural_system_key=layout_key,
    )
    layout = resolve_system_type(layout_key)
    fracs = tuple(float(f) for f in layout.web_chord_fracs)
    n_webs = len(fracs)

    chord_root = float(stations[0]["chord"])
    # Root-anchored webs keep the Y_B offset they were given at the root, which is
    # what "anchored" means here — not a fixed distance from the leading edge.
    y_root = [(f - PIVOT_FRACTION) * chord_root for f in fracs]

    z_root = float(tip_radius_m) * ROOT_RADIUS_FRACTION
    blade_length = float(tip_radius_m) - z_root

    rows: list[dict[str, Any]] = []
    for station in stations:
        chord = float(station["chord"])
        airfoil = airfoil_sdf_from_station(station["coords"], station["metadata"])

        # --- station-relative: fracs against the *local* chord ---
        xs_requested = [f * float(airfoil.chord) for f in fracs]
        try:
            anchors = resolve_web_anchors(airfoil, T_SKIN, xs_requested)
            xs_station = [float(a[0]) for a in anchors]
        except Exception:
            xs_station = []

        # --- root-anchored: the same Y_B offsets, tested against this cavity ---
        xs_root: list[float | None] = []
        for offset in y_root:
            x = offset + PIVOT_FRACTION * float(airfoil.chord)
            if not (0.0 < x < float(airfoil.chord)):
                xs_root.append(None)
                continue
            try:
                inner_y_at_x(airfoil, T_SKIN, float(x))
            except Exception:
                xs_root.append(None)
                continue
            xs_root.append(float(x))

        z_m = z_root + float(station["zL"]) * blade_length
        rows.append(
            {
                "idx": int(station["idx"]),
                "zL": float(station["zL"]),
                "z_m": z_m,
                "r_over_R": z_m / float(tip_radius_m),
                "chord_m": chord,
                "twist_deg": float(station["metadata"]["twist_deg"]),
                "airfoil": airfoil,
                "xs_station_relative_m": xs_station,
                "xs_root_anchored_m": xs_root,
            }
        )

    # First station at which each root-anchored web can no longer be placed.
    drops: list[dict[str, Any]] = []
    for web in range(n_webs):
        for row in rows:
            if row["xs_root_anchored_m"][web] is None:
                drops.append(
                    {
                        "web_index": web,
                        "z_m": row["z_m"],
                        "r_over_R": row["r_over_R"],
                        "x_over_c": float(
                            np.clip(
                                (y_root[web] + PIVOT_FRACTION * row["chord_m"]) / row["chord_m"],
                                0.0,
                                1.0,
                            )
                        ),
                        "y_b_m": y_root[web],
                    }
                )
                break

    return {
        "rows": rows,
        "layout": layout,
        "fracs": fracs,
        "n_webs": n_webs,
        "chord_root_m": chord_root,
        "y_root_m": y_root,
        "z_root_m": z_root,
        "tip_radius_m": float(tip_radius_m),
        "drops": drops,
    }


def _build_section(airfoil, xs: list[float], twist_rad: float):
    """A family-A multicell at the given web stations, or None if it cannot build."""
    if not xs:
        return None
    try:
        return MultiCellSection(
            airfoil,
            web_x_positions=list(xs),
            web_thickness=T_WEB,
            skin_thickness=T_SKIN,
            twist_angle=float(twist_rad),
            web_alignment="flapwise",
            structural_family="A",
            core_enabled=True,
            pitch_fraction_of_chord_from_le=PIVOT_FRACTION,
        )
    except Exception:
        return None


def _draw_section(ax, section, title: str) -> None:
    if section is None:
        ax.text(0.5, 0.5, "no anchorable web", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="#B71C1C")
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    plot_section_hybrid(
        section,
        residual_mesh=None,
        ax=ax,
        alpha=0.55,
        title=title,
        show_legend=False,
    )
    # plot_section_hybrid already fixes the aspect; re-asserting it here only makes
    # matplotlib discard one of the limit pairs and warn.
    apply_structural_section_axes(ax)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)


def _draw_plan(ax, span: dict[str, Any], *, anchored: bool, current) -> None:
    """Plan view looking down on the blade: web tracks inside the planform."""
    rows = span["rows"]
    z = np.array([r["z_m"] for r in rows], dtype=float)
    chord = np.array([r["chord_m"] for r in rows], dtype=float)

    # Leading and trailing edge in Y_B: the chord runs from -pivot*c to (1-pivot)*c.
    y_le = -PIVOT_FRACTION * chord
    y_te = (1.0 - PIVOT_FRACTION) * chord
    ax.fill_between(z, y_le, y_te, color=PLANFORM_FILL, alpha=0.55, lw=0.0, zorder=1)
    ax.plot(z, y_le, color=PLANFORM_EDGE, lw=1.1, zorder=2)
    ax.plot(z, y_te, color=PLANFORM_EDGE, lw=1.1, zorder=2)
    ax.axhline(0.0, color="#90A4AE", lw=0.8, ls=":", zorder=1)

    ax.axvspan(0.0, span["z_root_m"], color=ROOT_BAND, zorder=0)

    for web in range(span["n_webs"]):
        colour = _web_colour(web)
        if anchored:
            track = np.full_like(z, span["y_root_m"][web])
        else:
            track = (span["fracs"][web] - PIVOT_FRACTION) * chord
        ax.plot(z, track, color=colour, lw=1.7, zorder=3)

    if anchored:
        for drop in span["drops"]:
            ax.plot(
                drop["z_m"],
                drop["y_b_m"],
                marker="X",
                ms=9,
                color=DROP_MARKER,
                markeredgecolor="#1B5E20",
                zorder=6,
            )

    ax.axvline(current["z_m"], color=STATION_MARKER, lw=1.8, zorder=7)
    ax.set_xlim(0.0, span["tip_radius_m"])
    ax.set_xlabel(r"$Z_R$  (radial) [m]", fontsize=8)
    ax.set_ylabel(r"$Y_B$ (edgewise) [m]", fontsize=8)
    ax.text(0.995, 0.92, "TE", transform=ax.transAxes, ha="right", fontsize=7, color="#37474F")
    ax.text(0.995, 0.06, "LE", transform=ax.transAxes, ha="right", fontsize=7, color="#37474F")
    ax.set_title(
        "plan ($Z_R$-$Y_B$) — anchored webs" if anchored else "plan ($Z_R$-$Y_B$) — webs hold x/c",
        fontsize=9,
    )
    ax.grid(True, alpha=0.25)
    ax.tick_params(labelsize=7)


def _draw_chord_fraction(ax, span: dict[str, Any], current) -> None:
    """Where each web sits in x/c: the target, and what anchoring actually gives."""
    rows = span["rows"]
    z = np.array([r["z_m"] for r in rows], dtype=float)
    chord = np.array([r["chord_m"] for r in rows], dtype=float)

    ax.axvspan(0.0, span["z_root_m"], color=ROOT_BAND, zorder=0)
    for level, label in ((1.0, "trailing edge  (x/c = 1)"), (0.0, "leading edge  (x/c = 0)")):
        ax.axhline(level, color="#607D8B", lw=0.9, ls=":", zorder=1)
        ax.text(0.012, level + (-0.09 if level else 0.04), label, transform=ax.get_yaxis_transform(),
                fontsize=7, color="#455A64")

    drop_z = {d["web_index"]: d["z_m"] for d in span["drops"]}
    for web in range(span["n_webs"]):
        colour = _web_colour(web)
        # Dashed: the station-relative target, constant by construction. Drawn over
        # the blade only — there is no structure inboard of the root cut.
        ax.plot(
            [span["z_root_m"], span["tip_radius_m"]],
            [span["fracs"][web]] * 2,
            color=colour, lw=1.2, ls="--", alpha=0.85, zorder=2,
        )
        # Solid: root anchoring, which drifts as the chord changes. Solid while the
        # web can still be anchored, faint past its drop — the track is still a
        # well-defined position, but no web is built there any more.
        track = np.clip((span["y_root_m"][web] + PIVOT_FRACTION * chord) / chord, 0.0, 1.0)
        cut = drop_z.get(web)
        alive = z <= cut if cut is not None else np.ones_like(z, dtype=bool)
        ax.plot(z[alive], track[alive], color=colour, lw=1.7, zorder=3)
        if cut is not None:
            ax.plot(z[~alive], track[~alive], color=colour, lw=1.0, ls=":", alpha=0.45, zorder=3)

    for drop in span["drops"]:
        ax.plot(
            drop["z_m"],
            drop["x_over_c"],
            marker="X",
            ms=9,
            color=DROP_MARKER,
            markeredgecolor="#1B5E20",
            zorder=6,
        )

    ax.axvline(current["z_m"], color=STATION_MARKER, lw=1.8, zorder=7)
    ax.set_xlim(0.0, span["tip_radius_m"])
    ax.set_ylim(-0.05, 1.08)
    ax.set_xlabel(r"$Z_R$  (radial) [m]", fontsize=8)
    ax.set_ylabel("web station  x / c  [-]", fontsize=8)
    ax.grid(True, alpha=0.25)
    ax.tick_params(labelsize=7)


def _footer(fig, span: dict[str, Any], current) -> None:
    n_station = len(current["xs_station_relative_m"])
    n_anchored = sum(1 for x in current["xs_root_anchored_m"] if x is not None)
    nominal = span["n_webs"]
    family = str(span["layout"].structural_family)
    suffix = span["layout"].key.split("-")[-1]

    fig.text(
        0.02,
        0.055,
        f"r/R = {current['r_over_R']:.3f}    chord = {current['chord_m'] * 1e3:7.1f} mm"
        f"    twist = {current['twist_deg']:6.2f}°",
        fontsize=9,
        family="monospace",
    )
    fig.text(
        0.02,
        0.018,
        f"station-relative:  {n_station}{family}-{suffix}, {n_station} webs,"
        f" {nominal - n_station} dropped"
        f"      root-anchored:  {n_anchored}{family}-{suffix}, {n_anchored} webs,"
        f" {nominal - n_anchored} dropped",
        fontsize=9,
        family="monospace",
    )


def render_span_still(
    span: dict[str, Any], current, out_path: Path, *, presentation: bool
) -> None:
    """Clean two-panel plan + chord-fraction figure for the PDF export still."""
    fig = plt.figure(figsize=(13.0, 7.0))
    gs = fig.add_gridspec(
        2, 2, hspace=0.45, wspace=0.12, left=0.08, right=0.97, top=0.95, bottom=0.08
    )
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1])
    ax_chord = fig.add_subplot(gs[1, :])

    _draw_plan(ax_left, span, anchored=False, current=current)
    _draw_plan(ax_right, span, anchored=True, current=current)
    _draw_chord_fraction(ax_chord, span, current)

    ax_left.set_title("chord-relative — all webs survive", fontsize=12)
    ax_right.set_title("root-anchored — webs drop at ×", fontsize=12)
    for ax in (ax_left, ax_right):
        ax.set_xlabel("Span [m]", fontsize=12)
        ax.set_ylabel("Edgewise offset [m]", fontsize=12)
        ax.tick_params(labelsize=10)
    ax_chord.set_xlabel("Span [m]", fontsize=12)
    ax_chord.set_ylabel("Web position [x/c]", fontsize=12)
    ax_chord.tick_params(labelsize=10)

    save_figure(fig, out_path, presentation=presentation)
    plt.close(fig)
    print(f"Wrote {out_path}")


def render_frame(span: dict[str, Any], current, out_path: Path, *, dpi: int) -> Path:
    fig = plt.figure(figsize=(12.8, 7.91))
    gs = fig.add_gridspec(
        3, 2, height_ratios=(1.0, 1.05, 1.05), hspace=0.62, wspace=0.18,
        left=0.07, right=0.98, top=0.96, bottom=0.13,
    )

    twist_rad = float(np.radians(current["twist_deg"]))
    _draw_section(
        fig.add_subplot(gs[0, 0]),
        _build_section(current["airfoil"], current["xs_station_relative_m"], twist_rad),
        "station-relative — webs at fixed x/c",
    )
    _draw_section(
        fig.add_subplot(gs[0, 1]),
        _build_section(
            current["airfoil"],
            [x for x in current["xs_root_anchored_m"] if x is not None],
            twist_rad,
        ),
        r"root-anchored — webs at fixed $Y_B$ [m]",
    )

    _draw_plan(fig.add_subplot(gs[1, 0]), span, anchored=False, current=current)
    _draw_plan(fig.add_subplot(gs[1, 1]), span, anchored=True, current=current)
    _draw_chord_fraction(fig.add_subplot(gs[2, :]), span, current)

    _footer(fig, span, current)
    fig.savefig(out_path, dpi=dpi, facecolor="white")
    plt.close(fig)
    return out_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rotor", default=ROTOR)
    parser.add_argument("--layout", default=LAYOUT_KEY)
    parser.add_argument(
        "--tip-radius-m",
        type=float,
        default=1.6,
        help="Rotor tip radius, for the radial axis (default: 1.6, the 3.2 m rotor)",
    )
    parser.add_argument("--profile", default=DEFAULT_PROFILE, choices=sorted(PROFILES))
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--no-gif", action="store_true", help="Render frames only.")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )
    args = parser.parse_args(argv)

    apply_presentation_style(bool(args.presentation))
    profile = PROFILES[args.profile]
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    span = build_span(
        rotor=args.rotor,
        layout_key=args.layout,
        tip_radius_m=float(args.tip_radius_m),
    )
    rows = span["rows"]
    print(f"Resolved {len(rows)} stations; {len(span['drops'])} root-anchored web drop(s)")

    paths = _render.prepare_scenario_paths(
        out_dir, "spanwise_yz_plan_flyover", "spanwise_yz_plan_flyover.gif"
    )
    _render.reset_frame_dir(paths.frame_dir)

    indices = _render.select_frame_indices(len(rows), args.max_frames)
    frames: list[Path] = []
    for n, idx in enumerate(indices):
        frames.append(
            render_frame(
                span,
                rows[idx],
                _render.frame_path(paths.frame_dir, n),
                dpi=profile.dpi,
            )
        )
        if (n + 1) % 10 == 0 or n + 1 == len(indices):
            print(f"  frame {n + 1} / {len(indices)}")

    payload: dict[str, Any] = {
        "rotor": args.rotor,
        "layout": args.layout,
        "tip_radius_m": float(args.tip_radius_m),
        "root_radius_m": span["z_root_m"],
        "chord_root_m": span["chord_root_m"],
        "n_stations": len(rows),
        "n_frames": len(frames),
        "web_chord_fracs": list(span["fracs"]),
        "y_root_m": span["y_root_m"],
        "drops": span["drops"],
        "profile": profile.name,
    }

    if not args.no_gif:
        held = _render.apply_holds(frames, head=profile.hold_head, tail=profile.hold_tail)
        _render.encode_gif(held, paths.gif_path, profile)
        payload["gif"] = _render.validate_gif(paths.gif_path)
        print(f"Wrote {paths.gif_path}")

    paths.metadata_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {paths.metadata_path}")

    still_row = rows[int(0.8 * len(rows))]
    render_span_still(
        span,
        still_row,
        out_dir / "spanwise_yz_plan_flyover_still.png",
        presentation=bool(args.presentation),
    )
    for drop in span["drops"]:
        print(
            f"  web {drop['web_index']} dropped at r/R={drop['r_over_R']:.3f}"
            f"  (x/c={drop['x_over_c']:.3f})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
