"""
visualisation/demos/design_vector_morph.py
==========================================
One design vector in; solids, centrelines and a conforming shell mesh out.

Every other demo holds the design fixed and takes the pipeline apart. This one
does the opposite: it holds the pipeline fixed and moves the design vector,
frame by frame, so the whole chain can be watched re-deriving. Three panels,
all sharing one window:

1. **CSG solids** - the named implicit regions.
2. **Junction-prepared centrelines** - after ``topology.prepare_shell_midline_junctions``.
3. **Q4 shell mesh** - after ``section_mesh.SectionMeshRunner``.

and beneath them the vector actually driving the frame, with the component
currently moving picked out.

The path walks three kinds of change, because they are the three that a sizing
loop makes and they are usually claimed to be different in difficulty:

* a **continuous gauge sweep** - skin thickness, the dominant design variable;
* **discrete topology changes** - ``2D-F`` to ``3B-F`` to ``4C-F`` to ``5A-F``,
  which is a different web count *and* a different spar-cap family each time;
* a **rigid-body change** - station twist.

Nothing is re-modelled between frames and no mesh is repaired. Each frame is a
fresh evaluation of the same expression tree with different parameters, which
is the entire argument for carrying geometry as a program.

Writes under ``visualisation/output/design_vector_morph/``:

* ``design_vector_morph.gif`` - the animation.
* ``design_vector_morph_stages.png`` - a static contact sheet of the keyframes,
  for the PDF/PPTX export where the GIF cannot play.
* ``design_vector_morph.json`` - the sampled vector and per-frame mesh counts.

Requires ImageMagick (``magick``) on PATH for GIF encoding; ``--no-gif`` renders
frames and the contact sheet only.

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/design_vector_morph.py
"""

from __future__ import annotations

import argparse
import json
import sys
from argparse import Namespace
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

_demos_dir = Path(__file__).resolve().parent
if str(_demos_dir) not in sys.path:
    sys.path.insert(0, str(_demos_dir))

# geometry_to_mesh_flythrough already owns the "build one station and draw its
# stages" machinery, including the repo-path bootstrap. Reusing it keeps the two
# animations showing literally the same pipeline rather than two look-alikes.
import feature_flythrough_render as _render  # noqa: E402
import feature_flythrough_specs as _specs  # noqa: E402
import geometry_to_mesh_flythrough as _pipeline  # noqa: E402

from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (  # noqa: E402
    apply_presentation_style,
    demo_output_dir,
    resolve_save_dpi,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (  # noqa: E402
    INK,
    MUTED,
    NAVY,
    RED,
    SURFACE,
    TEAL,
)

#: Order the vector is displayed in, with the label used on the readout.
_VECTOR_FIELDS = (
    ("key", "layout"),
    ("skin_thickness", "t_skin"),
    ("spar_cap_thickness", "t_cap"),
    ("shear_web_thickness", "t_web"),
    ("twist_deg", "twist"),
)


@dataclass(frozen=True)
class Keyframe:
    """One corner of the design-vector path."""

    key: str
    skin_thickness: float
    spar_cap_thickness: float
    shear_web_thickness: float
    twist_deg: float
    #: Which vector component this leg is moving, for the readout highlight.
    moving: str
    #: What the leg is demonstrating, shown as the frame's caption.
    caption: str
    #: Frames spent travelling from the previous keyframe (0 for the first).
    steps: int = 8


def default_path() -> list[Keyframe]:
    """The design-vector path: gauge sweep, then topology, then twist.

    Gauges are pulled back to a mid value before the topology legs so that every
    layout in the cascade is comfortably feasible; a web that fails to anchor
    raises rather than rendering, and that behaviour has its own demo.
    """
    return [
        Keyframe("2D-F", 0.004, 0.012, 0.006, 0.0, "skin_thickness",
                 "one continuous gauge, one unchanged tree", steps=0),
        Keyframe("2D-F", 0.014, 0.012, 0.006, 0.0, "skin_thickness",
                 "skin gauge 4 -> 14 mm: no remesh, no repair", steps=10),
        Keyframe("2D-F", 0.006, 0.014, 0.006, 0.0, "spar_cap_thickness",
                 "cap gauge moves independently", steps=6),
        Keyframe("3B-F", 0.006, 0.014, 0.006, 0.0, "key",
                 "topology change: 2 webs -> 3, cap family D -> B", steps=4),
        Keyframe("4C-F", 0.006, 0.014, 0.006, 0.0, "key",
                 "topology change: 3 webs -> 4, per-web caps", steps=4),
        Keyframe("5A-F", 0.006, 0.014, 0.006, 0.0, "key",
                 "topology change: 5 webs, no caps", steps=4),
        Keyframe("5A-F", 0.006, 0.014, 0.008, 15.0, "twist_deg",
                 "and the whole station twists", steps=8),
    ]


def sample_path(path: list[Keyframe]) -> list[tuple[dict, str, str]]:
    """Expand keyframes into per-frame ``(vector, moving, caption)`` samples.

    Continuous components are interpolated with a smoothstep so the animation
    eases rather than ticking; ``key`` is discrete and snaps at the leg start,
    because there is no such thing as half a shear web.
    """
    samples: list[tuple[dict, str, str]] = []
    first = path[0]
    samples.append(
        (
            {
                "key": first.key,
                "skin_thickness": first.skin_thickness,
                "spar_cap_thickness": first.spar_cap_thickness,
                "shear_web_thickness": first.shear_web_thickness,
                "twist_deg": first.twist_deg,
            },
            first.moving,
            first.caption,
        )
    )

    for previous, current in zip(path, path[1:]):
        for step in range(1, int(current.steps) + 1):
            raw = step / float(current.steps)
            eased = raw * raw * (3.0 - 2.0 * raw)
            samples.append(
                (
                    {
                        # The new layout applies from the leg's first frame.
                        "key": current.key,
                        "skin_thickness": previous.skin_thickness
                        + eased * (current.skin_thickness - previous.skin_thickness),
                        "spar_cap_thickness": previous.spar_cap_thickness
                        + eased
                        * (current.spar_cap_thickness - previous.spar_cap_thickness),
                        "shear_web_thickness": previous.shear_web_thickness
                        + eased
                        * (current.shear_web_thickness - previous.shear_web_thickness),
                        "twist_deg": previous.twist_deg
                        + eased * (current.twist_deg - previous.twist_deg),
                    },
                    current.moving,
                    current.caption,
                )
            )
    return samples


def _context_for(vector: dict, args: argparse.Namespace):
    """Build the full geometry -> midlines -> mesh chain for one design vector."""
    return _pipeline.build_context(
        Namespace(
            key=str(vector["key"]),
            naca=str(args.naca),
            chord=float(args.chord),
            twist_deg=float(vector["twist_deg"]),
            nx=int(args.nx),
            ny=int(args.ny),
            skin_thickness=float(vector["skin_thickness"]),
            spar_cap_thickness=float(vector["spar_cap_thickness"]),
            shear_web_thickness=float(vector["shear_web_thickness"]),
            element_length_frac=float(args.element_length_frac),
            zoom_half_width_frac=0.035,
            frame_aspect=float(args.panel_aspect),
        )
    )


def _draw_readout(figure, vector: dict, moving: str) -> None:
    """The design vector along the bottom, with the moving component picked out."""
    left = 0.055
    width = (1.0 - 2.0 * left) / len(_VECTOR_FIELDS)
    for index, (field, label) in enumerate(_VECTOR_FIELDS):
        value = vector[field]
        if field == "key":
            text = str(value)
        elif field == "twist_deg":
            text = f"{float(value):.1f}°"
        else:
            text = f"{float(value) * 1e3:.1f} mm"
        active = field == moving
        figure.text(
            left + (index + 0.5) * width,
            0.175,
            label,
            ha="center",
            va="center",
            fontsize=10.5,
            color=RED if active else MUTED,
            fontweight="bold" if active else "normal",
        )
        figure.text(
            left + (index + 0.5) * width,
            0.098,
            text,
            ha="center",
            va="center",
            fontsize=16 if active else 14,
            color=RED if active else INK,
            fontweight="bold",
            bbox={
                "facecolor": "#FDE8ED" if active else SURFACE,
                "edgecolor": RED if active else "none",
                "linewidth": 1.4 if active else 0.0,
                "boxstyle": "round,pad=0.42",
            },
        )


_PANEL_TITLES = (
    "1 · implicit solids",
    "2 · junction-prepared centrelines",
    "3 · conforming Q4 shell mesh",
)


def render_frame(
    context, vector: dict, moving: str, caption: str, path: Path, dpi: int, window
) -> None:
    """Three panels on one shared window, plus the vector driving them."""
    # Wide and short: a station is a thin lens, and a taller frame just adds
    # dead space above and below it once the axes are set to equal aspect.
    figure, axes = plt.subplots(1, 3, figsize=(16.0, 5.4))

    _pipeline.draw_csg(axes[0], context, 1.0)
    _pipeline.draw_strips(
        axes[1], context.prepared_strips, 1.0, show_endpoints=True, linewidth=2.0, zorder=4
    )
    _pipeline.draw_mesh(axes[2], context, 1.0)

    for axis, title in zip(axes, _PANEL_TITLES):
        axis.set_xlim(window[0], window[1])
        axis.set_ylim(window[2], window[3])
        axis.set_aspect("equal", adjustable="box")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(title, fontsize=12.5, color=NAVY, pad=8)
        for spine in axis.spines.values():
            spine.set_color("#d7dce3")

    summary = context.mesh.summary()
    axes[2].text(
        0.98,
        0.04,
        f"{int(summary['n_total_elements'])} elements · "
        f"{int(summary['n_total_nodes'])} nodes · "
        f"{int(summary['n_clusters'])} shared junctions",
        transform=axes[2].transAxes,
        ha="right",
        fontsize=10,
        color=TEAL,
    )

    figure.suptitle(
        "One design vector, three products — nothing re-modelled between frames",
        fontsize=15,
        y=0.985,
    )
    # Clear of the panel titles, which sit just above the axes band.
    figure.text(0.5, 0.912, caption, ha="center", fontsize=11.5, color=MUTED, style="italic")
    _draw_readout(figure, vector, moving)
    figure.subplots_adjust(left=0.02, right=0.98, top=0.835, bottom=0.285, wspace=0.05)

    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, facecolor="white")
    plt.close(figure)


def render_contact_sheet(
    entries: list[tuple], path: Path, *, presentation: bool, window
) -> None:
    """Keyframes as a static grid, for the export where the GIF cannot play."""
    columns = len(entries)
    # Sized to land in a full-width slide panel of roughly 2.5:1. One column per
    # leg made seven narrow columns whose sections were too small to read.
    figure, axes = plt.subplots(
        2, columns, figsize=(3.0 * columns, 4.2), gridspec_kw={"hspace": 0.05, "wspace": 0.04}
    )
    axes = np.atleast_2d(axes)

    for column, (context, vector, caption, _moving) in enumerate(entries):
        _pipeline.draw_csg(axes[0, column], context, 1.0)
        _pipeline.draw_mesh(axes[1, column], context, 1.0)
        for row in (0, 1):
            axis = axes[row, column]
            axis.set_xlim(window[0], window[1])
            axis.set_ylim(window[2], window[3])
            axis.set_aspect("equal", adjustable="box")
            axis.set_xticks([])
            axis.set_yticks([])
            for spine in axis.spines.values():
                spine.set_color("#d7dce3")
        axes[0, column].set_title(
            f"{vector['key']}\n"
            f"t_skin {vector['skin_thickness'] * 1e3:.0f} mm · "
            f"twist {vector['twist_deg']:.0f}°",
            fontsize=12.5,
            color=NAVY,
        )
        summary = context.mesh.summary()
        axes[1, column].set_xlabel(
            f"{int(summary['n_total_elements'])} elements", fontsize=12, color=TEAL
        )

    axes[0, 0].set_ylabel("implicit solids", fontsize=13, color=NAVY)
    axes[1, 0].set_ylabel("Q4 shell mesh", fontsize=13, color=NAVY)
    save_figure(figure, path, presentation=presentation)
    plt.close(figure)
    print(f"Wrote {path}")


def render_still_3plus2(
    entries: list[tuple], path: Path, *, presentation: bool, window
) -> None:
    """3+2 two-row contact sheet for PDF/PPTX export — aspect fits a full-width slide panel.

    Shows CSG solids only (no mesh row) so the sections fill more of each cell.
    The last slot in row 2 carries a text annotation when there are fewer than
    3 entries in that row.
    """
    n_per_row = 3
    figure = plt.figure(figsize=(14.0, 6.2))
    gs = figure.add_gridspec(
        2, n_per_row,
        hspace=0.06, wspace=0.04,
        left=0.04, right=0.98, top=0.87, bottom=0.04,
    )

    def _draw_cell(row: int, col: int, context, vector) -> None:
        axis = figure.add_subplot(gs[row, col])
        _pipeline.draw_csg(axis, context, 1.0)
        axis.set_xlim(window[0], window[1])
        axis.set_ylim(window[2], window[3])
        axis.set_aspect("equal", adjustable="box")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_color("#d7dce3")
        skin_mm = vector["skin_thickness"] * 1e3
        cap_mm = vector["spar_cap_thickness"] * 1e3
        web_mm = vector["shear_web_thickness"] * 1e3
        twist = vector["twist_deg"]
        title = (
            f"{vector['key']}\n"
            f"t_skin {skin_mm:.0f}  t_cap {cap_mm:.0f}  t_web {web_mm:.0f} mm"
        )
        if twist != 0.0:
            title += f"  twist {twist:.0f}°"
        axis.set_title(title, fontsize=10.5, color=NAVY)

    row1_entries = entries[:n_per_row]
    row2_entries = entries[n_per_row:]

    for col, (context, vector, _caption, _moving) in enumerate(row1_entries):
        _draw_cell(0, col, context, vector)

    for col, (context, vector, _caption, _moving) in enumerate(row2_entries):
        _draw_cell(1, col, context, vector)

    for col in range(len(row2_entries), n_per_row):
        axis = figure.add_subplot(gs[1, col])
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
        axis.text(
            0.5, 0.5,
            "same constructor\n·\ndifferent numbers",
            transform=axis.transAxes,
            ha="center", va="center",
            fontsize=13, color=NAVY, fontweight="bold",
            linespacing=1.8,
        )

    figure.suptitle(
        "Design vector: gauge sweep → topology change → twist"
        " — one call each, no remesh",
        fontsize=12.5, color=NAVY, y=0.96,
    )
    save_figure(figure, path, presentation=presentation)
    plt.close(figure)
    print(f"Wrote {path}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--naca", default="2412")
    parser.add_argument("--chord", type=float, default=1.0)
    # Deliberately coarser than the single-station fly-through: this animation
    # rebuilds the entire chain per frame, so the field grid is the cost driver.
    parser.add_argument("--nx", type=int, default=460)
    parser.add_argument("--ny", type=int, default=210)
    parser.add_argument("--element-length-frac", type=float, default=0.016)
    parser.add_argument("--panel-aspect", type=float, default=1.70)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument(
        "--still-columns",
        type=int,
        default=5,
        help="Keyframe columns kept in the static contact sheet.",
    )
    parser.add_argument("--no-gif", action="store_true", help="Render frames only.")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--presentation",
        action="store_true",
        help="Larger fonts and higher DPI for slides/decks.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)
    dpi = resolve_save_dpi(presentation=presentation)

    out_dir = Path(args.out_dir or demo_output_dir("design_vector_morph"))
    out_dir.mkdir(parents=True, exist_ok=True)

    path = default_path()
    samples = sample_path(path)
    paths = _render.prepare_scenario_paths(
        out_dir, "design_vector_morph", "design_vector_morph.gif"
    )
    _render.reset_frame_dir(paths.frame_dir)
    selected = _render.select_frame_indices(len(samples), args.max_frames)

    records = []
    frame_files: list[Path] = []
    keyframe_entries: list[tuple] = []
    keyframe_captions = {frame.caption for frame in path}
    seen_captions: set[str] = set()

    # One window for the whole sweep. build_context sizes its window from the
    # station it just built, so using that per frame would make the camera
    # breathe as the section changes — the one thing a morph must not do. The
    # centre comes from the first station (plot coordinates put the pitch axis
    # at the origin, so the chord is *not* symmetric about it) and the extent is
    # then held fixed for every later frame.
    window = None

    for output_index, sample_index in enumerate(selected):
        vector, moving, caption = samples[sample_index]
        context = _context_for(vector, args)
        if window is None:
            half_width = 0.62 * float(args.chord)
            half_height = half_width / float(args.panel_aspect)
            centre = context.wide_center
            window = (
                centre[0] - half_width,
                centre[0] + half_width,
                centre[1] - half_height,
                centre[1] + half_height,
            )
        frame_file = _render.frame_path(paths.frame_dir, output_index)
        render_frame(context, vector, moving, caption, frame_file, dpi, window)
        frame_files.append(frame_file)

        summary = context.mesh.summary()
        records.append(
            {
                "frame": output_index,
                "vector": {k: (v if isinstance(v, str) else float(v)) for k, v in vector.items()},
                "elements": int(summary["n_total_elements"]),
                "nodes": int(summary["n_total_nodes"]),
                "junctions": int(summary["n_clusters"]),
            }
        )

        # One contact-sheet column per leg, taken at the frame that completes it.
        if caption in keyframe_captions and caption not in seen_captions:
            is_last_of_leg = (
                sample_index + 1 >= len(samples)
                or samples[sample_index + 1][2] != caption
            )
            if is_last_of_leg:
                seen_captions.add(caption)
                keyframe_entries.append((context, vector, caption, moving))

        print(f"  frame {output_index + 1} / {len(selected)}  {vector['key']}")

    # One column per leg is too many to read on a slide. Rather than subsample by
    # index — which silently dropped the 14 mm gauge extreme, the most visible
    # frame in the sweep — keep the baseline plus where each *kind* of change
    # ends up: gauge, cap gauge, topology, twist.
    if len(keyframe_entries) > int(args.still_columns):
        picks = [0]
        for index in range(1, len(keyframe_entries)):
            is_last_of_kind = (
                index == len(keyframe_entries) - 1
                or keyframe_entries[index + 1][3] != keyframe_entries[index][3]
            )
            if is_last_of_kind:
                picks.append(index)
        keyframe_entries = [keyframe_entries[index] for index in picks]

    if keyframe_entries:
        render_still_3plus2(
            keyframe_entries,
            out_dir / "design_vector_morph_stages.png",
            presentation=presentation,
            window=window,
        )

    payload = {
        "naca": str(args.naca),
        "chord_m": float(args.chord),
        "frames": records,
    }
    if not args.no_gif:
        profile = _specs.PROFILES["balanced"]
        # Hold the first and last frames so the loop reads as a sweep with a
        # beginning and an end rather than a continuous churn.
        held = _render.apply_holds(
            frame_files, head=profile.hold_head, tail=profile.hold_tail
        )
        _render.encode_gif(held, paths.gif_path, profile)
        payload["gif"] = _render.validate_gif(paths.gif_path)
        print(f"Wrote {paths.gif_path}")

    out_json = out_dir / "design_vector_morph.json"
    _render.write_metadata(out_json, payload)
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
