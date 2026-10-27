"""
visualisation/demos/bug_log_sg009_demo.py
========================================
SG-009 schematic: making a junction exact destroys the evidence used to find it.

Three panels:

1. A raw web midline crosses the skin midline (strict-interior detector finds it).
2. The skin is split and the web trimmed at the crossing: one shared node.
3. The same detector returns nothing (touching is not an interior crossing).

Writes under ``visualisation/output/bug_log_sg009/``:

* ``bug_log_sg009_schematic.png``
* ``bug_log_sg009.json``

Run from the repository root::

    python src/blade_geometry_model/section_geometry/visualisation/demos/bug_log_sg009_demo.py --presentation
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_geometry_model.section_geometry.midlines.ring_geometry import (
    intersect_open_polylines_2d,
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


def _crossing_polylines() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A skin midline and a web midline that overshoots through it.

    Schematic geometry in the shape of a real web-skin junction: a gently
    curved skin sampled at a handful of nodes, and a web whose raw end runs
    past the skin, so the two strictly cross.
    """
    skin_x = np.linspace(-1.0, 1.0, 6)
    skin = np.column_stack([skin_x, 0.55 - 0.12 * skin_x**2])
    web = np.array([[0.08, -0.95], [0.16, 0.85]], dtype=float)
    hits = intersect_open_polylines_2d(skin, web)
    if len(hits) != 1:
        raise RuntimeError(f"Schematic must cross exactly once; got {len(hits)} hits.")
    shared = np.asarray(hits[0][4], dtype=float)
    return skin, web, shared


def _insert_shared_node(poly: np.ndarray, shared: np.ndarray) -> np.ndarray:
    """Insert ``shared`` into ``poly`` at its position along the polyline."""
    projections = [
        np.dot(shared - poly[i], poly[i + 1] - poly[i]) / np.dot(poly[i + 1] - poly[i], poly[i + 1] - poly[i])
        for i in range(len(poly) - 1)
    ]
    for index, fraction in enumerate(projections):
        if 0.0 <= fraction <= 1.0:
            return np.vstack([poly[: index + 1], shared, poly[index + 1 :]])
    raise RuntimeError("Shared node does not lie on the polyline.")


def _trim_at(poly: np.ndarray, shared: np.ndarray) -> np.ndarray:
    """The web ends at the junction: keep the part below the shared node."""
    return np.vstack([poly[0], shared])


def _draw_panel(
    axis,
    skin: np.ndarray,
    web: np.ndarray,
    *,
    title: str,
    hits: list,
    shared: np.ndarray | None = None,
    verdict: str,
    verdict_colour: str,
) -> None:
    axis.plot(skin[:, 0], skin[:, 1], color=NAVY, lw=3.2, solid_capstyle="round", zorder=2)
    axis.plot(web[:, 0], web[:, 1], color=TEAL, lw=3.2, solid_capstyle="round", zorder=2)
    axis.scatter(skin[:, 0], skin[:, 1], s=34, c=NAVY, zorder=3, edgecolors="white", linewidths=0.8)
    axis.scatter(web[:, 0], web[:, 1], s=34, c=TEAL, zorder=3, edgecolors="white", linewidths=0.8)
    if shared is not None:
        axis.scatter([shared[0]], [shared[1]], s=320, facecolors="none", edgecolors=TEAL,
                     linewidths=2.4, zorder=4)
    for hit in hits:
        point = np.asarray(hit[4], dtype=float)
        axis.scatter([point[0]], [point[1]], s=320, facecolors="none", edgecolors=RED,
                     linewidths=2.4, zorder=4)
    axis.text(0.5, 1.02, title, transform=axis.transAxes, ha="center", va="bottom",
              color=INK, fontsize=14, fontweight="bold")
    axis.text(0.5, -0.04, verdict, transform=axis.transAxes, ha="center", va="top",
              color=verdict_colour, fontsize=13, fontweight="bold")
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-1.15, 1.15)
    axis.set_ylim(-1.05, 1.0)
    strip_axes_chrome(axis)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="bug_log_sg009")
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    skin, web, shared = _crossing_polylines()
    hits_before = intersect_open_polylines_2d(skin, web)
    skin_split = _insert_shared_node(skin, shared)
    web_trimmed = _trim_at(web, shared)
    hits_after = intersect_open_polylines_2d(skin_split, web_trimmed)

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), gridspec_kw={"wspace": 0.16})
    _draw_panel(
        axes[0], skin, web, title="1 · Raw: web crosses skin", hits=hits_before,
        verdict=f"detector: {len(hits_before)} crossing", verdict_colour=NAVY,
    )
    _draw_panel(
        axes[1], skin_split, web_trimmed, title="2 · Fixed: one shared node", hits=[],
        shared=shared, verdict="junction now exact", verdict_colour=TEAL,
    )
    _draw_panel(
        axes[2], skin_split, web_trimmed, title="3 · Re-check", hits=hits_after,
        shared=shared,
        verdict=f"detector: {len(hits_after)} — touching ≠ crossing",
        verdict_colour=RED if len(hits_after) == 0 else NAVY,
    )
    for left, right in ((axes[0], axes[1]), (axes[1], axes[2])):
        fig.add_artist(
            FancyArrowPatch(
                (left.get_position().x1 - 0.005, 0.5), (right.get_position().x0 + 0.005, 0.5),
                transform=fig.transFigure, arrowstyle="-|>", mutation_scale=22, color=MUTED, lw=1.6,
            )
        )

    png_path = out_dir / "bug_log_sg009_schematic.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    payload = {
        "hits_before": len(hits_before),
        "hits_after_shared_node": len(hits_after),
        "oracle": "intersect_web_axis_open_polyline_2d (post-resample)",
        "detector": "intersect_open_polylines_2d (strict interior; tools only)",
        "invariant": (
            "A component must emit its mid-surface at its true extent, "
            "so junctions coincide by construction."
        ),
    }
    json_path = out_dir / "bug_log_sg009.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
