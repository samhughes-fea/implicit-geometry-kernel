"""
visualisation/demos/feasibility_cascade_demo.py
===============================================
What the package does when a layout cannot fit: drop → step down → raise.

Sweeps shrinking chord against fixed root-anchored web positions and plots the
resolve ladder for a multi-web catalogue key (default ``5D-F``):

* requested key stays fixed
* surviving web count steps down
* each dropped web is marked with its ``reason``
* the final chord where resolve raises ``InfeasibleSectionGeometry``

Uses ``SectionGeometryRunner(...).resolve().get_layout_resolution()`` — the
cheap path documented in ``docs/guides/feasibility_and_drops.md``.

Writes under ``visualisation/output/feasibility_cascade/``:

* ``feasibility_cascade.png``
* ``feasibility_cascade.json``

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/feasibility_cascade_demo.py --presentation
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

from blade_structure_utilities.geometry_utilities.section_contract import (
    InfeasibleSectionGeometry,
)
from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SectionGeometryRunner,
    resolve_system_type,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    INK,
    MUTED,
    NAVY,
    RED,
    SURFACE,
    TEAL,
    strip_axes_chrome,
)


def _root_web_positions_m(layout_key: str, root_chord_m: float) -> tuple[float, ...]:
    layout = resolve_system_type(layout_key)
    return tuple(float(frac) * float(root_chord_m) for frac in layout.web_chord_fracs)


def _sweep_cascade(
    *,
    layout_key: str,
    root_chord_m: float,
    chord_values_m: np.ndarray,
    skin_thickness: float,
) -> dict:
    layout = resolve_system_type(layout_key)
    web_x = _root_web_positions_m(layout_key, root_chord_m)
    rows: list[dict] = []
    raise_chord: float | None = None
    raise_message: str | None = None

    for chord in chord_values_m:
        airfoil = AirfoilSDF.from_naca("2412", chord=float(chord), n_points=160)
        runner = SectionGeometryRunner(
            layout=layout,
            airfoil=airfoil,
            skin_thickness=float(skin_thickness),
            spar_cap_thickness=0.0,
            shear_web_thickness=0.004,
            twist_rad=0.0,
            web_x_positions_m=web_x,
            outputs=(),
        )
        try:
            resolution = runner.resolve().get_layout_resolution()
        except InfeasibleSectionGeometry as exc:
            raise_chord = float(chord)
            raise_message = str(exc)
            rows.append(
                {
                    "chord_m": float(chord),
                    "status": "raise",
                    "nominal_key": layout_key,
                    "effective_key": None,
                    "n_webs_surviving": 0,
                    "stepped_down": True,
                    "dropped": [],
                    "message": raise_message,
                }
            )
            break

        dropped = [
            {"x_m": float(item.x_m), "reason": str(item.reason)}
            for item in (resolution.dropped or ())
        ]
        rows.append(
            {
                "chord_m": float(chord),
                "status": "ok",
                "nominal_key": str(resolution.nominal_layout.key),
                "effective_key": str(resolution.layout.key),
                "n_webs_surviving": int(len(resolution.web_x_positions_m)),
                "stepped_down": bool(resolution.stepped_down),
                "dropped": dropped,
                "message": None,
            }
        )

    return {
        "layout_key": layout_key,
        "root_chord_m": float(root_chord_m),
        "web_x_positions_m": list(web_x),
        "skin_thickness_m": float(skin_thickness),
        "raise_chord_m": raise_chord,
        "raise_message": raise_message,
        "rows": rows,
    }


# Human wording for the resolver's drop reasons. An unknown reason is shown
# verbatim rather than guessed at.
_DROP_REASON_TEXT = {
    "no_inner_skin_intercept": "web lands in the trailing-edge wedge",
}


def _plateaus(rows: list[dict]) -> list[dict]:
    """First row of every surviving-web level, plus the raising row."""
    plateaus: list[dict] = []
    previous_count: int | None = None
    for row in rows:
        count = int(row["n_webs_surviving"])
        if row["status"] == "raise" or count != previous_count:
            plateaus.append(row)
        previous_count = count
    return plateaus


def _vertical_extent(vertices: np.ndarray, x: float) -> tuple[float, float] | None:
    """Lowest and highest crossing of the closed outline with the line ``X = x``."""
    start = vertices
    end = np.roll(vertices, -1, axis=0)
    spans = (start[:, 0] - x) * (end[:, 0] - x) <= 0.0
    spans &= start[:, 0] != end[:, 0]
    if not np.any(spans):
        return None
    fraction = (x - start[spans, 0]) / (end[spans, 0] - start[spans, 0])
    crossings = start[spans, 1] + fraction * (end[spans, 1] - start[spans, 1])
    return float(crossings.min()), float(crossings.max())


def _draw_station_inset(axis: plt.Axes, row: dict, payload: dict, *, height_limit: float) -> None:
    """The station itself: outline, surviving webs solid, dropped webs ghosted in red."""
    chord = float(row["chord_m"])
    vertices = np.asarray(AirfoilSDF.from_naca("2412", chord=chord, n_points=160).vertices, dtype=float)
    axis.fill(vertices[:, 0], vertices[:, 1], color=SURFACE, zorder=1)
    axis.plot(
        np.append(vertices[:, 0], vertices[0, 0]),
        np.append(vertices[:, 1], vertices[0, 1]),
        color=NAVY, lw=1.4, zorder=2,
    )

    dropped_x = {round(float(item["x_m"]), 9) for item in row["dropped"]}
    for web_x in payload["web_x_positions_m"]:
        extent = _vertical_extent(vertices, float(web_x))
        is_dropped = row["status"] == "raise" or round(float(web_x), 9) in dropped_x
        if is_dropped:
            if extent is None:
                axis.scatter([web_x], [0.0], marker="x", s=36, color=RED, lw=1.8, zorder=3)
            else:
                axis.plot([web_x, web_x], extent, color=RED, lw=1.4, ls=(0, (2, 2)), zorder=3)
        elif extent is not None:
            axis.plot([web_x, web_x], extent, color=TEAL, lw=2.6, zorder=3)

    axis.set_xlim(-0.02, float(payload["root_chord_m"]) + 0.02)
    axis.set_ylim(-height_limit, height_limit)
    axis.set_aspect("equal")
    axis.invert_xaxis()
    strip_axes_chrome(axis)


def _plot_cascade(payload: dict, *, presentation: bool) -> plt.Figure:
    rows = payload["rows"]
    chords = np.array([row["chord_m"] for row in rows], dtype=float)
    n_webs = np.array(
        [row["n_webs_surviving"] if row["status"] == "ok" else 0 for row in rows],
        dtype=float,
    )
    plateaus = _plateaus(rows)

    fig = plt.figure(figsize=(14.0, 6.0))
    grid = fig.add_gridspec(2, len(plateaus), height_ratios=(0.75, 2.6), hspace=0.02, wspace=0.06)
    step_axis = fig.add_subplot(grid[1, :])

    # Top strip: every station is drawn at the root's scale, so the chord visibly
    # shrinks under webs that stay where the root put them.
    height_limit = 0.09 * float(payload["root_chord_m"])
    for index, row in enumerate(plateaus):
        inset = fig.add_subplot(grid[0, index])
        _draw_station_inset(inset, row, payload, height_limit=height_limit)
        label = row["effective_key"] if row["status"] == "ok" else "raise"
        inset.set_title(
            f"{label}\nchord {float(row['chord_m']):.2f} m",
            fontsize=12, color=RED if row["status"] == "raise" else INK, pad=2,
        )

    step_axis.step(chords, n_webs, where="post", color=NAVY, lw=2.6, label="surviving webs")
    step_axis.plot(chords[:-1] if rows[-1]["status"] == "raise" else chords,
                   n_webs[:-1] if rows[-1]["status"] == "raise" else n_webs,
                   "o", color=NAVY, ms=4)

    previous_key: str | None = None
    drop_reasons: set[str] = set()
    for row in plateaus:
        if row["status"] != "ok":
            continue
        key = str(row["effective_key"])
        drop_reasons.update(drop["reason"] for drop in row["dropped"])
        if previous_key is not None:
            # Key grammar is <webs><family>-<axis>: a changed family letter is a
            # hand-over (a box cap needs two webs to close), not just a drop.
            handover = f" ({previous_key[1]}→{key[1]})" if key[1] != previous_key[1] else ""
            count = int(row["n_webs_surviving"])
            plateau_chords = chords[(n_webs == count) & np.array([r["status"] == "ok" for r in rows])]
            step_axis.annotate(
                f"−1 web · {key}{handover}",
                xy=(float(plateau_chords.mean()), count),
                xytext=(0, 10), textcoords="offset points", ha="center",
                fontsize=12, color=TEAL, fontweight="bold",
            )
        previous_key = key

    raise_chord = payload.get("raise_chord_m")
    if raise_chord is not None:
        step_axis.axvline(raise_chord, color=RED, ls="--", lw=1.6)
        step_axis.scatter([raise_chord], [0.0], s=90, c=RED, zorder=4)
        step_axis.annotate(
            "InfeasibleSectionGeometry\nnot one web can be anchored:\nraise, never return a guess",
            xy=(raise_chord, 0.0), xytext=(-14, 4), textcoords="offset points",
            ha="right", color=RED, fontsize=11.5, fontweight="bold",
        )

    step_axis.set_xlabel("local chord [m]  →  outboard", color=INK)
    step_axis.set_ylabel("surviving webs", color=INK)
    step_axis.set_yticks(range(0, int(max(n_webs.max(), 1.0)) + 1))
    step_axis.set_ylim(-0.5, max(n_webs.max(), 1.0) + 0.9)
    step_axis.grid(True, axis="y", alpha=0.5)
    step_axis.invert_xaxis()

    legend_handles = [
        plt.Line2D([], [], color=TEAL, lw=2.6, label="web anchored"),
        plt.Line2D(
            [], [], color=RED, lw=1.4, ls=(0, (2, 2)),
            label="web dropped, never slid: "
            + "; ".join(sorted(_DROP_REASON_TEXT.get(reason, reason) for reason in drop_reasons)),
        ),
    ]
    step_axis.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize=11)
    return fig


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="feasibility_cascade")
    parser.add_argument("--layout", default="5D-F")
    parser.add_argument("--root-chord", type=float, default=1.0)
    parser.add_argument("--skin-thickness", type=float, default=0.006)
    parser.add_argument("--chord-min", type=float, default=0.08)
    parser.add_argument("--chord-max", type=float, default=1.0)
    parser.add_argument("--n-chords", type=int, default=28)
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    chord_values = np.linspace(float(args.chord_max), float(args.chord_min), int(args.n_chords))
    payload = _sweep_cascade(
        layout_key=str(args.layout),
        root_chord_m=float(args.root_chord),
        chord_values_m=chord_values,
        skin_thickness=float(args.skin_thickness),
    )
    fig = _plot_cascade(payload, presentation=bool(args.presentation))
    png_path = out_dir / "feasibility_cascade.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    json_path = out_dir / "feasibility_cascade.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
