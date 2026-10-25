"""
visualisation/demos/brep_extend_trim_demo.py
============================================
Extend and trim: how many code paths does placing one web cost?

Companion to ``brep_vs_implicit_offset_demo.py``, which covers *offset*. Here a
single chord-normal web of fixed thickness is swept across the chord of a
NACA 2412 with an inner skin face and a fixed box spar cap, and placed two ways.

**Boundary representation.** Each web face is a line that must be *extended*
to meet the inner skin face and *trimmed* there. Which branch that takes
depends on where the web is, and every branch is code somebody maintains:

* ``extend + trim``  — each face line meets the inner face exactly twice, cleanly;
* ``near-tangent``   — an intersection exists but the line meets the face at a
  grazing angle, so the trim point is ill-conditioned (the nose);
* ``split at cap``   — a cap end falls between the two web faces, so the web
  end meets two different faces and must be split into a T;
* ``repair fold``    — the curve-offset inner face has folded through itself
  (the trailing edge, as on the offset slide), so a face line meets it more
  than twice and the loop must be repaired before the trim is defined;
* ``drop``           — the face line meets the inner face nowhere.

The inner face is the naive mitred curve offset from the offset demo, reused
so the two slides describe the same geometry.

**Implicit.** The web is one expression at every position:
``Intersect(Capsule(axis, t_web/2), Offset(airfoil, -t_skin))``. The axis is
the IML-to-IML segment, lengthened so the stadium's round ends sit outside
the cavity. The cap is a **separate** field: overlap is legal. The figure
maps the B-rep branches along the chord, then zooms into the three places
the sweep leaves the clean branch, with that same capsule ∩ cavity drawn
in teal.

Sizes are illustrative (5 mm skin, 12 mm cap over 0.25–0.55 c, 6 mm web), not a
design.

Writes under ``visualisation/output/brep_extend_trim/``:

* ``brep_extend_trim.png``
* ``brep_extend_trim.json``

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/brep_extend_trim_demo.py --presentation
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_geometry import AirfoilSDF
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
    TEAL,
    strip_axes_chrome,
)

# The naive curve offset lives with the offset demo; reuse it so both slides
# describe the same inner face.
_offset_spec = importlib.util.spec_from_file_location(
    "brep_vs_implicit_offset_demo", Path(__file__).resolve().parent / "brep_vs_implicit_offset_demo.py"
)
_offset_mod = importlib.util.module_from_spec(_offset_spec)
assert _offset_spec.loader is not None
_offset_spec.loader.exec_module(_offset_mod)
inward_mitre_offset = _offset_mod._inward_mitre_offset

_SKIN_M = 0.005
_CAP_M = 0.012
_CAP_SPAN = (0.25, 0.55)
_WEB_M = 0.006
_GRAZING_DEG = 15.0

#: Branch name → colour, in the order the legend lists them.
_BRANCHES: dict[str, str] = {
    "extend + trim": "#d9dde3",
    "near-tangent": "#f3c77a",
    "split at cap": TEAL,
    "repair fold": RED,
    "drop": NAVY,
}


def _vertical_crossings(polygon: np.ndarray, x: float) -> list[tuple[float, float]]:
    """``(y, angle_deg)`` for every crossing of the line ``X = x`` with the closed polyline.

    ``angle_deg`` is between the line and the crossed edge: 90° is square-on,
    near 0° is grazing.
    """
    starts = polygon
    ends = np.roll(polygon, -1, axis=0)
    crossings: list[tuple[float, float]] = []
    for start, end in zip(starts, ends):
        if (start[0] - x) * (end[0] - x) > 0.0 or start[0] == end[0]:
            continue
        fraction = (x - start[0]) / (end[0] - start[0])
        if not 0.0 <= fraction < 1.0:
            continue
        edge = end - start
        cosine = abs(edge[1]) / max(float(np.linalg.norm(edge)), 1e-300)
        crossings.append((float(start[1] + fraction * edge[1]), float(np.degrees(np.arccos(min(cosine, 1.0))))))
    return crossings


def _brep_branch(inner_face: np.ndarray, x_web: float) -> str:
    """Which extend/trim branch a B-rep web placement at ``x_web`` has to take."""
    faces = (x_web - 0.5 * _WEB_M, x_web + 0.5 * _WEB_M)
    counts = []
    angles = []
    for x_face in faces:
        crossings = _vertical_crossings(inner_face, x_face)
        counts.append(len(crossings))
        angles.extend(angle for _, angle in crossings)
    if min(counts) == 0:
        return "drop"
    if max(counts) > 2:
        return "repair fold"
    if any(faces[0] < cap_end < faces[1] for cap_end in _CAP_SPAN):
        return "split at cap"
    if min(angles) < _GRAZING_DEG:
        return "near-tangent"
    return "extend + trim"


def _sdf_capsule(x, y, ax: float, ay: float, bx: float, by: float, r: float) -> np.ndarray:
    """Metric SDF of a stadium: distance to the segment minus radius."""
    abx = bx - ax
    aby = by - ay
    apx = x - ax
    apy = y - ay
    ab2 = abx * abx + aby * aby + 1e-30
    t = np.clip((apx * abx + apy * aby) / ab2, 0.0, 1.0)
    return np.hypot(apx - t * abx, apy - t * aby) - r


def _iml_anchors(airfoil, x_web: float) -> tuple[float, float] | None:
    """Upper and lower inner-mould-line intercepts at ``x_web``, or None if empty."""
    y_vals = np.linspace(-0.20, 0.20, 801)
    phi = np.asarray(airfoil(np.full_like(y_vals, x_web), y_vals), dtype=float) + _SKIN_M
    hits = np.where(np.diff(np.sign(phi)))[0]
    if len(hits) < 2:
        return None
    ys: list[float] = []
    for i in hits:
        y0, y1 = float(y_vals[i]), float(y_vals[i + 1])
        p0, p1 = float(phi[i]), float(phi[i + 1])
        ys.append(y0 if p1 == p0 else y0 - p0 * (y1 - y0) / (p1 - p0))
    ys.sort()
    return ys[-1], ys[0]


def _implicit_web_field(airfoil, x_grid: np.ndarray, y_grid: np.ndarray, x_web: float) -> np.ndarray:
    """Capsule ∩ cavity. Cavity is Offset(airfoil, −t_skin); the cap is not subtracted."""
    cavity = np.asarray(airfoil(x_grid, y_grid), dtype=float) + _SKIN_M
    anchors = _iml_anchors(airfoil, x_web)
    if anchors is None:
        return np.ones_like(x_grid)
    y_top, y_bot = anchors
    radius = 0.5 * _WEB_M
    extend = 2.0 * _WEB_M + _SKIN_M
    return np.maximum(
        _sdf_capsule(x_grid, y_grid, x_web, y_top + extend, x_web, y_bot - extend, radius),
        cavity,
    )


def _cap_field(phi: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    # The band from the inner skin face (inner = 0) to _CAP_M inside it
    # (inner = -_CAP_M), limited to the cap span.
    inner = phi + _SKIN_M
    cap_slab = np.maximum(_CAP_SPAN[0] - x_grid, x_grid - _CAP_SPAN[1])
    return np.maximum(np.maximum(inner, -inner - _CAP_M), cap_slab)


#: Three stations, one expression. Captions lead with the kernel; B-rep is the aside.
_ZOOMS = (
    ("near-tangent", 0.0080, (-0.004, 0.040), (-0.022, 0.030),
     "Nose: same capsule intersect interior.\nB-rep: grazing trim, ill-conditioned"),
    ("split at cap", 0.2500, (0.222, 0.278), (0.040, 0.082),
     "Cap end: capsule overlaps the cap field;\nwinner names the laminate. B-rep would split a T"),
    ("repair fold", 0.9850, (0.955, 1.004), (-0.017, 0.019),
     "Trailing edge: cavity closed, implicit web empty. B-rep: folded offset, repair first"),
)

_ZOOM_TITLES = {
    "near-tangent": "Implicit — nose",
    "split at cap": "Implicit — cap end",
    "repair fold": "Implicit — trailing edge",
}


def _draw_zoom(axis, airfoil, outline, inner_face, branch: str, x_web: float, x_window, y_window,
               caption: str) -> None:
    x_axis = np.linspace(*x_window, 700)
    y_axis = np.linspace(*y_window, 500)
    x_grid, y_grid = np.meshgrid(x_axis, y_axis)
    phi = np.asarray(airfoil(x_grid, y_grid), dtype=float)

    axis.contourf(x_grid, y_grid, phi, levels=[-np.inf, 0.0], colors=["#f4f6f9"], zorder=0)
    axis.contourf(x_grid, y_grid, _cap_field(phi, x_grid), levels=[-np.inf, 0.0], colors=[MUTED],
                  alpha=0.55, zorder=1)
    axis.contourf(
        x_grid, y_grid, _implicit_web_field(airfoil, x_grid, y_grid, x_web),
        levels=[-np.inf, 0.0], colors=[TEAL], alpha=0.9, zorder=2,
    )

    closed_outline = np.vstack([outline, outline[:1]])
    closed_inner = np.vstack([inner_face, inner_face[:1]])
    axis.plot(closed_outline[:, 0], closed_outline[:, 1], color=INK, lw=2.0, zorder=3)
    axis.plot(closed_inner[:, 0], closed_inner[:, 1], color=NAVY, lw=1.4, ls="--", zorder=3)

    for x_face in (x_web - 0.5 * _WEB_M, x_web + 0.5 * _WEB_M):
        axis.plot([x_face, x_face], list(y_window), color=RED, lw=1.1, ls=":", zorder=4, alpha=0.85)
        crossings = _vertical_crossings(inner_face, x_face)
        if crossings:
            axis.scatter([x_face] * len(crossings), [y for y, _ in crossings], s=55,
                         facecolors="none", edgecolors=RED, linewidths=1.6, zorder=5)

    anchors = _iml_anchors(airfoil, x_web)
    if anchors is not None:
        y_top, y_bot = anchors
        extend = 2.0 * _WEB_M + _SKIN_M
        axis.plot(
            [x_web, x_web], [y_bot - extend, y_top + extend],
            color="#0d3b4c", lw=2.2, zorder=6, solid_capstyle="round",
        )

    axis.set_xlim(*x_window)
    axis.set_ylim(*y_window)
    axis.set_aspect("equal", adjustable="box")
    strip_axes_chrome(axis, keep_frame=True)
    for spine in axis.spines.values():
        spine.set_color(TEAL)
        spine.set_linewidth(2.4)
    axis.set_title(_ZOOM_TITLES[branch], loc="left", fontsize=13, color=TEAL, fontweight="bold")
    axis.text(0.0, -0.04, caption, transform=axis.transAxes, ha="left", va="top", fontsize=10.5, color=INK)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_output_and_presentation_args(parser, default_stem="brep_extend_trim")
    parser.add_argument("--positions", type=int, default=400)
    args = parser.parse_args(argv)
    apply_presentation_style(bool(args.presentation))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    airfoil = AirfoilSDF.from_naca("2412", chord=1.0, n_points=220)
    outline = np.asarray(airfoil.vertices, dtype=float)
    inner_face = inward_mitre_offset(outline, _SKIN_M)

    positions = np.linspace(0.003, 0.997, int(args.positions))
    branches = [_brep_branch(inner_face, float(x)) for x in positions]
    for branch, x_web, *_ in _ZOOMS:
        found = _brep_branch(inner_face, x_web)
        if found != branch:
            raise RuntimeError(f"zoom at x = {x_web} was meant to show {branch!r} but classifies as {found!r}")

    runs: list[tuple[float, float, str]] = []
    run_start = 0
    half_step = 0.5 * (positions[1] - positions[0])
    for index in range(1, len(positions) + 1):
        if index == len(positions) or branches[index] != branches[run_start]:
            runs.append((positions[run_start] - half_step, positions[index - 1] + half_step, branches[run_start]))
            run_start = index
    used = [name for name in _BRANCHES if name in branches]
    n_switches = len(runs) - 1

    fig = plt.figure(figsize=(14.0, 6.2))
    outer = fig.add_gridspec(2, 1, height_ratios=(0.75, 1.25), hspace=0.3)
    ax_section = fig.add_subplot(outer[0])
    zoom_grid = outer[1].subgridspec(1, 3, wspace=0.12)

    # Top: where along the chord each branch is taken. Bands are widened to a
    # visible minimum so a 3 mm band is still findable; the zooms show them true.
    minimum_band = 0.006
    for start, end, name in runs:
        if name == "extend + trim":
            continue
        centre = 0.5 * (start + end)
        half = max(0.5 * (end - start), 0.5 * minimum_band)
        ax_section.axvspan(centre - half, centre + half, color=_BRANCHES[name], alpha=0.85, lw=0, zorder=0)
    closed_outline = np.vstack([outline, outline[:1]])
    closed_inner = np.vstack([inner_face, inner_face[:1]])
    ax_section.fill(closed_outline[:, 0], closed_outline[:, 1], color="#f4f6f9", zorder=1)
    x_axis = np.linspace(-0.01, 1.01, 1400)
    y_axis = np.linspace(-0.08, 0.10, 300)
    x_grid, y_grid = np.meshgrid(x_axis, y_axis)
    phi = np.asarray(airfoil(x_grid, y_grid), dtype=float)
    ax_section.contourf(x_grid, y_grid, _cap_field(phi, x_grid), levels=[-np.inf, 0.0], colors=[MUTED],
                        alpha=0.55, zorder=2)
    ax_section.plot(closed_outline[:, 0], closed_outline[:, 1], color=INK, lw=1.8, zorder=3)
    ax_section.plot(closed_inner[:, 0], closed_inner[:, 1], color=NAVY, lw=1.0, ls="--", zorder=3)
    ax_section.set_xlim(-0.02, 1.02)
    ax_section.set_aspect("equal", adjustable="datalim")
    strip_axes_chrome(ax_section)
    ax_section.set_title(
        f"One expression along the chord: Capsule intersect Offset(airfoil, -t_skin). "
        f"Colour bands: where B-rep would change branch ({n_switches} switches)",
        loc="left", fontsize=12.5, color=INK,
    )
    ax_section.legend(
        handles=[Patch(facecolor=_BRANCHES[name], label=f"B-rep: {name}") for name in used],
        loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=len(used), frameon=False, fontsize=10.5,
    )

    for column, (branch, x_web, x_window, y_window, caption) in enumerate(_ZOOMS):
        _draw_zoom(fig.add_subplot(zoom_grid[0, column]), airfoil, outline, inner_face, branch, x_web,
                   x_window, y_window, caption)

    fig.text(
        0.5, 0.0,
        "Teal fill + dark axis: the kernel web, same Capsule intersect interior in every zoom. "
        "Red dotted lines and rings: B-rep face-meets-offset — not the recipe.",
        ha="center", va="top", fontsize=11.5, color=TEAL, fontweight="bold",
    )

    png_path = out_dir / "brep_extend_trim.png"
    save_figure(fig, png_path, presentation=bool(args.presentation))
    plt.close(fig)
    print(f"Wrote {png_path}")

    payload = {
        "airfoil": "NACA 2412, chord 1 m",
        "sizes_m": {"skin": _SKIN_M, "cap": _CAP_M, "cap_span_chord": list(_CAP_SPAN), "web": _WEB_M},
        "grazing_threshold_deg": _GRAZING_DEG,
        "positions": int(args.positions),
        "brep_branches_used": used,
        "brep_branch_switches": n_switches,
        "runs": [{"from": float(a), "to": float(b), "branch": name} for a, b, name in runs],
        "implicit_web": "Intersect(Capsule, Offset(airfoil, -t_skin))",
        "implicit_expressions": 1,
    }
    json_path = out_dir / "brep_extend_trim.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
