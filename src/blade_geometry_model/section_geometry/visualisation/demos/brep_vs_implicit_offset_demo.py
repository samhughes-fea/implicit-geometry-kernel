"""
visualisation/demos/brep_vs_implicit_offset_demo.py
===================================================
Why the representation was changed: what an inward offset costs each way.

Wall thickness is the design variable the sizing loop actually turns, and
applying it means offsetting the mould line inward. This demo runs that one
operation two ways on the same aerofoil:

* **Curve offset** (the boundary-representation move) - displace every vertex
  along its mitred normal. Cheap, exact where the boundary is locally convex,
  and *self-intersecting* wherever the offset distance exceeds the local radius
  of curvature. On an aerofoil that happens first at the trailing-edge wedge,
  where the two surfaces run toward each other: past a few millimetres the two
  offset surfaces cross and the curve folds into a swallowtail loop. The
  enclosed area computed from that curve is then wrong - the shoelace integral
  counts the folded lobe with the opposite sign.

* **Level-set offset** (what the package does) - evaluate ``phi = -t``. The
  boundary that comes back is the set of points at distance ``t``, which cannot
  self-intersect: where the wedge closes, the region simply stops existing.
  There is no fold to detect and nothing to trim.

**On fairness.** A production B-rep kernel does not ship the naive curve: it
detects these loops and trims them, and for an aerofoil at modest offsets it
will get the right answer. The claim here is narrower and is the one the talk
makes - the *repair step exists at all*, it is the part that gets hard as
sections get thinner and topologies change, and the implicit path never runs
it. The mitre is capped at the usual limit of 4, exactly as a real
implementation would, so the failure shown is folding and not a spike artifact.

Writes under ``visualisation/output/brep_vs_implicit_offset/``:

* ``brep_vs_implicit_offset.png`` - two-panel comparison: curve offset (left)
  vs level-set offset (right). Sized for full-width slide display (~2.5:1).

Run from the repository root::

    python blade_structure_utilities/geometry_utilities/section_geometry/visualisation/demos/brep_vs_implicit_offset_demo.py
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

_erp = Path(__file__).resolve().parent / "_ensure_repo_path.py"
_spec = importlib.util.spec_from_file_location("_ensure_repo_path", _erp)
_ensure_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_ensure_mod)
_ensure_mod.ensure_repo_path()

from blade_structure_utilities.geometry_utilities.section_geometry import (
    AirfoilSDF,
    SDFGrid,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.figure_io import (
    add_output_and_presentation_args,
    apply_presentation_style,
    save_figure,
)
from blade_structure_utilities.geometry_utilities.section_geometry.visualisation.palette import (
    INK,
    TEAL,
)

#: The standard mitre limit used by stroke/offset implementations. Without it a
#: near-cusp vertex throws its offset point to infinity, which would make this
#: comparison a strawman about spikes rather than about folding.
_MITRE_LIMIT = 4.0


def _signed_area(polygon: np.ndarray) -> float:
    """Shoelace area. Negative for a clockwise loop, and *cancelling* over a fold."""
    x = polygon[:, 0]
    y = polygon[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _inward_mitre_offset(polygon: np.ndarray, distance: float) -> np.ndarray:
    """Displace every vertex along its mitred inward normal by ``distance``.

    This is the naive curve offset: a local operation with no global awareness,
    which is precisely why it folds where the offset outruns the curvature.
    """
    points = np.asarray(polygon, dtype=float)
    # Orient counter-clockwise so "left of travel" is reliably inward.
    if _signed_area(points) < 0.0:
        points = points[::-1]

    following = np.roll(points, -1, axis=0)
    preceding = np.roll(points, 1, axis=0)

    def _left_normals(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        edges = ends - starts
        lengths = np.linalg.norm(edges, axis=1, keepdims=True)
        lengths[lengths == 0.0] = 1.0
        edges = edges / lengths
        return np.column_stack((-edges[:, 1], edges[:, 0]))

    incoming_normal = _left_normals(preceding, points)
    outgoing_normal = _left_normals(points, following)

    # Mitre vector: bisector scaled so the offset edges still meet at a point.
    bisector = incoming_normal + outgoing_normal
    scale = 1.0 + np.sum(incoming_normal * outgoing_normal, axis=1)
    scale[np.abs(scale) < 1e-12] = 1e-12
    mitre = bisector / scale[:, None]

    lengths = np.linalg.norm(mitre, axis=1)
    over_limit = lengths > _MITRE_LIMIT
    if np.any(over_limit):
        mitre[over_limit] *= (_MITRE_LIMIT / lengths[over_limit])[:, None]

    return points + float(distance) * mitre


def _count_self_intersections(polygon: np.ndarray) -> int:
    """Proper crossings between non-adjacent segments of a closed polyline.

    O(N^2) and deliberately simple: at a couple of hundred vertices it costs
    nothing, and a sweep line would only obscure what is being counted.
    """
    starts = np.asarray(polygon, dtype=float)
    ends = np.roll(starts, -1, axis=0)
    count = int(starts.shape[0])

    def _orientation(a, b, c):
        return np.sign(
            (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1])
            - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])
        )

    crossings = 0
    for i in range(count):
        # Skip the two neighbours (they share an endpoint) and, for segment 0,
        # the closing segment as well.
        j_lo = i + 2
        j_hi = count - 1 if i == 0 else count
        if j_lo >= j_hi:
            continue
        p1 = starts[i]
        p2 = ends[i]
        q1 = starts[j_lo:j_hi]
        q2 = ends[j_lo:j_hi]
        d1 = _orientation(q1, q2, p1)
        d2 = _orientation(q1, q2, p2)
        d3 = _orientation(p1, p2, q1)
        d4 = _orientation(p1, p2, q2)
        crossings += int(np.sum((d1 * d2 < 0) & (d3 * d4 < 0)))
    return crossings




def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--naca", default="2412")
    parser.add_argument("--chord", type=float, default=1.0)
    parser.add_argument("--n-points", type=int, default=220)
    parser.add_argument("--nx", type=int, default=1100)
    parser.add_argument("--ny", type=int, default=480)
    add_output_and_presentation_args(parser, default_stem="brep_vs_implicit_offset")
    args = parser.parse_args(argv)
    presentation = bool(args.presentation)
    apply_presentation_style(presentation)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    airfoil = AirfoilSDF.from_naca(
        str(args.naca), n_points=int(args.n_points), chord=float(args.chord)
    )
    polygon = np.asarray(airfoil.vertices, dtype=float)
    grid = SDFGrid.from_airfoil(airfoil, nx=int(args.nx), ny=int(args.ny))
    phi = np.asarray(airfoil(grid.X, grid.Y), dtype=float)

    # Four gauges spanning "fine" to "wedge closed at the TE".
    drawn = [0.004, 0.010, 0.020, 0.030]

    # Two-panel figure sized for a full-width slide panel (~2.5:1).
    figure, (axis_curve, axis_level) = plt.subplots(1, 2, figsize=(14.0, 5.5))
    shades = plt.cm.viridis(np.linspace(0.15, 0.85, len(drawn)))

    # --- left: curve offset, folding ------------------------------------------
    axis_curve.plot(
        polygon[:, 0], polygon[:, 1], color=INK, linewidth=1.8, zorder=5, label="outer skin"
    )
    for colour, distance in zip(shades, drawn):
        offset_polygon = _inward_mitre_offset(polygon, distance)
        closed = np.vstack([offset_polygon, offset_polygon[:1]])
        crossings = _count_self_intersections(offset_polygon)
        axis_curve.plot(
            closed[:, 0],
            closed[:, 1],
            color=colour,
            linewidth=1.4,
            label=f"t = {distance * 1e3:.0f} mm  ({crossings} crossings)",
        )
    axis_curve.set_title("Curve offset: displace each vertex along its normal", fontsize=12)
    axis_curve.legend(fontsize=9.0, loc="upper center", bbox_to_anchor=(0.5, -0.02),
                      ncol=3, frameon=False)

    # --- right: level-set offset, clean ---------------------------------------
    axis_level.plot(polygon[:, 0], polygon[:, 1], color=INK, linewidth=1.8, zorder=5)
    for colour, distance in zip(shades, drawn):
        # linestyles must be forced: matplotlib dashes negative contour levels by
        # default, which would make the two panels look different for a reason
        # that has nothing to do with the comparison.
        axis_level.contour(
            grid.X, grid.Y, phi,
            levels=[-float(distance)],
            colors=[colour],
            linewidths=1.4,
            linestyles="solid",
        )
    axis_level.set_title(r"Level-set offset: evaluate $\varphi = -t$", fontsize=12)
    axis_level.text(
        0.98, 0.08,
        "0 crossings, every gauge — by construction",
        transform=axis_level.transAxes,
        ha="right", fontsize=10.5, style="italic", color=TEAL,
    )

    vertical_pad = 0.035 * float(args.chord)
    for axis in (axis_curve, axis_level):
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(polygon[:, 0].min() - 0.02, polygon[:, 0].max() + 0.02)
        axis.set_ylim(
            float(polygon[:, 1].min()) - vertical_pad,
            float(polygon[:, 1].max()) + vertical_pad,
        )
        axis.set_xticks([])
        axis.set_yticks([])

    figure.suptitle(
        "One inward offset, two representations — the fold is the work the implicit path never does",
        fontsize=13.5,
    )
    figure.tight_layout()
    out_png = out_dir / "brep_vs_implicit_offset.png"
    save_figure(figure, out_png, presentation=presentation)
    plt.close(figure)
    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
