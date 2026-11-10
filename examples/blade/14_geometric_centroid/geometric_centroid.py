"""Geometric centroid of the walls: boundary integral, thin-wall strip, area integration.

The kernel has no materials or loads. The only centroid it can form is geometric.
This example compares the three realisations on a made-up two-web loft:

1. Boundary integral — Green / shoelace on constructive wall-solid polygons.
2. Thin-wall strip — midline length times thickness (t^3 terms dropped).
3. Area integration — adaptive quadtree quadrature of wall SDF solids.

    python examples/blade/14_geometric_centroid/geometric_centroid.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.blade_viz import _map_host_xy  # noqa: E402
from _shared.plotting import NAVY, RED, SERIES, TEAL, draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.contract import InfeasibleSectionGeometry  # noqa: E402
from implicit_geometry_kernel.section.boundary.cores import emit_constructive_polygons  # noqa: E402
from implicit_geometry_kernel.section.boundary.labels import is_wall_solid_label  # noqa: E402
from implicit_geometry_kernel.section.properties.backend import quadtree_spec_for_airfoil  # noqa: E402
from implicit_geometry_kernel.section.properties.compute import compute_section_properties  # noqa: E402
from implicit_geometry_kernel.section.properties.polygons import (  # noqa: E402
    section_properties_from_polygons,
)

N_SPAN = 13
STATION_Z = 5.0
COLOUR_SHOE = SERIES[0]
COLOUR_STRIP = SERIES[1]
COLOUR_QT = SERIES[2]


def _loft() -> BladeGeometry:
    foil = Airfoil.naca4("0018")
    stations = [
        Station(0.0, 2.0, foil, twist=8.0, pitch_axis=0.3, offset=(0.0, 0.0)),
        Station(10.0, 1.0, foil, twist=0.0, pitch_axis=0.3, offset=(0.0, 0.0)),
    ]
    return BladeGeometry(
        stations, layout="2D-F-LC-0R", skin_thickness=0.004,
        spar_cap_thickness=0.012, shear_web_thickness=0.006, pivot_fraction=0.3,
    )


def _wall_props(props: dict) -> dict:
    return {k: v for k, v in props.items() if is_wall_solid_label(k)}


def _aggregate(props: dict) -> tuple[float, float, float] | None:
    walls = _wall_props(props)
    area = sx = sy = 0.0
    for p in walls.values():
        a = float(p["area"])
        cx, cy = float(p["cx"]), float(p["cy"])
        if a <= 1e-18 or not np.isfinite(cx) or not np.isfinite(cy):
            continue
        area += a
        sx += a * cx
        sy += a * cy
    if area <= 1e-18:
        return None
    return area, sx / area, sy / area


def _to_blade(cut, cx: float, cy: float) -> tuple[float, float]:
    xy = _map_host_xy(cut, np.array([[cx, cy]], float))
    return float(xy[0, 0]), float(xy[0, 1])


def _shoelace(cut) -> tuple[float, float, float] | None:
    host = getattr(cut, "host", None)
    if host is None:
        return None
    polys = emit_constructive_polygons(host)
    tot = _aggregate(section_properties_from_polygons(polys))
    if tot is None:
        return None
    area, cx, cy = tot
    k = float(getattr(cut, "scale", 1.0) or 1.0)
    bx, by = _to_blade(cut, cx, cy)
    return area / (k * k), bx, by


def _strip(cut) -> tuple[float, float, float] | None:
    area = sx = sy = 0.0
    for wall in cut.walls:
        p, t = wall.xy, float(wall.thickness)
        if len(p) < 2 or t <= 0.0:
            continue
        a, b = p[:-1], p[1:]
        length = np.linalg.norm(b - a, axis=1)
        area += t * float(length.sum())
        sx += t * float(np.sum(length * 0.5 * (a[:, 0] + b[:, 0])))
        sy += t * float(np.sum(length * 0.5 * (a[:, 1] + b[:, 1])))
    if area <= 1e-18:
        return None
    return area, sx / area, sy / area


def _quadtree(cut) -> tuple[float, float, float] | None:
    host = getattr(cut, "host", None)
    if host is None:
        return None
    af = getattr(host, "airfoil", None) or getattr(host, "_af", None)
    if af is None:
        return None
    k = float(getattr(cut, "scale", 1.0) or 1.0)
    chord_h = float(cut.chord) * k
    spec = quadtree_spec_for_airfoil(af, chord=chord_h, max_depth=7, thin_band=True)
    props, _meta = compute_section_properties(
        host, backend="quadtree", quadtree_spec=spec, chord=chord_h,
        twist_rad=0.0, airfoil_for_spec=af,
    )
    tot = _aggregate(props)
    if tot is None:
        return None
    area, cx, cy = tot
    bx, by = _to_blade(cut, cx, cy)
    return area / (k * k), bx, by


def sample(geo: BladeGeometry, n: int = N_SPAN) -> dict:
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = np.linspace(z0, z1, n)
    le, te, pitch = [], [], []
    shoe, strip, qt = [], [], []
    for zi in zs:
        ring = geo.oml_polyline(float(zi))
        le.append(ring[np.argmin(ring[:, 0])])
        te.append(ring[np.argmax(ring[:, 0])])
        _, _, _, offset, _ = geo.interpolate(float(zi))
        pitch.append((offset[0], float(zi)))
        try:
            cut = geo.section_at(float(zi))
        except InfeasibleSectionGeometry:
            continue
        s = _shoelace(cut)
        if s is not None:
            shoe.append((s[1], float(zi), s[0]))
        m = _strip(cut)
        if m is not None:
            strip.append((m[1], float(zi), m[0]))
        q = _quadtree(cut)
        if q is not None:
            qt.append((q[1], float(zi), q[0]))
    return {
        "le": np.asarray(le),
        "te": np.asarray(te),
        "pitch": np.asarray(pitch),
        "shoelace": np.asarray(shoe) if shoe else np.zeros((0, 3)),
        "strip": np.asarray(strip) if strip else np.zeros((0, 3)),
        "quadtree": np.asarray(qt) if qt else np.zeros((0, 3)),
    }


def report(geo: BladeGeometry, z: float) -> None:
    cut = geo.section_at(z)
    s, m, q = _shoelace(cut), _strip(cut), _quadtree(cut)
    print(f"z = {z:.1f} m  chord = {cut.chord:.3f} m  layout = {cut.layout}")
    print(f"{'route':<28} {'area m^2':>12} {'x_B m':>12} {'y_B m':>12} {'|dx| vs shoe':>14}")
    rows = (
        ("boundary integral", s),
        ("thin-wall strip", m),
        ("area integration", q),
    )
    xref = s[1] if s is not None else float("nan")
    yref = s[2] if s is not None else float("nan")
    for name, row in rows:
        if row is None:
            print(f"{name:<28} {'—':>12}")
            continue
        a, cx, cy = row
        dx = float(np.hypot(cx - xref, cy - yref))
        print(f"  {name:<26} {a:12.6e} {cx:12.5f} {cy:12.5f} {dx:14.3e}")


def figure_planform(data: dict) -> None:
    le, te = data["le"], data["te"]
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    ax.plot(le[:, 2], le[:, 0], color=NAVY, lw=1.2, label="LE")
    ax.plot(te[:, 2], te[:, 0], color=TEAL, lw=1.2, label="TE")
    ax.fill_between(le[:, 2], le[:, 0], te[:, 0], color=NAVY, alpha=0.08, lw=0)
    ax.plot(data["pitch"][:, 1], data["pitch"][:, 0], color=RED, lw=1.0, ls="--",
            label="pitch axis")
    if len(data["shoelace"]):
        ax.plot(data["shoelace"][:, 1], data["shoelace"][:, 0], color=COLOUR_SHOE, lw=1.1,
                label="boundary integral")
    if len(data["strip"]):
        ax.plot(data["strip"][:, 1], data["strip"][:, 0], color=COLOUR_STRIP, lw=1.1, ls="--",
                label="thin-wall strip")
    if len(data["quadtree"]):
        ax.plot(data["quadtree"][:, 1], data["quadtree"][:, 0], color=COLOUR_QT, lw=1.1, ls=":",
                label="area integration")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_title("Geometric centroid of the walls (blade-frame x–z)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=8, frameon=False)
    note(fig, "Uniform density. No materials or loads — not a shear centre or elastic axis.")
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    save(fig, __file__, "planform")


def figure_station(geo: BladeGeometry, z: float) -> None:
    cut = geo.section_at(z)
    s, m, q = _shoelace(cut), _strip(cut), _quadtree(cut)
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    draw_cut(ax, cut)
    marks = (
        ("boundary integral", s, COLOUR_SHOE, "o"),
        ("thin-wall strip", m, COLOUR_STRIP, "s"),
        ("area integration", q, COLOUR_QT, "D"),
    )
    for label, row, colour, marker in marks:
        if row is None:
            continue
        ax.plot(row[1], row[2], marker, color=colour, ms=6, label=f"{label} ({row[1]:.3f}, {row[2]:.4f}) m")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"z = {z:.1f} m — geometric centroids of the walls")
    ax.legend(loc="upper right", fontsize=7)
    note(fig, "Skins navy, caps red, webs teal. Markers are wall-solid centroids, not enclosed-bay centroids.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, __file__, "station")


def main() -> None:
    geo = _loft()
    report(geo, STATION_Z)
    data = sample(geo, n=N_SPAN)
    figure_planform(data)
    figure_station(geo, STATION_Z)


if __name__ == "__main__":
    main()
