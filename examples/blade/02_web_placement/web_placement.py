"""Placing shear webs along a tapered blade: chord fractions against fixed distances from the LE.

``BladeGeometry`` places webs in one of two ways:

``web_positions``  chord fractions, the same at every station. The webs taper with the
                   chord, so in planform they converge towards the tip.
``web_offsets``    metres from the leading edge, the same at every station, as on a blade
                   whose webs are straight. On a tapering blade a rear web runs into the
                   thin trailing-edge region; where the skins' inner faces close up behind
                   it the web no longer fits, and the cut drops it and says so in
                   ``SectionCut.dropped_webs``.

The blade is 20 m long, from a NACA 0030 at 3 m chord and 12 degrees of twist to a NACA
0015 at 1 m chord and no twist, with a 2D-F layout (flapwise webs, box caps). The script
measures each web where it crosses the chord line, in the untwisted chord frame, and
checks it against the requested fraction or offset. It then finds the span position
where the rear web (1.5 m from the LE) is dropped, and checks it against the point
where the interpolated airfoil's thickness at the web's rear face falls to two skin
thicknesses.

    python examples/blade/02_web_placement/web_placement.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from _shared.plotting import REFERENCE, SERIES, TEXT_SECONDARY, draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402

LENGTH, ROOT_CHORD, TIP_CHORD = 20.0, 3.0, 1.0
SKIN, CAP, WEB = 0.01, 0.03, 0.015
FRACTIONS = (0.2, 0.5)
OFFSETS = (0.6, 1.5)
Z = np.linspace(0.0, LENGTH, 21)
ROOT, TIP = Airfoil.naca4("0030"), Airfoil.naca4("0015")


def stations():
    return [Station(0.0, ROOT_CHORD, ROOT, twist=12.0, pitch_axis=0.3),
            Station(LENGTH, TIP_CHORD, TIP, twist=0.0, pitch_axis=0.3)]


def blade(**webs):
    return BladeGeometry(
        stations(), "2D-F-LC-0R", SKIN, pivot_fraction=0.3,
        spar_cap_thickness=CAP, shear_web_thickness=WEB, **webs,
    )


def from_le(cut, wall):
    """Distance from the LE where the web crosses the chord line, in the untwisted chord frame."""
    t = np.radians(cut.twist)
    X, Y = wall.xy[:, 0] - cut.offset[0], wall.xy[:, 1] - cut.offset[1]
    x, y = np.cos(t) * X - np.sin(t) * Y, np.sin(t) * X + np.cos(t) * Y
    order = np.argsort(y)
    return float(np.interp(0.0, y[order], x[order])) + cut.pitch_axis * cut.chord


def chord_at(z):
    return ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * z / LENGTH


def survey(geo):
    rows = []
    for z in Z:
        cut = geo.section_at(z)
        rows.append((z, cut.chord, [from_le(cut, w) for w in cut.walls_of("web")], cut.dropped_webs))
    return rows


def first_drop(geo, lo=0.0, hi=LENGTH, tol=1e-3):
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if geo.section_at(mid).dropped_webs else (mid, hi)
    return 0.5 * (lo + hi)


def predicted_drop(x_web):
    """Where the airfoil thickness at the web's rear face falls to two skin thicknesses."""
    def gap(z):
        w = z / LENGTH
        c = chord_at(z)
        xc = (x_web + 0.5 * WEB) / c
        return ((1 - w) * ROOT.thickness(xc) + w * TIP.thickness(xc)) * c - 2 * SKIN
    lo, hi = 0.0, LENGTH
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if gap(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def check(by_fraction, by_offset):
    print(f"{'z [m]':>6} {'chord':>6} | {'fractions: x/c of each web':>27} | {'offsets: x from LE [m]':>24}  dropped")
    worst_f = worst_o = 0.0
    for (z, c, xf, _), (_, _, xo, dropped) in zip(by_fraction, by_offset):
        worst_f = max(worst_f, *(abs(x / c - f) for x, f in zip(xf, FRACTIONS)))
        worst_o = max(worst_o, *(abs(x - o) for x, o in zip(xo, OFFSETS[: len(xo)])))
        print(f"{z:6.1f} {c:6.2f} | {'  '.join(f'{x / c:7.4f}' for x in xf):>27} | "
              f"{'  '.join(f'{x:7.4f}' for x in xo):>24}  {'; '.join(dropped) or '-'}")
    print(f"\nlargest error: fractions {worst_f:.1e} of the chord, offsets {worst_o * 1e3:.2f} mm")
    assert worst_f < 1e-3 and worst_o < 1e-3


def figure(by_fraction, by_offset, z_drop, z_pred):
    fig = plt.figure(figsize=(14, 6))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.3, 1.0])
    ax = fig.add_subplot(grid[0, :])
    z = np.array([r[0] for r in by_fraction])
    c = np.array([r[1] for r in by_fraction])
    ax.plot(z, np.zeros_like(z), color=REFERENCE, lw=1.0)
    ax.plot(z, c, color=REFERENCE, lw=1.0)
    ax.text(LENGTH + 0.2, 0.0, "LE", va="center", fontsize=8)
    ax.text(LENGTH + 0.2, TIP_CHORD, "TE", va="center", fontsize=8)
    for k in range(2):
        ax.plot(z, [r[2][k] for r in by_fraction], "-", color=SERIES[0], lw=1.5,
                label="web_positions = (0.2, 0.5): taper with the chord" if k == 0 else None)
        zs = [r[0] for r in by_offset if len(r[2]) > k]
        xs = [r[2][k] for r in by_offset if len(r[2]) > k]
        ax.plot(zs, xs, "-", color=SERIES[1], lw=1.5,
                label="web_offsets = (0.6, 1.5) m: straight, the rear one dropped" if k == 0 else None)
    ax.axvline(z_drop, color=SERIES[1], ls=":", lw=1.0)
    ax.annotate(f"rear web dropped from z = {z_drop:.2f} m\n(thickness check predicts {z_pred:.2f} m)",
                (z_drop, OFFSETS[1]), xytext=(z_drop + 0.6, 2.2), fontsize=8, color=TEXT_SECONDARY,
                arrowprops=dict(arrowstyle="-", color=TEXT_SECONDARY, lw=0.6))
    ax.set_xlim(0, LENGTH + 1)
    ax.set_ylim(3.2, -0.2)
    ax.set_xlabel("span z [m]")
    ax.set_ylabel("distance from the LE [m]")
    ax.set_title("Planform: where the webs cross the chord line")
    ax.legend(loc="lower left")

    geo = blade(web_offsets=OFFSETS)
    for i, zc in enumerate((0.0, 10.0, 17.0)):
        ax = fig.add_subplot(grid[1, i])
        cut = geo.section_at(zc)
        draw_cut(ax, cut)
        n = len(cut.walls_of("web"))
        ax.set_title(f"web_offsets at z = {zc:g} m: {n} web{'s' * (n != 1)}, {len(cut.cells)} cells, "
                     f"layout {cut.layout}", fontsize=8.5)
        ax.set_xlim(-1.0, 2.2)
        ax.set_ylim(-0.6, 0.6)
    note(fig, "20 m blade, NACA 0030 at 3 m chord and 12° twist to NACA 0015 at 1 m and 0°; 2D-F, skin 10 mm, "
              "caps 30 mm, webs 15 mm. Web positions match the request to under 1 mm.")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, __file__, "web_placement")


if __name__ == "__main__":
    by_fraction = survey(blade(web_positions=FRACTIONS))
    offset_blade = blade(web_offsets=OFFSETS)
    by_offset = survey(offset_blade)
    check(by_fraction, by_offset)
    z_drop, z_pred = first_drop(offset_blade), predicted_drop(OFFSETS[1])
    print(f"rear web dropped from z = {z_drop:.3f} m; thickness at its rear face equals two skins at z = {z_pred:.3f} m")
    assert abs(z_drop - z_pred) < 0.1
    figure(by_fraction, by_offset, z_drop, z_pred)
