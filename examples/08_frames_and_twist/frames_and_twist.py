"""Coordinate frames and the twist convention: points, second moments and web directions.

Three frames carry a section's geometry (``implicit_geometry_kernel.frames``):

NACA   x along the chord from the LE, y through the thickness (the airfoil's own axes).
S      the zero-twist station frame on the pitch axis, with the components swapped:
       X_S = y (flapwise), Y_S = x - p c (chordwise). ``naca_to_S``.
B      the twisted blade frame: S rotated by the twist,
       X_B = cos t X_S + sin t Y_S, Y_B = -sin t X_S + cos t Y_S. ``rotate_S_to_B``.

``BladeGeometry`` hands its cuts over in its own blade frame (x along the chord towards
the TE, y towards the suction side, origin on the pitch axis, positive twist nose up),
so most users never touch S and B. They matter when section properties or loads move
between the station and the blade.

The script checks

1. the point transforms: NACA -> S -> B -> NACA returns every point to 1e-15 m;
2. the second-moment transform: ``transform_second_moments_S_to_B`` against Ixx, Iyy and
   Ixy integrated afresh from the rotated outline, at twists from -90 to 90 degrees,
   and that Ixx + Iyy is invariant;
3. BladeGeometry's twist sign: a chord-normal (CN) web lies along (sin t, cos t), normal
   to a chord that runs along (cos t, -sin t), so the LE rises for positive twist,
   while a flapwise (F) web stays along (0, 1) whatever the twist.

    python examples/08_frames_and_twist/frames_and_twist.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402

from _shared.plotting import REFERENCE, SERIES, draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.frames import S_to_naca, naca_to_S, rotate_B_to_S, rotate_S_to_B  # noqa: E402
from implicit_geometry_kernel.frames.section_frames.tensors import transform_second_moments_S_to_B  # noqa: E402
from implicit_geometry_kernel.integration import section_properties_S  # noqa: E402

CHORD, PITCH = 2.0, 0.3
TWISTS = np.linspace(-90.0, 90.0, 13)


def outline_naca():
    a = Airfoil.naca4("2418")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 201)))
    return np.vstack([np.column_stack([x[::-1], a.surface("upper", x)[::-1]]),
                      np.column_stack([x[1:], a.surface("lower", x)[1:]])])[:-1] * CHORD


def check_points(naca):
    worst = 0.0
    for t in TWISTS:
        S = naca_to_S(naca, CHORD, PITCH)
        back = S_to_naca(rotate_B_to_S(rotate_S_to_B(S, np.radians(t)), np.radians(t)), CHORD, PITCH)
        worst = max(worst, float(np.abs(back - naca).max()))
    print(f"1. NACA -> S -> B -> NACA round trip: largest error {worst:.1e} m")
    assert worst < 1e-12
    return worst


def check_second_moments(naca):
    S = naca_to_S(naca, CHORD, PITCH)
    _, _, ixx, iyy, ixy = section_properties_S(Polygon(S))
    rows = []
    print("\n2. Second moments of the outline about the pitch axis, S -> B")
    print(f"{'twist':>6} {'Ixx_B':>10} {'integrated':>10} {'Iyy_B':>10} {'integrated':>10} "
          f"{'Ixy_B':>10} {'integrated':>10}")
    worst = 0.0
    for t in TWISTS:
        got = transform_second_moments_S_to_B(ixx, iyy, ixy, twist_deg=t)
        _, _, bxx, byy, bxy = section_properties_S(Polygon(rotate_S_to_B(S, np.radians(t))))
        ref = (bxx, byy, bxy)
        worst = max(worst, max(abs(g - r) for g, r in zip(got, ref)) / (ixx + iyy))
        rows.append((t, got, ref))
        print(f"{t:6.1f} {got[0]:10.5f} {bxx:10.5f} {got[1]:10.5f} {byy:10.5f} {got[2]:10.5f} {bxy:10.5f}")
        assert abs(got[0] + got[1] - (ixx + iyy)) < 1e-12 * (ixx + iyy)
    print(f"largest difference {worst:.1e} of Ixx + Iyy; Ixx + Iyy = {ixx + iyy:.5f} m^4 at every twist")
    assert worst < 1e-10
    return rows, worst


def web_direction(wall):
    d = wall.xy[-1] - wall.xy[0]
    d = d / np.linalg.norm(d)
    return d if d[1] > 0 else -d


def check_twist_sign():
    airfoil = Airfoil.naca4("2418")
    print("\n3. BladeGeometry twist: web directions in the blade frame")
    print(f"{'twist':>6} {'CN web':>18} {'(sin t, cos t)':>18} {'F web':>18}")
    worst = 0.0
    cuts = {}
    for t in (-20.0, -10.0, 0.0, 10.0, 20.0):
        st = [Station(0.0, CHORD, airfoil, twist=t, pitch_axis=PITCH), Station(1.0, CHORD, airfoil, twist=t, pitch_axis=PITCH)]
        cn = BladeGeometry(st, "1A-CN", 0.006, shear_web_thickness=0.01, web_positions=(0.4,)).section_at(0.5)
        f = BladeGeometry(st, "1A-F", 0.006, shear_web_thickness=0.01, web_positions=(0.4,)).section_at(0.5)
        dcn, df = web_direction(cn.walls_of("web")[0]), web_direction(f.walls_of("web")[0])
        expect = np.array([np.sin(np.radians(t)), np.cos(np.radians(t))])
        worst = max(worst, float(np.linalg.norm(dcn - expect)), float(np.linalg.norm(df - [0.0, 1.0])))
        cuts[t] = (cn, f)
        print(f"{t:6.1f} ({dcn[0]:+.5f}, {dcn[1]:.5f}) ({expect[0]:+.5f}, {expect[1]:.5f}) ({df[0]:+.5f}, {df[1]:.5f})")
    print(f"largest direction error {worst:.1e}")
    assert worst < 1e-6
    return cuts, worst


def figure(naca, rows, cuts, errs):
    fig = plt.figure(figsize=(14, 7.0))
    grid = fig.add_gridspec(2, 3)
    t = 20.0
    S = naca_to_S(naca, CHORD, PITCH)
    B = rotate_S_to_B(S, np.radians(t))
    for i, (pts, title, xl, yl) in enumerate((
            (naca, "NACA: x from the LE, y thickness", "x [m]", "y [m]"),
            (S, f"S: X_S = y, Y_S = x − {PITCH:g} c", "X_S (flapwise) [m]", "Y_S (chordwise) [m]"),
            (B, f"B: S rotated by the twist, here {t:g}°", "X_B [m]", "Y_B [m]"))):
        ax = fig.add_subplot(grid[0, i])
        ax.fill(pts[:, 0], pts[:, 1], color=SERIES[0], alpha=0.25, lw=1.0, ec=SERIES[0])
        ax.plot(pts[0, 0], pts[0, 1], "o", color=SERIES[1], ms=4, label="TE")
        le = pts[np.argmin(naca[:, 0])]
        ax.plot(le[0], le[1], "s", color=SERIES[2], ms=4, label="LE")
        ax.axhline(0.0, color="0.75", lw=0.5)
        ax.axvline(0.0, color="0.75", lw=0.5)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
    fig.axes[0].legend(loc="upper right", fontsize=7)

    ax = fig.add_subplot(grid[1, :2])
    tw = np.array([r[0] for r in rows])
    for k, name in enumerate(("Ixx_B", "Iyy_B", "Ixy_B")):
        ax.plot(tw, [r[1][k] for r in rows], "-", color=SERIES[k], label=f"{name}: transform_second_moments_S_to_B")
        ax.plot(tw, [r[2][k] for r in rows], "o", mfc="none", color=REFERENCE, ms=5,
                label="integrated from the rotated outline" if k == 0 else None)
    ax.set_xlabel("twist [°]")
    ax.set_ylabel("second moment about the pitch axis [m⁴]")
    ax.set_title("Second moments rotate with the outline")
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)

    ax = fig.add_subplot(grid[1, 2])
    for tt, colour in ((-20.0, SERIES[0]), (20.0, SERIES[1])):
        cn, _ = cuts[tt]
        draw_cut(ax, cn, alpha=0.25)
        w = cn.walls_of("web")[0]
        ax.plot(w.xy[:, 0], w.xy[:, 1], color=colour, lw=2.0, label=f"CN web at {tt:+g}°")
    ax.plot(0.0, 0.0, "k+")
    ax.set_title("Positive twist is nose up (1A-CN, ±20°)")
    ax.legend(fontsize=7.5, loc="lower left")
    note(fig, f"NACA 2418, 2 m chord, pitch axis at {PITCH:g} c. Round trips {errs[0]:.0e} m, second moments "
              f"{errs[1]:.0e} of Ixx + Iyy, web directions {errs[2]:.0e}.")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, __file__, "frames_and_twist")


if __name__ == "__main__":
    naca = outline_naca()
    e1 = check_points(naca)
    rows, e2 = check_second_moments(naca)
    cuts, e3 = check_twist_sign()
    figure(naca, rows, cuts, (e1, e2, e3))
