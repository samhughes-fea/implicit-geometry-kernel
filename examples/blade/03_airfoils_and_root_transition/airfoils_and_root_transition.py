"""Airfoils from codes and coordinates, and a blade root that turns a circle into an airfoil.

``Airfoil`` is a chord-normalised profile with ``surface``, ``thickness`` and
``max_thickness`` queries. It comes from a NACA 4-digit code (``Airfoil.naca4``), an
ellipse or circle (``Airfoil.ellipse``, as at a blade root) or any coordinates in Selig
order (``Airfoil(points)``), with a sharp or a blunt (flatback) trailing edge.

``BladeGeometry.interpolate`` morphs a root circle into an airfoil with a smoothstep
blend (not a raw y(x) lerp), then thickness-blends foil to foil.

The script checks

- the NACA 4-digit generator: maximum thickness t and maximum camber m of 0012, 0018,
  0030, 2412 and 4415 against the digits;
- a flatback built from coordinates (a NACA 0030 opened by 5 % of the chord at the
  TE): its TE thickness;
- the root transition, a circle at z = 0 to a NACA 0030 at 4 m and a NACA 0021 at 12 m:
  t/c is monotone in the root interval, the TE gap never exceeds the blunt parent, and
  every cut along it keeps the 2D-F layout's three cells.

    python examples/blade/03_airfoils_and_root_transition/airfoils_and_root_transition.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from _shared.plotting import REFERENCE, SERIES, draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402

X = np.linspace(0.0, 1.0, 2001)
FLATBACK_TE = 0.05


def flatback():
    """A NACA 0030 opened linearly to a blunt TE of 5 % chord, given as coordinates."""
    base = Airfoil.naca4("0030")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 81)))
    upper = base.surface("upper", x) + 0.5 * FLATBACK_TE * x
    lower = base.surface("lower", x) - 0.5 * FLATBACK_TE * x
    points = np.vstack([np.column_stack([x[::-1], upper[::-1]]), np.column_stack([x[1:], lower[1:]])])
    return Airfoil(points, "flatback 0030, 5 % TE")


def check_profiles():
    print(f"{'airfoil':<24} {'t/c':>7} {'expected':>8} {'max camber':>10} {'expected':>8}")
    worst = 0.0
    for code in ("0012", "0018", "0030", "2412", "4415"):
        a = Airfoil.naca4(code)
        camber = 0.5 * (a.surface("upper", X) + a.surface("lower", X))
        t, m = int(code[2:]) / 100, int(code[0]) / 100
        worst = max(worst, abs(a.max_thickness - t), abs(camber.max() - m))
        print(f"{a.name:<24} {a.max_thickness:7.4f} {t:8.4f} {camber.max():10.4f} {m:8.4f}")
    fb = flatback()
    te = float(fb.thickness(1.0))
    print(f"{fb.name:<24} TE thickness {te:.4f} (built with {FLATBACK_TE})")
    worst = max(worst, abs(te - FLATBACK_TE))
    assert worst < 1e-3
    return worst


def root_blade():
    return BladeGeometry([Station(0.0, 3.0, Airfoil.ellipse(1.0), pitch_axis=0.5),
                          Station(4.0, 3.4, Airfoil.naca4("0030"), twist=10.0, pitch_axis=0.35),
                          Station(12.0, 2.6, Airfoil.naca4("0021"), twist=5.0, pitch_axis=0.3)],
                         "2D-F-LC-0R", 0.015, pivot_fraction=0.35,
                         spar_cap_thickness=0.04, shear_web_thickness=0.02, web_positions=(0.2, 0.5))


def check_transition(geo):
    z = np.linspace(0.0, 12.0, 49)
    tc, te_gap = [], []
    st = geo.stations
    for zi in z:
        _chord, *_rest, (_x, upper, lower) = geo.interpolate(zi)
        tc.append(float(np.max(upper - lower)))
        te_gap.append(abs(float(upper[-1]) - float(lower[-1])))
        k = 0 if zi <= st[1].z else 1
        gap_lo = abs(float(st[k].airfoil.surface("upper", 1.0) - st[k].airfoil.surface("lower", 1.0)))
        gap_hi = abs(float(st[k + 1].airfoil.surface("upper", 1.0) - st[k + 1].airfoil.surface("lower", 1.0)))
        assert te_gap[-1] <= max(gap_lo, gap_hi) + 1e-6
    tc = np.array(tc)
    root = z <= st[1].z + 1e-12
    assert np.all(np.diff(tc[root]) <= 1e-6)
    print(f"\nroot transition: t/c {tc[0]:.3f} → {tc[np.where(root)[0][-1]]:.3f} (monotone), "
          f"max TE gap {max(te_gap):.2e}")
    cells = {zi: len(geo.section_at(zi).cells) for zi in (0.0, 1.0, 2.0, 3.0, 4.0, 8.0, 12.0)}
    print("cells along the span:", cells)
    assert set(cells.values()) == {3}
    linear = np.interp(z, [s.z for s in st], [s.airfoil.max_thickness for s in st])
    return z, tc, linear


def figure(geo, z, tc, linear, worst):
    fig = plt.figure(figsize=(14, 6.2))
    grid = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.0])
    ax = fig.add_subplot(grid[0, :2])
    for i, a in enumerate([Airfoil.naca4("0012"), Airfoil.naca4("2412"), Airfoil.naca4("4415"), flatback(),
                           Airfoil.ellipse(0.6)]):
        for side, ls in (("upper", "-"), ("lower", "-")):
            ax.plot(X, a.surface(side, X), ls, color=SERIES[i], lw=1.2, label=a.name if side == "upper" else None)
    ax.set_aspect("equal")
    ax.set_title("Airfoil.naca4, Airfoil(points) and Airfoil.ellipse")
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))

    ax = fig.add_subplot(grid[0, 2:])
    ax.plot(z, tc, color=SERIES[0], lw=1.5, label="interpolated surfaces")
    ax.plot(z, linear, "--", color=REFERENCE, lw=0.9, label="linear blend of the stations' t/c")
    for s in geo.stations:
        ax.axvline(s.z, color="0.8", lw=0.6)
        ax.text(s.z + 0.1, 0.95, s.airfoil.name, fontsize=7.5, va="top")
    ax.set_xlabel("span z [m]")
    ax.set_ylabel("maximum thickness / chord")
    ax.set_title("Root transition: circle to NACA 0030 to NACA 0021")
    ax.legend(loc="center right")

    for i, zc in enumerate((0.0, 1.5, 4.0, 12.0)):
        ax = fig.add_subplot(grid[1, i])
        draw_cut(ax, geo.section_at(zc))
        ax.set_xlim(-1.8, 2.6)
        ax.set_ylim(-1.7, 1.7)
        ax.set_title(f"z = {zc:g} m", fontsize=8.5)
    note(fig, f"NACA t/c and camber, and the flatback's TE, match their definitions to {worst:.0e}. "
              "2D-F, skin 15 mm, caps 40 mm, webs 20 mm at 20 % and 50 % chord.")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, __file__, "airfoils_and_root_transition")


if __name__ == "__main__":
    worst = check_profiles()
    geo = root_blade()
    z, tc, linear = check_transition(geo)
    figure(geo, z, tc, linear, worst)
