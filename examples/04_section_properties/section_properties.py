"""Areas and second moments of a section's walls: exact solids against closed forms and midlines.

``SectionGeometry(..., outputs=("properties",))`` integrates every wall solid of
a section (skins, caps, webs, and the cavities between them, labelled
``sandwich_core_*``) and returns its area, centroid and second moments about its own
centroid in ``properties_map``. Downstream models often use the thin-walled idealisation
instead: each wall reduced to its midline (``BladeGeometry.section_at(z).walls``) times
its thickness. This example checks both.

1. A circular tube, D = 1 m with a 10 mm skin (layout 0A, then 2A-F with two 10 mm
   webs): the skin's area and Ixx against pi (R^2 - r^2) and pi (R^4 - r^4) / 4, and each
   web against the strip of the inner circle it fills.
2. A NACA 2418 at 2 m chord with box spar caps (2D-F), at wall thicknesses scaled from
   1 to 8 times a baseline of 4 mm skin, 12 mm caps and 6 mm webs: the whole section's
   area, centroid, Ixx, Iyy and principal angle from the solids and from the midlines.
   The midline model counts each junction twice (a web runs to the skin midline, the
   solid web stops at the skin's inner face) and drops the t^3 terms, so its error
   grows linearly with the wall thickness.

    python examples/04_section_properties/section_properties.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from _shared.plotting import REFERENCE, SERIES, TEXT_SECONDARY, draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.sdf.airfoil import AirfoilSDF
from implicit_geometry_kernel.section import (
    SectionGeometry,
    resolve_system_type,
)

CHORD = 2.0
BASE = {"skin": 0.004, "cap": 0.012, "web": 0.006}
SCALES = (1, 2, 4, 8)  # below 1 the solids would hit section_geometry's minimum wall thicknesses


# --------------------------------------------------------------- integrals
def solid_props(layout, vertices, chord, skin, cap=0.0, web=0.0):
    """Per-label properties of the wall solids (cavities left out), in the section frame (x from the LE)."""
    result = SectionGeometry(layout=resolve_system_type(layout), airfoil=AirfoilSDF(vertices, chord=chord),
                                   skin_thickness=skin, spar_cap_thickness=cap, shear_web_thickness=web,
                                   outputs=("properties",)).execute().get_results()
    return {k: v for k, v in result.properties_map.items() if not k.startswith("sandwich_core")}


def totals_from_solids(props):
    """Area, first moments and second moments about the origin of all walls together."""
    A = Sx = Sy = Ixx = Iyy = Ixy = 0.0
    for p in props.values():
        a, cx, cy = p["area"], p["cx"], p["cy"]
        A += a
        Sx += a * cx
        Sy += a * cy
        Ixx += p["Ixx"] + a * cy * cy
        Iyy += p["Iyy"] + a * cx * cx
        Ixy += p["Ixy"] + a * cx * cy
    return A, Sx, Sy, Ixx, Iyy, Ixy


def totals_from_midlines(cut):
    """The thin-walled model: each wall a line of its midline with thickness t (t^3 terms dropped)."""
    A = Sx = Sy = Ixx = Iyy = Ixy = 0.0
    for w in cut.walls:
        p, t = w.xy, w.thickness
        a, b = p[:-1], p[1:]
        L = np.linalg.norm(b - a, axis=1)
        A += t * L.sum()
        Sx += t * np.sum(L * 0.5 * (a[:, 0] + b[:, 0]))
        Sy += t * np.sum(L * 0.5 * (a[:, 1] + b[:, 1]))
        Ixx += t * np.sum(L * (a[:, 1] ** 2 + a[:, 1] * b[:, 1] + b[:, 1] ** 2) / 3)
        Iyy += t * np.sum(L * (a[:, 0] ** 2 + a[:, 0] * b[:, 0] + b[:, 0] ** 2) / 3)
        Ixy += t * np.sum(L * (2 * a[:, 0] * a[:, 1] + a[:, 0] * b[:, 1] + b[:, 0] * a[:, 1] + 2 * b[:, 0] * b[:, 1]) / 6)
    return A, Sx, Sy, Ixx, Iyy, Ixy


def centroidal(tot):
    """Area, centroid, centroidal Ixx, Iyy, Ixy and the principal angle [deg]."""
    A, Sx, Sy, Ixx, Iyy, Ixy = tot
    cx, cy = Sx / A, Sy / A
    ixx, iyy, ixy = Ixx - A * cy * cy, Iyy - A * cx * cx, Ixy - A * cx * cy
    # the principal axis nearest the chord (the flapwise bending axis), as an angle from the chord
    vals, vecs = np.linalg.eigh(np.array([[iyy, ixy], [ixy, ixx]]))
    v = vecs[:, int(np.argmax(vals))]
    theta = float(np.degrees(np.arctan2(v[1], v[0])))
    theta = (theta + 90.0) % 180.0 - 90.0
    return A, cx, cy, ixx, iyy, ixy, theta


# --------------------------------------------------------------- 1. circular tube
def strip_area(r, d, w):
    """Area of the vertical strip |x - d| < w/2 inside a circle of radius r (the web between the skin's inner faces)."""
    def F(x):  # antiderivative of 2 sqrt(r^2 - x^2)
        return x * np.sqrt(r * r - x * x) + r * r * np.arcsin(x / r)
    return F(d + 0.5 * w) - F(d - 0.5 * w)


def tube():
    D, t, tw = 1.0, 0.010, 0.010
    R, r = 0.5 * D, 0.5 * D - t
    phi = np.linspace(0.0, 2 * np.pi, 1441)[:-1]
    verts = np.column_stack([0.5 + 0.5 * np.cos(phi), 0.5 * np.sin(phi)]) * D
    print("1. Circular tube, D = 1 m, 10 mm skin")
    print(f"{'quantity':<44} {'solid':>12} {'closed form':>12} {'error':>9}")
    rows = []
    p = solid_props("0A", verts, D, t)
    skin = [v for k, v in p.items() if k.startswith("skin")]
    A = sum(v["area"] for v in skin)
    Ixx = sum(v["Ixx"] + v["area"] * v["cy"] ** 2 for v in skin)
    rows += [("0A skin area [m^2]", A, np.pi * (R**2 - r**2)), ("0A skin Ixx [m^4]", Ixx, 0.25 * np.pi * (R**4 - r**4))]
    p = solid_props("2A-F", verts, D, t, web=tw)
    for k in sorted(k for k in p if k.startswith("shear_web")):
        d = p[k]["cx"] - 0.5 * D  # web position from the centre
        rows.append((f"2A-F {k} at x/c = {p[k]['cx']:.2f}, area [m^2]", p[k]["area"], strip_area(r, d, tw)))
    worst = 0.0
    for name, got, ref in rows:
        worst = max(worst, abs(got / ref - 1))
        print(f"{name:<44} {got:12.6e} {ref:12.6e} {got / ref - 1:+9.1e}")
    assert worst < 1e-4
    return worst


# --------------------------------------------------------------- 2. NACA 2418, 2D-F
def airfoil_section():
    airfoil = Airfoil.naca4("2418")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 161)))
    verts = np.vstack([np.column_stack([x[::-1], airfoil.surface("upper", x)[::-1]]),
                       np.column_stack([x[1:], airfoil.surface("lower", x)[1:]])])[:-1] * CHORD
    # pitch axis at the LE and no twist: the cut's blade frame is the solid's section frame
    stations = [Station(0.0, CHORD, airfoil, pitch_axis=0.0), Station(1.0, CHORD, airfoil, pitch_axis=0.0)]
    rows, cuts = [], {}
    print("\n2. NACA 2418, 2 m chord, 2D-F: whole section, solids against midlines")
    print(f"{'scale':>5} {'skin [mm]':>9} {'A solid':>10} {'A error':>8} {'Ixx error':>9} {'Iyy error':>9} "
          f"{'theta solid':>11} {'theta mid':>9}")
    for s in SCALES:
        t = {k: v * s for k, v in BASE.items()}
        solid = centroidal(totals_from_solids(solid_props("2D-F", verts, CHORD, t["skin"], t["cap"], t["web"])))
        geo = BladeGeometry(stations, "2D-F", t["skin"], spar_cap_thickness=t["cap"], shear_web_thickness=t["web"])
        cut = geo.section_at(0.5)
        cuts[s] = cut
        mid = centroidal(totals_from_midlines(cut))
        err = [mid[i] / solid[i] - 1 for i in (0, 3, 4)]
        rows.append((s, t["skin"], solid, mid, err))
        print(f"{s:5g} {t['skin'] * 1e3:9.1f} {solid[0]:10.4e} {err[0]:+8.2%} {err[1]:+9.2%} {err[2]:+9.2%} "
              f"{solid[6]:10.2f}° {mid[6]:8.2f}°")
    # thin-wall error is first order in t: halving the walls halves it
    ratio = rows[1][4][0] / rows[0][4][0]
    print(f"area error ratio, scale 2 to 1: {ratio:.2f} (first order: 2)")
    assert 1.6 < ratio < 2.4
    return rows, cuts


def figure(tube_err, rows, cuts):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw={"width_ratios": [1.5, 1.0]})
    ax = axes[0]
    cut = cuts[4]
    draw_cut(ax, cut, alpha=0.5)
    s, _, solid, mid, _ = rows[SCALES.index(4)]
    for (A, cx, cy, ixx, iyy, ixy, th), colour, label, ls in ((solid, REFERENCE, "solids", "-"),
                                                              (mid, SERIES[1], "midlines", "--")):
        d = np.array([np.cos(np.radians(th)), np.sin(np.radians(th))])
        for v, length in ((d, 0.6), (np.array([-d[1], d[0]]), 0.2)):
            ax.plot([cx - length * v[0], cx + length * v[0]], [cy - length * v[1], cy + length * v[1]],
                    ls, color=colour, lw=1.0)
        ax.plot(cx, cy, "o", ms=4, color=colour, label=f"{label}: centroid ({cx:.3f}, {cy:.4f}) m, θ = {th:.2f}°")
    ax.set_title("NACA 2418, 2 m chord, 2D-F at 4x the baseline walls: centroid and principal axes")
    ax.legend(loc="lower right")
    ax.set_xlabel("x from the leading edge [m]")
    ax.set_ylabel("y [m]")

    ax = axes[1]
    t = np.array([r[1] for r in rows]) * 1e3
    for i, name in enumerate(("area", "Ixx", "Iyy")):
        ax.plot(t, [100 * r[4][i] for r in rows], "o-", color=SERIES[i], label=name)
    ax.plot([0, t[-1]], [0, 100 * rows[-1][4][0]], "--", color=REFERENCE, lw=0.8, label="proportional to t")
    ax.axhline(0.0, color=TEXT_SECONDARY, lw=0.5)
    ax.set_xlim(0, None)
    ax.set_xlabel("skin thickness [mm] (caps 3x, webs 1.5x)")
    ax.set_ylabel("midline / solid − 1 [%]")
    ax.set_title("The midline model's area error grows in proportion to t")
    ax.legend()
    note(fig, f"Solid wall properties match the circular tube's closed forms to {tube_err:.0e}. "
              "The midline model counts junctions twice and drops t³ terms.")
    fig.text(0.99, 0.005, "skin blue, caps red, webs green", fontsize=7.5, color=TEXT_SECONDARY, ha="right")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, __file__, "section_properties")


if __name__ == "__main__":
    err = tube()
    rows, cuts = airfoil_section()
    figure(err, rows, cuts)
