"""All 40 structural layouts on one twisted airfoil, with a consistency check of each.

A layout key reads ``{webs}{cap family}-{web orientation}``:

    webs          0 to 5 shear webs, at the layout's default chord fractions
    cap family    A no spar caps, B one cap band across the webs, C a discrete cap
                  on each web, D a continuous box cap between the outer webs (two
                  or more webs); 0B has caps but no webs
    orientation   F webs stay flapwise (normal to the rotor plane) as the section
                  twists, CN webs stay normal to the chord and turn with it

Every layout is cut from the same NACA 2418 at 2 m chord and 15 degrees of twist, so
the F and CN variants differ visibly. For each cut the script checks the counts of
webs, caps and cells against the layout's spec, and checks that the cells tile the
region inside the skin midline: their areas add up to the area of the outer mould line
offset inward by half the skin thickness, whatever the webs.

    python examples/section/01_layout_catalogue/layout_catalogue.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import SYSTEM_TYPE_LAYOUT_KEYS, Airfoil, BladeGeometry, Station, resolve_system_type  # noqa: E402
from implicit_geometry_kernel.section.boundary import OPS, Polygon  # noqa: E402

CHORD, TWIST, PITCH = 2.0, 15.0, 0.3
SKIN, CAP, WEB = 0.008, 0.03, 0.015
DISPLAY = 4
FAMILIES = ("A", "B", "C", "D")
ORIENTATIONS = ("F", "CN")


def blade():
    airfoil = Airfoil.naca4("2418")
    return [Station(0.0, CHORD, airfoil, twist=TWIST, pitch_axis=PITCH),
            Station(1.0, CHORD, airfoil, twist=TWIST, pitch_axis=PITCH)]


def skin_midline_area(airfoil: Airfoil) -> float:
    """Reference: the outer mould line offset inward by half the skin thickness."""
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 400)))
    upper, lower = airfoil.surface("upper", x), airfoil.surface("lower", x)
    outline = np.vstack([np.column_stack([x[::-1], upper[::-1]]), np.column_stack([x[1:], lower[1:]])]) * CHORD
    return OPS.offset_inward(Polygon(outline), 0.5 * SKIN).area


def cut_all(stations, skin=SKIN, cap=CAP, web=WEB):
    cuts = {}
    for key in SYSTEM_TYPE_LAYOUT_KEYS:
        geo = BladeGeometry(stations, key, skin, spar_cap_thickness=cap, shear_web_thickness=web, pivot_fraction=0.3)
        cuts[key] = geo.section_at(0.5)
    return cuts


def check(cuts, reference_area):
    print(f"{'layout':<7} {'webs':>4} {'caps':>4} {'cells':>5} {'cell area [m^2]':>16} {'vs offset OML':>13}  description")
    worst = 0.0
    for key, cut in cuts.items():
        spec = resolve_system_type(key)
        webs, caps = len(cut.walls_of("web")), len(cut.walls_of("cap"))
        assert webs == spec.n_webs, key
        assert len(cut.cells) == spec.n_cells == spec.n_webs + 1, key
        assert (caps > 0) == spec.has_spar_caps, key
        area = sum(c.area for c in cut.cells)
        err = area / reference_area - 1.0
        worst = max(worst, abs(err))
        print(f"{key:<7} {webs:>4} {caps:>4} {len(cut.cells):>5} {area:16.5f} {err:+13.3%}  {spec.description}")
    print(f"\nreference area inside the skin midline: {reference_area:.5f} m^2; worst difference {worst:.3%}")
    assert worst < 0.005
    return worst


def figure(cuts, worst):
    fig, axes = plt.subplots(6, 8, figsize=(16, 8.4))
    for ax in axes.ravel():
        ax.axis("off")
    # row 0: the two web-less layouts; rows 1-5: 1 to 5 webs, columns A-F, A-CN, B-F, ...
    placed = {}
    for key in cuts:
        if key.startswith("0A"):
            placed[key] = (0, 0)
        elif key.startswith("0B"):
            placed[key] = (0, 2)
        else:
            n = int(key[0])
            fam = key[1]
            ori = "CN" if "-CN-" in key else "F"
            placed[key] = (n, 2 * FAMILIES.index(fam) + ORIENTATIONS.index(ori))
    for key, (r, c) in placed.items():
        ax = axes[r, c]
        draw_cut(ax, cuts[key], midlines=True, alpha=0.55)
        ax.set_xlim(-0.75, 1.45)
        ax.set_ylim(-0.42, 0.42)
        ax.set_title(key, fontsize=8, pad=1)
    for i, fam in enumerate(FAMILIES):
        label = {"A": "A: no caps", "B": "B: one cap band", "C": "C: a cap per web", "D": "D: box caps"}[fam]
        fig.text((2 * i + 1) / 8, 0.955, label, ha="center", fontsize=9)
    fig.suptitle("The 40 layouts on a NACA 2418, 2 m chord, 15° twist: flapwise (F) webs stay vertical, "
                 "chord-normal (CN) webs turn with the chord", fontsize=10, y=0.995)
    note(fig, f"Display thickness {DISPLAY}× (skin {DISPLAY * SKIN * 1e3:g} mm, cap {DISPLAY * CAP * 1e3:g} mm, "
              f"web {DISPLAY * WEB * 1e3:g} mm). Printed web counts use the real gauges. "
              f"Cells tile the skin-midline area to within {worst:.2%}.")
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    save(fig, __file__, "layout_catalogue")


if __name__ == "__main__":
    stations = blade()
    cuts = cut_all(stations)
    worst = check(cuts, skin_midline_area(stations[0].airfoil))
    figure(cut_all(stations, DISPLAY * SKIN, DISPLAY * CAP, DISPLAY * WEB), worst)
