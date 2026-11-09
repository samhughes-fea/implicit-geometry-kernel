"""Loft a section layout along a tapered, twisted wing and cut it at any span position.

Writes figures/layouts.png (eight of the 40 layouts on one airfoil) and
figures/span_cuts.png (a 2D-F wing cut at five span positions, in the blade frame).

    python examples/blade/01_loft_and_cut/loft_and_cut.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import draw_cut, draw_shell, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402


def layouts():
    airfoil = Airfoil.naca4("2418")
    stations = [Station(0.0, 2.0, airfoil, pitch_axis=0.3), Station(1.0, 2.0, airfoil, pitch_axis=0.3)]
    keys = ["0A-LC-0R", "1A-F-LC-0R", "2A-F-LC-0R", "2B-F-LC-0R", "2D-F-LC-0R", "2D-CN-LC-0R", "3D-F-LC-0R", "4C-F-LC-0R"]
    fig, axes = plt.subplots(4, 2, figsize=(9, 6.4))
    for ax, key in zip(axes.ravel(), keys):
        cut = BladeGeometry(stations, key, 0.032, spar_cap_thickness=0.12, shear_web_thickness=0.06, pivot_fraction=0.3).section_at(0.5)
        n_web, n_cell = len(cut.walls_of("web")), len(cut.cells)
        webs = "web" if n_web == 1 else "webs"
        cells = "cell" if n_cell == 1 else "cells"
        draw_cut(ax, cut)
        ax.set_title(f"{key}: {n_web} {webs}, {n_cell} {cells}", fontsize=9)
        ax.axis("off")
    fig.suptitle("Eight of the 40 layouts on a NACA 2418, 2 m chord (skin navy, caps red, webs teal)", fontsize=10)
    note(fig, "Display thickness 4× (skin 32 mm, cap 120 mm, web 60 mm on a 2 m chord).")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, __file__, "layouts")


def span_cuts():
    stations = [Station(0.0, 3.0, Airfoil.naca4("0030"), twist=12.0, pitch_axis=0.3),
                Station(20.0, 1.2, Airfoil.naca4("0018"), twist=0.0, pitch_axis=0.3, offset=(0.0, 0.4))]
    wing = BladeGeometry(stations, "2D-F-LC-0R", skin_thickness=0.01, spar_cap_thickness=lambda z: 0.04 - 0.0015 * z,
                         shear_web_thickness=0.015, web_positions=(0.2, 0.5), pivot_fraction=0.3)
    fig, axes = plt.subplots(1, 5, figsize=(13, 2.5), sharey=True)
    for ax, z in zip(axes, (0.0, 5.0, 10.0, 15.0, 20.0)):
        cut = wing.section_at(z)
        draw_shell(ax, cut)
        ax.plot(0.0, 0.0, "k+")
        ax.set_title(f"z = {z:g} m, chord {cut.chord:.2f} m, twist {cut.twist:.1f}°", fontsize=8.5)
        ax.set_xlim(-1.0, 2.2)
        ax.set_ylim(-0.6, 0.6)
        ax.tick_params(labelsize=7)
        if ax is axes[0]:
            ax.set_xlabel("x / m")
            ax.set_ylabel("y / m")
    fig.suptitle("2D-F wing cut along the span: flapwise webs stay vertical as the section twists", fontsize=10)
    note(fig, "Shared axes: the tip is small because the planform tapers. Manufacturing gauges, not display thickness.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, __file__, "span_cuts")
    return wing


if __name__ == "__main__":
    layouts()
    wing = span_cuts()
    for z in (0.0, 10.0, 20.0):
        cut = wing.section_at(z)
        print(f"z = {z:4.1f} m: {len(cut.walls)} walls, {len(cut.cells)} cells, "
              f"enclosed area {sum(c.area for c in cut.cells):.3f} m^2")
