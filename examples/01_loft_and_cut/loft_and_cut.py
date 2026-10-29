"""Loft a section layout along a tapered, twisted wing and cut it at any span position.

Writes figures/layouts.png (eight of the 40 layouts on one airfoil) and
figures/span_cuts.png (a 2D-F wing cut at five span positions, in the blade frame).

    python examples/01_loft_and_cut/loft_and_cut.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from shapely.geometry import LineString  # noqa: E402

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402

FIGURES = Path(__file__).resolve().parent / "figures"
COLOURS = {"skin": "#2a6fdb", "cap": "#d1495b", "web": "#2e9e6a"}


def draw(ax, cut):
    """Each wall as its midline and a band of its thickness."""
    for wall in cut.walls:
        band = LineString(wall.xy).buffer(0.5 * wall.thickness, cap_style=2, join_style=2)
        for poly in getattr(band, "geoms", [band]):
            ax.fill(*poly.exterior.xy, color=COLOURS[wall.kind], alpha=0.35, lw=0)
        ax.plot(wall.xy[:, 0], wall.xy[:, 1], color=COLOURS[wall.kind], lw=1.0)
    ax.set_aspect("equal")


def layouts():
    airfoil = Airfoil.naca4("2418")
    stations = [Station(0.0, 2.0, airfoil, pitch_axis=0.3), Station(1.0, 2.0, airfoil, pitch_axis=0.3)]
    keys = ["0A", "1A-F", "2A-F", "2B-F", "2D-F", "2D-CN", "3D-F", "4C-F"]
    fig, axes = plt.subplots(4, 2, figsize=(9, 6.4))
    for ax, key in zip(axes.ravel(), keys):
        cut = BladeGeometry(stations, key, 0.008, spar_cap_thickness=0.03, shear_web_thickness=0.015).section_at(0.5)
        draw(ax, cut)
        ax.set_title(f"{key}: {len(cut.walls_of('web'))} webs, {len(cut.cells)} cells", fontsize=9)
        ax.axis("off")
    fig.suptitle("Eight of the 40 layouts on a NACA 2418, 2 m chord (skin blue, caps red, webs green)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "layouts.png", dpi=150)
    plt.close(fig)


def span_cuts():
    stations = [Station(0.0, 3.0, Airfoil.naca4("0030"), twist=12.0, pitch_axis=0.3),
                Station(20.0, 1.2, Airfoil.naca4("0018"), twist=0.0, pitch_axis=0.3, offset=(0.0, 0.4))]
    wing = BladeGeometry(stations, "2D-F", skin_thickness=0.01, spar_cap_thickness=lambda z: 0.04 - 0.0015 * z,
                         shear_web_thickness=0.015, web_positions=(0.2, 0.5))
    fig, axes = plt.subplots(1, 5, figsize=(13, 2.5), sharey=True)
    for ax, z in zip(axes, (0.0, 5.0, 10.0, 15.0, 20.0)):
        cut = wing.section_at(z)
        draw(ax, cut)
        ax.plot(0.0, 0.0, "k+")
        ax.set_title(f"z = {z:g} m, chord {cut.chord:.2f} m, twist {cut.twist:.1f}°", fontsize=8.5)
        ax.set_xlim(-1.0, 2.2)
        ax.set_ylim(-0.6, 0.6)
        ax.tick_params(labelsize=7)
    fig.suptitle("2D-F wing cut along the span: flapwise webs stay vertical as the section twists", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "span_cuts.png", dpi=150)
    plt.close(fig)
    return wing


if __name__ == "__main__":
    FIGURES.mkdir(exist_ok=True)
    layouts()
    wing = span_cuts()
    for z in (0.0, 10.0, 20.0):
        cut = wing.section_at(z)
        print(f"z = {z:4.1f} m: {len(cut.walls)} walls, {len(cut.cells)} cells, "
              f"enclosed area {sum(c.area for c in cut.cells):.3f} m^2")
