"""Plot NACA, CST, t/c family and a flatback on one chord.

    python examples/blade/04_airfoil_families/airfoil_families.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import SERIES, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil  # noqa: E402


def _draw(ax, foils, colours):
    x = np.linspace(0.0, 1.0, 200)
    for foil, colour in zip(foils, colours):
        yu, yl = foil.surface("upper", x), foil.surface("lower", x)
        ax.plot(x, yu, color=colour, lw=1.2, label=foil.name)
        ax.plot(x, yl, color=colour, lw=1.2)
        print(f"{foil.name}: t/c = {foil.max_thickness:.4f}")
    ax.set_aspect("equal")
    ax.set_xlabel("x/c")
    ax.set_ylabel("y/c")
    legend_outside(ax)


def main():
    families = [
        Airfoil.naca4("2412"),
        Airfoil.naca5("23012"),
        Airfoil.cst([0.18, 0.16, 0.14], [0.16, 0.14, 0.12], name="CST"),
    ]
    thickness = [
        Airfoil.naca4("0012"),
        Airfoil.naca4("0012").with_thickness(0.18),
        Airfoil.naca4("0012").with_te_thickness(0.04),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 4.8), sharex=True)
    _draw(axes[0], families, SERIES[:3])
    axes[0].set_title("families: NACA 4, NACA 5, CST")
    _draw(axes[1], thickness, SERIES[3:6])
    axes[1].set_title("thickness / TE: NACA 0012, t/c = 0.18, flatback TE = 0.04c")
    note(fig, "CST is a different family, not a NACA 2412 overlay. Lower row holds thickness and TE changes.")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, __file__, "families")


if __name__ == "__main__":
    main()
