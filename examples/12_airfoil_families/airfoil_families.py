"""Plot NACA, CST, t/c family and a flatback on one chord.

    python examples/12_airfoil_families/airfoil_families.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _shared.plotting import plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil  # noqa: E402


def main():
    foils = [
        Airfoil.naca4("2412"),
        Airfoil.naca5("23012"),
        Airfoil.cst([0.18, 0.16, 0.14], [0.16, 0.14, 0.12], name="CST"),
        Airfoil.naca4("0012").with_thickness(0.18),
        Airfoil.naca4("0012").with_te_thickness(0.04),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    x = np.linspace(0.0, 1.0, 200)
    for foil in foils:
        yu, yl = foil.surface("upper", x), foil.surface("lower", x)
        ax.plot(x, yu, lw=1.0, label=foil.name)
        ax.plot(x, yl, lw=1.0, color=ax.lines[-1].get_color())
        print(f"{foil.name}: t/c = {foil.max_thickness:.4f}")
    ax.set_aspect("equal")
    ax.set_xlabel("x/c")
    ax.set_ylabel("y/c")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_title("NACA, CST, thickness family and flatback")
    fig.tight_layout()
    save(fig, __file__, "families")


if __name__ == "__main__":
    main()
