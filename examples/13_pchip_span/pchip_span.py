"""Compare linear and PCHIP chord laws; print a root-cylinder blend.

    python examples/13_pchip_span/pchip_span.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _shared.plotting import plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402


def main():
    foil = Airfoil.naca4("0018")
    stations = [Station(0.0, 2.0, foil, pitch_axis=0.3),
                Station(5.0, 1.2, foil, pitch_axis=0.3),
                Station(10.0, 1.0, foil, pitch_axis=0.3)]
    linear = BladeGeometry(stations, "2D-F", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004)
    pchip = BladeGeometry(stations, "2D-F", 0.004, spar_cap_thickness=0.008, shear_web_thickness=0.004,
                          span_kind="pchip")
    z = np.linspace(0.0, 10.0, 81)
    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    ax.plot(z, [linear.interpolate(zz)[0] for zz in z], label="linear")
    ax.plot(z, [pchip.interpolate(zz)[0] for zz in z], label="PCHIP")
    ax.plot([s.z for s in stations], [s.chord for s in stations], "ko", label="stations")
    ax.set_xlabel("z / m")
    ax.set_ylabel("chord / m")
    ax.legend()
    ax.set_title("Planform chord: linear vs monotone PCHIP")
    fig.tight_layout()
    save(fig, __file__, "chord")
    for s in (0.0, 0.5, 1.0):
        blended = Airfoil.cylinder_to_airfoil(foil, s)
        print(f"s = {s:g}: {blended.name}, t/c = {blended.max_thickness:.4f}")


if __name__ == "__main__":
    main()
