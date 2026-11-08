"""Write lofted mid-surface patches next to the OML.

    python examples/blade/07_mid_surfaces/mid_surfaces.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402


def main():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    blade = BladeGeometry(
        stations,
        "2D-F-LC-0R",
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48), pivot_fraction=0.3)
    out = Path(__file__).resolve().parent / "figures"
    out.mkdir(exist_ok=True)
    oml = blade.write_oml(out / "oml.stl", n_span=8, n_chord=24)
    mids = blade.write_mid_surfaces(out / "mid_surfaces.obj", n_span=8, n_chord=24)
    print("wrote", Path(oml).name, "and", Path(mids).name)


if __name__ == "__main__":
    main()
