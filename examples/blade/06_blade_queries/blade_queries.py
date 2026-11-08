"""Ply-drop thickness, a point query, OML STL and validate().

    python examples/blade/06_blade_queries/blade_queries.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from implicit_geometry_kernel import Airfoil, BladeGeometry, PlyDrop, Station  # noqa: E402


def main():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    drop = PlyDrop("spar_buildup", s0=0.2, s1=0.5, z0=2.0, z1=8.0, thickness=0.01, taper=0.2)
    blade = BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                          shear_web_thickness=0.002, web_positions=(0.18, 0.48), ply_drops=(drop,), pivot_fraction=0.3)
    print("T(z=5, s=0.35) =", float(blade.thickness_at(5.0, 0.35)))
    print("key points:", sorted(blade.key_points(5.0)))
    p = blade.section_at(5.0).walls[0].points[4]
    hit = blade.query(p.reshape(1, 3))
    print("query member", hit["member"][0], "distance", float(hit["distance"][0]))
    out = Path(__file__).resolve().parent / "figures" / "oml.stl"
    out.parent.mkdir(exist_ok=True)
    blade.write_oml(out, n_span=12, n_chord=24)
    print("wrote", out.relative_to(Path(__file__).resolve().parents[2]), "issues", blade.validate())


if __name__ == "__main__":
    main()
