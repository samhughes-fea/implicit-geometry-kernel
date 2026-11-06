"""A two-station blade with one midspan rib: closed faces on the shell graph.

Writes figures/rib_face.png (OML, walls and rib cell faces at z = 5 m).

    python examples/10_transverse_rib/transverse_rib.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _shared.plotting import draw_cut, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, TransverseMember  # noqa: E402
from implicit_geometry_kernel.blade.shell_graph import topology_diff  # noqa: E402


def main():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    rib = TransverseMember("rib_mid", z=5.0, thickness=0.012)
    blade = BladeGeometry(
        stations,
        "2D-F",
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48),
        transverse=(rib,),
    )
    at_rib = blade.section_at(5.0)
    other = blade.section_at(2.0)
    faces = at_rib.shell_graph().faces
    print(f"forced stations: {blade.forced_stations()}")
    print(f"z = 5 m: {len(faces)} rib faces, {len(at_rib.cells)} cells")
    print(f"z = 2 m: {len(other.shell_graph().faces)} rib faces")
    events = topology_diff(at_rib.shell_graph(), other.shell_graph()).events
    print("topology events:", ", ".join(f"{e.kind}:{e.name}" for e in events if "rib" in e.kind))

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    draw_cut(ax, at_rib)
    oml = blade.oml_polyline(5.0)
    ax.plot(oml[:, 0], oml[:, 1], color="#0b0b0b", lw=0.7, ls="--", label="OML")
    for face in faces:
        ax.plot(face.points[:, 0], face.points[:, 1], color="#eda100", lw=1.4, label=face.name)
    ax.set_aspect("equal")
    ax.set_title("2D-F cut at the rib: closed cell faces (amber) inside the OML")
    handles, labels = ax.get_legend_handles_labels()
    by = dict(zip(labels, handles))
    ax.legend(by.values(), by.keys(), loc="upper right", fontsize=7)
    fig.tight_layout()
    save(fig, __file__, "rib_face")


if __name__ == "__main__":
    main()
