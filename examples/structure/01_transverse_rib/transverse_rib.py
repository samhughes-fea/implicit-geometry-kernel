"""A two-station blade with one midspan rib: closed faces on the shell graph.

Writes figures/rib_face.png (OML, walls and rib cell faces at z = 5 m).

    python examples/structure/01_transverse_rib/transverse_rib.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import SERIES, draw_cut, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import (  # noqa: E402
    Airfoil,
    BladeConcept,
    BladeGeometry,
    SpanwiseScheme,
    Station,
    TransverseMember,
)
from implicit_geometry_kernel.blade.shell_graph import topology_diff  # noqa: E402

CELL_LABELS = ("LE cell", "box", "TE cell")


def main():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    rib = TransverseMember("rib_mid", z=5.0, thickness=0.012)
    concept = BladeConcept("2D-F-LC-0R", SpanwiseScheme(transverse=(rib,)))
    blade = BladeGeometry.from_concept(
        concept,
        stations,
        0.003,
        spar_cap_thickness=0.008,
        shear_web_thickness=0.002,
        web_positions=(0.18, 0.48), pivot_fraction=0.3)
    at_rib = blade.section_at(5.0)
    other = blade.section_at(2.0)
    faces = at_rib.shell_graph().faces
    print(f"forced stations: {blade.forced_stations()}")
    print(f"z = 5 m: {len(faces)} rib faces, {len(at_rib.cells)} cells")
    print(f"z = 2 m: {len(other.shell_graph().faces)} rib faces")
    events = topology_diff(at_rib.shell_graph(), other.shell_graph()).events
    print("topology events:", ", ".join(f"{e.kind}:{e.name}" for e in events if "rib" in e.kind))

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    oml = blade.oml_polyline(5.0)
    ax.plot(oml[:, 0], oml[:, 1], color="#0b0b0b", lw=0.7, ls="--", label="OML")
    for i, face in enumerate(faces):
        xy = face.points[:, :2]
        colour = SERIES[i % len(SERIES)]
        label = CELL_LABELS[i] if i < len(CELL_LABELS) else face.name
        ax.fill(xy[:, 0], xy[:, 1], color=colour, alpha=0.35, lw=0, label=label)
        ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=0.8)
    draw_cut(ax, at_rib)
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("2D-F cut at the rib: closed cell faces inside the OML")
    legend_outside(ax)
    note(fig, "Filled LE / box / TE cells; walls are manufacturing gauges (skin 3 mm).")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "rib_face")


if __name__ == "__main__":
    main()
