"""rib_face_polygons and partition_edge_curve on a loft with a midspan rib.

    python examples/structure/04_partition_rib_faces/partition_rib_faces.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from _shared.plotting import SERIES, draw_cut, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import (  # noqa: E402
    Airfoil,
    BladeGeometry,
    Station,
    TransverseMember,
    partition_edge_curve,
    rib_face_polygons,
)

CELL_LABELS = ("LE cell", "box", "TE cell")


def main() -> None:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
        shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3,
        transverse=(TransverseMember("rib_mid", z=5.0, thickness=0.012),),
    )
    faces = rib_face_polygons(geo, n_chord=8)
    print(f"rib_face_polygons: {len(faces)}")
    for name, ring in faces:
        print(f"  {name}: n={len(ring)} z={ring[0, 2]:.2f}")
    cut = geo.section_at(5.0)
    edges = cut.shell_graph().partition_edges()
    print(f"partition edges at rib: {len(edges)}")
    if edges:
        e = edges[0]
        zs = np.array([2.0, 5.0, 8.0])
        curve = partition_edge_curve(geo, e.name_a, e.name_b, zs)
        print(f"  {e.name_a} ∩ {e.name_b}: {len(curve)} span samples")

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    for i, (name, ring) in enumerate(faces):
        colour = SERIES[i % len(SERIES)]
        label = CELL_LABELS[i] if i < len(CELL_LABELS) else name
        ax.fill(ring[:, 0], ring[:, 1], color=colour, alpha=0.35, lw=0, label=label)
        ax.plot(ring[:, 0], ring[:, 1], color=colour, lw=0.8)
    draw_cut(ax, cut)
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("rib n-gons at z = 5 m")
    legend_outside(ax)
    note(fig, "Filled rib n-gons (LE / box / TE); walls at manufacturing gauges.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "partition_rib_faces")


if __name__ == "__main__":
    main()
