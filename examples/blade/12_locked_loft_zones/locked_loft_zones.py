"""Locked loft plan: constant-signature zones and shared partition vertices.

    python examples/blade/12_locked_loft_zones/locked_loft_zones.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from _shared.plotting import WALL_COLOURS, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import (  # noqa: E402
    Airfoil,
    BladeGeometry,
    LockedLoftPlan,
    Station,
    resample_locked,
    shell_loft_zones,
)


def _kind(name: str, plan: LockedLoftPlan) -> str:
    if name in plan.web_names:
        return "web"
    if "cap" in name:
        return "cap"
    return "skin"


def main() -> None:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    dropped = BladeGeometry(
        stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
        shear_web_thickness=0.004, web_offsets=(0.4, 1.2), pivot_fraction=0.3,
    )
    zs = np.linspace(0.0, 10.0, 11)
    graphs = [dropped.section_at(float(z)).shell_graph() for z in zs]
    zones = shell_loft_zones(graphs)
    print(f"web_offsets loft: {len(zones)} signature zones")
    for sl in zones:
        names = tuple(r.name for r in graphs[sl][0].webs)
        print(f"  z in [{zs[sl.start]:.1f}, {zs[sl.stop - 1]:.1f}] webs={names}")

    geo = BladeGeometry(
        stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
        shear_web_thickness=0.002, web_positions=(0.18, 0.48), pivot_fraction=0.3,
    )
    g = geo.section_at(4.0).shell_graph()
    plan = LockedLoftPlan.from_graph(g, n_chord=12)
    samples = resample_locked(g, plan)
    print(f"locked plan names={plan.names}")
    print(f"resampled {len(samples)} runs")

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    for name, pts in samples.items():
        ax.plot(pts[:, 0], pts[:, 1], color=WALL_COLOURS[_kind(name, plan)], lw=1.2)
        ax.plot(pts[[0, -1], 0], pts[[0, -1], 1], "k.", ms=5)
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("resample_locked at z = 4 m")
    legend_outside(ax, handles=[
        Line2D([0], [0], color=WALL_COLOURS["skin"], label="skin"),
        Line2D([0], [0], color=WALL_COLOURS["cap"], label="cap"),
        Line2D([0], [0], color=WALL_COLOURS["web"], label="web"),
        Line2D([0], [0], color="k", marker=".", ls="none", label="shared ends"),
    ])
    note(fig, "Coloured runs; black dots are shared partition vertices.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "locked_loft_zones")


if __name__ == "__main__":
    main()
