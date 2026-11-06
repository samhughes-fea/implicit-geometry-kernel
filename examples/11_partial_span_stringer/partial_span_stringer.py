"""Partial-span stringer: a free-end branch on the shell graph between two stations.

Writes figures/stringer.png (cut at z = 4 m with the stringer drawn).

    python examples/11_partial_span_stringer/partial_span_stringer.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _shared.plotting import draw_cut, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, Stringer  # noqa: E402
from implicit_geometry_kernel.blade.shell_graph import topology_diff  # noqa: E402


def main():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    stringer = Stringer("str_le", position=0.12, thickness=0.008, depth=0.04,
                        extent=(2.0, 6.0), side="pressure")
    blade = BladeGeometry(stations, "2D-F", 0.003, spar_cap_thickness=0.008,
                          shear_web_thickness=0.002, web_positions=(0.18, 0.48),
                          stringers=(stringer,))
    inside = blade.section_at(4.0)
    outside = blade.section_at(8.0)
    run = inside.shell_graph().branches[0]
    print(f"forced stations: {blade.forced_stations()}")
    print(f"z = 4 m: {len(inside.shell_graph().branches)} stringer branch, length "
          f"{float(np.linalg.norm(run.points[1, :2] - run.points[0, :2])):.3f} m")
    print(f"z = 8 m: {len(outside.shell_graph().branches)} stringer branches")
    events = topology_diff(inside.shell_graph(), outside.shell_graph()).events
    print("topology events:", ", ".join(f"{e.kind}:{e.name}" for e in events if "stringer" in e.kind))

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    draw_cut(ax, inside)
    ax.plot(run.points[:, 0], run.points[:, 1], color="#6250d6", lw=1.8, marker="o", ms=3, label=run.name)
    ax.set_aspect("equal")
    ax.set_title("2D-F cut inside the stringer span: free-end branch off the pressure skin")
    ax.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    save(fig, __file__, "stringer")


if __name__ == "__main__":
    main()
