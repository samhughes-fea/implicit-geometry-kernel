"""Partial-span stringer: a free-end branch on the shell graph between two stations.

Writes figures/stringer.png (cut at z = 4 m with the stringer drawn).

    python examples/structure/02_partial_span_stringer/partial_span_stringer.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import SERIES, draw_cut, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, Stringer  # noqa: E402
from implicit_geometry_kernel.blade.shell_graph import topology_diff  # noqa: E402


def main():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    stringer = Stringer("str_le", position=0.12, thickness=0.008, depth=0.04,
                        extent=(2.0, 6.0), side="pressure")
    blade = BladeGeometry(stations, "2D-F-LC-0R", 0.003, spar_cap_thickness=0.008,
                          shear_web_thickness=0.002, web_positions=(0.18, 0.48),
                          stringers=(stringer,), pivot_fraction=0.3)
    inside = blade.section_at(4.0)
    outside = blade.section_at(8.0)
    run = inside.shell_graph().branches[0]
    print(f"forced stations: {blade.forced_stations()}")
    print(f"z = 4 m: {len(inside.shell_graph().branches)} stringer branch, length "
          f"{float(np.linalg.norm(run.points[1, :2] - run.points[0, :2])):.3f} m")
    print(f"z = 8 m: {len(outside.shell_graph().branches)} stringer branches")
    events = topology_diff(inside.shell_graph(), outside.shell_graph()).events
    print("topology events:", ", ".join(f"{e.kind}:{e.name}" for e in events if "stringer" in e.kind))

    foot, tip = run.points[0], run.points[-1]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))
    ax = axes[0]
    draw_cut(ax, inside)
    ax.plot(run.points[:, 0], run.points[:, 1], color=SERIES[5], lw=2.0, marker="o", ms=4, label="stringer")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("2D-F cut inside the stringer span")
    legend_outside(ax)
    ax = axes[1]
    draw_cut(ax, inside)
    ax.plot(run.points[:, 0], run.points[:, 1], color=SERIES[5], lw=2.4, marker="o", ms=6)
    ax.annotate("free end, 40 mm depth", xy=(tip[0], tip[1]),
                xytext=(tip[0] - 0.02, tip[1] + 0.04),
                fontsize=7, ha="right",
                arrowprops={"arrowstyle": "->", "lw": 0.6})
    pad = 0.08
    ax.set_xlim(min(foot[0], tip[0]) - pad, max(foot[0], tip[0]) + pad)
    ax.set_ylim(min(foot[1], tip[1]) - pad, max(foot[1], tip[1]) + pad)
    ax.set_xlabel("x / m")
    ax.set_title("zoom on the free end")
    note(fig, "Stringer on the pressure skin at 12% chord; 8 mm thick, 40 mm stand-off.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "stringer")


if __name__ == "__main__":
    main()
