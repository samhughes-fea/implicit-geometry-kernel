"""Station network and FeatureTable ids from a loft cut (no Gmsh).

    python examples/blade/11_feature_station_network/feature_station_network.py
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
    Station,
    station_network,
    station_network_from_meta,
    wall_meta_from_cut,
)


def blade() -> BladeGeometry:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
        shear_web_thickness=0.004, web_positions=(0.18, 0.48), pivot_fraction=0.3,
    )


def main() -> None:
    geo = blade()
    z = 4.0
    cut = geo.section_at(z)
    net = station_network(cut)
    field, table = wall_meta_from_cut(cut)
    print(f"z={z:g}: {len(net.walls)} network walls, loop_features={net.loop_features}")
    for w in net.walls:
        print(f"  id={w.feature_id} kind={w.kind!r} label={w.label!r} n={len(w.points)} T={w.thickness:.4f}")
    xs = np.concatenate([w.xy[:, 0] for w in cut.walls])
    ys = np.concatenate([w.xy[:, 1] for w in cut.walls])
    pad = 0.05
    bounds = (float(xs.min()) - pad, float(xs.max()) + pad, float(ys.min()) - pad, float(ys.max()) + pad)
    from_meta = station_network_from_meta(geo.wall_mid_sdf(), field, z, bounds, table, n=36)
    print(f"from wall_meta: {len(from_meta.walls)} branches, features={from_meta.loop_features}")
    assert set(net.loop_features) <= set(range(len(table)))

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    for w in net.walls:
        ax.plot(w.points[:, 0], w.points[:, 1], color=WALL_COLOURS.get(w.kind, "#888"), lw=1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("station_network: labelled shell runs")
    legend_outside(ax, handles=[
        Line2D([0], [0], color=WALL_COLOURS["skin"], label="skin"),
        Line2D([0], [0], color=WALL_COLOURS["cap"], label="cap"),
        Line2D([0], [0], color=WALL_COLOURS["web"], label="web"),
    ])
    note(fig, "One colour per wall kind. Feature ids are printed, not stacked on the trailing edge.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "feature_station_network")


if __name__ == "__main__":
    main()
