"""OML polyline versus oml_sdf / wall_mid_sdf / wall_meta at one station.

    python examples/blade/10_loft_oracles/loft_oracles.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from _shared.plotting import draw_shell, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.sdf import GeometrySolver  # noqa: E402


def blade() -> BladeGeometry:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-0R", 0.004, spar_cap_thickness=0.008,
        shear_web_thickness=0.004, pivot_fraction=0.3,
    )


def main() -> None:
    geo = blade()
    z = 4.0
    cut = geo.section_at(z)
    ring = geo.oml_polyline(z)
    oml = geo.oml_sdf()
    mid = geo.wall_mid_sdf()
    meta = geo.wall_meta()
    solver = GeometrySolver(tolerance=1e-8, max_step=0.02)
    le = ring[int(np.argmin(ring[:, 0]))]
    outside = le + np.array([-0.05, 0.04, 0.0])
    snapped = solver.project(oml, outside.reshape(1, 3), plane=True)
    print(f"z={z:g}: oml n={len(ring)}  walls={len(cut.walls)}")
    print(f"oml_sdf outside phi={oml.evaluate(outside.reshape(1, 3)).phi[0]:.4f}")
    print(f"snap |phi|={abs(oml.evaluate(snapped.points).phi[0]):.2e}  dz={snapped.points[0, 2] - z:.1e}")
    web = next(w for w in cut.walls if w.kind == "web")
    mid_pt = np.append(web.xy[len(web.xy) // 2], z)
    sample = meta.evaluate(mid_pt.reshape(1, 3))
    spec = meta.table.by_name(web.label)
    print(f"wall_meta at {web.label}: phi={sample.phi[0]:.4f} meta={int(sample.meta[0])} (table id {spec.id})")
    mid_phi = mid.evaluate(mid_pt.reshape(1, 3)).phi[0]
    print(f"wall_mid_sdf at web midpoint phi={mid_phi:.4f} (web is a branch, not the skin loop)")

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))
    ax = axes[0]
    draw_shell(ax, cut)
    ax.plot(ring[:, 0], ring[:, 1], color="#0b0b0b", lw=0.7, ls="--", label="oml_polyline")
    ax.plot([outside[0], snapped.points[0, 0]], [outside[1], snapped.points[0, 1]], color="k", lw=0.7)
    ax.plot(outside[0], outside[1], "o", ms=4, color="#eb6834", label="guess")
    ax.plot(snapped.points[0, 0], snapped.points[0, 1], "k+", ms=8, label="snap oml_sdf")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("cut + OML oracle")
    legend_outside(ax)
    ax = axes[1]
    loop = cut.shell_graph().loop_polygon()
    xs = np.linspace(loop[:, 0].min() - 0.05, loop[:, 0].max() + 0.05, 80)
    ys = np.linspace(loop[:, 1].min() - 0.05, loop[:, 1].max() + 0.05, 50)
    X, Y = np.meshgrid(xs, ys)
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, z)])
    phi = mid.evaluate(pts).phi.reshape(X.shape)
    cf = ax.contourf(X, Y, phi, levels=16, cmap="RdBu_r")
    ax.contour(X, Y, phi, levels=[0.0], colors="k", linewidths=0.7)
    ax.plot(loop[:, 0], loop[:, 1], color="k", lw=0.6, ls="--")
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_title("wall_mid_sdf φ at z = 4 m")
    fig.colorbar(cf, ax=ax, shrink=0.8, label="φ / m")
    note(fig, "Left: GeometrySolver snap onto oml_sdf. Right: zero contour of the wall-mid oracle.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "loft_oracles")


if __name__ == "__main__":
    main()
