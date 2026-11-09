"""LoftedPolygonSDF, TagNode union, FeatureTable, GeometrySolver.

    python examples/sdf/02_feature_field_loft/feature_field_loft.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from matplotlib.patches import Patch  # noqa: E402

from _shared.plotting import SERIES, legend_outside, meta_cmap, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, FeatureTable, Station  # noqa: E402
from implicit_geometry_kernel.sdf import GeometrySolver, LoftedPolygonSDF, PriorityUnionNode, TagNode  # noqa: E402
from implicit_geometry_kernel.sdf.nodes.tagged import CircleSDF  # noqa: E402


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
    table = FeatureTable()
    oml_spec = table.add("oml", "skin")
    probe_spec = table.add("probe", "web")
    lofted = LoftedPolygonSDF(lambda z: geo.oml_polyline(float(z))[:, :2], meta=oml_spec.id)
    tagged = TagNode(lofted, oml_spec.id)
    probe = TagNode(CircleSDF(0.0, 0.0, 0.05), probe_spec.id)
    union = PriorityUnionNode([tagged, probe], [oml_spec.priority, probe_spec.priority])
    z = 4.0
    ring = geo.oml_polyline(z)
    le = ring[int(np.argmin(ring[:, 0]))]
    guess = le + np.array([-0.06, 0.05, 0.0])
    hit = GeometrySolver(tolerance=1e-8, max_step=0.02).project(union, guess.reshape(1, 3), plane=True)
    sample = union.evaluate(hit.points)
    print(f"LoftedPolygonSDF snap |phi|={abs(sample.phi[0]):.2e} meta={int(sample.meta[0])} ({table[int(sample.meta[0])].name})")
    print(f"table: {[(s.id, s.name, s.kind) for s in table]}")
    assert type(geo.oml_sdf()).__name__ == "LoftedPolygonSDF"

    xs = np.linspace(float(ring[:, 0].min()) - 0.15, float(ring[:, 0].max()) + 0.08, 160)
    ys = np.linspace(float(ring[:, 1].min()) - 0.08, float(ring[:, 1].max()) + 0.08, 90)
    X, Y = np.meshgrid(xs, ys)
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, z)])
    phi = lofted.evaluate(pts).phi.reshape(X.shape)
    meta = union.evaluate(pts).meta.reshape(X.shape)
    cmap, norm = meta_cmap(len(table), [SERIES[0], SERIES[1]])

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))
    ax = axes[0]
    cf = ax.contourf(X, Y, phi, levels=16, cmap="RdBu_r")
    ax.contour(X, Y, phi, levels=[0.0], colors="k", linewidths=0.8)
    ax.plot([guess[0], hit.points[0, 0]], [guess[1], hit.points[0, 1]], color="k", lw=0.8)
    ax.plot(guess[0], guess[1], "o", color=SERIES[1], ms=5, label="guess")
    ax.plot(hit.points[0, 0], hit.points[0, 1], "+", color="k", ms=10, label="project")
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title("LoftedPolygonSDF φ at z = 4 m")
    ax.legend(loc="lower right", fontsize=7, frameon=True)
    fig.colorbar(cf, ax=ax, shrink=0.8, label="φ / m")
    ax = axes[1]
    ax.pcolormesh(X, Y, meta, cmap=cmap, norm=norm, shading="nearest")
    ax.plot(ring[:, 0], ring[:, 1], color="k", lw=0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_title("PriorityUnionNode meta (OML ∪ probe)")
    legend_outside(ax, handles=[
        Patch(facecolor=cmap(0), label="0:oml"),
        Patch(facecolor=cmap(1), label="1:probe"),
    ])
    note(fig, "Guess is off the leading edge. Probe is a 50 mm CircleSDF at the pitch axis, tagged as web.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "feature_field_loft")


if __name__ == "__main__":
    main()
