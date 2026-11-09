"""Extract phi=0 loops from an OML oracle and a spanwise TE feature curve.

Station cuts already carry wall midlines. This package traces signed-distance
fields: marching squares plus Newton (``GeometrySolver``), and
``feature_curve`` for a named shell vertex along the span.

    python examples/extract/01_zero_loops_and_features/zero_loops_and_features.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from _shared.plotting import SERIES, legend_outside, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.extract import extract_zero_loops, feature_curve  # noqa: E402
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


def bounds_at(geo: BladeGeometry, z: float, pad: float = 0.08) -> tuple[float, float, float, float]:
    ring = geo.oml_polyline(z)
    return (
        float(ring[:, 0].min()) - pad, float(ring[:, 0].max()) + pad,
        float(ring[:, 1].min()) - pad, float(ring[:, 1].max()) + pad,
    )


def main() -> None:
    geo = blade()
    z = 4.0
    oml = geo.oml_sdf()
    solver = GeometrySolver(tolerance=1e-8, max_step=0.05)
    loops = extract_zero_loops(oml, z, bounds_at(geo, z), n=48, solver=solver)
    cut = geo.section_at(z)
    zs = np.linspace(0.0, 10.0, 6)
    te = feature_curve(geo, ("TE", "suction"), zs)
    print(f"z={z:g}: {len(loops)} zero loop(s); cut has {len(cut.walls)} walls (station product)")
    print(f"TE suction curve: {len(te)} stations, z in [{te[0, 2]:.1f}, {te[-1, 2]:.1f}]")
    assert loops, "OML oracle should close at midspan"
    r = np.linalg.norm(loops[0][:, :2] - geo.oml_polyline(z)[:, :2].mean(axis=0), axis=1)
    print(f"extracted loop n={len(loops[0])}  radial span {r.min():.3f}–{r.max():.3f} m")

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))
    ax = axes[0]
    ring = geo.oml_polyline(z)
    ax.plot(ring[:, 0], ring[:, 1], color="#888", lw=0.8, ls="--", label="oml_polyline")
    ax.plot(loops[0][:, 0], loops[0][:, 1], color=SERIES[0], lw=1.2, label="extract_zero_loops")
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title(f"OML zero set at z = {z:g} m")
    legend_outside(ax)
    ax = axes[1]
    for zi in (0.0, 5.0, 10.0):
        r = geo.oml_polyline(zi)
        ax.plot([zi, zi], [float(r[:, 0].min()), float(r[:, 0].max())],
                color="#ccc", lw=6, solid_capstyle="butt", label="OML x-extent" if zi == 0.0 else None)
    ax.plot(te[:, 2], te[:, 0], "o-", color=SERIES[1], ms=4, lw=1.4, label="TE suction")
    ax.set_xlabel("z / m")
    ax.set_ylabel("x / m")
    ax.set_title("TE suction locus on the planform")
    legend_outside(ax)
    note(fig, "Left: marching-squares loop vs authored OML. Right: TE vertex x(z) on station OML extents.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, __file__, "zero_loops_and_features")


if __name__ == "__main__":
    main()
