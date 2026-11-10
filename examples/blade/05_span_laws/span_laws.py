"""Piecewise-linear span tables and ply-drop ramps.

Writes figures/table_law.png, ply_drop.png and loft_scalars.png.

    python examples/blade/05_span_laws/span_laws.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import SERIES, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, PlyDrop, Station  # noqa: E402
from implicit_geometry_kernel.structure import TableLaw  # noqa: E402

LINE, KNOT = SERIES[0], "k"


def _stations():
    foil = Airfoil.naca4("0018")
    return [
        Station(0.0, 2.0, foil, twist=12.0, pitch_axis=0.25, offset=(0.0, 0.0)),
        Station(5.0, 1.2, foil, twist=4.0, pitch_axis=0.35, offset=(0.20, 0.15)),
        Station(10.0, 1.0, foil, twist=0.0, pitch_axis=0.40, offset=(0.10, 0.40)),
    ]


def table_law():
    zk = np.array([0.0, 2.0, 3.0, 6.0])
    vk = np.array([0.030, 0.030, 0.008, 0.008])
    law = TableLaw(zk, vk)
    z = np.linspace(-0.5, 6.5, 401)
    fig, ax = plt.subplots(figsize=(6.8, 3.0))
    ax.plot(z, law(z), color=LINE, label="TableLaw")
    ax.plot(zk, vk, "o", color=KNOT, ms=4.5, label="knots")
    ax.set_xlabel("z / m")
    ax.set_ylabel("thickness / m")
    ax.set_title("TableLaw: plateau, linear ramp, plateau (clamp past the end knots)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, __file__, "table_law")
    return law


def ply_drop():
    drop = PlyDrop("spar_buildup", s0=0.2, s1=0.5, z0=2.0, z1=8.0, thickness=0.010, taper=0.5)
    z = np.linspace(0.0, 10.0, 401)
    extra = [drop.extra(zz, 0.35) for zz in z]
    fig, ax = plt.subplots(figsize=(6.8, 3.0))
    ax.plot(z, extra, color=LINE, label="PlyDrop.extra at s = 0.35")
    ax.axvline(drop.z0, color="0.7", ls=":", lw=0.8)
    ax.axvline(drop.z1, color="0.7", ls=":", lw=0.8)
    ax.set_xlabel("z / m")
    ax.set_ylabel("extra thickness / m")
    ax.set_title("PlyDrop: linear taper, constant window, linear taper")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, __file__, "ply_drop")
    return drop


def loft_scalars(stations):
    geo = BladeGeometry(
        stations, layout="2D-F-LC-0R", skin_thickness=0.004,
        spar_cap_thickness=0.008, shear_web_thickness=0.004, pivot_fraction=0.3,
    )
    z = np.linspace(0.0, 10.0, 201)
    rows = [geo.interpolate(zz) for zz in z]
    fig, axes = plt.subplots(3, 1, figsize=(6.8, 6.4), sharex=True)
    axes[0].plot(z, [r[0] for r in rows], color=LINE)
    axes[0].plot([s.z for s in stations], [s.chord for s in stations], "o", color=KNOT, ms=4.5)
    axes[0].set_ylabel("chord / m")
    axes[1].plot(z, [r[1] for r in rows], color=LINE)
    axes[1].plot([s.z for s in stations], [s.twist for s in stations], "o", color=KNOT, ms=4.5)
    axes[1].set_ylabel("twist / deg")
    axes[2].plot(z, [r[3][0] for r in rows], color=LINE)
    axes[2].plot([s.z for s in stations], [s.offset[0] for s in stations], "o", color=KNOT, ms=4.5)
    axes[2].set_ylabel("offset x / m")
    axes[2].set_xlabel("z / m")
    axes[0].set_title("Station chord, twist and offset: linear between knots")
    fig.tight_layout()
    save(fig, __file__, "loft_scalars")
    return geo


def report(law, drop, geo, z_mid: float = 2.5):
    print(f"TableLaw plateau then ramp: z=1 → {law.at(1.0):.4f} m, z=2.5 → {law.at(2.5):.4f} m")
    print(f"TableLaw clamp: z=-1 → {law.at(-1.0):.4f} m, z=9 → {law.at(9.0):.4f} m")
    print(f"PlyDrop extra at s=0.35: z=1 → {drop.extra(1.0, 0.35):.4f}, "
          f"z=5 → {drop.extra(5.0, 0.35):.4f}, z=8.25 → {drop.extra(8.25, 0.35):.4f}")
    print(f"station-chord at z={z_mid:g}: {geo.interpolate(z_mid)[0]:.4f} m "
          f"(linear 1.6 between 2.0 and 1.2)")


def main():
    stations = _stations()
    law = table_law()
    drop = ply_drop()
    geo = loft_scalars(stations)
    report(law, drop, geo)


if __name__ == "__main__":
    main()
