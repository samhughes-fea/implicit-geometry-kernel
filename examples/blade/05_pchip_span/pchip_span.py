"""Compare linear and PCHIP interpolation in every place this repo can choose.

Writes figures/table_law.png, loft_scalars.png and oml_delta.png.

    python examples/blade/05_pchip_span/pchip_span.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.plotting import SERIES, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.structure import TableLaw  # noqa: E402

LINEAR, PCHIP, KNOT = SERIES[0], SERIES[1], "k"


def _stations():
    foil = Airfoil.naca4("0018")
    return [
        Station(0.0, 2.0, foil, twist=12.0, pitch_axis=0.25, offset=(0.0, 0.0)),
        Station(5.0, 1.2, foil, twist=4.0, pitch_axis=0.35, offset=(0.20, 0.15)),
        Station(10.0, 1.0, foil, twist=0.0, pitch_axis=0.40, offset=(0.10, 0.40)),
    ]


def _blades(stations):
    kw = dict(layout="2D-F-LC-0R", skin_thickness=0.004, spar_cap_thickness=0.008,
              shear_web_thickness=0.004, pivot_fraction=0.3)
    linear = BladeGeometry(stations, span_kind="linear", **kw)
    pchip = BladeGeometry(stations, span_kind="pchip", **kw)
    return linear, pchip


def table_law():
    zk = np.array([0.0, 1.0, 2.0, 4.0, 6.0])
    vk = np.array([0.030, 0.030, 0.020, 0.008, 0.008])
    linear = TableLaw(zk, vk, kind="linear")
    pchip = TableLaw(zk, vk, kind="pchip")
    z = np.linspace(-0.5, 6.5, 401)
    fig, ax = plt.subplots(figsize=(6.8, 3.0))
    ax.plot(z, linear(z), color=LINEAR, label="linear")
    ax.plot(z, pchip(z), color=PCHIP, label="PCHIP")
    ax.plot(zk, vk, "o", color=KNOT, ms=4.5, label="knots")
    ax.axhline(vk.min(), color="0.7", ls=":", lw=0.8)
    ax.set_xlabel("z / m")
    ax.set_ylabel("thickness / m")
    ax.set_title("TableLaw: default PCHIP vs kind='linear' (clamp past the end knots)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, __file__, "table_law")
    return linear, pchip, zk, vk


def loft_scalars(linear, pchip, stations):
    z = np.linspace(0.0, 10.0, 201)
    fig, axes = plt.subplots(3, 1, figsize=(6.8, 6.4), sharex=True)

    def series(geo):
        rows = [geo.interpolate(zz) for zz in z]
        chord = [r[0] for r in rows]
        twist = [r[1] for r in rows]
        pitch = [r[2] for r in rows]
        off_x = [r[3][0] for r in rows]
        return chord, twist, pitch, off_x

    lc, lt, lp, lox = series(linear)
    pc, pt, pp, pox = series(pchip)
    axes[0].plot(z, lc, color=LINEAR, label="linear")
    axes[0].plot(z, pc, color=PCHIP, label="PCHIP")
    axes[0].plot([s.z for s in stations], [s.chord for s in stations], "o", color=KNOT, ms=4.5)
    axes[0].set_ylabel("chord / m")
    axes[1].plot(z, lt, color=LINEAR)
    axes[1].plot(z, pt, color=PCHIP)
    axes[1].plot([s.z for s in stations], [s.twist for s in stations], "o", color=KNOT, ms=4.5)
    axes[1].set_ylabel("twist / deg")
    axes[2].plot(z, lox, color=LINEAR)
    axes[2].plot(z, pox, color=PCHIP)
    axes[2].plot([s.z for s in stations], [s.offset[0] for s in stations], "o", color=KNOT, ms=4.5)
    axes[2].set_ylabel("offset x / m")
    axes[2].set_xlabel("z / m")
    axes[0].legend(loc="upper right")
    axes[0].set_title("Station planform via span_kind (pitch is not in this set)")
    fig.tight_layout()
    save(fig, __file__, "loft_scalars")
    return lp, pp


def oml_delta(linear, pchip, z_cut: float = 2.5):
    fig, ax = plt.subplots(figsize=(6.8, 3.2))
    for geo, colour, label in ((linear, LINEAR, "linear"), (pchip, PCHIP, "PCHIP")):
        xy = geo.oml_polyline(z_cut)[:, :2]
        ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=1.2, label=label)
    ax.set_aspect("equal")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title(f"OML at z = {z_cut:g} m (mid-bay): chord/twist/offset differ, airfoil blend does not")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, __file__, "oml_delta")


def report(linear, pchip, stations, lp, pp, z_mid: float = 2.5):
    zs = [s.z for s in stations]
    override_z = (0.0, 5.0, 10.0)
    override_c = (2.0, 1.6, 1.0)
    kw = dict(layout="2D-F-LC-0R", skin_thickness=0.004, spar_cap_thickness=0.008,
              shear_web_thickness=0.004, pivot_fraction=0.3, chord=(override_z, override_c))
    over_lin = BladeGeometry(stations, span_kind="linear", **kw)
    over_pchip = BladeGeometry(stations, span_kind="pchip", **kw)
    c_lin, c_pchip = over_lin.interpolate(z_mid)[0], over_pchip.interpolate(z_mid)[0]
    expected = float(np.interp(z_mid, override_z, override_c))
    print(f"pitch at sample zs identical: {np.allclose(lp, pp)}")
    print(f"station pitch_axis unused by interpolate (pivot_fraction=0.3): "
          f"{linear.interpolate(z_mid)[2]:.3f} vs station values {[s.pitch_axis for s in stations]}")
    print(f"override SpanTable chord at z={z_mid:g}: linear={c_lin:.6f}, pchip={c_pchip:.6f}, "
          f"np.interp={expected:.6f} (span_kind ignored)")
    print(f"airfoil morph is shared: n_airfoil_points={linear.n_airfoil_points}, "
          f"same x-grid at z={z_mid:g}: "
          f"{np.allclose(linear.interpolate(z_mid)[4][0], pchip.interpolate(z_mid)[4][0])}")
    cl, cp = linear.interpolate(z_mid)[0], pchip.interpolate(z_mid)[0]
    print(f"station-chord at z={z_mid:g}: linear={cl:.4f} m, PCHIP={cp:.4f} m "
          f"(knots at {zs})")


def main():
    stations = _stations()
    linear, pchip = _blades(stations)
    table_law()
    lp, pp = loft_scalars(linear, pchip, stations)
    oml_delta(linear, pchip)
    report(linear, pchip, stations, lp, pp)


if __name__ == "__main__":
    main()
