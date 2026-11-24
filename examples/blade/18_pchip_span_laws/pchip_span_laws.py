"""Span laws along the blade: linear tables, the PCHIP cubic, its sampling check and the airfoil morph.

Draws Figure R3.1 of the R manual (docs/r/R3_span_laws.md) and checks results
derived there:

1. the kernel's PCHIP equals SciPy's ``PchipInterpolator`` (value and slope), (R3.4)-(R3.9);
2. collinear knots give the straight line, (R3.10);
3. monotone knots give a monotone law with no overshoot, and the interior slope
   ratios lie in (0, 3), (R3.8) and (R3.11);
4. the law is C1 at the knots and its second derivative may jump, (R3.6);
5. the leave-one-out error at a smooth maximum is kappa h^2 / 2, so the check
   passes only for h <= sqrt(2 tol / kappa), (R3.12)-(R3.13);
6. the thickness morph: the blend weight equals the span weight and the rescale
   factor is at least 1; the trailing-edge gap follows each rule's own weight,
   (R3.15)-(R3.18).

    python examples/blade/18_pchip_span_laws/pchip_span_laws.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import NAVY, REFERENCE, SERIES, note, plt, save  # noqa: E402
from scipy.interpolate import CubicSpline, PchipInterpolator  # noqa: E402

from implicit_geometry_kernel import Airfoil  # noqa: E402
from implicit_geometry_kernel.blade.loft.loft_tables import _morph_surfaces, _smoothstep  # noqa: E402
from implicit_geometry_kernel.structure import PchipLaw, TableLaw, UnderSampledSpanLaw, check_pchip_sampling  # noqa: E402
from implicit_geometry_kernel.structure.laws import pchip_slopes  # noqa: E402

RNG = np.random.default_rng(20261124)
# A chord law: a root cylinder, a ramp to the maximum chord, a plateau, then taper.
Z = np.array([0.0, 1.0, 3.0, 4.5, 6.0, 8.0, 11.0, 14.0, 17.0, 20.0])
C = np.array([1.2, 1.2, 1.9, 2.1, 2.1, 1.9, 1.6, 1.3, 1.0, 0.6])


def random_knots(n, monotone=False):
    z = np.sort(RNG.uniform(0.0, 10.0, n))
    while np.min(np.diff(z)) < 0.05:
        z = np.sort(RNG.uniform(0.0, 10.0, n))
    v = np.cumsum(RNG.uniform(0.0, 1.0, n)) if monotone else RNG.normal(size=n)
    if monotone:
        k = RNG.integers(1, n - 1)
        v[k] = v[k - 1]  # a plateau
    return z, v


def check_scipy():
    worst = 0.0
    for _ in range(300):
        z, v = random_knots(int(RNG.integers(3, 12)))
        law = PchipLaw(z, v, sampling="flag")
        zs = np.linspace(z[0], z[-1], 400)
        ref = PchipInterpolator(z, v)
        worst = max(worst, float(np.abs(law(zs) - ref(zs)).max()),
                    float(np.abs(law.derivative(zs) - ref.derivative()(zs)).max()))
    print(f"1. PchipLaw against scipy PchipInterpolator, 300 random knot sets: largest difference {worst:.1e}")
    assert worst < 1e-10
    return worst


def check_line():
    z = np.array([0.0, 0.7, 2.0, 2.5, 6.0])
    v = 0.3 - 0.12 * z
    law = PchipLaw(z, v)
    zs = np.linspace(0.0, 6.0, 601)
    err = float(np.abs(law(zs) - (0.3 - 0.12 * zs)).max())
    print(f"\n2. Collinear knots, unequal spacing: largest distance from the line {err:.1e}")
    assert err < 1e-15
    return err


def check_monotone():
    worst_ratio, worst_overshoot, worst_dip = 0.0, 0.0, 0.0
    for _ in range(300):
        z, v = random_knots(int(RNG.integers(4, 12)), monotone=True)
        d = pchip_slopes(z, v)
        delta = np.diff(v) / np.diff(z)
        live = delta > 0.0
        alpha = d[:-1][live] / delta[live]
        beta = d[1:][live] / delta[live]
        worst_ratio = max(worst_ratio, float(np.max(np.concatenate([alpha, beta]))))
        assert np.all(alpha >= 0.0) and np.all(beta >= 0.0)
        zs = np.linspace(z[0], z[-1], 2000)
        y = PchipLaw(z, v, sampling="flag")(zs)
        worst_dip = max(worst_dip, float(-np.min(np.diff(y))))
        k = np.clip(np.searchsorted(z, zs, side="right") - 1, 0, len(z) - 2)
        lo, hi = np.minimum(v[k], v[k + 1]), np.maximum(v[k], v[k + 1])
        worst_overshoot = max(worst_overshoot, float(np.max(np.maximum(y - hi, lo - y))))
    print(f"\n3. Monotone knots (300 sets, each with a plateau): largest slope ratio {worst_ratio:.3f} (< 3); "
          f"largest decrease {worst_dip:.1e}; largest overshoot of a bay {worst_overshoot:.1e}")
    assert worst_ratio <= 3.0 and worst_dip <= 1e-12 and worst_overshoot <= 1e-12
    return worst_ratio


def check_c1():
    law = PchipLaw(Z, C, sampling="flag")
    eps = 1e-9
    zi = Z[1:-1]
    slope_jump = float(np.abs(law.derivative(zi + eps) - law.derivative(zi - eps)).max())
    curv_jump = float(np.abs(law.derivative(zi + eps, 2) - law.derivative(zi - eps, 2)).max())
    print(f"\n4. Chord law at its interior knots: largest slope jump {slope_jump:.1e} 1/m; "
          f"largest curvature jump {curv_jump:.3f} 1/m^2 (C1, not C2)")
    assert slope_jump < 1e-7
    return slope_jump, curv_jump


def check_sampling():
    """A smooth chord law with an interior maximum: the leave-one-out error at the peak is kappa h^2 / 2."""
    amp, width, z_peak = 1.5, 6.0, 5.0
    kappa = 2.0 * amp / width**2  # -c''(z_peak)

    def chord(z):
        return 0.6 + amp * np.exp(-(((z - z_peak) / width) ** 2))

    rows = []
    for n in (21, 41, 61, 81):
        z = np.linspace(0.0, 20.0, n)
        rep = check_pchip_sampling(z, chord(z), name="chord")
        h = z[1] - z[0]
        rows.append((n, h, rep.loo_error[np.argmin(abs(z - z_peak))], 0.5 * kappa * h * h, rep.tolerance, rep.ok))
    print("\n5. Leave-one-out check on a smooth chord law with a maximum at z = 5 m:")
    worst = 0.0
    for n, h, got, pred, tol, ok in rows:
        worst = max(worst, abs(got - pred) / pred)
        print(f"   {n:3d} stations (h = {h:.3f} m): error at the peak {got:.5f} m, kappa h^2 / 2 = {pred:.5f} m, "
              f"tolerance {tol:.5f} m, ok = {ok}")
    h_max = np.sqrt(2.0 * rows[0][4] / kappa)
    print(f"   (R3.13): the peak passes for h <= sqrt(2 tol / kappa) = {h_max:.3f} m")
    assert worst < 0.02 and not rows[1][5] and rows[2][5]
    try:
        z = np.linspace(0.0, 20.0, 41)
        PchipLaw(z, chord(z), name="chord")
        raise AssertionError("41 stations should be refused")
    except UnderSampledSpanLaw as exc:
        print(f"   PchipLaw with 41 stations raises: {exc}")
    ramp = check_pchip_sampling(Z, C, name="chord")
    print(f"   the ten-station chord law of the figure: ok = {ramp.ok} (worst {ramp.loo_error.max():.3f} m), "
          "so span_interp='auto' builds it linear")
    return ramp


def check_morph():
    fa, fb = Airfoil.naca4("0030"), Airfoil.naca4("4412")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 161)))
    dense = np.linspace(0.0, 1.0, 401)
    ta = float(np.max(fa.surface("upper", dense) - fa.surface("lower", dense)))
    tb = float(np.max(fb.surface("upper", dense) - fb.surface("lower", dense)))
    worst_t, worst_scale = 0.0, np.inf
    for w in np.linspace(0.05, 0.95, 19):
        t_c = (1 - w) * ta + w * tb
        s = (t_c - ta) / (tb - ta)
        assert abs(s - w) < 1e-15
        yu = (1 - w) * fa.surface("upper", dense) + w * fb.surface("upper", dense)
        yl = (1 - w) * fa.surface("lower", dense) + w * fb.surface("lower", dense)
        worst_scale = min(worst_scale, t_c / float(np.max(yu - yl)))
        u, lo = _morph_surfaces(fa, fb, w, dense, "thickness")
        worst_t = max(worst_t, abs(float(np.max(u - lo)) - t_c))
    print(f"\n6. Thickness morph NACA 0030 -> 4412: s = w exactly; smallest rescale factor {worst_scale:.6f} (>= 1); "
          f"peak t/c against the lerped t/c {worst_t:.1e}")
    assert worst_scale >= 1.0 and worst_t < 1e-12
    flat = Airfoil.naca4("2418").with_te_thickness(0.04)
    gaps = []
    for rule in ("linear", "smoothstep", "thickness"):
        for w in (0.25, 0.5, 0.75):
            u, lo = _morph_surfaces(fb, flat, w, x, rule)
            weight = _smoothstep(w) if rule == "smoothstep" else w
            gaps.append(abs((u[-1] - lo[-1]) - weight * 0.04))
    print(f"   trailing-edge gap against s * 0.04, s the rule's own weight, all three rules: {max(gaps):.1e}")
    assert max(gaps) < 1e-12
    assert _smoothstep(0.0) == 0.0 and _smoothstep(1.0) == 1.0 and _smoothstep(0.5) == 0.5
    return worst_scale


def figure(ramp, jumps):
    fig = plt.figure(figsize=(9.0, 7.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0])
    ax = fig.add_subplot(grid[0, :])
    zs = np.linspace(Z[0], Z[-1], 800)
    pch = PchipLaw(Z, C, sampling="flag")
    ax.plot(zs, CubicSpline(Z, C, bc_type="natural")(zs), color=SERIES[3], lw=1.0, ls="--",
            label="natural cubic spline (overshoots)")
    ax.plot(zs, TableLaw(Z, C)(zs), color=REFERENCE, lw=0.9, ls=":", label="linear TableLaw (R3.3)")
    ax.plot(zs, pch(zs), color=NAVY, lw=1.6, label="PchipLaw (R3.4), (R3.9)")
    ax.plot(Z, C, "o", color=NAVY, ms=4, label="stations")
    ax.set_xlabel("span $z$ (m)")
    ax.set_ylabel("chord $c$ (m)")
    ax.set_title("A chord law through ten stations: plateaus stay flat and no bay overshoots its knots")
    ax.legend(loc="upper right", frameon=False)

    ax = fig.add_subplot(grid[1, 0])
    ax.plot(zs, pch.derivative(zs), color=NAVY, lw=1.4, label="$c'(z)$, PCHIP (continuous)")
    lin = np.gradient(TableLaw(Z, C)(zs), zs)
    ax.plot(zs, lin, color=REFERENCE, lw=0.9, ls=":", label="$c'(z)$, linear (jumps)")
    ax.set_xlabel("span $z$ (m)")
    ax.set_ylabel("$\\mathrm{d}c/\\mathrm{d}z$")
    ax.set_title("Slope: C1 at every station (R3.6)")
    ax.legend(loc="lower left", frameon=False)

    ax = fig.add_subplot(grid[1, 1])
    ellipse, foil = Airfoil.ellipse(1.0), Airfoil.naca4("2418").with_te_thickness(0.03)
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 161)))
    for i, w in enumerate((0.0, 0.25, 0.5, 0.75, 1.0)):
        u, lo = _morph_surfaces(ellipse, foil, w, x, "auto")
        colour = SERIES[i % len(SERIES)]
        ax.plot(x, u, color=colour, lw=1.1, label=f"$w$ = {w:g}")
        ax.plot(x, lo, color=colour, lw=1.1)
    ax.set_aspect("equal")
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$y/c$")
    ax.set_title("Circle to flatback: smoothstep morph (R3.14)")
    ax.legend(loc="center left", frameon=False, fontsize=7, bbox_to_anchor=(1.0, 0.5))
    note(fig, f"These stations fail the leave-one-out check ({ramp.loo_error.max():.2f} m against "
              f"{ramp.tolerance:.4f} m), so span_interp='auto' lofts them linearly; PCHIP drawn with "
              f"sampling='flag'. Curvature jumps up to {jumps[1]:.2f} 1/m².")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, __file__, "pchip_span_laws")


if __name__ == "__main__":
    check_scipy()
    check_line()
    check_monotone()
    jumps = check_c1()
    ramp = check_sampling()
    check_morph()
    figure(ramp, jumps)
