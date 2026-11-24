"""How the kernel builds an aerofoil: the NACA 4-digit construction and the other families.

Draws Figure R2.1 of the R manual (docs/r/R2_aerofoils.md) and checks results
derived there:

1. the closed trailing-edge coefficient a4 = -0.1036 makes h(1) = 0, (R2.6);
2. NACA 2412: maximum camber m = 0.02 at x = p = 0.4 and maximum thickness 0.12, (R2.7);
3. the surfaces are the camber line offset by h along its normal, (R2.9);
4. the kernel's CST surfaces against (R2.16), the TE gap, and symmetry for equal
   weights;
5. the 5-digit join r solves 3p^2 - 6rp + r^2(3 - r) = 0, (R2.12), and the
   kernel's NACA 23012 against Abbott and von Doenhoff's mean line, (R2.11).

    python examples/blade/17_aerofoil_construction/aerofoil_construction.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import NAVY, REFERENCE, SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil  # noqa: E402
from implicit_geometry_kernel.sdf.shapes.naca import naca_four_digit_vertices  # noqa: E402

A = (0.2969, -0.1260, -0.3516, 0.2843)
X = np.linspace(0.0, 1.0, 4001)


def half_thickness(x, t, a4):
    return (t / 0.2) * (A[0] * np.sqrt(x) + A[1] * x + A[2] * x**2 + A[3] * x**3 + a4 * x**4)


def camber4(x, m, p):
    yc = np.where(x <= p, m / p**2 * (2 * p * x - x**2), m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x**2))
    dyc = np.where(x <= p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    return yc, dyc


def check_closed_te():
    a4 = -sum(A)
    print(f"1. a4 that closes the trailing edge: {a4:+.4f} (kernel -0.1036); "
          f"open-TE thickness 2 h(1) = {2 * half_thickness(1.0, 1.0, -0.1015):.4f} t")
    assert abs(a4 - -0.1036) < 5e-5
    return abs(a4 + 0.1036)


def check_naca2412():
    foil = Airfoil.naca4("2412", n=201)
    yu, yl = foil.surface("upper", X), foil.surface("lower", X)
    cam = 0.5 * (yu + yl)
    print(f"\n2. NACA 2412: max camber {cam.max():.5f} at x = {X[cam.argmax()]:.3f} (0.02 at 0.4); "
          f"Airfoil.max_thickness {foil.max_thickness:.5f}")
    assert abs(cam.max() - 0.02) < 2e-5 and abs(X[cam.argmax()] - 0.4) < 5e-3
    assert abs(foil.max_thickness - 0.12) < 5e-5
    return abs(cam.max() - 0.02)


def check_normal_offset():
    v = naca_four_digit_vertices(0.02, 0.4, 0.12, 2 * 201 - 1, 1.0)
    n = len(v) // 2
    up, lo = v[: n + 1][::-1], v[n:]
    xc = 0.5 * (up[:, 0] + lo[:, 0])
    yc, dyc = camber4(xc, 0.02, 0.4)
    mid = 0.5 * (up + lo)
    half = 0.5 * np.linalg.norm(up - lo, axis=1)
    err_mid = float(np.abs(mid[:, 1] - yc).max())
    err_t = float(np.abs(half[1:-1] - half_thickness(xc[1:-1], 0.12, -0.1036)).max())
    normal = np.column_stack([-dyc, np.ones_like(dyc)]) / np.sqrt(1 + dyc**2)[:, None]
    d = (up - lo)[1:-1]
    err_dir = float(np.abs(d[:, 0] * normal[1:-1, 1] - d[:, 1] * normal[1:-1, 0]).max())
    worst = max(err_mid, err_t, err_dir)
    print(f"\n3. Surfaces = camber line +/- y_t along its normal: largest error {worst:.1e}")
    assert worst < 1e-12
    return worst


def check_cst():
    from math import comb

    au, al, dte = [0.17, 0.15, 0.16], [0.12, 0.05, 0.10], 0.01
    foil = Airfoil.cst(au, al, te_thickness=dte, n=161)

    def surf(w, sign):
        n = len(w) - 1
        shape = sum(w[i] * comb(n, i) * X**i * (1 - X) ** (n - i) for i in range(n + 1))
        return sign * (np.sqrt(X) * (1 - X) * shape + X * 0.5 * dte)

    err = max(float(np.abs(foil.surface("upper", X) - surf(au, 1)).max()),
              float(np.abs(foil.surface("lower", X) - surf(al, -1)).max()))
    sym = Airfoil.cst(au, au)
    sym_err = float(np.abs(sym.surface("upper", X) + sym.surface("lower", X)).max())
    gap = float(foil.surface("upper", 1.0) - foil.surface("lower", 1.0))
    print(f"\n4. CST: kernel surfaces against (R2.16) {err:.1e} of chord; TE gap {gap:.4f} (given {dte}); "
          f"equal weights, |y_u + y_l| <= {sym_err:.1e}")
    assert err < 2e-4 and sym_err < 1e-12 and abs(gap - dte) < 1e-12
    return err


def check_five_digit():
    rows = []
    for p, r_abbott in ((0.05, 0.0580), (0.10, 0.1260), (0.15, 0.2025), (0.20, 0.2900), (0.25, 0.3910)):
        roots = np.roots([-1.0, 3.0, -6.0 * p, 3.0 * p * p])
        r = float(min(z.real for z in roots if abs(z.imag) < 1e-12 and p < z.real < 1.0))
        rows.append((p, r, r_abbott))
    worst = max(abs(r - ra) for _, r, ra in rows)
    print("\n5. 5-digit join r from (R2.12) against Abbott and von Doenhoff:")
    for p, r, ra in rows:
        print(f"   p = {p:.2f}: r = {r:.4f} (tabulated {ra:.4f})")
    assert worst < 5e-4
    r, k1 = 0.2025, 15.957
    yc = np.where(X < r, k1 / 6 * (X**3 - 3 * r * X**2 + r * r * (3 - r) * X), k1 * r**3 / 6 * (1 - X))
    foil = Airfoil.naca5("23012", n=401)
    kc = 0.5 * (foil.surface("upper", X) + foil.surface("lower", X))
    print(f"   NACA 23012 max camber: Abbott {yc.max():.4f} at x = {X[yc.argmax()]:.3f}; "
          f"kernel section mid-line {kc.max():.4f} at x = {X[kc.argmax()]:.3f}")
    assert abs(kc.max() - yc.max()) < 2e-4
    return X, yc, kc


def figure(five):
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.8), gridspec_kw={"height_ratios": [1.25, 1.0], "hspace": 0.5})
    ax = axes[0]
    m, p, t = 0.04, 0.4, 0.18
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 41)))
    yc, dyc = camber4(x, m, p)
    th = half_thickness(x, t, -0.1036)
    nrm = np.column_stack([-dyc, np.ones_like(dyc)]) / np.sqrt(1 + dyc**2)[:, None]
    up = np.column_stack([x, yc]) + th[:, None] * nrm
    lo = np.column_stack([x, yc]) - th[:, None] * nrm
    xf = np.linspace(0, 1, 400)
    ycf, _ = camber4(xf, m, p)
    ax.plot(xf, ycf, color=SERIES[1], lw=1.2, label="camber line $y_c$ (R2.7)")
    ax.plot(xf, half_thickness(xf, t, -0.1036), color=REFERENCE, lw=0.9, ls="--", label="half-thickness $h$ about the chord (R2.5)")
    for i in range(2, len(x) - 1, 3):
        ax.plot([lo[i, 0], up[i, 0]], [lo[i, 1], up[i, 1]], color="0.6", lw=0.6)
    ax.plot(up[:, 0], up[:, 1], color=NAVY, lw=1.5, label="surfaces: $h$ along the camber normal (R2.9)")
    ax.plot(lo[:, 0], lo[:, 1], color=NAVY, lw=1.5)
    ax.plot(p, m, "o", color=SERIES[1], ms=4)
    ax.annotate("$(p, m)$", (p, m), xytext=(6, 6), textcoords="offset points", fontsize=8, color=SERIES[1])
    ax.set_aspect("equal")
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$y/c$")
    ax.set_title("NACA 4418: camber line, thickness, and the normal offset that builds the surfaces")
    ax.legend(loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.62))

    ax = axes[1]
    fams = [(Airfoil.naca4("2412"), "NACA 2412 (R2.3)"), (Airfoil.naca6(63, 15), "NACA 63-415 (R2.13)"),
            (Airfoil.cst([0.2, 0.25, 0.2], [0.12, 0.05, 0.1], te_thickness=0.01), "CST, 3 weights a side (R2.16)"),
            (Airfoil.naca4("2418").with_te_thickness(0.05), "NACA 2418 flatback, TE 0.05 (R2.21)"),
            (Airfoil.ellipse(0.4), "ellipse t/c = 0.4 (R2.18)")]
    for i, (foil, label) in enumerate(fams):
        pts = foil.contour_points(301)
        ax.plot(pts[:, 0], pts[:, 1], color=SERIES[i % len(SERIES)], lw=1.2, label=label)
    xs, yc5, kc5 = five
    ax.plot(xs, yc5, color=REFERENCE, lw=0.9, ls="--", label="23012 mean line, Abbott (R2.11)")
    ax.plot(xs, kc5, color=REFERENCE, lw=0.9, ls=":", label="23012 mid-line, kernel")
    ax.set_aspect("equal")
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$y/c$")
    ax.set_title("Families the kernel builds, chord-normalised")
    ax.legend(loc="upper left", ncol=2, frameon=False, bbox_to_anchor=(0.0, -0.22))
    note(fig, "Checks: a4 closes the TE; NACA 2412 camber and thickness; normal offset to 1e-12; "
              "CST surfaces and symmetry; 5-digit join r and NACA 23012 camber.")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.2)
    save(fig, __file__, "aerofoil_construction")


if __name__ == "__main__":
    check_closed_te()
    check_naca2412()
    check_normal_offset()
    check_cst()
    five = check_five_digit()
    figure(five)
