"""The section frames of R1 on one twisted station: B, S, S' and C, and the twist angle.

Draws Figure R1.1 of the R manual (docs/r/R1_frames.md) and checks four results
derived there:

1. the swap identity P Q(theta) P = Q(-theta), (R1.10);
2. that ``BladeGeometry`` cuts rotate by Q(tau) in their own (x, y) order with
   tau = -theta, (R1.12): the trailing edge of a twisted cut against the rotated
   trailing edge of the untwisted cut;
3. the parallel-axis relations between S and S', (R1.24);
4. the invariants of the second-moment tensor under S -> B, (R1.22).

The figure uses the section-plot convention of ``implicit_geometry_kernel.plotting``:
horizontal Y_B, vertical X_B, horizontal axis inverted so that the leading edge is
on the right. The data are not negated.

    python examples/frames/02_frame_triad/frame_triad.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import NAVY, REFERENCE, SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.frames import naca_to_S, rotate_S_to_B  # noqa: E402
from implicit_geometry_kernel.frames.tensors.inertia import transform_second_moments_S_to_B  # noqa: E402
from implicit_geometry_kernel.geos.geom import Polygon  # noqa: E402
from implicit_geometry_kernel.integration import section_properties_S  # noqa: E402
from implicit_geometry_kernel.plotting.frames import section_plot_axes, setup_section_axes  # noqa: E402

CHORD, PITCH, THETA_DEG = 2.0, 0.3, 20.0


def Q(theta):
    """S -> B rotation on (X, Y) components, (R1.4)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s], [-s, c]])


P = np.array([[0.0, 1.0], [1.0, 0.0]])


def outline_naca():
    a = Airfoil.naca4("2418")
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 201)))
    return np.vstack([np.column_stack([x[::-1], a.surface("upper", x)[::-1]]),
                      np.column_stack([x[1:], a.surface("lower", x)[1:]])])[:-1] * CHORD


def check_swap():
    worst = max(float(np.abs(P @ Q(t) @ P - Q(-t)).max()) for t in np.radians(np.linspace(-90, 90, 19)))
    print(f"1. P Q(theta) P = Q(-theta): largest entry error {worst:.1e}")
    assert worst < 1e-15
    return worst


def trailing_edge(tau_deg):
    airfoil = Airfoil.naca4("2418")
    st = [Station(0.0, CHORD, airfoil, twist=tau_deg, pitch_axis=PITCH),
          Station(1.0, CHORD, airfoil, twist=tau_deg, pitch_axis=PITCH)]
    cut = BladeGeometry(st, "1A-CN-LC-0R", 0.006, shear_web_thickness=0.01, web_positions=(0.4,),
                        pivot_fraction=PITCH).section_at(0.5)
    pts = np.vstack([w.xy[:, :2] for w in cut.walls_of("skin")])
    return pts[np.argmax(np.hypot(pts[:, 0], pts[:, 1]))]


def check_cut_rotation():
    te0 = trailing_edge(0.0)
    print("\n2. BladeGeometry cut frame: trailing edge at twist tau against Q(tau) applied at zero twist")
    worst = 0.0
    for tau in (-20.0, -10.0, 10.0, 20.0):
        got, expect = trailing_edge(tau), Q(np.radians(tau)) @ te0
        err = float(np.linalg.norm(got - expect))
        worst = max(worst, err)
        print(f"   tau {tau:+5.1f} deg: cut ({got[0]:+.6f}, {got[1]:+.6f})  Q(tau) ({expect[0]:+.6f}, {expect[1]:+.6f})")
    print(f"   largest distance {worst:.1e} m ({worst / CHORD:.1e} of chord)")
    assert worst < 1e-6 * CHORD
    return worst


def check_parallel_axis(S):
    area, (xb, yb), ixx, iyy, ixy = section_properties_S(Polygon(S))
    _, _, ixx_c, iyy_c, ixy_c = section_properties_S(Polygon(S - [xb, yb]))
    errs = (ixx - (ixx_c + area * yb**2), iyy - (iyy_c + area * xb**2), ixy - (ixy_c + area * xb * yb))
    worst = max(abs(e) for e in errs) / (ixx + iyy)
    print(f"\n3. Parallel axis S -> S': largest residual {worst:.1e} of Ixx + Iyy (A = {area:.5f} m^2)")
    assert worst < 1e-12
    return (xb, yb), worst


def check_invariants(S):
    _, _, ixx, iyy, ixy = section_properties_S(Polygon(S))
    worst = 0.0
    for t in np.linspace(-90.0, 90.0, 13):
        bxx, byy, bxy = transform_second_moments_S_to_B(ixx, iyy, ixy, twist_deg=t)
        worst = max(worst, abs((bxx + byy) - (ixx + iyy)) / (ixx + iyy),
                    abs((bxx * byy - bxy**2) - (ixx * iyy - ixy**2)) / (ixx * iyy - ixy**2))
    print(f"\n4. Invariants Ixx + Iyy and Ixx Iyy - Ixy^2 under S -> B: largest relative change {worst:.1e}")
    assert worst < 1e-12
    return worst


def arrow(ax, origin_B, vec_B, colour, label, *, ls="-", scale=0.55, offset=1.12):
    """Draw a B-frame vector; plot coordinates are (Y_B, X_B)."""
    o = np.asarray(origin_B, float)
    v = np.asarray(vec_B, float) * scale
    ax.annotate("", xy=(o[1] + v[1], o[0] + v[0]), xytext=(o[1], o[0]),
                arrowprops=dict(arrowstyle="-|>", color=colour, lw=1.6, ls=ls, shrinkA=0, shrinkB=0))
    ax.text(o[1] + offset * v[1], o[0] + offset * v[0], label, color=colour, fontsize=10,
            ha="center", va="center")


def figure(naca, centroid_S, errs):
    theta = np.radians(THETA_DEG)
    S = naca_to_S(naca, CHORD, PITCH)
    B = rotate_S_to_B(S, theta)
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    h, v = section_plot_axes(B)
    ax.fill(h, v, color=SERIES[0], alpha=0.15, ec=NAVY, lw=1.0)

    # Basis vectors of S in B components, (R1.3): e_XS = (c, -s), e_YS = (s, c) as (X_B, Y_B).
    eXS, eYS = Q(theta) @ [1.0, 0.0], Q(theta) @ [0.0, 1.0]
    o = np.zeros(2)
    arrow(ax, o, [1.0, 0.0], REFERENCE, r"$\mathbf{e}_{X_B}$")
    arrow(ax, o, [0.0, 1.0], REFERENCE, r"$\mathbf{e}_{Y_B}$")
    arrow(ax, o, eXS, SERIES[0], r"$\mathbf{e}_{X_S}$")
    arrow(ax, o, eYS, SERIES[0], r"$\mathbf{e}_{Y_S}$")

    # Twist angle between e_YB and e_YS: in plot coordinates (Y_B, X_B), e_YS = (cos, sin) of theta.
    phi = np.linspace(0.0, theta, 40)
    ax.plot(0.36 * np.cos(phi), 0.36 * np.sin(phi), color=SERIES[1], lw=1.2)
    ax.text(0.44 * np.cos(theta / 2), 0.44 * np.sin(theta / 2), r"$\theta$", color=SERIES[1], fontsize=11,
            ha="center", va="center")

    # S' at the centroid: same orientation as S, dashed.
    cB = Q(theta) @ np.asarray(centroid_S)
    arrow(ax, cB, eXS, SERIES[2], r"$\mathbf{e}_{X_{S'}}$", ls="--", scale=0.4, offset=1.25)
    arrow(ax, cB, eYS, SERIES[2], r"$\mathbf{e}_{Y_{S'}}$", ls="--", scale=0.4, offset=1.25)
    ax.plot(cB[1], cB[0], "o", color=SERIES[2], ms=5)

    # C at a point on the suction side: tangent along the stored point order, normal inward.
    i = 60
    p, q = B[i], B[i + 1]
    t = (q - p) / np.linalg.norm(q - p)
    n = np.array([-t[1], t[0]])
    if np.dot(n, cB - p) < 0:
        n = -n
    arrow(ax, p, t, SERIES[3], r"$\mathbf{e}_s$", scale=0.3, offset=1.35)
    arrow(ax, p, n, SERIES[3], r"$\mathbf{e}_n$", scale=0.3, offset=1.35)
    ax.plot(p[1], p[0], "o", color=SERIES[3], ms=4)

    le, te = B[np.argmin(naca[:, 0])], B[0]
    ax.plot(0.0, 0.0, "k+", ms=10)
    ax.text(0.02, -0.06, "pitch axis", fontsize=8, ha="right", va="top")
    ax.text(le[1] - 0.05, le[0], "LE", fontsize=9, ha="right", va="center")
    ax.text(te[1] + 0.05, te[0], "TE", fontsize=9, ha="left", va="center")

    setup_section_axes(ax)
    ax.axhline(0.0, color="0.85", lw=0.5, zorder=0)
    ax.axvline(0.0, color="0.85", lw=0.5, zorder=0)
    ax.set_title(f"Frames B, S, S′ and C on a station twisted by θ = {THETA_DEG:g}° (frames convention)")
    note(fig, f"NACA 2418, chord {CHORD:g} m, pitch axis at {PITCH:g} c. Checks: swap identity {errs[0]:.0e}, "
              f"cut rotation {errs[1]:.0e} m, parallel axis {errs[2]:.0e}, invariants {errs[3]:.0e}.")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, __file__, "frame_triad")


if __name__ == "__main__":
    naca = outline_naca()
    S = naca_to_S(naca, CHORD, PITCH)
    e1 = check_swap()
    e2 = check_cut_rotation()
    centroid, e3 = check_parallel_axis(S)
    e4 = check_invariants(S)
    figure(naca, centroid, (e1, e2, e3, e4))
