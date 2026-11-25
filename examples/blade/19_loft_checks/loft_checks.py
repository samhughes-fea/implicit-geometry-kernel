"""The loft: placement, scaling, anchored members, locked resampling and patches.

Draws Figure R8.1 of the R manual (docs/r/R8_loft_members.md) and checks:

1. placement of a unit section in the blade frame and its inverse, (R8.1)-(R8.2);
2. a thin-walled cut is built k times larger and scaled back exactly, (R8.3),
   with the defect BG-058 printed;
3. tilted webs: both feet on the line through the pivot at the tilt angle, (R8.4)-(R8.5);
4. metre anchors: an LE-frozen web keeps its distance from the leading edge, (R8.6);
5. cap edges given on the outer surface move half a cap and one skin inward, (R8.7),
   with the defect BG-059 printed;
6. forced stations collect member ends, (R8.8);
7. locked resampling: equal counts, web feet kept as skin vertices, (R8.9)-(R8.11);
8. the partition patch: not-a-knot spline exact on cubics, curvature of a cylinder,
   Newton projection, (R8.12)-(R8.16);
9. the vertex-normal skin offset, its corner error and the defect BG-060, (R8.17)-(R8.18).

    python examples/blade/19_loft_checks/loft_checks.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import REFERENCE, SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel import (  # noqa: E402
    Airfoil,
    BladeGeometry,
    PlyDrop,
    Station,
    Stringer,
    TransverseMember,
)
from implicit_geometry_kernel.blade.loft.locked_loft import (  # noqa: E402
    LockedLoftPlan,
    _resample_with_pins,
    resample_locked,
)
from implicit_geometry_kernel.blade.loft.loft_cut import _tilted_feet  # noqa: E402
from implicit_geometry_kernel.blade.loft.span_place import (  # noqa: E402
    offset_rings,
    place_chord_normal,
    place_oml,
    to_chord_normal,
)
from implicit_geometry_kernel.blade.patch import (  # noqa: E402
    PartitionPatch,
    _spline_eval,
    _spline_second,
)
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule  # noqa: E402
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout  # noqa: E402

X = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 150)))


def surfaces(code):
    af = Airfoil.naca4(code)
    return np.asarray(af.surface("upper", X)), np.asarray(af.surface("lower", X))


def ring_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def blade(af, chords=(1.0, 1.0), length=2.0, **kw):
    st = [Station(z=0.0, chord=chords[0], airfoil=af), Station(z=length, chord=chords[1], airfoil=af)]
    kw.setdefault("pivot_fraction", 0.25)
    return BladeGeometry(st, kw.pop("layout", "2D-F-LC-0R"), kw.pop("skin_thickness", 0.004), **kw)


# --- 1. placement ---------------------------------------------------------------------------


def check_placement():
    u, l = surfaces("2418")
    c, tau, p, o = 1.7, 12.0, 0.3, (0.05, -0.02)
    ring = place_oml(
        np.array([c]), np.array([tau]), np.array([p]), np.array([o]), X, u[None], l[None], np.array([3.0])
    )[0]
    unit = np.vstack([np.column_stack([X[::-1], u[::-1]]), np.column_stack([X[1:], l[1:]])])
    if np.allclose(unit[0], unit[-1]):
        unit = unit[:-1]
    t = np.radians(tau)
    q = np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])
    hand = (np.column_stack([c * (unit[:, 0] - p), c * unit[:, 1]]) @ q.T) + o
    e1 = float(np.max(np.abs(ring[:-1, :2] - hand)))
    back = to_chord_normal(ring[:-1, :2], c, tau, p, o)
    e2 = float(np.max(np.abs(back - unit)))
    again = place_chord_normal(back, c, tau, p, o, 3.0)[0]
    e3 = float(np.max(np.abs(again[:, :2] - ring[:-1, :2])))
    ea = abs(ring_area(ring[:-1, :2]) / (c * c * ring_area(unit)) - 1)
    assert max(e1, e2, e3) < 1e-14 and ea < 1e-13
    print(
        f"1. Placement (R8.1): OML against the hand map {e1:.1e} m; inverse (R8.2) {e2:.1e}; "
        f"place(inverse) {e3:.1e} m; area / (c^2 * unit area) - 1 = {ea:.1e}"
    )
    return unit, ring, (c, tau, p, o)


# --- 2. similar scaling ---------------------------------------------------------------------


def check_scaling():
    af = Airfoil.naca4("2418")
    a = blade(af, skin_thickness=0.0025, shear_web_thickness=0.004, spar_cap_thickness=0.005).section_at(1.0)
    b = blade(
        af, chords=(0.4, 0.4), skin_thickness=0.001, shear_web_thickness=0.0016, spar_cap_thickness=0.002
    ).section_at(1.0)
    worst = 0.0
    for wa, wb in zip(a.walls, b.walls):
        assert wa.label == wb.label
        worst = max(worst, float(np.max(np.abs(wb.points[:, :2] - 0.4 * wa.points[:, :2]))) / 0.4)
    ca = sum(cell.area for cell in a.cells)
    cb = sum(cell.area for cell in b.cells)
    assert a.scale == 1.0 and abs(b.scale - 2.5) < 1e-6 and worst < 1e-12 and abs(cb / (0.16 * ca) - 1) < 1e-12
    print(
        f"2. Scaling (R8.3): chord 0.4 m with a 1 mm skin is built at k = {b.scale:.4f}; walls equal 0.4 x the "
        f"1 m cut to {worst:.1e} (relative); cell area ratio / 0.16 - 1 = {cb / (0.16 * ca) - 1:.1e}"
    )
    rows = []
    for t in (0.004, 0.001):
        lay = SectionLayout(
            web_chord_fracs=(0.3, 0.6),
            orientation="vertical",
            cap_placement=CapPlacementRule.none(),
            web_stations_m=(0.3, 0.6),
        )
        cut = blade(Airfoil.naca4("0018"), layout=lay, skin_thickness=t, shear_web_thickness=0.005).section_at(1.0)
        xs = [float(np.mean(w.points[:, 0])) + 0.25 for w in cut.walls if w.kind == "web"]
        rows.append((t, cut.scale, xs))
    print(
        "   web_stations_m = (0.3, 0.6) m on a 1 m chord: "
        + "; ".join(f"skin {t * 1e3:.0f} mm, k = {k:.1f}: webs at {xs[0]:.3f}, {xs[1]:.3f} m" for t, k, xs in rows)
        + "  <- BG-058"
    )


# --- 3. tilted webs -------------------------------------------------------------------------


def check_tilt():
    u, l = surfaces("0018")
    worst = 0.0
    out = []
    for f, deg in ((0.3, 10.0), (0.5, -15.0), (0.6, 25.0)):
        ((su, pr),) = _tilted_feet((f,), (deg,), X, u, l)
        top = np.array([su, np.interp(su, X, u)])
        bot = np.array([pr, np.interp(pr, X, l)])
        th = np.radians(deg)
        d = np.array([-np.sin(th), np.cos(th)])  # suction end leans to the LE for deg > 0
        for p in (top, bot):
            r = p - np.array([f, 0.0])
            worst = max(worst, abs(r[0] * d[1] - r[1] * d[0]))
        out.append((f, deg, su, pr))
    ((su_a, pr_a),) = _tilted_feet((0.4,), (12.0,), X, u, l, anchor="suction")
    assert worst < 1e-12 and su_a == 0.4
    print(
        f"3. Tilted webs (R8.4)-(R8.5): both feet on the line through (f, 0) at the tilt, largest residual {worst:.1e}; "
        + "; ".join(f"f = {f}, {deg:+.0f} deg: feet {su:.4f} / {pr:.4f}" for f, deg, su, pr in out)
        + f"; suction anchor keeps 0.4 and moves the pressure foot to {pr_a:.4f}"
    )
    return out


# --- 4. metre anchors -----------------------------------------------------------------------


def check_anchors():
    af = Airfoil.naca4("2421")
    geo = blade(
        af, chords=(2.0, 1.0), length=4.0, skin_thickness=0.005, shear_web_thickness=0.006, web_offsets=(0.3, 0.7)
    )
    rows = []
    for z in (0.0, 2.0, 3.5):
        c = float(geo.interpolate(z)[0])
        cut = geo.section_at(z)
        xs = sorted(
            float(np.mean(to_chord_normal(w.points[:, :2], c, 0.0, 0.25, (0.0, 0.0))[:, 0])) * c
            for w in cut.walls
            if w.kind == "web"
        )
        rows.append((z, c, xs))
    worst = max(abs(x - e) for _, _, xs in rows for x, e in zip(xs, (0.3, 0.7)))
    assert worst < 1e-12
    print(
        "4. LE freeze (R8.6): web distance from the LE on a chord tapering 2 m to 1 m: "
        + "; ".join(f"z = {z:g} (c = {c:.2f}): {xs[0]:.4f}, {xs[1]:.4f} m" for z, c, xs in rows)
        + f" (largest error {worst:.1e} m)"
    )


# --- 5. cap edges on the outer surface ------------------------------------------------------


def check_cap_edges():
    af = Airfoil.naca4("4406")
    geo = blade(
        af,
        chords=(2.0, 2.0),
        skin_thickness=0.01,
        spar_cap_thickness=0.04,
        shear_web_thickness=0.005,
        cap_edge_reference="outer",
        cap_windows=lambda z: (0.2, 0.5),
    )
    c, _, _, _, (x, up, lo) = geo.interpolate(1.0)
    inset = 0.01 + 0.5 * 0.04
    rows = []
    for side, y in (("suction", up), ("pressure", lo)):
        for f in (0.2, 0.5, 0.7):
            got = geo._outer_fraction_to_midplane(1.0, f, side)
            slope = float(np.interp(f, x, np.gradient(y, x)))
            tangent = np.array([1.0, slope]) / np.hypot(1.0, slope)
            inward = np.array([-tangent[1], tangent[0]]) * (-1.0 if side == "suction" else 1.0)
            rows.append((side, f, float(np.interp(f, x, y)), got, f + inward[0] * inset / c))
    suction = max(abs(r[3] - r[4]) for r in rows if r[0] == "suction")
    assert suction < 1e-4
    print(
        f"5. Cap edges from the outer surface (R8.7), NACA 4406, c = 2 m, inset {inset * 1e3:.0f} mm: suction side "
        f"within {suction:.1e}; pressure side "
        + ", ".join(f"f = {f}: {g:.5f} against {e:.5f}" for s, f, _, g, e in rows if s == "pressure")
        + "  <- BG-059 (moved outward where the pressure surface is above the chord line)"
    )


# --- 6. forced stations ---------------------------------------------------------------------


def check_forced():
    geo = blade(
        Airfoil.naca4("2418"),
        length=10.0,
        skin_thickness=0.005,
        transverse=(TransverseMember("rib_1", 3.0, 0.01),),
        stringers=(Stringer("s0", 0.5, 0.004, 0.02, extent=(1.5, 6.0)),),
        ply_drops=(PlyDrop("pd", 0.2, 0.4, 2.5, 7.5, 0.003, taper=0.5),),
    )
    got = geo.forced_stations()
    assert got == (1.5, 2.5, 3.0, 6.0, 7.5)
    print(
        f"6. Forced stations (R8.8): rib at 3, stringer over [1.5, 6], ply drop over [2.5, 7.5] with a 0.5 m taper: "
        f"{got}; the taper ends 2.0 and 8.0 are not included"
    )


# --- 7. locked resampling -------------------------------------------------------------------


def check_locked():
    geo = blade(
        Airfoil.naca4("2418"),
        chords=(2.0, 2.0),
        skin_thickness=0.01,
        shear_web_thickness=0.01,
        spar_cap_thickness=0.02,
        layout="3D-F-LC-0R",
    )
    graph = geo.section_at(1.0).shell_graph()
    plan = LockedLoftPlan.from_graph(graph, n_chord=41)
    out = resample_locked(graph, plan)
    counts = {len(v) for v in out.values()}
    loop = np.vstack([out[n] for n in plan.loop_names])
    gap = 0.0
    for name in plan.web_names:
        for end in (out[name][0], out[name][-1]):
            gap = max(gap, float(np.min(np.linalg.norm(loop - end, axis=1))))
    shared = 0.0
    for edge in graph.partition_edges():
        if edge.name_a in out and edge.name_b in out:
            pa, pb = out[edge.name_a], out[edge.name_b]
            da = min(np.linalg.norm(pa[0] - edge.point), np.linalg.norm(pa[-1] - edge.point))
            db = min(np.linalg.norm(pb[0] - edge.point), np.linalg.norm(pb[-1] - edge.point))
            shared = max(shared, float(da), float(db))
    pts = np.column_stack([np.linspace(0.0, 1.0, 23), np.zeros(23), np.zeros(23)])
    pts[:, 1] = np.sin(3 * pts[:, 0])
    pins = [pts[7], pts[15]]
    res = _resample_with_pins(pts, 12, pins, 1e-12)
    kept = [int(np.argmin(np.linalg.norm(res - pin, axis=1))) for pin in pins]
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    expect = [int(round(s[j] / s[-1] * 11)) for j in (7, 15)]
    assert counts == {41} and gap == 0.0 and shared == 0.0 and kept == expect
    print(
        f"7. Locked resampling (R8.9)-(R8.11), 3D-F: {len(out)} runs at 41 points each; web ends on skin vertices "
        f"(gap {gap:.0e} m); partition points shared exactly ({shared:.0e}); pins of a 23-point run kept at "
        f"output indices {kept} = round(s_j / L * (n - 1)) = {expect}"
    )


# --- 8. partition patch ---------------------------------------------------------------------


def check_patch():
    t = np.array([0.0, 0.1, 0.35, 0.5, 0.8, 1.0])
    y = 2 * t**3 - t**2 + 0.5 * t - 0.2
    m = _spline_second(t, y)
    tq = np.linspace(0, 1, 41)
    err = max(abs(float(_spline_eval(t, y, m, q, 0)) - (2 * q**3 - q**2 + 0.5 * q - 0.2)) for q in tq)
    err2 = max(abs(float(_spline_eval(t, y, m, q, 2)) - (12 * q - 2)) for q in tq)
    assert err < 1e-13 and err2 < 1e-11
    radius, ts = 0.5, 0.01
    geo = blade(
        Airfoil.ellipse(1.0, n=241),
        chords=(2 * radius, 2 * radius),
        length=2.0,
        layout="0A-LC-0R",
        skin_thickness=ts,
        span_interp="linear",
    )
    graph = geo.section_at(0.5).shell_graph()
    name = graph.loop[0].name
    patch = PartitionPatch(geo, name, (0.0, 0.5, 1.0, 1.5, 2.0), n_t=33)
    r_mid = radius - 0.5 * ts
    ks, tt = [], np.linspace(0.1, 0.9, 17)
    for q in tt:
        sp = patch.at(float(q), 1.0)
        ks.append(max(abs(sp.k1), abs(sp.k2)))
    ks = np.array(ks)
    kerr = float(np.max(np.abs(ks * r_mid - 1)))
    probe = patch(0.4, 1.2)
    centre = np.array([(0.5 - 0.25) * 2 * radius, 0.0])
    radial = (probe[:2] - centre) / np.linalg.norm(probe[:2] - centre)
    p = probe + 0.03 * np.concatenate([radial, [0.0]])
    tq_, zq_, dist = patch.closest(p, t0=0.5, z0=1.0)
    assert kerr < 2e-3 and abs(abs(dist) - 0.03) < 1e-4 and abs(zq_ - 1.2) < 1e-6
    print(
        f"8. Partition patch (R8.12)-(R8.16): not-a-knot spline reproduces a cubic to {err:.1e} (second derivative "
        f"{err2:.1e}); cylinder r_mid = {r_mid} m: |k| r_mid - 1 within {kerr:.1e} over t in [0.1, 0.9]; a point "
        f"0.03 m off the patch projects back to z = {zq_:.6f} at distance {abs(dist):.5f} m"
    )
    return tt, ks, r_mid


# --- 9. vertex-normal offset ----------------------------------------------------------------


def check_offset():
    d = 0.01
    rows = []
    for n in (6, 12, 24):
        a = np.linspace(0, 2 * np.pi, n + 1)
        ring = np.column_stack([np.cos(a), np.sin(a), np.zeros(n + 1)])[None]
        moved = offset_rings(ring, d)[0, :-1, :2]
        edge_dist = float(np.cos(np.pi / n)) - float(np.linalg.norm((moved[0] + moved[1]) / 2))
        rows.append((n, edge_dist, d * np.cos(np.pi / n)))
    worst = max(abs(e - x) for _, e, x in rows)
    assert worst < 1e-14
    flips = {}
    for code in ("0012", "4412", "4406", "6409"):
        up, lo = surfaces(code)
        ring = place_oml(np.ones(1), np.zeros(1), np.zeros(1), np.zeros((1, 2)), X, up[None], lo[None], np.zeros(1))
        xy = ring[0, :-1, :2]
        mv = offset_rings(ring, 0.002)[0, :-1, :2] - xy
        tang = np.roll(xy, -1, 0) - np.roll(xy, 1, 0)
        inward = np.column_stack([-tang[:, 1], tang[:, 0]])
        flips[code] = (int(np.sum(np.sum(inward * mv, axis=1) < 0)), len(xy))
    print(
        f"9. Vertex-normal offset (R8.17)-(R8.18): on a regular N-gon the edge moves d cos(pi/N), to {worst:.0e} "
        f"(N = 6: {rows[0][1] / d:.4f} d). Vertices moved outward: "
        + ", ".join(f"NACA {k} {v[0]} of {v[1]}" for k, v in flips.items())
        + "  <- BG-060"
    )
    return flips


# --- figure ---------------------------------------------------------------------------------


def figure(unit, ring, params, tilts, curv):
    c, tau, p, o = params
    fig, axes = plt.subplots(
        2, 2, figsize=(10.5, 6.6), gridspec_kw={"wspace": 0.22, "hspace": 0.55, "height_ratios": [1.0, 1.25]}
    )

    ax = axes[0, 0]
    ax.plot(unit[:, 0], unit[:, 1], color=SERIES[0], lw=1, label="unit section $(\\xi, \\eta)$")
    scaled = np.column_stack([c * (unit[:, 0] - p), c * unit[:, 1]])
    ax.plot(scaled[:, 0], scaled[:, 1], color=SERIES[2], lw=1, ls="--", label="scaled, pitch axis at 0")
    ax.plot(ring[:, 0], ring[:, 1], color=SERIES[1], lw=1.4, label=f"twisted {tau:g} deg and offset (R8.1)")
    ax.plot([0], [0], "+", color=REFERENCE, ms=9)
    ax.plot([o[0]], [o[1]], "x", color=SERIES[1], ms=6)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$y$ (m)")
    ax.set_title("Placing a section in the blade frame", fontsize=9)
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.3))

    ax = axes[0, 1]
    up, lo = surfaces("0018")
    ax.plot(X, up, color=REFERENCE, lw=0.8)
    ax.plot(X, lo, color=REFERENCE, lw=0.8)
    ax.axhline(0, color="0.7", lw=0.5)
    for k, (f, deg, su, pr) in enumerate(tilts):
        col = SERIES[k]
        ax.plot(
            [su, pr],
            [np.interp(su, X, up), np.interp(pr, X, lo)],
            color=col,
            lw=1.4,
            label=f"$f$ = {f}, {deg:+.0f} deg",
        )
        ax.plot([f], [0], "o", color=col, ms=3)
    ax.set_aspect("equal")
    ax.set_xlabel("$\\xi$")
    ax.set_ylabel("$\\eta$")
    ax.set_title("Tilted webs pivot about $(f, 0)$ (R8.4)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.42))

    ax = axes[1, 0]
    tt, ks, r_mid = curv
    ax.plot(tt, ks * r_mid, "o-", color=SERIES[0], ms=3, label="$|\\kappa|\\, r_m$ from (R8.15)")
    ax.axhline(1.0, color=REFERENCE, ls="--", lw=0.8, label="exact")
    ax.set_ylim(0.995, 1.005)
    ax.set_xlabel("run coordinate $t$")
    ax.set_ylabel("$|\\kappa|\\, r_m$")
    ax.set_title("Patch curvature on a cylinder", fontsize=9)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    up, lo = surfaces("4406")
    ring = place_oml(np.ones(1), np.zeros(1), np.zeros(1), np.zeros((1, 2)), X, up[None], lo[None], np.zeros(1))
    xy = ring[0, :-1, :2]
    d = 0.004
    moved = offset_rings(ring, d)[0, :-1, :2]
    tang = np.roll(xy, -1, 0) - np.roll(xy, 1, 0)
    inward = np.column_stack([-tang[:, 1], tang[:, 0]])
    bad = np.sum(inward * (moved - xy), axis=1) < 0
    ax.plot(xy[:, 0], xy[:, 1], color=REFERENCE, lw=0.8, label="outer mould line, NACA 4406")
    ax.plot(moved[~bad, 0], moved[~bad, 1], ".", color=SERIES[0], ms=2, label="moved inward")
    ax.plot(moved[bad, 0], moved[bad, 1], ".", color=SERIES[1], ms=2, label="moved outward (BG-060)")
    ax.plot([xy[:, 0].mean()], [xy[:, 1].mean()], "+", color=REFERENCE, ms=8, label="vertex mean")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.03, 0.09)
    ax.set_xlabel("$\\xi$")
    ax.set_ylabel("$\\eta$")
    ax.set_title(f"Vertex-normal offset by d = {d} (R8.17)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    note(
        fig,
        "Checks: placement; scaling (BG-058); tilted webs; LE freeze; cap edges (BG-059); forced stations; "
        "locked resampling; partition patch; vertex-normal offset (BG-060).",
    )
    fig.subplots_adjust(left=0.07, right=0.98, top=0.96, bottom=0.1)
    save(fig, __file__, "loft_checks")


if __name__ == "__main__":
    unit, ring, params = check_placement()
    check_scaling()
    tilts = check_tilt()
    check_anchors()
    check_cap_edges()
    check_forced()
    check_locked()
    curv = check_patch()
    check_offset()
    figure(unit, ring, params, tilts, curv)
