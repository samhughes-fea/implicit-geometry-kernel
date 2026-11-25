"""Mid-surfaces, cells and the thin-wall idealisation.

Draws Figure R6.1 of the R manual (docs/r/R6_midsurfaces_cells.md) and checks:

1. the mid-surface depth of each reference face, (R6.1);
2. the skin mid ring is the exact offset locus at that depth, (R6.2);
3. the offset length and area laws and the band identity A_band = t L_mid, (R6.3)-(R6.6);
4. the trailing-edge retreat of the mid ring, (R6.7);
5. the thin-wall second moment and its O(t^2) error, exact on a circle, (R6.8)-(R6.10);
6. the cap mid-surface depth and the merged-laminate eccentricity h/2, (R6.11)-(R6.12);
7. web feet on the skin mid ring and their reach past the inner skin, (R6.13)-(R6.14);
8. the wall graph: n + 1 cells, Euler's count and the area identity, (R6.15)-(R6.17).

    python examples/section/05_midsurfaces_cells/midsurfaces_cells.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import NAVY, RED, REFERENCE, SERIES, TEAL, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, resolve_system_type  # noqa: E402
from implicit_geometry_kernel.frames import naca_to_twisted_pitch_axis_blade_frame  # noqa: E402
from implicit_geometry_kernel.section.boundary import airfoil_offset_ring  # noqa: E402
from implicit_geometry_kernel.section.cells.geom import shoelace_area  # noqa: E402
from implicit_geometry_kernel.section.layout.anchoring import airfoil_inward_depth_for_midsurface  # noqa: E402
from implicit_geometry_kernel.section.spec import SectionGeometry, SectionOptions, SectionSpec  # noqa: E402

CHORD, SKIN, CAP, WEB = 2.0, 0.008, 0.03, 0.015
PIVOT = 0.35
TWIST = np.radians(15.0)
AF = Airfoil.naca4("2418").to_sdf(chord=CHORD)
OML = np.asarray(AF.vertices, dtype=float)


def edge_distance(points, ring):
    """Unsigned distance from each point to a closed polygon (exact, per edge)."""
    a = ring
    b = np.roll(ring, -1, axis=0)
    e = b - a
    t = np.clip(np.einsum("mnk,nk->mn", points[:, None, :] - a[None], e) / np.sum(e * e, axis=1), 0.0, 1.0)
    q = a[None] + t[..., None] * e[None]
    return np.min(np.linalg.norm(points[:, None, :] - q, axis=2), axis=1)


def perimeter(ring):
    return float(np.sum(np.linalg.norm(np.roll(ring, -1, axis=0) - ring, axis=1)))


def second_moment_x(ring):
    """|I_xx| about y = 0 of a closed polygon by the shoelace rule."""
    x, y = ring[:, 0], ring[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    c = x * y1 - x1 * y
    return abs(float(np.sum((y * y + y * y1 + y1 * y1) * c)) / 12.0)


def thin_wall_second_moment(mid, t):
    """t * integral of y^2 ds along the closed mid ring (exact per straight segment)."""
    p = np.vstack([mid, mid[:1]])
    ds = np.linalg.norm(np.diff(p, axis=0), axis=1)
    y0, y1 = p[:-1, 1], p[1:, 1]
    return t * float(np.sum(ds * (y0 * y0 + y0 * y1 + y1 * y1) / 3.0))


def turning_angles(ring):
    e = np.roll(ring, -1, axis=0) - ring
    a = np.arctan2(e[:, 1], e[:, 0])
    return np.angle(np.exp(1j * (a - np.roll(a, 1))))


def run(key, twist=TWIST):
    spec = SectionSpec.from_kwargs(
        resolve_system_type(key),
        AF,
        skin_thickness=SKIN,
        spar_cap_thickness=CAP,
        shear_web_thickness=WEB,
        twist_rad=twist,
        pivot_fraction=PIVOT,
    )
    return SectionGeometry(spec, SectionOptions(outputs=("midlines", "cell_wall_loops"))).execute().get_results()


def to_b(xy, twist=TWIST):
    return naca_to_twisted_pitch_axis_blade_frame(np.asarray(xy, float), float(twist), CHORD, PIVOT)


# --- 1-2. mid-surface depth and the offset locus ------------------------------------------


def check_mid_depth():
    t = SKIN
    depths = {
        ref: airfoil_inward_depth_for_midsurface(
            ref, t, reference_iml_depth_m=0.02 if ref != "outer_mould_line" else 0.0
        )
        for ref in ("outer_mould_line", "inner_mould_line", "midsurface", "skin_bond_line")
    }
    assert depths == {
        "outer_mould_line": t / 2,
        "inner_mould_line": 0.02 - t / 2,
        "midsurface": 0.02,
        "skin_bond_line": 0.02 + t / 2,
    }
    print(
        f"1. Mid depth d_m for t = {t * 1e3:g} mm, IML depth 20 mm: "
        + ", ".join(f"{k} {v * 1e3:g} mm" for k, v in depths.items())
    )
    mid = airfoil_offset_ring(AF, t / 2)
    err = float(np.max(np.abs(edge_distance(mid, OML) - t / 2)))
    assert err < 1e-12 and np.all(turning_angles(OML) > 0.0)
    print(
        f"2. Skin mid ring: {len(mid)} vertices, max |dist(OML) - t/2| = {err:.1e} m; the aerofoil polygon "
        "has no reflex vertex, so the mitred offset is the exact locus"
    )


# --- 3-5. offset laws, band identity, TE retreat, second moment ---------------------------


def check_offset_laws():
    tau = turning_angles(OML)
    big_t = float(np.sum(np.tan(tau / 2.0)))
    length, area = perimeter(OML), abs(shoelace_area(OML))
    edge = np.linalg.norm(np.roll(OML, -1, axis=0) - OML, axis=1)
    vanish = edge / (np.tan(tau / 2.0) + np.tan(np.roll(tau, -1) / 2.0))
    d_star = float(vanish.min())
    print(
        f"3. Turning angles sum to {tau.sum() / np.pi:.12f} pi; T = sum tan(tau/2) = {big_t:.4f} "
        f"(pi + {big_t - np.pi:.4f} from the TE wedge); first edge vanishes at d* = {d_star * 1e3:.3f} mm"
    )
    rows = []
    for d in (0.5e-3, 1e-3, 2e-3, 2.5e-3):
        ring = airfoil_offset_ring(AF, d)
        e_len = perimeter(ring) - (length - 2.0 * d * big_t)
        e_area = abs(shoelace_area(ring)) - (area - d * length + d * d * big_t)
        assert abs(e_len) < 1e-12 and abs(e_area) < 1e-12, (d, e_len, e_area)
        rows.append(f"d = {d * 1e3:g} mm: dL = {e_len:.0e}, dA = {e_area:.0e}")
    print("   (R6.3)-(R6.4) below d*: " + "; ".join(rows))
    band = []
    for t in (1e-3, 2e-3, 4e-3, 8e-3):
        exact = area - abs(shoelace_area(airfoil_offset_ring(AF, t)))
        thin = t * perimeter(airfoil_offset_ring(AF, t / 2))
        band.append((t, exact / thin - 1.0))
    assert abs(band[0][1]) < 1e-10 and abs(band[1][1]) < 1e-10 and abs(band[3][1]) > 1e-6
    print(
        "   (R6.6) band area / (t L_mid) - 1: "
        + ", ".join(f"t = {t * 1e3:g} mm {r:+.1e}" for t, r in band)
        + f"; exact while t <= d* = {d_star * 1e3:.2f} mm"
    )

    beta = np.pi - tau[0]
    rows = []
    for d in (0.5e-3, 1e-3, 2e-3):
        ring = airfoil_offset_ring(AF, d)
        tip = ring[np.argmax(ring[:, 0])]
        r = float(np.linalg.norm(tip - OML[0]))
        assert abs(r / (d / np.sin(beta / 2.0)) - 1.0) < 1e-9
        rows.append(f"d = {d * 1e3:g} mm: {r * 1e3:.4f} mm")
    print(
        f"4. TE wedge beta = {np.degrees(beta):.3f} deg; mid-ring tip lies d / sin(beta/2) from the TE (R6.7): "
        + ", ".join(rows)
    )
    return d_star


def check_second_moment():
    radius = 0.5
    rows = []
    for t in (0.01, 0.05, 0.2):
        exact = np.pi / 4.0 * ((radius + t / 2) ** 4 - (radius - t / 2) ** 4)
        formula = np.pi * radius**3 * t + np.pi * radius * t**3 / 4.0
        assert abs(formula / exact - 1.0) < 1e-14
        rows.append(f"t/R = {t / radius:g}: thin-wall {(np.pi * radius**3 * t) / exact - 1:+.2e}")
    print("5. Circle R = 0.5 m: (R6.9) with the t^3 term equals the annulus I exactly; " + ", ".join(rows))
    errors = []
    for t in (1e-3, 2e-3, 4e-3, 8e-3, 16e-3):
        exact = second_moment_x(OML) - second_moment_x(airfoil_offset_ring(AF, t))
        thin = thin_wall_second_moment(airfoil_offset_ring(AF, t / 2), t)
        errors.append((t, thin / exact - 1.0))
    ts, es = np.array(errors).T
    slope = np.polyfit(np.log(ts[:3]), np.log(np.abs(es[:3])), 1)[0]
    assert 1.9 < slope < 2.1
    print(
        "   Aerofoil, I_xx about the chord: thin-wall error "
        + ", ".join(f"t = {t * 1e3:g} mm {e:+.1e}" for t, e in errors)
        + f"; log-log slope {slope:.2f} (R6.10)"
    )
    return ts, es


# --- 6. cap mid-surface and eccentricity --------------------------------------------------


def check_cap_eccentricity():
    res = run("2D-F-LC-0R", twist=0.0)
    caps = {str(getattr(k, "value", k)): c for k, c in res.section.spar_cap_components_unrotated.items()}
    worst = 0.0
    for name, cap in caps.items():
        mid = np.asarray(cap.midline_polyline(n=120), float)
        worst = max(worst, float(np.max(np.abs(edge_distance(mid, OML) - (SKIN + CAP / 2)))))
    assert worst < 1e-9
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("2418"), pitch_axis=PIVOT),
        Station(10.0, 1.0, Airfoil.naca4("2418"), pitch_axis=PIVOT),
    ]
    geo = BladeGeometry(
        stations,
        "2D-F-LC-0R",
        SKIN,
        spar_cap_thickness=CAP,
        shear_web_thickness=WEB,
        web_positions=(0.2, 0.5),
        pivot_fraction=PIVOT,
    )
    cut = geo.section_at(0.0)
    graph = cut.shell_graph()
    runs = [r for r in graph.loop if r.kind == "cap"]
    assert runs and all(r.thickness == SKIN and r.eccentricity_m == CAP / 2 and r.requires_eccentricity for r in runs)
    p = np.asarray(runs[0].points[len(runs[0].points) // 2]).reshape(1, 3)
    oml, cavity = float(geo.oml_sdf().evaluate(p).phi[0]), float(geo.cavity_sdf().evaluate(p).phi[0])
    assert abs(oml + SKIN / 2) < 1e-4 and abs(cavity - (SKIN / 2 + CAP)) < 1e-3
    print(
        f"6. Cap solid mid-surface at depth t_s + h/2 = {(SKIN + CAP / 2) * 1e3:g} mm to {worst:.0e} m (R6.11); "
        f"{len(runs)} cap runs on the skin mid ring carry t = t_s and e = h/2 = {CAP / 2 * 1e3:g} mm (R6.12); "
        f"at a cap node phi_OML = {oml * 1e3:.3f} mm, phi_cavity = {cavity * 1e3:.3f} mm"
    )


# --- 7. web feet -------------------------------------------------------------------------


def line_ring_hit(p, u, ring):
    """Intersections p + s u with the closed ring; returns the s values."""
    a, b = ring, np.roll(ring, -1, axis=0)
    e = b - a
    out = []
    for ai, ei in zip(a, e):
        den = u[0] * ei[1] - u[1] * ei[0]
        if abs(den) < 1e-15:
            continue
        w = ai - p
        s = (w[0] * ei[1] - w[1] * ei[0]) / den
        r = (w[0] * u[1] - w[1] * u[0]) / den
        if -1e-12 <= r <= 1.0 + 1e-12:
            out.append(s)
    return np.array(out)


def check_web_feet():
    res = run("3D-F-LC-0R")
    mid_b, iml_b, oml_b = to_b(airfoil_offset_ring(AF, SKIN / 2)), to_b(airfoil_offset_ring(AF, SKIN)), to_b(OML)
    feet = []
    for strip in (s for s in res.midline_strips if s.kind == "web"):
        pts = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, float)
        u = (pts[-1] - pts[0]) / np.linalg.norm(pts[-1] - pts[0])
        for end, sign in ((pts[0], -1.0), (pts[-1], 1.0)):
            on_mid = float(edge_distance(end.reshape(1, 2), mid_b)[0])
            d_oml = float(edge_distance(end.reshape(1, 2), oml_b)[0])
            s_iml = line_ring_hit(end, u, iml_b)
            reach = float(np.min(np.abs(s_iml)))
            # skin normal at the foot from the nearest mid-ring edge
            k = int(np.argmin(np.linalg.norm(mid_b - end, axis=1)))
            tang = mid_b[(k + 1) % len(mid_b)] - mid_b[k - 1]
            n = np.array([-tang[1], tang[0]]) / np.linalg.norm(tang)
            cos_g = abs(float(u @ n))
            feet.append((strip.label, on_mid, d_oml, reach, SKIN / 2 / cos_g, np.degrees(np.arccos(cos_g))))
    worst_mid = max(f[1] for f in feet)
    worst_depth = max(abs(f[2] - SKIN / 2) for f in feet)
    worst_reach = max(abs(f[3] / f[4] - 1.0) for f in feet)
    assert worst_mid < 1e-9 and worst_depth < 1e-9 and worst_reach < 1e-3
    gam = [f[5] for f in feet]
    print(
        f"7. {len(feet)} web feet (3D-F, 15 deg twist): on the skin mid ring to {worst_mid:.0e} m, depth t_s/2 to "
        f"{worst_depth:.0e} m (R6.13); reach past the inner skin = (t_s/2)/cos(gamma) to {worst_reach:.1e} "
        f"for gamma in [{min(gam):.1f}, {max(gam):.1f}] deg (R6.14)"
    )
    return res


# --- 8. cells ----------------------------------------------------------------------------


def check_cells():
    rows = []
    for n in range(0, 6):
        key = "0A-LC-0R" if n == 0 else f"{n}A-CN-LC-0R"
        res = run(key)
        strips = res.midline_strips
        skins = [s for s in strips if s.kind == "skin"]
        webs = [s for s in strips if s.kind == "web"]
        loops = res.cell_wall_loops
        ends = [np.asarray(s.midline_twisted_pitch_axis_blade_frame)[[0, -1]] for s in skins]
        nodes = {tuple(np.round(p, 9)) for e in ends for p in e}
        v, e = len(nodes), len(skins) + len(webs)
        faces = e - v + 2  # Euler, V - E + F = 2 for the connected planar wall graph
        ring = (
            np.vstack([np.asarray(s.midline_twisted_pitch_axis_blade_frame)[:-1] for s in skins])
            if len(skins) > 1
            else np.asarray(skins[0].midline_twisted_pitch_axis_blade_frame)
        )
        ring_area = abs(shoelace_area(ring))
        total = sum(c.area_m2 for c in loops)
        uses = {s.label: sum(s.label in c.wall_labels for c in loops) for s in strips if s.kind in ("skin", "web")}
        assert (v, e) == (2 * n + 2, 3 * n + 2) and len(webs) == n
        assert len(loops) == n + 1 and faces - 1 == n + 1
        assert abs(total / ring_area - 1.0) < 1e-12
        assert all(uses[w.label] == 2 for w in webs) and all(uses[s.label] == 1 for s in skins)
        rows.append(f"n = {n}: V = {v}, E = {e}, cells {len(loops)}, sum A / A_ring - 1 = {total / ring_area - 1:+.0e}")
    print("8. Wall graph (R6.15)-(R6.17): " + "; ".join(rows) + ". Every web bounds two cells, every skin arc one")


# --- figure ------------------------------------------------------------------------------


def figure(res, d_star, ts, es):
    fig = plt.figure(figsize=(10.5, 7.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.8, 1.05], wspace=0.38, hspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    oml_b = to_b(OML)
    loops = res.cell_wall_loops
    for i, c in enumerate(loops):
        p = np.asarray(c.polyline_twisted_pitch_axis_blade_frame)
        ax.fill(p[:, 1], p[:, 0], color=SERIES[i % len(SERIES)], alpha=0.16, lw=0)
        cy, cx = p[:, 1].mean(), p[:, 0].mean()
        ax.annotate(f"cell {c.id}\n{c.area_m2:.3f} m²", (cy, cx), ha="center", va="center", fontsize=7)
    ax.plot(*np.vstack([oml_b, oml_b[:1]])[:, ::-1].T, color=REFERENCE, lw=0.6, label="OML")
    for s in res.midline_strips:
        p = np.asarray(s.midline_twisted_pitch_axis_blade_frame)
        col = TEAL if s.kind == "web" else NAVY
        ax.plot(p[:, 1], p[:, 0], color=col, lw=1.1)
        ax.plot(p[[0, -1], 1], p[[0, -1], 0], "o", color=col, ms=2.5)
    ax.plot([], [], color=NAVY, lw=1.2, label="skin mid ring, split into labelled arcs")
    ax.plot([], [], color=TEAL, lw=1.2, label="web mid-surfaces, feet on the skin mid ring")
    ax.set_aspect("equal")
    ax.set_xlabel("$Y_B$ (m), chordwise")
    ax.set_ylabel("$X_B$ (m), flapwise")
    ax.set_title("3D-F, NACA 2418, 2 m chord, 15° twist: n = 3 webs, n + 1 = 4 cells", fontsize=9)
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)

    ax = fig.add_subplot(gs[1, 0])
    for d, col, lab in (
        (0.0, REFERENCE, "OML"),
        (SKIN / 2, NAVY, "mid ring, $d = t_s/2$"),
        (SKIN, SERIES[1], "IML, $d = t_s$"),
    ):
        ring = OML if d == 0 else airfoil_offset_ring(AF, d)
        r = np.vstack([ring, ring[:1]])
        ax.plot(r[:, 0], r[:, 1], color=col, lw=1.0, label=lab, marker="." if d == 0 else None, ms=3)
    beta = np.pi - turning_angles(OML)[0]
    ax.annotate(f"tip retreat $d/\\sin(\\beta/2)$\n$\\beta$ = {np.degrees(beta):.1f}°", (1.957, -0.03), fontsize=7)
    ax.set_xlim(1.95, 2.004)
    ax.set_ylim(-0.035, 0.03)
    ax.set_aspect("equal")
    ax.set_title("TE: the ring ends short (R6.7)", fontsize=9)
    ax.set_xlabel("$x$ (m), chord frame", fontsize=8)
    ax.legend(frameon=False, fontsize=7, loc="upper left")

    ax = fig.add_subplot(gs[1, 1])
    ax.loglog(ts * 1e3, np.abs(es), "o-", color=SERIES[0], ms=3.5, label="$|I_{thin}/I - 1|$, aerofoil")
    ax.loglog(ts * 1e3, np.abs(es[0]) * (ts / ts[0]) ** 2, "--", color=REFERENCE, lw=0.8, label="slope 2 (R6.10)")
    ax.axvline(d_star * 1e3, color="0.5", lw=0.6, ls=":")
    ax.annotate("$t = d^*$", (d_star * 1e3 * 1.05, np.abs(es).max() / 2), fontsize=7, color="0.4")
    ax.set_xlabel("skin thickness $t$ (mm)", fontsize=8)
    ax.set_title("Thin-wall $I_{xx}$ error (R6.10)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper left")

    ax = fig.add_subplot(gs[1, 2])
    web = next(s for s in res.midline_strips if s.label == "shear_web_1")
    w = np.asarray(web.midline_twisted_pitch_axis_blade_frame)
    foot = w[-1]
    for d, col, lab in ((0.0, REFERENCE, "OML"), (SKIN / 2, NAVY, "skin mid ring"), (SKIN, SERIES[1], "IML")):
        ring = to_b(OML if d == 0 else airfoil_offset_ring(AF, d))
        r = np.vstack([ring, ring[:1]])
        ax.plot(r[:, 1], r[:, 0], color=col, lw=1.0, label=lab)
    u = (w[-1] - w[0]) / np.linalg.norm(w[-1] - w[0])
    nrm = np.array([-u[1], u[0]])
    for side in (-1, 1):
        edge = w + side * WEB / 2 * nrm
        ax.plot(edge[:, 1], edge[:, 0], color=TEAL, lw=0.6, ls="--")
    ax.plot(w[:, 1], w[:, 0], color=TEAL, lw=1.3, label="web mid-surface")
    ax.plot(foot[1], foot[0], "o", color=RED, ms=4, label="foot: axis $\\cap$ mid ring (R6.13)")
    ax.set_xlim(foot[1] - 0.022, foot[1] + 0.022)
    ax.set_ylim(foot[0] - 0.035, foot[0] + 0.02)
    ax.set_aspect("equal")
    ax.set_xlabel("$Y_B$ (m)", fontsize=8)
    ax.set_ylabel("$X_B$ (m)", fontsize=8)
    ax.set_title("Web foot (R6.13)-(R6.14)", fontsize=9)
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    note(
        fig,
        "Checks: mid depths; exact offset locus; offset laws and band identity; TE retreat; thin-wall I; "
        "cap eccentricity; web feet; n + 1 cells and the area identity.",
    )
    fig.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.09)
    save(fig, __file__, "midsurfaces_cells")


if __name__ == "__main__":
    check_mid_depth()
    d_star = check_offset_laws()
    ts, es = check_second_moment()
    check_cap_eccentricity()
    res = check_web_feet()
    check_cells()
    figure(res, d_star, ts, es)
