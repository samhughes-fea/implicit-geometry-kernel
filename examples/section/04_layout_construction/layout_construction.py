"""How a layout key becomes walls: the catalogue, the wall fields and the web anchoring.

Draws Figure R5.1 of the R manual (docs/r/R5_section_layout.md) and checks:

1. the catalogue grammar and its counts, 40 cores and 120 keys, (R5.1)-(R5.2);
2. the step-down ladders and the cap run-out, (R5.3);
3. the default web stations, (R5.4);
4. the skin, cap and web fields: wall depths and widths exact, (R5.5)-(R5.9);
5. web axes: chord-normal and flapwise directions and inner-skin feet, (R5.10)-(R5.12);
6. the cap bands of families B, C and D, (R5.13)-(R5.15);
7. web stations along the span for LC, LE and PA, (R5.16)-(R5.18).

    python examples/section/04_layout_construction/layout_construction.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import NAVY, RED, REFERENCE, SERIES, TEAL, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, resolve_system_type  # noqa: E402
from implicit_geometry_kernel.blade.concepts.web_anchoring import (  # noqa: E402
    blade_x_le_te,
    chord_fraction_from_pitch,
    le_web_offsets_m,
    pa_web_x_from_pitch,
)
from implicit_geometry_kernel.catalogue.system_layout import (  # noqa: E402
    SYSTEM_TYPE_KEYS,
    SYSTEM_TYPE_LAYOUT_KEYS,
    layout_without_spar_caps,
    step_down_cascade,
)
from implicit_geometry_kernel.sdf import eval_expr  # noqa: E402
from implicit_geometry_kernel.section.assembly.from_layout import multicell_from_layout  # noqa: E402
from implicit_geometry_kernel.section.layout.section_axes import max_thickness_chord_x  # noqa: E402

CHORD, SKIN, CAP, WEB = 2.0, 0.008, 0.03, 0.015
TWIST = np.radians(15.0)


def section(key, twist=TWIST, chord=CHORD, code="2418"):
    af = Airfoil.naca4(code).to_sdf(chord=chord)
    sec = multicell_from_layout(af, resolve_system_type(key), twist_rad=twist, skin_thickness=SKIN,
                                spar_cap_thickness=CAP, shear_web_thickness=WEB)
    return af, sec


def comp(sec, name):
    for k, c in sec.components_unrotated.items():
        if str(getattr(k, "value", k)) == name:
            return c
    raise KeyError(name)


# --- 1-3. catalogue ---------------------------------------------------------------------

def check_catalogue():
    cores = {k.rsplit("-", 2)[0] for k in SYSTEM_TYPE_LAYOUT_KEYS}
    expect = {"0A", "0B"} | {f"{w}{x}-{y}" for w in range(1, 6) for x in "ABCD" for y in ("CN", "F")
                             if not (w == 1 and x == "D")}
    assert cores == expect and len(SYSTEM_TYPE_LAYOUT_KEYS) == 40 and len(SYSTEM_TYPE_KEYS) == 120
    refused = []
    for bad in ("2D-F", "2D-F-LC", "2D-F-LC-R", "1D-F-LC-0R", "0C-LC-0R", "2D-F-XX-0R"):
        try:
            resolve_system_type(bad)
        except KeyError:
            refused.append(bad)
    assert len(refused) == 6
    spec = resolve_system_type("4D-F-PA-10R")
    assert (spec.n_webs, spec.n_cells, spec.structural_family, spec.web_anchoring, spec.n_ribs) == \
        (4, 5, "D", "pitch_axis", 10)
    print(f"1. Catalogue: {len(cores)} cores = 2 + 5 x 4 x 2 - 2 (no 1D); {len(SYSTEM_TYPE_KEYS)} keys = 40 x 3 "
          f"anchorings at 0R; refused {', '.join(refused)}")
    ladders = {f: step_down_cascade(resolve_system_type(f"5{f}-CN-LC-0R")) for f in "ABCD"}
    for f, lad in ladders.items():
        cores_l = [k.rsplit("-", 2)[0] for k in lad]
        want = [f"{n}{f}-CN" for n in range(5, 0, -1)] if f != "D" else \
            ["5D-CN", "4D-CN", "3D-CN", "2D-CN", "1C-CN"]
        assert cores_l == want, (f, cores_l)
        print(f"2. Step-down from 5{f}-CN: {' -> '.join(cores_l)}")
    runout = layout_without_spar_caps(resolve_system_type("3D-F-LE-4R")).key
    assert runout == "3A-F-LE-4R"
    print(f"   cap run-out of 3D-F-LE-4R: {runout} (webs, Z and R kept)")
    fr = {n: resolve_system_type(f"{n}A-CN-LC-0R").web_chord_fracs for n in range(1, 6)}
    assert fr[1] == (0.5,) and fr[2] == (0.15, 0.5)
    for n in range(3, 6):
        assert np.allclose(fr[n], [k / (n + 1) for k in range(1, n + 1)])
    print("3. Default web stations: " + "; ".join(f"W={n}: {', '.join(f'{v:.3f}' for v in fr[n])}" for n in fr))


# --- 4. wall fields ---------------------------------------------------------------------

def edge_probes(verts):
    """Edge mid-points and unit outward edge normals: along these the distance is exact."""
    a, b = verts, np.roll(verts, -1, axis=0)
    t = b - a
    n = np.column_stack([t[:, 1], -t[:, 0]])
    n /= np.linalg.norm(n, axis=1)[:, None]
    area = 0.5 * np.sum(a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1])
    return 0.5 * (a + b), (n if area > 0 else -n)


def check_walls():
    af, sec = section("2D-CN-LC-0R", twist=0.0)
    v, n = edge_probes(af.vertices)
    # stations away from the nose and tail, where the wall is thin against the local radius
    keep = (v[:, 0] > 0.1 * CHORD) & (v[:, 0] < 0.8 * CHORD)
    depths = np.linspace(-0.004, 0.024, 15)
    skin = comp(sec, "skin_suction_surface").expr
    errs = []
    for d in depths:
        p = v[keep] - d * n[keep]
        up = p[:, 1] > np.interp(p[:, 0], *af.camber_line(400))
        phi = eval_expr(skin, p[up, 0], p[up, 1])
        errs.append(np.abs(phi - np.maximum(-d, d - SKIN)).max())
    err_skin = float(max(errs))
    cap = comp(sec, "spar_cap_suction_surface")
    p_mid = v[keep] - (SKIN + 0.5 * CAP) * n[keep]
    sel = (p_mid[:, 1] > 0.0) & (p_mid[:, 0] > cap.x_start + 0.05) & (p_mid[:, 0] < cap.x_end - 0.05)
    err_cap = 0.0
    for d in np.linspace(SKIN - 0.004, SKIN + CAP + 0.004, 15):
        p = v[keep][sel] - d * n[keep][sel]
        phi = eval_expr(cap.expr, p[:, 0], p[:, 1])
        err_cap = max(err_cap, float(np.abs(phi - np.maximum(SKIN - d, d - SKIN - CAP)).max()))
    web = comp(sec, "shear_web_0")
    a, b = np.array([web.x_top, web.y_top]), np.array([web.x_bot, web.y_bot])
    u = (a - b) / np.linalg.norm(a - b)
    w = np.array([-u[1], u[0]])
    s = np.linspace(-0.02, 0.02, 41)
    err_web = 0.0
    for frac in (0.25, 0.5, 0.75):
        c = b + frac * (a - b)
        p = c + s[:, None] * w
        phi = eval_expr(web.expr, p[:, 0], p[:, 1])
        err_web = max(err_web, float(np.abs(phi - (np.abs(s) - 0.5 * WEB)).max()))
    phi_i = af(np.array([a[0], b[0]]), np.array([a[1], b[1]])) + SKIN
    err_feet = float(np.abs(phi_i).max())
    print(f"\n4. Wall fields on a NACA 2418, 2 m chord: skin slab against max(-d, d - t_s) {err_skin:.1e}; "
          f"cap slab {err_cap:.1e}; web strip against |s| - t_w/2 {err_web:.1e}")
    print(f"   web feet on the inner mould line: |phi_A + t_s| <= {err_feet:.1e}")
    assert err_skin < 1e-12 and err_cap < 1e-12 and err_web < 1e-12 and err_feet < 1e-6
    return err_skin, err_cap, err_web


# --- 5. web axes ------------------------------------------------------------------------

def check_web_axes():
    worst_f = worst_cn = 0.0
    for deg in (-10.0, 0.0, 8.0, 15.0, 30.0):
        th = np.radians(deg)
        for key, worst in (("2A-F-LC-0R", "f"), ("2A-CN-LC-0R", "cn")):
            _, sec = section(key, twist=th)
            for i in range(2):
                web = comp(sec, f"shear_web_{i}")
                d = np.array([web.x_top - web.x_bot, web.y_top - web.y_bot])
                ang = np.arctan2(d[0], d[1])  # angle of the axis from the chord normal
                if worst == "f":
                    worst_f = max(worst_f, abs(ang - th))
                else:
                    worst_cn = max(worst_cn, abs(ang))
    print(f"\n5. Web axes for twists -10 to 30 deg: flapwise axis at the twist angle from the chord normal "
          f"to {worst_f:.1e} rad; chord-normal axis vertical to {worst_cn:.1e} rad")
    assert worst_f < 1e-9 and worst_cn < 1e-12


# --- 6. cap bands -----------------------------------------------------------------------

def check_cap_bands():
    rows = []
    af, sec = section("2D-CN-LC-0R", twist=0.0)
    c = comp(sec, "spar_cap_suction_surface")
    rows.append(("D", c.x_start, c.x_end, 0.15 * CHORD, 0.5 * CHORD))
    af, sec = section("3C-CN-LC-0R", twist=0.0)
    caps = [w for k, w in sec.spar_cap_wall_recipes.items() if k[-1].isdigit()]
    xc = sorted({round(0.5 * (cc.x_start + cc.x_end), 9) for cc in caps})
    hw = sorted({round(0.5 * (cc.x_end - cc.x_start), 9) for cc in caps})
    rows.append(("C", min(xc), max(xc), 0.25 * CHORD, 0.75 * CHORD))
    af, sec = section("1B-CN-LC-0R", twist=0.0)
    cb = comp(sec, "spar_cap_suction_surface")
    x_t = max_thickness_chord_x(af)
    x = np.linspace(0.0, 1.0, 20001)
    foil = Airfoil.naca4("2418")
    x_ref = CHORD * x[np.argmax(foil.surface("upper", x) - foil.surface("lower", x))]
    print("\n6. Cap bands:")
    print(f"   D (2 webs): x from {rows[0][1]:.4f} to {rows[0][2]:.4f} m = the outer webs {rows[0][3]:.4f}, {rows[0][4]:.4f}")
    print(f"   C (3 webs): pad centres {', '.join(f'{v:.4f}' for v in xc)} m, half-width {hw[0]:.4f} m = 0.04 c")
    print(f"   B (1 web, max-thickness anchor): band {cb.x_start:.4f} to {cb.x_end:.4f} m, centre "
          f"{0.5 * (cb.x_start + cb.x_end):.4f}; sampled max-thickness station {x_t:.4f}, "
          f"fine-grid vertical thickness peak {x_ref:.4f}")
    assert abs(rows[0][1] - rows[0][3]) < 1e-12 and abs(rows[0][2] - rows[0][4]) < 1e-12
    assert np.allclose(xc, [0.25 * CHORD, 0.5 * CHORD, 0.75 * CHORD], atol=1e-9) and abs(hw[0] - 0.04 * CHORD) < 1e-9
    assert abs(0.5 * (cb.x_end - cb.x_start) - 0.05 * CHORD) < 1e-9
    return x_t, x_ref


# --- 7. anchoring along the span --------------------------------------------------------

Z = np.linspace(0.0, 10.0, 101)
C0, C1, P = 2.0, 0.8, 0.3


def chord_at(z):
    return C0 + (C1 - C0) * z / Z[-1]


class _St:
    def __init__(self, z):
        self.z, self.chord = float(z), float(chord_at(z))


def check_anchoring():
    fr = (0.15, 0.5)
    stations = [_St(z) for z in Z]
    lc = np.array([[f * chord_at(z) for f in fr] for z in Z])                       # from the LE
    le = np.array([le_web_offsets_m(fr, C0) for _ in Z])
    pa_x = pa_web_x_from_pitch(stations, lambda z: P, 2)
    pa = np.array([[chord_fraction_from_pitch(x, chord_at(z), P) * chord_at(z) for x in pa_x] for z in Z])
    x_le0, x_te0 = blade_x_le_te(C0, P)
    lattice = [x_le0 + k * (x_te0 - x_le0) / 3 for k in (1, 2)]
    pa4 = pa_web_x_from_pitch(stations, lambda z: P, 4)
    dx = (x_te0 - x_le0) / 4
    err_lc = float(np.abs(lc / chord_at(Z)[:, None] - np.array(fr)).max())
    err_le = float(np.abs(le - np.array(fr) * C0).max())
    err_pa = float(np.abs((pa - P * chord_at(Z)[:, None]) - np.array(pa_x)).max())
    print("\n7. Web stations along a 2.0 -> 0.8 m taper, pitch axis at 0.3 c:")
    print(f"   LC: x/c held at {fr} to {err_lc:.1e}; LE: x held at {tuple(f * C0 for f in fr)} m to {err_le:.1e}; "
          f"PA: x from the pitch axis held at ({pa_x[0]:.4f}, {pa_x[1]:.4f}) m to {err_pa:.1e}")
    print(f"   PA lattice, 2 webs: {np.round(lattice, 6).tolist()} m (thirds of the root band); "
          f"4 webs: {np.round(pa4, 6).tolist()} m, the fourth at the root TE {x_te0:.4f} m")
    assert np.allclose(pa_x, lattice) and abs(pa4[3] - x_te0) < 1e-12 and np.allclose(np.diff(pa4[:3]), dx)
    assert err_lc < 1e-15 and err_le < 1e-15 and err_pa < 1e-12
    return lc, le, pa


# --- figure -----------------------------------------------------------------------------

def figure(span):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1.45, 1.0], "wspace": 0.25})
    ax = axes[0]
    af, sec = section("2D-F-LC-0R", twist=TWIST)
    xs, ys = np.meshgrid(np.linspace(-0.1, 2.1, 1101), np.linspace(-0.3, 0.35, 326))
    colours = {"skin": NAVY, "spar_cap": RED, "shear_web": TEAL}
    for k, c in sec.components_unrotated.items():
        name = str(getattr(k, "value", k))
        col = next((v for kk, v in colours.items() if name.startswith(kk)), None)
        if col is None:
            continue
        phi = eval_expr(c.expr, xs, ys)
        ax.contourf(xs, ys, phi, levels=[-1.0, 0.0], colors=[col], alpha=0.85)
    v = np.vstack([af.vertices, af.vertices[:1]])
    ax.plot(v[:, 0], v[:, 1], color=REFERENCE, lw=0.6)
    xc, yc = af.camber_line(200)
    ax.plot(xc, yc, color="0.5", lw=0.7, ls="--", label="camber line: suction / pressure split (R5.6)")
    for i in range(2):
        w = comp(sec, f"shear_web_{i}")
        ax.plot([w.x_bot, w.x_top], [w.y_bot, w.y_top], color="k", lw=0.6, ls=":")
        ax.annotate(f"$x_{i + 1}$", (0.5 * (w.x_bot + w.x_top), w.y_bot - 0.03), ha="center", fontsize=8)
    ax.plot([], [], color=NAVY, lw=4, label="skin: slab of depth $t_s$ (R5.5)")
    ax.plot([], [], color=RED, lw=4, label="caps: slab of $\\phi_I$, band $[x_1, x_2]$ (R5.7)")
    ax.plot([], [], color=TEAL, lw=4, label="webs: capsule $\\cap$ cavity (R5.9), flapwise (R5.11)")
    ax.set_aspect("equal")
    ax.set_xlim(-0.1, 2.1)
    ax.set_ylim(-0.3, 0.35)
    ax.set_xlabel("$x$ (m), chord frame")
    ax.set_ylabel("$y$ (m)")
    ax.set_title("2D-F on a NACA 2418, 2 m chord, 15° twist, chord frame", fontsize=9)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False, fontsize=7.5)

    ax = axes[1]
    lc, le, pa = span
    x_le, x_te = -P * chord_at(Z), (1 - P) * chord_at(Z)
    ax.fill_betweenx(Z, x_le, x_te, color="0.94")
    ax.plot(x_le, Z, color=REFERENCE, lw=0.9)
    ax.plot(x_te, Z, color=REFERENCE, lw=0.9)
    ax.axvline(0.0, color="0.5", lw=0.6, ls="-.")
    for j in range(2):
        ax.plot(lc[:, j] + x_le, Z, color=SERIES[0], lw=1.3, label="LC (R5.16)" if j == 0 else None)
        ax.plot(le[:, j] + x_le, Z, color=SERIES[1], lw=1.3, ls="--", label="LE (R5.17)" if j == 0 else None)
        ax.plot(pa[:, j] + x_le, Z, color=SERIES[2], lw=1.3, ls=":", label="PA (R5.18)" if j == 0 else None)
    ax.set_xlabel("$x$ from the pitch axis (m)")
    ax.set_ylabel("$z$ (m)")
    ax.set_title("Two webs on a 2.0 → 0.8 m taper", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    note(fig, "Checks: catalogue counts and ladders; skin, cap and web fields exact; web axes; "
              "cap bands B, C, D; LC, LE and PA stations.")
    fig.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.12)
    save(fig, __file__, "layout_construction")


if __name__ == "__main__":
    check_catalogue()
    check_walls()
    check_web_axes()
    check_cap_bands()
    span = check_anchoring()
    figure(span)
