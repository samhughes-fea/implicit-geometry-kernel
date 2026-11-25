"""Feasibility, drops and failure modes.

Draws Figure R10.1 of the R manual (docs/r/R10_feasibility.md) and checks:

1. the wall thickness floors, (R10.1), and the resolve/assemble skin mismatch (SG-040);
2. the intercept probe grid and the cavities it can miss, (R10.2)-(R10.3) (SG-038);
3. the accuracy of a refined crossing, (R10.4);
4. the three drop reasons in their order, (R10.5);
5. the step-down of the layout after a drop, (R10.6);
6. the cap allowance and the symmetric pad extent, (R10.7) and (R10.9);
7. the clip sample count, (R10.8);
8. the two height scalings and the pad run-out, (R10.10)-(R10.12);
9. the aft half-depth rule of family D, (R10.13);
10. the run-out gate against the clamp, (R10.14) (SG-039, SG-041);
11. the failure modes that raise.

    python examples/section/06_feasibility_checks/feasibility_checks.py
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import NAVY, RED, REFERENCE, SERIES, TEAL, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, resolve_system_type  # noqa: E402
from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry  # noqa: E402
from implicit_geometry_kernel.section.assembly.cap_feasibility import (  # noqa: E402
    MIN_CAP_VOID_ABS_M,
    MIN_CAP_VOID_FRAC,
    cap_feasibility_clamp,
    default_camber_clip_n_samples,
)
from implicit_geometry_kernel.section.assembly.web_layout import (  # noqa: E402
    inner_y_at_x,
    inner_y_at_xs,
    resolve_web_anchors,
)
from implicit_geometry_kernel.section.layout.section_axes import max_thickness_chord_x  # noqa: E402
from implicit_geometry_kernel.section.layout.spar_cap_runout import max_spar_from_cavity  # noqa: E402
from implicit_geometry_kernel.section.layout.thickness_limits import (  # noqa: E402
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M as T_SKIN_MIN,
    MIN_REALISTIC_SPAR_WALL_THICKNESS_M as T_CAP_MIN,
    MIN_REALISTIC_WEB_WALL_THICKNESS_M as T_WEB_MIN,
    clamp_skin_thickness_m,
    clamp_spar_thickness_m,
    clamp_web_thickness_m,
)
from implicit_geometry_kernel.section.spec import SectionGeometry, SectionSpec  # noqa: E402
from implicit_geometry_kernel.section.spec.layout_resolution import resolve_station_layout  # noqa: E402

CHORD, SKIN = 1.0, 0.006
AF = Airfoil.naca4("2412").to_sdf(chord=CHORD)
OML = np.asarray(AF.vertices, dtype=float)
THIN = Airfoil.naca4("0006").to_sdf(chord=CHORD)
N_PROBE = 240


def fine_depth(af, t, x, n=20_001):
    """Cavity depth y_top - y_bot from a dense probe (the reference for every check)."""
    try:
        yt, yb = inner_y_at_x(af, t, x, n_samples=n, y_search=(-0.2, 0.2))
    except InfeasibleSectionGeometry:
        return 0.0, None, None
    return yt - yb, yt, yb


def allowed(depth):
    """(R10.7): room for the two caps together."""
    return depth - np.maximum(MIN_CAP_VOID_ABS_M, MIN_CAP_VOID_FRAC * depth)


def section(t_skin, cap=0.0, web=0.001, key="2D-F-LC-0R", twist=0.0):
    spec = SectionSpec.from_kwargs(
        resolve_system_type(key),
        AF,
        skin_thickness=t_skin,
        spar_cap_thickness=cap,
        shear_web_thickness=web,
        twist_rad=twist,
        pivot_fraction=0.35,
    )
    return SectionGeometry(spec).execute().get_results()


# --- 1. floors -----------------------------------------------------------------------------


def check_floors():
    for f, floor in (
        (clamp_skin_thickness_m, T_SKIN_MIN),
        (clamp_web_thickness_m, T_WEB_MIN),
        (clamp_spar_thickness_m, T_CAP_MIN),
    ):
        assert f(0.0) == floor and f(0.001) == floor and f(0.010) == 0.010
    walls = section(0.002, cap=0.0, web=0.001).section.wall_recipes
    t_s = walls["skin_suction_surface"].thickness
    t_w = walls["shear_web_0"].thickness
    h_c = walls["spar_cap_suction_surface"].cap_height
    assert (t_s, t_w, h_c) == (T_SKIN_MIN, T_WEB_MIN, T_CAP_MIN)
    print(
        f"1. Floors (R10.1): asked skin 2 mm, web 1 mm, cap 0 mm; built {t_s * 1e3:g}, {t_w * 1e3:g}, {h_c * 1e3:g} mm"
    )
    rows = []
    for t in (0.0010, 0.0012, 0.0013):
        try:
            section(t)
            rows.append(f"t = {t * 1e3:g} mm builds")
            assert t >= T_SKIN_MIN / 2
        except InfeasibleSectionGeometry as exc:
            gap = float(str(exc).split(" lies ")[1].split(" m ")[0])
            assert t < T_SKIN_MIN / 2 and abs(gap - (T_SKIN_MIN / 2 - t)) < 1e-9
            rows.append(f"t = {t * 1e3:g} mm raises, anchor {gap * 1e3:.2f} mm outboard")
    print(
        f"   SG-040: resolve probes the unclamped skin, assembly builds at {T_SKIN_MIN * 1e3:g} mm, so a skin "
        f"under t_min/2 = {T_SKIN_MIN / 2 * 1e3:g} mm raises: " + "; ".join(rows)
    )


# --- 2-3. probe grid and crossings ------------------------------------------------------------


def probe_window(af, pad_scale):
    v = np.asarray(af.vertices, dtype=float)
    pad = max(0.05 * pad_scale, 0.01)
    return v[:, 1].min() - pad, v[:, 1].max() + pad


def check_probe_grid():
    x = 0.8
    v = np.asarray(THIN.vertices, dtype=float)
    thick = float(v[:, 1].max() - v[:, 1].min())
    d_outer = fine_depth(THIN, 0.0, x)[0]
    win_chord = probe_window(THIN, CHORD)
    win_thick = probe_window(THIN, thick)
    dy_chord = (win_chord[1] - win_chord[0]) / (N_PROBE - 1)
    dy_thick = (win_thick[1] - win_thick[0]) / (N_PROBE - 1)
    rows = []
    for target in (0.2e-3, 0.4e-3, 0.6e-3, 0.8e-3, 1.0e-3, 2.0e-3):
        lo, hi = 0.0, d_outer / 2
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if fine_depth(THIN, mid, x)[0] > target else (lo, mid)
        d, yt, yb = fine_depth(THIN, lo, x)
        kept, dropped = resolve_web_anchors(THIN, lo, [x], return_dropped=True)
        ys = np.linspace(*win_chord, N_PROBE)
        predicted = bool(np.any((ys > yb) & (ys < yt)))
        assert predicted == bool(kept)
        if d > dy_chord:
            assert kept
        ys2 = np.linspace(*win_thick, N_PROBE)
        thick_pad = bool(np.any((ys2 > yb) & (ys2 < yt)))
        rows.append((d, bool(kept), thick_pad))
    print(
        f"2. Probe grid (R10.2)-(R10.3), NACA 0006, x = 0.8c: n = {N_PROBE} samples, dy = {dy_chord * 1e3:.3f} mm "
        f"with the chord pad of resolve_web_anchors, {dy_thick * 1e3:.3f} mm with the thickness pad; "
        + ", ".join(f"d = {d * 1e3:.1f} mm {'kept' if k else 'dropped'}" for d, k, _ in rows)
        + "; every outcome matches 'a sample lies in the cavity'"
    )
    missed = [d for d, k, tp in rows if not k and tp]
    print(
        f"   SG-038: {len(missed)} real cavities ({', '.join(f'{d * 1e3:.1f}' for d in missed)} mm) are dropped as "
        "no_inner_skin_intercept; the thickness pad of _default_y_search would see them"
    )
    return rows, dy_chord, dy_thick


def check_crossings():
    win = probe_window(AF, CHORD)
    dy = (win[1] - win[0]) / (N_PROBE - 1)
    worst, worst0 = 0.0, 0.0
    xs = np.linspace(0.05, 0.9, 18)
    coarse0 = inner_y_at_xs(AF, SKIN, xs, bisection_iters=0)
    for x, c0 in zip(xs, coarse0):
        _, yt, yb = fine_depth(AF, SKIN, x)
        yt10, yb10 = inner_y_at_x(AF, SKIN, x)
        worst = max(worst, abs(yt10 - yt), abs(yb10 - yb))
        worst0 = max(worst0, abs(c0[0] - yt), abs(c0[1] - yb))
    assert worst < dy / 2**10 and worst0 < dy
    print(
        f"3. Crossings (R10.4), NACA 2412, 6 mm skin, 18 columns: 10 bisections + secant {worst:.1e} m "
        f"(bound dy/2^10 = {dy / 2**10:.1e}); secant only (cap clamp) {worst0:.1e} m (bound dy = {dy:.1e}); "
        "phi is linear in y across each crossing, so the secant is exact to round-off"
    )


# --- 4-5. drops and step-down -----------------------------------------------------------------


def check_drops():
    xs = [0.3, 0.31, 0.7, 0.95, 0.985]
    d_min = 0.010
    kept, dropped = resolve_web_anchors(AF, SKIN, xs, min_cavity_depth_m=d_min, return_dropped=True)
    reasons = {d.x_m: d.reason for d in dropped}
    assert [k[0] for k in kept] == [0.3, 0.7]
    assert reasons == {0.31: "merged_into_previous_web", 0.95: "cavity_too_shallow", 0.985: "no_inner_skin_intercept"}
    assert fine_depth(AF, SKIN, 0.985)[0] == 0.0 and 0.0 < fine_depth(AF, SKIN, 0.95)[0] < d_min
    assert abs(0.31 - 0.3) < 0.02 * CHORD
    print(
        f"4. Drops (R10.5), d_min = {d_min * 1e3:g} mm: kept {[k[0] for k in kept]}; "
        + ", ".join(f"{x:g} {r}" for x, r in reasons.items())
        + f" (depth at 0.95c = {fine_depth(AF, SKIN, 0.95)[0] * 1e3:.2f} mm)"
    )


def check_step_down():
    rows = []
    for key, xs, want in (("2D-F-LC-0R", [0.3, 0.985], "1C-F-LC-0R"), ("3D-F-LC-0R", [0.2, 0.5, 0.985], "2D-F-LC-0R")):
        r = resolve_station_layout(resolve_system_type(key), AF, skin_thickness=SKIN, web_x_positions_m=xs)
        assert r.layout.key == want and r.stepped_down
        rows.append(f"{key} -> {r.layout.key} ({r.layout.cap_placement.kind})")
    try:
        resolve_station_layout(
            resolve_system_type("2D-F-LC-0R"), AF, skin_thickness=SKIN, web_x_positions_m=[0.985, 0.99]
        )
        raise AssertionError("expected a raise")
    except InfeasibleSectionGeometry:
        rows.append("2D-F with no web anchoring raises")
    print("5. Step-down (R10.6): " + "; ".join(rows))


# --- 6-9. the cap clamp -----------------------------------------------------------------------


def band_depth(x_lo, x_hi):
    """The clamp's own depth samples, (R10.8) count, secant crossings."""
    n = default_camber_clip_n_samples(x_lo, x_hi)
    xs = np.linspace(x_lo, x_hi, n)
    d = np.array([-np.inf if p is None else p[0] - p[1] for p in inner_y_at_xs(AF, SKIN, xs, bisection_iters=0)])
    return xs, d


def pad(x_c, h, dhw=0.04):
    return cap_feasibility_clamp(AF, SKIN, x_c - dhw, x_c + dhw, h, h, y_search=None, phi_inner_ref=None, center_x=x_c)


def check_pads():
    d_cross = MIN_CAP_VOID_ABS_M / MIN_CAP_VOID_FRAC
    assert abs(allowed(np.array([d_cross]))[0] - 0.85 * d_cross) < 1e-15
    out = {}
    for x_c, h in ((0.3, 0.02), (0.8, 0.012), (0.84, 0.01), (0.9, 0.01)):
        lo, hi, hu, hl, info = pad(x_c, h)
        xs, d = band_depth(x_c - 0.04, x_c + 0.04)
        a = allowed(d)
        inside = (xs >= lo - 1e-15) & (xs <= hi + 1e-15)
        assert abs((lo + hi) / 2 - x_c) < 1e-12
        assert np.all(hu + hl <= a[inside] + 1e-12)
        out[x_c] = (lo, hi, hu, hl, info, xs, a)
    for x_c in (0.3, 0.8):
        assert abs(out[x_c][0] - (x_c - 0.04)) < 1e-12 and abs(out[x_c][1] - (x_c + 0.04)) < 1e-12
        assert not out[x_c][4]["scaled_heights"]
    lo, hi, hu, _, info, xs, a = out[0.84]
    beyond = (xs > hi + 1e-12) & (xs < hi + (xs[1] - xs[0]) + 1e-12)
    assert info["clamped_extent"] and not info["scaled_heights"] and np.all(2 * hu > a[beyond])
    lo9, hi9, hu9, hl9, info9, xs9, a9 = out[0.9]
    assert info9["scaled_heights"] and abs((hu9 + hl9) - float(np.min(a9))) < 1e-12
    assert abs(lo9 - 0.86) < 1e-12 and abs(hi9 - 0.94) < 1e-12
    print(
        f"6. Allowance (R10.7): void max(1 mm, 0.15 d), crossover d = {d_cross * 1e3:.2f} mm. Pads of half-width 0.04c "
        f"(R10.9): 0.3c and 0.8c full; 0.84c trimmed symmetrically to [{lo:.4f}, {hi:.4f}] (next sample infeasible); "
        f"0.9c scaled to min a over the band, h = {hu9 * 1e3:.3f} mm each (R10.11)"
    )
    return out


def check_samples():
    rows = []
    for band in (0.08, 0.5, 1.0, 2.0):
        n = default_camber_clip_n_samples(0.0, band)
        c = 18.56e-6 * 33**2 / 0.50
        assert n == max(33, math.ceil(math.sqrt(c * band / 1e-6)))
        assert c * band / n**2 <= 1e-6 < c * band / (n - 1) ** 2 or n == 33
        rows.append(f"B = {band:g} m: n = {n}")
    print("7. Clip samples (R10.8), C B / n^2 <= 1 um: " + ", ".join(rows))


def check_scaling():
    x_lo, x_hi = 0.2, 0.6
    xs, d = band_depth(x_lo, x_hi)
    a = allowed(d)
    i0, i1 = len(a) // 4, (3 * len(a)) // 4 + 1
    target = max(float(np.min(a[i0:i1])), 0.5 * float(np.max(a)))
    lo, hi, hu, hl, info = cap_feasibility_clamp(AF, SKIN, x_lo, x_hi, 0.06, 0.06, y_search=None, phi_inner_ref=None)
    assert info["scaled_heights"] and abs(hu + hl - target) < 1e-12 and abs(info["scale"] - target / 0.12) < 1e-12
    print(
        f"8. Global scale (R10.10): 2 x 60 mm caps in [0.2c, 0.6c] exceed max a = {np.max(a) * 1e3:.2f} mm; "
        f"scaled by {info['scale']:.4f} to h = {hu * 1e3:.3f} mm each over [{lo:.4f}, {hi:.4f}]"
    )
    xs_r, d_r = band_depth(0.95 - 0.04, 0.95 + 0.04)
    a_r = allowed(d_r)
    lo, hi, hu, hl, info = pad(0.95, 0.01)
    i_c = int(np.argmin(np.abs(xs_r - 0.95)))
    assert np.min(a_r[np.isfinite(a_r)]) <= 0.0 < a_r[i_c] < 0.02
    assert info.get("pad_run_out") and (hu, hl) == (0.0, 0.0) and lo == hi == 0.95
    print(
        f"   Run-out (R10.12): pad at 0.95c, min a over the band {np.min(a_r[np.isfinite(a_r)]) * 1e3:.2f} mm <= 0 "
        f"so no rescale; centre a = {a_r[i_c] * 1e3:.2f} mm < 2h = 20 mm: heights set to 0"
    )


def check_family_d():
    x_lo, x_hi = 0.2, 0.6
    xs, d = band_depth(x_lo, x_hi)
    lo, hi, hu, hl, info = cap_feasibility_clamp(
        AF, SKIN, x_lo, x_hi, 0.05, 0.03, y_search=None, phi_inner_ref=None, aft_half_depth_rule=True
    )
    h_allow = 0.5 * d[-1]
    assert (lo, hi) == (x_lo, x_hi) and hu == h_allow and hl == 0.03 and info["scaled_heights"]
    print(
        f"9. Family D (R10.13): h_allow = d(x_aft)/2 = {h_allow * 1e3:.2f} mm; suction 50 -> {hu * 1e3:.2f} mm, "
        f"pressure 30 mm kept; extent stays [{lo:g}, {hi:g}]"
    )
    return h_allow


# --- 10-11. run-out gate and raises -----------------------------------------------------------


def check_runout_gate():
    t, g = 0.004, 0.002
    rows, curves = [], []
    for d in (0.005, 0.010, 0.020, 0.100):
        gate = max_spar_from_cavity(d + t, t, gap_min=g)
        clamp = 0.5 * float(allowed(np.array([d]))[0])
        assert abs(gate - 0.5 * (d - g)) < 1e-15
        rows.append(f"d = {d * 1e3:g} mm: gate {gate * 1e3:.2f}, clamp {clamp * 1e3:.2f}")
        curves.append((d, gate, clamp))
    x_mt = max_thickness_chord_x(AF)
    print(
        "10. Run-out gate (R10.14) with h_local the mid-to-mid depth d + t_s, per-cap ceiling (mm): "
        + "; ".join(rows)
        + ". SG-039: the gate is stricter below d = g/0.15 = 13.3 mm and more lenient above it"
    )
    print(f"    SG-041: the gate centres a fixed band at 0.5c; the build centres it at max thickness, x = {x_mt:.3f}c")
    assert abs(x_mt - 0.5) > 0.1
    assert abs(0.5 * (0.002 / 0.15 - 0.002) - 0.5 * float(allowed(np.array([0.002 / 0.15]))[0])) < 1e-15


def check_raises():
    rows = []
    try:
        section(SKIN, cap=0.01, web=0.008, twist=np.pi / 2)
        raise AssertionError("expected a raise")
    except InfeasibleSectionGeometry as exc:
        assert "SG-033" in str(exc)
        rows.append("flapwise web at 90 deg twist (SG-033)")
    try:
        resolve_web_anchors(AF, SKIN, [0.2, 0.985], return_dropped=True)
        from implicit_geometry_kernel.section.assembly.plan import _anchor_webs

        _anchor_webs(
            AF, SKIN, [0.2, 0.985], explicit_endpoints=None, web_y_coords=None, y_search=None, phi_inner_ref=None
        )
        raise AssertionError("expected a raise")
    except InfeasibleSectionGeometry as exc:
        assert "resolve_station_layout" in str(exc)
        rows.append("assembly refuses a web set that would drop")
    print("11. Raises: " + "; ".join(rows) + "; no-web stations raise in check 5")


# --- figure -----------------------------------------------------------------------------------


def figure(pads, probe_rows, dy_chord, dy_thick, h_allow):
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), gridspec_kw={"width_ratios": [1.5, 1.0, 1.0], "wspace": 0.34})
    ax = axes[0]
    xs = np.linspace(0.01, 0.99, 300)
    d = np.array([fine_depth(AF, SKIN, x, n=20_001)[0] for x in xs])
    ax.plot(xs, d * 1e3, color=NAVY, lw=1.3, label="cavity depth $d(x)$")
    ax.plot(
        xs,
        np.where(d > 0, allowed(d), np.nan) * 1e3,
        color=REFERENCE,
        lw=0.9,
        ls="--",
        label="allowance $a(x)$ (R10.7)",
    )
    for i, (x_c, (lo, hi, hu, hl, info, _, _)) in enumerate(pads.items()):
        ax.plot([x_c - 0.04, x_c + 0.04], [(hu + hl) * 1e3] * 2, color=SERIES[i], lw=0.8, ls=":")
        ax.plot(
            [lo, hi],
            [(hu + hl) * 1e3] * 2,
            color=SERIES[i],
            lw=3.0,
            solid_capstyle="butt",
            label=f"pad at {x_c:g}c, $h_u + h_l$ = {(hu + hl) * 1e3:.1f} mm",
        )
    ax.axhline(MIN_CAP_VOID_ABS_M / MIN_CAP_VOID_FRAC * 1e3, color=TEAL, lw=0.6)
    ax.annotate("void floor crossover 6.67 mm", (0.03, 7.3), fontsize=7, color=TEAL)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 165)
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("height (mm)")
    ax.set_title("(a) Cap allowance and pad extents, NACA 2412, 6 mm skin")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")

    ax = axes[1]
    ds = np.array([r[0] for r in probe_rows]) * 1e3
    kept = np.array([r[1] for r in probe_rows])
    ax.plot(ds[kept], np.ones(kept.sum()), "o", color=TEAL, label="kept")
    ax.plot(ds[~kept], np.ones((~kept).sum()), "x", color=RED, ms=7, label="dropped: no intercept")
    ax.axvline(dy_chord * 1e3, color=RED, lw=0.9, ls="--", label=f"$\\Delta y$ chord pad {dy_chord * 1e3:.2f} mm")
    ax.axvline(dy_thick * 1e3, color=TEAL, lw=0.9, ls="--", label=f"$\\Delta y$ thickness pad {dy_thick * 1e3:.2f} mm")
    ax.set_yticks([])
    ax.set_ylim(0.9, 1.25)
    ax.set_xlabel("cavity depth $d$ at 0.8c (mm)")
    ax.set_title("(b) Probe grid misses shallow cavities (SG-038)")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")

    ax = axes[2]
    dd = np.linspace(0.0, 0.03, 300)
    ax.plot(dd * 1e3, 0.5 * np.maximum(allowed(dd), 0) * 1e3, color=REFERENCE, lw=1.0, label="clamp $a/2$ (R10.7)")
    ax.plot(
        dd * 1e3,
        np.maximum(0.5 * (dd - 0.002), 0) * 1e3,
        color=RED,
        lw=1.0,
        ls="--",
        label="run-out gate $(d - g)/2$ (R10.14)",
    )
    ax.axhline(T_CAP_MIN * 1e3, color=TEAL, lw=0.6, label="cap floor 5 mm")
    ax.set_xlabel("cavity depth $d$ (mm)")
    ax.set_ylabel("per-cap ceiling (mm)")
    ax.set_title("(c) Two void rules (SG-039)")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    note(
        fig,
        "Checks: floors; probe grid and crossings; drop reasons; step-down; allowance and symmetric extent; "
        "sample count; scalings and run-out; family D; run-out gate; raises.",
    )
    fig.subplots_adjust(left=0.06, right=0.99, top=0.9, bottom=0.2)
    save(fig, __file__, "feasibility_checks")


if __name__ == "__main__":
    check_floors()
    probe_rows, dy_chord, dy_thick = check_probe_grid()
    check_crossings()
    check_drops()
    check_step_down()
    pads = check_pads()
    check_samples()
    check_scaling()
    h_allow = check_family_d()
    check_runout_gate()
    check_raises()
    figure(pads, probe_rows, dy_chord, dy_thick, h_allow)
