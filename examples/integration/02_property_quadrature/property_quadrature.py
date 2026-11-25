"""Section properties by Green's theorem and by quadtree quadrature.

Draws Figure R7.1 of the R manual (docs/r/R7_section_properties.md) and checks:

1. the shoelace sums are exact on polygons, including the product moment, (R7.1)-(R7.5);
2. holes subtract and the parts of a multipolygon add about a common origin, (R7.6);
3. the parallel-axis reduction to the centroid, (R7.7);
4. a sampled ring converges to a circle as 1/N^2, (R7.8);
5. the rectangle moments of an interior leaf, (R7.9);
6. the cut polygon of a boundary leaf and the leaf classes, (R7.10)-(R7.12);
7. quadtree area and second moment converge as h^2 on a circle, (R7.13);
8. the hybrid backend: wall solids by shoelace on their polygons, (R7.14).

    python examples/integration/02_property_quadrature/property_quadrature.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import REFERENCE, SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel.geos.geom import MultiPolygon, Polygon  # noqa: E402
from implicit_geometry_kernel.integration import (  # noqa: E402
    QuadtreeSpec,
    build_shared_quadtree,
    section_properties_from_quadtree,
)
from implicit_geometry_kernel.integration.quadtree.moments import (  # noqa: E402
    _inside_polygon_from_cell,
    _polygon_moments,
    _rectangle_moments,
)
from implicit_geometry_kernel.integration.quadtree.tree import QuadCell, leaf_bounds  # noqa: E402
from implicit_geometry_kernel.integration.shoelace import (  # noqa: E402
    section_properties_S,
    section_properties_Sprime,
)

R = 0.3


def circle(cx, cy, r):
    return lambda p: np.hypot(*(np.asarray(p, float).reshape(-1, 2) - [cx, cy]).T) - r


def rect(x0, x1, y0, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], float)


# --- 1-3. shoelace ------------------------------------------------------------------------


def check_shoelace():
    x0, x1, y0, y1 = 0.4, 1.3, -0.2, 0.5
    a, (cx, cy), ixx, iyy, ixy = section_properties_S(Polygon(rect(x0, x1, y0, y1)))
    exact = (
        (x1 - x0) * (y1 - y0),
        0.5 * (x0 + x1),
        0.5 * (y0 + y1),
        (x1 - x0) * (y1**3 - y0**3) / 3,
        (y1 - y0) * (x1**3 - x0**3) / 3,
        (x1**2 - x0**2) * (y1**2 - y0**2) / 4,
    )
    got = (a, cx, cy, ixx, iyy, ixy)
    err_rect = max(abs(g - e) for g, e in zip(got, exact))
    tri = np.array([[0.2, 0.1], [1.1, 0.3], [0.5, 0.9]])
    at, (tcx, tcy), tixx, tiyy, tixy = section_properties_S(Polygon(tri))
    # Triangle: I about the origin = A/6 (sum of squares + sum of pairwise products) of the vertices.
    xs, ys = tri[:, 0], tri[:, 1]
    a_ref = 0.5 * abs((xs[1] - xs[0]) * (ys[2] - ys[0]) - (xs[2] - xs[0]) * (ys[1] - ys[0]))
    ixx_ref = a_ref / 6 * (np.sum(ys**2) + ys[0] * ys[1] + ys[1] * ys[2] + ys[2] * ys[0])
    ixy_ref = (
        a_ref
        / 12
        * (
            np.sum(xs * ys) * 2
            + xs[0] * ys[1]
            + xs[1] * ys[0]
            + xs[1] * ys[2]
            + xs[2] * ys[1]
            + xs[2] * ys[0]
            + xs[0] * ys[2]
        )
    )
    err_tri = max(abs(at - a_ref), abs(tcx - xs.mean()), abs(tcy - ys.mean()), abs(tixx - ixx_ref), abs(tixy - ixy_ref))
    assert err_rect < 1e-15 and err_tri < 1e-15
    print(
        f"1. Shoelace: rectangle off the origin (A, centroid, Ixx, Iyy, Ixy) to {err_rect:.0e}; "
        f"triangle against the vertex formulas to {err_tri:.0e}"
    )

    outer, hole = rect(0.0, 2.0, 0.0, 1.0), rect(0.5, 1.5, 0.25, 0.75)
    ph = section_properties_S(Polygon(outer, [hole[::-1]]))
    po, pi_ = section_properties_S(Polygon(outer)), section_properties_S(Polygon(hole))
    err_hole = max(abs(ph[0] - (po[0] - pi_[0])), abs(ph[2] - (po[2] - pi_[2])), abs(ph[4] - (po[4] - pi_[4])))
    pa, pb = Polygon(rect(0.0, 1.0, 0.0, 1.0)), Polygon(rect(3.0, 4.0, 2.0, 3.0))
    pm = section_properties_S(MultiPolygon([pa, pb]))
    sa, sb = section_properties_S(pa), section_properties_S(pb)
    err_multi = max(abs(pm[2] - (sa[2] + sb[2])), abs(pm[3] - (sa[3] + sb[3])), abs(pm[4] - (sa[4] + sb[4])))
    assert err_hole < 1e-14 and err_multi < 1e-13
    print(
        f"2. A hole subtracts to {err_hole:.0e}; two disjoint squares add about the origin to {err_multi:.0e} (SG-036)"
    )

    p = Polygon(tri)
    s, sp = section_properties_S(p), section_properties_Sprime(p)
    a, (cx, cy) = s[0], s[1]
    err_pa = max(abs(sp[2] - (s[2] - a * cy * cy)), abs(sp[4] - (s[4] - a * cx * cy)))
    centred = section_properties_S(Polygon(tri - [cx, cy]))
    err_pa2 = max(abs(sp[2] - centred[2]), abs(sp[3] - centred[3]), abs(sp[4] - centred[4]))
    assert err_pa2 < 1e-15
    print(
        f"3. Parallel axis: S' moments equal the moments of the polygon moved to its centroid to {err_pa2:.0e} "
        f"(R7.7, R1.24; {err_pa:.0e} against the formula)"
    )


def check_ring():
    rows = []
    for n in (16, 64, 256, 1024):
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        a = section_properties_S(Polygon(np.column_stack([R * np.cos(t), R * np.sin(t)])))[0]
        exact_n = 0.5 * n * R * R * np.sin(2 * np.pi / n)
        assert abs(a - exact_n) < 1e-15
        rows.append((n, (np.pi * R * R - a) / (np.pi * R * R)))
    ratio = rows[-2][1] / rows[-1][1]
    assert 15.9 < ratio < 16.1
    print(
        "4. Inscribed N-gon: area = (N/2) R^2 sin(2 pi/N) exactly; relative deficit "
        + ", ".join(f"N = {n} {e:.2e}" for n, e in rows)
        + f"; ratio per doubling squared {ratio:.2f} (R7.8)"
    )
    return rows


# --- 5-7. quadtree ------------------------------------------------------------------------


def check_rectangle_moments():
    x0, x1, y0, y1 = 1.0, 2.0, 1.0, 2.0
    m = _rectangle_moments(x0, x1, y0, y1)
    p = _polygon_moments(rect(x0, x1, y0, y1))
    exact_sxy = (x1**2 - x0**2) * (y1**2 - y0**2) / 4
    err = max(abs(m[i] - p[i]) for i in range(6))
    assert err < 1e-15 and abs(m[5] - exact_sxy) < 1e-15
    print(
        f"5. Leaf rectangle [1, 2]^2: area, first, second and product moments equal the shoelace values to "
        f"{err:.0e} (R7.9); integral of xy = {m[5]:.4f} against {exact_sxy:.4f} (IN-001 fixed)"
    )
    return m[5], exact_sxy


def check_cut_cells():
    cell = QuadCell(0.0, 1.0, 0.0, 1.0)
    cases = (  # straight boundaries: one, two and three inside corners, with exact areas
        ((1.0, 1.0), 0.8, 0.5 * 0.8**2),
        ((1.0, 0.0), 0.3, 0.3),
        ((1.0, 1.0), 1.5, 1.0 - 0.5 * 0.5**2),
    )
    worst = 0.0
    for (nx, ny), c, exact in cases:
        n = np.hypot(nx, ny)
        phi = (cell.corners() @ [nx, ny] - c) / n
        poly = _inside_polygon_from_cell(cell, phi)
        worst = max(worst, abs(abs(_polygon_moments(poly)[0]) - exact))
    assert worst < 1e-15
    print(
        f"6. Cut leaf: straight boundaries with one, two and three inside corners give the exact areas "
        f"0.32, 0.3, 0.875 to {worst:.0e} (R7.10)"
    )


def quadtree_props(sampler, h, box=(-0.5, 0.5, -0.5, 0.5)):
    spec = QuadtreeSpec(
        *box,
        max_depth=16,
        min_cell_size=h,
        refine_band=2 * h,
        max_cells=2_000_000,
        root_cell_size=0.125,
        refine_probes="full",
    )
    mesh = build_shared_quadtree(spec, {"a": sampler})
    return section_properties_from_quadtree(mesh, mesh.phi_by_label, "a"), mesh


def check_quadtree():
    hs = np.array([0.02, 0.01, 0.005, 0.0025])
    exact_a, exact_i = np.pi * R**2, np.pi * R**4 / 4
    ea, ei = [], []
    for h in hs:
        pr, mesh = quadtree_props(circle(0.0, 0.0, R), h)
        ea.append(abs(pr["area"] - exact_a) / exact_a)
        ei.append(abs(pr["Ixx"] - exact_i) / exact_i)
    slope_a = np.polyfit(np.log(hs), np.log(ea), 1)[0]
    slope_i = np.polyfit(np.log(hs), np.log(ei), 1)[0]
    assert slope_a > 1.8 and slope_i > 1.8 and ea[-1] < 1e-4
    _, mesh = quadtree_props(circle(0.0, 0.0, R), 0.01)
    b = leaf_bounds(mesh)
    phi = np.asarray(mesh.phi_by_label["a"])[:, :4]
    cut = (phi.min(axis=1) <= 0) & (phi.max(axis=1) > 0)
    sizes = b[cut, 1] - b[cut, 0]
    assert sizes.max() <= 0.01 + 1e-12
    print(
        "7. Quadtree on a circle R = 0.3: relative area error "
        + ", ".join(f"h = {h:g} {e:.1e}" for h, e in zip(hs, ea))
        + f"; slope {slope_a:.2f} (area), {slope_i:.2f} (Ixx) (R7.13); every cut leaf at h: "
        f"largest {sizes.max():.4f} m, {len(b)} leaves"
    )
    off, _ = quadtree_props(circle(0.2, 0.1, 0.2), 0.005, box=(-0.1, 0.5, -0.2, 0.4))
    assert abs(off["Ixy"]) < 1e-3 * off["Ixx"]
    print(f"   Ixy of a circle centred at (0.2, 0.1), exactly 0: {off['Ixy']:.1e} against Ixx {off['Ixx']:.2e}")
    return hs, np.array(ea), np.array(ei), off["Ixy"]


# --- 8. hybrid ----------------------------------------------------------------------------


def check_hybrid():
    from implicit_geometry_kernel import Airfoil, resolve_system_type
    from implicit_geometry_kernel.section.spec import SectionGeometry, SectionOptions, SectionSpec

    af = Airfoil.naca4("2418").to_sdf(chord=2.0)
    spec = SectionSpec.from_kwargs(
        resolve_system_type("2D-F-LC-0R"),
        af,
        skin_thickness=0.008,
        spar_cap_thickness=0.03,
        shear_web_thickness=0.015,
        twist_rad=0.0,
    )
    out = {}
    for backend in ("hybrid", "quadtree"):
        res = SectionGeometry(spec, SectionOptions(outputs=("properties",), props_backend=backend)).execute()
        out[backend] = res.get_results().properties_map
    rows = []
    for lbl in ("skin_suction_surface", "spar_cap_suction_surface", "shear_web_0"):
        a_h, a_q = out["hybrid"][lbl]["area"], out["quadtree"][lbl]["area"]
        rows.append((lbl, a_h, a_q, a_q / a_h - 1))
    worst = max(abs(r[3]) for r in rows)
    assert worst < 5e-3
    print(
        "8. Hybrid (shoelace on wall polygons) against quadtree on the wall fields, 2D-F on a 2 m NACA 2418: "
        + "; ".join(f"{r[0]} {r[1] * 1e4:.3f} cm^2 vs {r[2] * 1e4:.3f} ({r[3]:+.1e})" for r in rows)
    )


# --- figure ------------------------------------------------------------------------------


def figure(ring, hs, ea, ei):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.9), gridspec_kw={"wspace": 0.35})
    ax = axes[0]
    _, mesh = quadtree_props(circle(0.0, 0.0, R), 0.02)
    b = leaf_bounds(mesh)
    phi = np.asarray(mesh.phi_by_label["a"])[:, :4]
    interior = phi.max(axis=1) <= 0
    cut = (phi.min(axis=1) <= 0) & ~interior
    for k in range(len(b)):
        x0, x1, y0, y1 = b[k]
        if x1 < -0.02 or y1 < -0.02:
            continue
        col = SERIES[0] if interior[k] else (SERIES[1] if cut[k] else "white")
        ax.add_patch(
            plt.Rectangle(
                (x0, y0), x1 - x0, y1 - y0, facecolor=col, alpha=0.35 if col != "white" else 1, edgecolor="0.4", lw=0.3
            )
        )
        if cut[k]:
            poly = _inside_polygon_from_cell(QuadCell(x0, x1, y0, y1), phi[k])
            if poly is not None:
                ax.fill(poly[:, 0], poly[:, 1], color=SERIES[1], alpha=0.7, lw=0)
    t = np.linspace(0, np.pi / 2, 200)
    ax.plot(R * np.cos(t), R * np.sin(t), color=REFERENCE, lw=0.8)
    ax.set_xlim(-0.02, 0.4)
    ax.set_ylim(-0.02, 0.4)
    ax.set_aspect("equal")
    ax.plot([], [], "s", color=SERIES[0], alpha=0.5, label="interior leaf: rectangle (R7.9)")
    ax.plot([], [], "s", color=SERIES[1], alpha=0.8, label="cut leaf: inside polygon (R7.10)")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.set_title("Quadtree leaves on a circle, h = 0.02", fontsize=9)

    ax = axes[1]
    ax.loglog(hs, ea, "o-", color=SERIES[0], ms=4, label="area")
    ax.loglog(hs, ei, "s-", color=SERIES[2], ms=4, label="$I_{xx}$")
    ax.loglog(hs, ea[0] * (hs / hs[0]) ** 2, "--", color=REFERENCE, lw=0.8, label="slope 2 (R7.13)")
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel("smallest leaf $h$ (m)")
    ax.set_ylabel("relative error")
    ax.set_title("Quadtree quadrature converges as $h^2$", fontsize=9)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2]
    ns = np.array([r[0] for r in ring])
    es = np.array([r[1] for r in ring])
    ax.loglog(ns, es, "o-", color=SERIES[1], ms=4, label="shoelace on an inscribed $N$-gon")
    ax.loglog(ns, 2 * np.pi**2 / 3 / ns**2, "--", color=REFERENCE, lw=0.8, label="$2\\pi^2/(3N^2)$ (R7.8)")
    ax.set_xlabel("vertices $N$")
    ax.set_ylabel("relative area deficit")
    ax.set_title("Polygon sampling of a circle", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    note(
        fig,
        "Checks: shoelace exact; holes and parts; parallel axis; N-gon; leaf rectangles; "
        "cut leaves; h^2 convergence; hybrid against quadtree.",
    )
    fig.subplots_adjust(left=0.05, right=0.98, top=0.9, bottom=0.17)
    save(fig, __file__, "property_quadrature")


if __name__ == "__main__":
    check_shoelace()
    ring = check_ring()
    check_rectangle_moments()
    check_cut_cells()
    hs, ea, ei, _ = check_quadtree()
    check_hybrid()
    figure(ring, hs, ea, ei)
