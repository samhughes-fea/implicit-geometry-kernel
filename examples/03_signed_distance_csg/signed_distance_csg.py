"""Signed-distance solids and CSG: build a two-web box beam, then check the field and its integrals.

``implicit_geometry_kernel.sdf`` is a stand-alone 2-D toolbox: primitives (circle, box,
rounded box, polygon, ...), CSG operators (union = min, intersect = max, subtract,
offset, shell) and a structured ``SDFGrid``. ``integration`` turns any field into area
and second moments, on the grid by counting cells or on an adaptive quadtree that
reconstructs the boundary inside each cut cell.

The script

1. builds a thin-walled box beam with two webs from a rounded box, ``shell`` and two
   boxes, and plots its field;
2. checks the eikonal property |grad phi| = 1 that makes a field a distance: the exact
   operators (min, max, offset, shell) keep it away from the kinks they create, while
   ``smooth_union`` and ``blend`` trade it for smoothness;
3. integrates a thin circular tube (R = 0.30 m outside, 30 mm wall) on the grid and on
   the quadtree, at shrinking cell sizes, against the exact area and second moment;
4. checks the shell of a rounded box against Steiner's formula: a wall of thickness t
   centred on a convex midline of perimeter P has area exactly P t.

    python examples/03_signed_distance_csg/signed_distance_csg.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from _shared.plotting import REFERENCE, SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel.integration import (  # noqa: E402
    QuadtreeSpec,
    build_shared_quadtree,
    section_properties_from_quadtree,
)
from implicit_geometry_kernel.sdf import (  # noqa: E402
    SDFGrid,
    blend,
    intersect,
    offset,
    sdf_box,
    sdf_circle,
    sdf_rounded_box,
    shell,
    smooth_union,
    union,
)

# box beam: midline 1.2 m x 0.4 m with 0.08 m corner radius, 20 mm skin, 30 mm webs
A_HALF, B_HALF, R_CORNER, T_SKIN, T_WEB = 0.6, 0.2, 0.08, 0.02, 0.03
WEB_X = (-0.25, 0.25)
# thin tube for the convergence study
R_OUT, R_IN = 0.30, 0.27


def box_beam():
    """Skin: a shell on the rounded-box midline. Webs: boxes clipped to the skin's inside."""
    midline = lambda x, y: sdf_rounded_box(x, y, half_w=A_HALF, half_h=B_HALF, r=R_CORNER)  # noqa: E731
    skin = shell(midline, T_SKIN)
    inside = offset(midline, 0.0)  # the region inside the skin midline: webs end on it
    webs = [intersect(lambda x, y, c=c: sdf_box(x, y, cx=c, half_w=0.5 * T_WEB, half_h=1.0), inside) for c in WEB_X]
    beam = union(skin, union(*webs))
    smooth = smooth_union(skin, smooth_union(*webs, k=0.03), k=0.03)
    return midline, skin, webs, beam, smooth


def as_sampler(f):
    return lambda pts: f(np.asarray(pts, float)[:, 0], np.asarray(pts, float)[:, 1])


def quadtree_props(f, h, half=0.75):
    spec = QuadtreeSpec(-half, half, -half, half, max_depth=16, min_cell_size=h, refine_band=2 * h,
                        max_cells=2_000_000, root_cell_size=0.05, refine_probes="full")
    mesh = build_shared_quadtree(spec, {"f": as_sampler(f)})
    return section_properties_from_quadtree(mesh, mesh.phi_by_label, "f", sampler=as_sampler(f)), len(mesh.leaf_cells())


def eikonal_table(grid, fields, band):
    """Share of the grid within ``band`` of the boundary, away from kinks, where | |grad phi| - 1 | > 1 %."""
    print(f"\n{'field':<26} {'cells with | |grad phi| - 1 | > 1 %':>36}")
    out = {}
    for name, f in fields.items():
        phi = grid.eval(f)
        err = grid.eikonal_error(phi)
        near = np.abs(phi) < band
        out[name] = float(np.mean(err[near] > 0.01))
        print(f"{name:<26} {out[name]:36.1%}")
    return out


def tube_convergence():
    t = R_OUT - R_IN
    tube = shell(lambda x, y: sdf_circle(x, y, r=0.5 * (R_OUT + R_IN)), t)
    A, Ixx = np.pi * (R_OUT**2 - R_IN**2), 0.25 * np.pi * (R_OUT**4 - R_IN**4)
    rows = []
    print(f"\nThin tube: exact area {A:.6e} m^2, Ixx {Ixx:.6e} m^4")
    print(f"{'method':<9} {'cell [mm]':>9} {'cells':>8} {'area error':>11} {'Ixx error':>11}")
    for n in (100, 200, 400, 800, 1600):
        grid = SDFGrid.from_bbox(-0.4, 0.4, -0.4, 0.4, nx=n, ny=n)
        p = grid.section_properties(grid.eval(tube))
        rows.append(("grid", grid.dx, n * n, p["area"] / A - 1, p["Ixx"] / Ixx - 1))
    for h in (0.02, 0.01, 0.005, 0.0025, 0.00125):
        p, n = quadtree_props(tube, h, half=0.4)
        rows.append(("quadtree", h, n, p["area"] / A - 1, p["Ixx"] / Ixx - 1))
    for m, h, n, ea, ei in rows:
        print(f"{m:<9} {h * 1e3:9.2f} {n:8d} {ea:+11.2e} {ei:+11.2e}")
    return rows


def steiner_check():
    """Shell of a convex rounded box: area = P_mid * t exactly (Steiner)."""
    midline = lambda x, y: sdf_rounded_box(x, y, half_w=A_HALF, half_h=B_HALF, r=R_CORNER)  # noqa: E731
    perimeter = 4 * (A_HALF - R_CORNER) + 4 * (B_HALF - R_CORNER) + 2 * np.pi * R_CORNER
    exact = perimeter * T_SKIN
    p, _ = quadtree_props(shell(midline, T_SKIN), 0.001)
    err = p["area"] / exact - 1
    print(f"\nRounded-box shell: P t = {exact:.6e} m^2, quadtree {p['area']:.6e} m^2 ({err:+.2e})")
    assert abs(err) < 1e-4
    return err


def figure(grid, beam, smooth, rows, shares, steiner_err):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1.5, 1.5, 1.0]})
    phi = grid.eval(beam)
    ax = axes[0]
    lim = 0.15
    im = ax.contourf(grid.X, grid.Y, np.clip(phi, 0.0, lim), levels=np.linspace(0.0, lim, 16), cmap="Blues")
    ax.contour(grid.X, grid.Y, phi, levels=np.arange(0.02, lim, 0.02), colors="k", linewidths=0.3, alpha=0.4)
    ax.contourf(grid.X, grid.Y, phi, levels=[phi.min() - 1.0, 0.0], colors=[REFERENCE])
    fig.colorbar(im, ax=ax, shrink=0.75, label="distance to the solid, φ [m]")
    ax.set_title("Box beam: shell(rounded box) ∪ two clipped webs\n(solid φ < 0 in black, iso-distance lines every 20 mm)")
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

    ax = axes[1]
    err_exact = grid.eikonal_error(phi)
    err_smooth = grid.eikonal_error(grid.eval(smooth))
    show = np.where(np.abs(phi) < 0.1, err_smooth - err_exact, np.nan)
    im = ax.imshow(show, extent=(grid.X.min(), grid.X.max(), grid.Y.min(), grid.Y.max()), origin="lower",
                   cmap="magma_r", vmin=0.0, vmax=0.3)
    ax.contour(grid.X, grid.Y, phi, levels=[0.0], colors="k", linewidths=0.6)
    fig.colorbar(im, ax=ax, shrink=0.75, label="extra | |∇φ| − 1 | from smooth_union")
    ax.set_title(f"smooth_union is no longer a distance near the joints\n"
                 f"(|∇φ| off by > 1 % on {shares['smooth_union']:.0%} of the band, union {shares['union']:.0%})")
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")

    ax = axes[2]
    for i, (method, marker) in enumerate((("grid", "o"), ("quadtree", "s"))):
        h = np.array([r[1] for r in rows if r[0] == method]) * 1e3
        e = np.array([abs(r[4]) for r in rows if r[0] == method])
        ax.loglog(h, e, marker + "-", color=SERIES[i], label=f"{method}: Ixx")
    h = np.array([r[1] for r in rows if r[0] == "quadtree"]) * 1e3
    ax.loglog(h, 4e-4 * (h / h[0]) ** 2, "--", color=REFERENCE, lw=0.8, label="slope 2")
    ax.set_xlabel("cell size [mm]")
    ax.set_ylabel("|relative error| in Ixx")
    ax.set_title("Thin tube against π(R⁴ − r⁴)/4")
    ax.legend()
    ax.invert_xaxis()
    note(fig, "Box beam 1.2 m x 0.4 m midline, 80 mm corner radius, 20 mm skin, 30 mm webs. Tube 0.30 m outside, "
              f"30 mm wall. Rounded-box shell area matches Steiner's P·t to {abs(steiner_err):.0e}.")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, __file__, "signed_distance_csg")


if __name__ == "__main__":
    midline, skin, webs, beam, smooth = box_beam()
    grid = SDFGrid.from_bbox(-0.75, 0.75, -0.35, 0.35, nx=900, ny=420)
    shares = eikonal_table(grid, {"rounded box (exact)": midline, "union": beam, "smooth_union": smooth,
                                  "blend(skin, web, 0.5)": blend(skin, webs[0], 0.5)}, band=0.1)
    assert shares["union"] < shares["smooth_union"]
    rows = tube_convergence()
    assert abs(rows[-1][4]) < 1e-5  # quadtree at 1.25 mm
    err = steiner_check()
    figure(grid, beam, smooth, rows, shares, err)
