"""Signed distance, CSG and the Newton projection: the results derived in R4.

Draws Figure R4.1 of the R manual (docs/r/R4_sdf_csg.md) and checks:

1. the primitive fields against the distance to a densely sampled boundary, (R4.3)-(R4.10);
2. CSG by min and max: exact outside a union and inside an intersection or a
   difference, a lower bound on |d| elsewhere, (R4.11)-(R4.14);
3. offsets and shells: exact outside a dilation and inside an erosion, (R4.15)-(R4.17);
4. the smooth union: equal to min beyond the blend width, k/4 below it at a = b,
   1-Lipschitz and C1, (R4.19)-(R4.21);
5. transforms: rigid motions and uniform scaling keep the distance, anisotropic
   scaling keeps a bound, (R4.22)-(R4.23);
6. the Newton projection lands on a circle in one step and the two-field solve
   converges quadratically, (R4.27)-(R4.28);
7. the conservative bounding box and the zero-loop extraction, (R4.29)-(R4.31),
   including the cases fixed as SD-002 and SD-003.

    python examples/sdf/03_distance_and_csg/distance_and_csg.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from _shared.plotting import NAVY, REFERENCE, SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel.extract import extract_zero_loops  # noqa: E402
from implicit_geometry_kernel.sdf import (  # noqa: E402
    Box,
    Circle,
    GeometrySolver,
    HalfPlane,
    Intersect,
    Offset,
    Polygon,
    Shell,
    Union,
    eval_expr,
    rotate_field,
    scale_field,
    sdf_box,
    sdf_capsule,
    sdf_circle,
    sdf_ellipse,
    sdf_oriented_box,
    sdf_polygon,
    sdf_rounded_box,
    smooth_union,
    translate_field,
)

RNG = np.random.default_rng(4)
N_BOUNDARY = 400_000


# --- reference distances: nearest point of a densely sampled boundary ----------------


def _ring(t, cx, cy, rx, ry):
    return np.column_stack([cx + rx * np.cos(t), cy + ry * np.sin(t)])


def circle_boundary(cx, cy, r, n=N_BOUNDARY):
    return _ring(np.linspace(0.0, 2 * np.pi, n, endpoint=False), cx, cy, r, r)


def polyline_boundary(verts, n=N_BOUNDARY):
    v = np.vstack([verts, verts[:1]])
    seg = np.linalg.norm(np.diff(v, axis=0), axis=1)
    s = np.linspace(0.0, seg.sum(), n, endpoint=False)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    k = np.clip(np.searchsorted(cum, s, side="right") - 1, 0, len(seg) - 1)
    w = (s - cum[k]) / seg[k]
    return v[k] + w[:, None] * (v[k + 1] - v[k])


def box_vertices(cx, cy, hw, hh):
    return np.array([[cx - hw, cy - hh], [cx + hw, cy - hh], [cx + hw, cy + hh], [cx - hw, cy + hh]])


def rounded_box_boundary(cx, cy, hw, hh, r, n=N_BOUNDARY):
    inner = box_vertices(cx, cy, hw - r, hh - r)
    t = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    pts = []
    for (qx, qy), a0 in zip(inner[[2, 3, 0, 1]], (0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi)):
        tt = t[(t >= a0) & (t < a0 + 0.5 * np.pi)]
        pts.append(_ring(tt, qx, qy, r, r))
    sides = [
        np.column_stack([np.full(n // 4, cx + hw), np.linspace(cy - hh + r, cy + hh - r, n // 4)]),
        np.column_stack([np.full(n // 4, cx - hw), np.linspace(cy - hh + r, cy + hh - r, n // 4)]),
        np.column_stack([np.linspace(cx - hw + r, cx + hw - r, n // 4), np.full(n // 4, cy + hh)]),
        np.column_stack([np.linspace(cx - hw + r, cx + hw - r, n // 4), np.full(n // 4, cy - hh)]),
    ]
    return np.vstack(pts + sides)


def in_rounded_box(x, y, hw, hh, r):
    qx, qy = np.abs(x) - (hw - r), np.abs(y) - (hh - r)
    corner = (qx > 0) & (qy > 0) & (qx**2 + qy**2 > r * r)
    return (np.abs(x) < hw) & (np.abs(y) < hh) & ~corner


def capsule_boundary(a, b, r, n=N_BOUNDARY):
    a, b = np.asarray(a, float), np.asarray(b, float)
    u = (b - a) / np.linalg.norm(b - a)
    nrm = np.array([-u[1], u[0]])
    s = np.linspace(0.0, 1.0, n // 4)[:, None]
    sides = [a + s * (b - a) + r * nrm, a + s * (b - a) - r * nrm]
    th = np.arctan2(u[1], u[0])
    t = np.linspace(-0.5 * np.pi, 0.5 * np.pi, n // 4)
    caps = [_ring(th + t, *b, r, r), _ring(th + np.pi + t, *a, r, r)]
    return np.vstack(sides + caps)


def reference(points, boundary, inside):
    d, _ = cKDTree(boundary).query(points)
    return np.where(inside, -d, d)


def sample(n=6000, lo=-1.6, hi=1.6):
    return RNG.uniform(lo, hi, (n, 2))


def away(phi_ref, tol=2e-3):
    """Points not within tol of the boundary, where the sign is unambiguous."""
    return np.abs(phi_ref) > tol


# --- 1. primitives --------------------------------------------------------------------


def check_primitives():
    p = sample()
    x, y = p[:, 0], p[:, 1]
    star = np.array(
        [
            [np.cos(a) * r, np.sin(a) * r]
            for a, r in zip(np.linspace(0.0, 2 * np.pi, 10, endpoint=False), [1.0, 0.45] * 5)
        ]
    )
    ang = 0.6
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    obox = box_vertices(0.0, 0.0, 0.9, 0.35) @ rot.T + [0.1, -0.2]
    cases = {
        "circle (R4.3)": (
            sdf_circle(x, y, 0.2, -0.1, 0.8),
            reference(p, circle_boundary(0.2, -0.1, 0.8), np.hypot(x - 0.2, y + 0.1) < 0.8),
        ),
        "box (R4.4)": (
            sdf_box(x, y, 0.1, 0.0, 0.9, 0.4),
            reference(
                p, polyline_boundary(box_vertices(0.1, 0.0, 0.9, 0.4)), (np.abs(x - 0.1) < 0.9) & (np.abs(y) < 0.4)
            ),
        ),
        "rounded box (R4.5)": (
            sdf_rounded_box(x, y, 0.0, 0.0, 1.0, 0.5, 0.2),
            reference(p, rounded_box_boundary(0.0, 0.0, 1.0, 0.5, 0.2), in_rounded_box(x, y, 1.0, 0.5, 0.2)),
        ),
        "oriented box (R4.6)": (
            sdf_oriented_box(x, y, 0.1, -0.2, 0.9, 0.35, ang),
            reference(p, polyline_boundary(obox), sdf_polygon(x, y, obox) < 0),
        ),
        "capsule (R4.7)": (sdf_capsule(x, y, -0.7, -0.3, 0.6, 0.5, 0.25), None),
        "star polygon (R4.9)": (sdf_polygon(x, y, star), None),
    }
    a, b, r = np.array([-0.7, -0.3]), np.array([0.6, 0.5]), 0.25
    t = np.clip(((p - a) @ (b - a)) / ((b - a) @ (b - a)), 0.0, 1.0)
    dseg = np.linalg.norm(p - (a + t[:, None] * (b - a)), axis=1)
    cases["capsule (R4.7)"] = (cases["capsule (R4.7)"][0], reference(p, capsule_boundary(a, b, r), dseg < r))
    # winding by the crossing count of a horizontal ray, independent of the kernel
    vs = np.vstack([star, star[:1]])
    cross = np.zeros(len(p), dtype=int)
    for (x0, y0), (x1, y1) in zip(vs[:-1], vs[1:]):
        hit = ((y0 > y) != (y1 > y)) & (x < x0 + (y - y0) * (x1 - x0) / (y1 - y0 + 1e-300))
        cross += hit
    cases["star polygon (R4.9)"] = (
        cases["star polygon (R4.9)"][0],
        reference(p, polyline_boundary(star), cross % 2 == 1),
    )
    print("1. Primitive fields against the distance to a densely sampled boundary:")
    worst = 0.0
    for name, (phi, ref) in cases.items():
        err = float(np.abs(phi - ref).max())
        worst = max(worst, err)
        print(f"   {name:22s} {err:.1e}")
    assert worst < 1e-5
    hp = HalfPlane(nx=1.0, ny=2.0, d=-0.5)
    hp_err = float(np.abs(eval_expr(hp, x, y) - (x + 2 * y - 0.5) / np.sqrt(5.0)).max())
    print(f"   half-plane (R4.8)      {hp_err:.1e}")
    assert hp_err < 1e-15
    ell = sdf_ellipse(x, y, 0.0, 0.0, 1.0, 0.5)
    ell_ref = reference(p, circle_boundary(0, 0, 1) * [1.0, 0.5], x**2 + (y / 0.5) ** 2 < 1.0)
    ell_err = float(np.abs(ell - ell_ref).max())
    print(f"   ellipse (R4.10)        {ell_err:.1e}  (SD-001 fixed: root of (R4.10) by bisection)")
    assert ell_err < 1e-5
    return worst, ell_err


# --- 2. CSG ---------------------------------------------------------------------------

A = (0.25, 0.1, 0.7)  # disk centre and radius
B = (-0.35, -0.1, 0.6, 0.45)  # box centre and half sizes


def in_a(p):
    return np.hypot(p[:, 0] - A[0], p[:, 1] - A[1]) < A[2]


def in_b(p):
    return (np.abs(p[:, 0] - B[0]) < B[2]) & (np.abs(p[:, 1] - B[1]) < B[3])


def csg_reference(p, op, n=N_BOUNDARY):
    ba, bb = circle_boundary(*A, n=n), polyline_boundary(box_vertices(*B), n=n)
    if op == "union":
        bd = np.vstack([ba[~in_b(ba)], bb[~in_a(bb)]])
        inside = in_a(p) | in_b(p)
    elif op == "intersect":
        bd = np.vstack([ba[in_b(ba)], bb[in_a(bb)]])
        inside = in_a(p) & in_b(p)
    else:  # A minus B
        bd = np.vstack([ba[~in_b(ba)], bb[in_a(bb)]])
        inside = in_a(p) & ~in_b(p)
    return reference(p, bd, inside), inside


def check_csg():
    p = sample(8000)
    a, b = Circle(*A), Box(*B)
    x, y = p[:, 0], p[:, 1]
    da, db = eval_expr(a, x, y), eval_expr(b, x, y)
    phis = {"union": np.minimum(da, db), "intersect": np.maximum(da, db), "subtract": np.maximum(da, -db)}
    exprs = {"union": Union(a, b), "intersect": Intersect(a, b)}
    for k, e in exprs.items():
        assert np.array_equal(eval_expr(e, x, y), phis[k])
    print("\n2. CSG against the distance to the true boundary (exact side | bound side):")
    out = {}
    for op, phi in phis.items():
        ref, inside = csg_reference(p, op)
        ok = away(ref)
        exact_side = ~inside if op == "union" else inside
        err = float(np.abs(phi - ref)[ok & exact_side].max())
        ratio = float((np.abs(phi) / np.abs(ref))[ok & ~exact_side].max())
        worst_gap = float((np.abs(ref) - np.abs(phi))[ok & ~exact_side].max())
        sign_ok = bool(np.all(np.sign(phi[ok]) == np.sign(ref[ok])))
        print(
            f"   {op:9s} exact side error {err:.1e}; bound side |phi|/|d| <= {ratio:.6f}, "
            f"largest shortfall {worst_gap:.3f}; signs agree: {sign_ok}"
        )
        assert err < 1e-5 and ratio <= 1.0 + 1e-5 and sign_ok
        out[op] = (err, worst_gap)
    # 1-Lipschitz for a deeper tree, (R4.14)
    tree = Union(Intersect(a, b), Offset(Circle(-0.6, 0.6, 0.3), 0.05))
    q1, q2 = sample(20000), sample(20000)
    lip = np.abs(eval_expr(tree, q1[:, 0], q1[:, 1]) - eval_expr(tree, q2[:, 0], q2[:, 1]))
    lip_ratio = float((lip / np.linalg.norm(q1 - q2, axis=1)).max())
    print(f"   min/max/offset tree: largest |phi(p) - phi(q)| / |p - q| = {lip_ratio:.6f} (<= 1, (R4.14))")
    assert lip_ratio <= 1.0 + 1e-12
    return out


# --- 3. offsets and shells ------------------------------------------------------------

L_SHAPE = np.array([[-0.8, -0.6], [0.8, -0.6], [0.8, -0.1], [-0.2, -0.1], [-0.2, 0.7], [-0.8, 0.7]])


def inside_polygon(p, verts):
    """Even-odd crossing count of a horizontal ray, independent of the kernel."""
    x, y = p[:, 0], p[:, 1]
    vs = np.vstack([verts, verts[:1]])
    cross = np.zeros(len(p), dtype=int)
    for (x0, y0), (x1, y1) in zip(vs[:-1], vs[1:]):
        cross += ((y0 > y) != (y1 > y)) & (x < x0 + (y - y0) * (x1 - x0) / (y1 - y0 + 1e-300))
    return cross % 2 == 1


def level_set(verts, delta, outside, n=20000):
    """Points at distance delta from the polygon boundary on one side: edge offsets and vertex arcs."""
    vs = np.vstack([verts, verts[:1]])
    cand = []
    s = np.linspace(0.0, 1.0, n)[:, None]
    for a, b in zip(vs[:-1], vs[1:]):
        t = (b - a) / np.linalg.norm(b - a)
        nrm = np.array([t[1], -t[0]])
        for sgn in (1.0, -1.0):
            cand.append(a + s * (b - a) + sgn * delta * nrm)
    th = np.linspace(0.0, 2 * np.pi, 4 * n, endpoint=False)
    for v in verts:
        cand.append(_ring(th, v[0], v[1], delta, delta))
    c = np.vstack(cand)
    dist = np.full(len(c), np.inf)
    for a, b in zip(vs[:-1], vs[1:]):  # distance to each edge in NumPy, not the kernel's loop
        t = np.clip(((c - a) @ (b - a)) / ((b - a) @ (b - a)), 0.0, 1.0)
        dist = np.minimum(dist, np.linalg.norm(c - a - t[:, None] * (b - a), axis=1))
    keep = (dist > delta - 1e-9) & (inside_polygon(c, verts) != outside)
    return c[keep]


def check_offsets():
    p = sample(8000, -1.4, 1.4)
    x, y = p[:, 0], p[:, 1]
    poly = Polygon(tuple(map(tuple, L_SHAPE)))
    d = eval_expr(poly, x, y)
    delta = 0.15
    phi_dil = eval_expr(Offset(poly, delta), x, y)
    out_dil = d > delta + 2e-3
    dil_ref, _ = cKDTree(level_set(L_SHAPE, delta, outside=True)).query(p[out_dil])
    err_dil = float(np.abs(phi_dil[out_dil] - dil_ref).max())
    phi_ero = eval_expr(Offset(poly, -0.1), x, y)
    in_ero = d < -0.1 - 2e-3
    ero_ref, _ = cKDTree(level_set(L_SHAPE, 0.1, outside=False)).query(p[in_ero])
    err_ero = float(np.abs(-phi_ero[in_ero] - ero_ref).max())
    sh = eval_expr(Shell(poly, 0.1), x, y)
    shell_out = np.abs(d) > 0.05 + 2e-3
    shell_ref = np.abs(d[shell_out]) - 0.05
    err_sh = float(np.abs(sh[shell_out] - shell_ref).max())
    print("\n3. Offsets of the L-shaped polygon (the level sets are sampled independently):")
    print(f"   dilation by {delta}: outside, phi against distance to the level set {err_dil:.1e} (R4.15)")
    print(f"   erosion by 0.1: inside, -phi against distance to the level set {err_ero:.1e} (R4.16)")
    print(f"   shell of thickness 0.1: outside the band, phi = |d| - t/2 to {err_sh:.1e} (R4.17)")
    assert err_dil < 1e-5 and err_ero < 1e-5 and err_sh < 1e-15
    return err_dil, err_ero


# --- 4. smooth union ------------------------------------------------------------------


def smin_scalar(a, b, k):
    return smooth_union(np.asarray(a, float), np.asarray(b, float), k)(None, None)


def check_smooth_union():
    k = 0.2
    u = np.linspace(-0.6, 0.6, 120001)
    f = smin_scalar(np.zeros_like(u), u, k)
    m = np.minimum(0.0, u)
    far = np.abs(u) >= k
    dev_far = float(np.abs(f - m)[far].max())
    dev_mid = float((m - f).max())
    at = float(u[np.argmax(m - f)])
    du = u[1] - u[0]
    g = np.gradient(f, du)
    jump = max(float(abs(g[np.searchsorted(u, s * k) + 2] - g[np.searchsorted(u, s * k) - 2])) for s in (-1, 1))
    a, b = Circle(-0.35, 0.0, 0.4).as_xy(), Box(0.35, 0.0, 0.3, 0.3).as_xy()
    su = smooth_union(a, b, k)
    q1, q2 = sample(20000, -1, 1), sample(20000, -1, 1)
    lip = float((np.abs(su(q1[:, 0], q1[:, 1]) - su(q2[:, 0], q2[:, 1])) / np.linalg.norm(q1 - q2, axis=1)).max())
    print("\n4. Smooth union, k = 0.2 (R4.19)-(R4.21):")
    print(
        f"   |smin - min| beyond |a - b| >= k: {dev_far:.1e}; largest min - smin {dev_mid:.6f} "
        f"at a - b = {at:.1e} (k/4 = {k / 4})"
    )
    print(f"   slope jump across |a - b| = k: {jump:.1e} (C1); Lipschitz ratio {lip:.6f} (<= 1)")
    assert dev_far < 1e-15 and abs(dev_mid - k / 4) < 1e-12 and jump < 1e-4 and lip <= 1 + 1e-12
    return u, f, m, k


# --- 5. transforms --------------------------------------------------------------------


def check_transforms():
    p = sample()
    x, y = p[:, 0], p[:, 1]
    box = lambda xx, yy: sdf_box(xx, yy, 0.0, 0.0, 0.8, 0.3)  # noqa: E731
    moved = translate_field(rotate_field(box, 0.6), 0.1, -0.2)
    err_rig = float(np.abs(moved(x, y) - sdf_oriented_box(x, y, 0.1, -0.2, 0.8, 0.3, 0.6)).max())
    circ = lambda xx, yy: sdf_circle(xx, yy, 0.0, 0.0, 0.5)  # noqa: E731
    err_uni = float(np.abs(scale_field(circ, 2.0)(x, y) - sdf_circle(x, y, 0.0, 0.0, 1.0)).max())
    an = scale_field(circ, 2.0, 0.6)
    q1, q2 = sample(20000), sample(20000)
    lip = float((np.abs(an(q1[:, 0], q1[:, 1]) - an(q2[:, 0], q2[:, 1])) / np.linalg.norm(q1 - q2, axis=1)).max())
    inside = (x / 1.0) ** 2 + (y / 0.3) ** 2 < 1.0
    sign_ok = bool(np.all((an(x, y) < 0) == inside))
    print("\n5. Transforms (R4.22)-(R4.23):")
    print(f"   rotate then translate a box against the oriented box: {err_rig:.1e}")
    print(f"   uniform scale 2 of a circle of radius 0.5 against radius 1: {err_uni:.1e}")
    print(f"   anisotropic (2, 0.6): Lipschitz ratio {lip:.6f} (<= 1), signs correct: {sign_ok}")
    assert err_rig < 1e-15 and err_uni < 1e-15 and lip <= 1 + 1e-12 and sign_ok


# --- 6. Newton ------------------------------------------------------------------------


def check_newton():
    circle = Circle(0.1, -0.2, 0.6)
    solver = GeometrySolver(tolerance=0.0, iterations=1, max_step=10.0)
    p0 = np.column_stack([sample(500, -1.5, 1.5), np.zeros(500)])
    keep = np.hypot(p0[:, 0] - 0.1, p0[:, 1] + 0.2) > 0.05
    one = solver.project(circle, p0[keep])
    err_one = float(np.abs(np.hypot(one.points[:, 0] - 0.1, one.points[:, 1] + 0.2) - 0.6).max())
    line = HalfPlane(nx=1.0, ny=-0.4, d=-0.05)
    hist = []
    p = np.array([[0.9, 0.5, 0.0]])
    for _ in range(6):
        res = GeometrySolver(tolerance=0.0, iterations=1, max_step=10.0).intersect(circle, line, p)
        p = res.points
        hist.append(float(res.residual[0]))
    rates = [np.log(hist[i + 1]) / np.log(hist[i]) for i in range(len(hist) - 1) if 1e-14 < hist[i + 1]]
    print("\n6. Newton on fields (R4.27)-(R4.28):")
    print(f"   one projection step onto a circle from 500 points: |phi| <= {err_one:.1e}")
    print(
        f"   circle x line residuals: {', '.join(f'{h:.1e}' for h in hist)}; "
        f"log-ratio of successive residuals {', '.join(f'{r:.2f}' for r in rates)} (2 = quadratic)"
    )
    assert err_one < 1e-9 and hist[3] < 1e-12
    return hist


# --- 7. bounding boxes and zero loops -------------------------------------------------


def check_bbox_and_loops():
    e = Offset(Intersect(Box(0, 0, 1, 1), Box(2.5, 0, 1, 1)), 0.5)
    xs, ys = np.meshgrid(np.linspace(-3, 6, 901), np.linspace(-3, 3, 601))
    m = eval_expr(e, xs, ys) <= 0
    bb = e.bbox()
    print("\n7. Bounding boxes and zero loops (R4.29)-(R4.31):")
    print(
        f"   Offset(Intersect(two boxes 0.5 apart), 0.5): region y in [{ys[m].min():.2f}, {ys[m].max():.2f}], "
        f"expr_bbox y in [{bb[2]:.2f}, {bb[3]:.2f}] (SD-002 fixed: the box contains the region)"
    )
    assert bb[0] <= xs[m].min() and xs[m].max() <= bb[1] and bb[2] <= ys[m].min() and ys[m].max() <= bb[3]
    two = Union(Circle(-0.37, -0.37, 0.5), Circle(0.37, 0.37, 0.5))
    loops = extract_zero_loops(two, 0.0, (-1.5, 1.5, -1.5, 1.5), n=48)
    one = extract_zero_loops(Circle(0.0, 0.0, 0.5), 0.0, (-1.5, 1.5, -1.5, 1.5), n=48)
    snap = float(np.abs(np.hypot(one[0][:, 0], one[0][:, 1]) - 0.5).max())
    print(
        f"   one circle: {len(one)} loop, on phi = 0 to {snap:.1e}; two circles 0.05 apart: "
        f"{len(loops)} loops (SD-003 fixed: saddle cells by the asymptotic decider)"
    )
    assert len(one) == 1 and snap < 1e-8 and len(loops) == 2
    return len(loops)


# --- figure ---------------------------------------------------------------------------


def figure(smooth, hist):
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.6), gridspec_kw={"hspace": 0.42, "wspace": 0.28})
    ax = axes[0, 0]
    xs, ys = np.meshgrid(np.linspace(-1.3, 1.3, 261), np.linspace(-1.0, 1.0, 201))
    a, b = Circle(*A), Box(*B)
    phi = eval_expr(Union(a, b), xs, ys)
    pts = np.column_stack([xs.ravel(), ys.ravel()])
    ref, _ = csg_reference(pts, "union", n=40_000)
    gap = (np.abs(ref) - np.abs(phi.ravel())).reshape(xs.shape)
    im = ax.contourf(xs, ys, gap, levels=[1e-3, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16], cmap="Blues")
    ax.contour(xs, ys, phi, levels=np.arange(-0.4, 0.61, 0.1), colors="0.45", linewidths=0.5)
    ax.contour(xs, ys, phi, levels=[0.0], colors=[NAVY], linewidths=1.4)
    fig.colorbar(im, ax=ax, shrink=0.85, label="$|d| - |\\phi|$")
    ax.set_aspect("equal")
    ax.set_title("Union by min (R4.11): exact outside, a bound inside", fontsize=9)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

    ax = axes[0, 1]
    u, f, m, k = smooth
    ax.plot(u, m, color=REFERENCE, lw=0.9, ls="--", label="$\\min(0, u)$")
    ax.plot(u, f, color=SERIES[0], lw=1.4, label="smooth union (R4.19)")
    ax.axvspan(-k, k, color="0.93")
    ax.annotate("$k/4$", (0.0, -k / 8), xytext=(10, 0), textcoords="offset points", fontsize=8)
    ax.plot([0, 0], [0, -k / 4], color=SERIES[1], lw=1.0)
    ax.set_xlabel("$u = b - a$ with $a = 0$")
    ax.set_ylabel("field value")
    ax.set_title("Smooth union: blend width $k = 0.2$", fontsize=9)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    p = np.random.default_rng(0).uniform(-1.6, 1.6, (4000, 2))
    xx, yy = p[:, 0], p[:, 1]
    ell = sdf_ellipse(xx, yy, 0.0, 0.0, 1.0, 0.5)
    ref = reference(p, circle_boundary(0, 0, 1) * [1.0, 0.5], xx**2 + (yy / 0.5) ** 2 < 1.0)
    sc = ax.scatter(xx, yy, c=np.abs(ell - ref), s=3, cmap="Reds", vmin=0.0)
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), 0.5 * np.sin(t), color=NAVY, lw=1.2)
    fig.colorbar(sc, ax=ax, shrink=0.85, label="$|\\phi - d|$")
    ax.set_aspect("equal")
    ax.set_title("Ellipse field against the exact distance (R4.10)", fontsize=9)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

    ax = axes[1, 1]
    it = np.arange(1, len(hist) + 1)
    ax.semilogy(it, np.maximum(hist, 1e-17), "o-", color=SERIES[2], ms=4, label="circle $\\times$ line, (R4.28)")
    ax.axhline(1e-16, color="0.6", lw=0.6, ls=":")
    ax.set_xlabel("Newton iteration")
    ax.set_ylabel("$\\max(|\\phi_a|, |\\phi_b|)$")
    ax.set_title("Two-field Newton solve converges quadratically", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    note(
        fig,
        "Checks: primitives against a sampled boundary; CSG exact and bound sides; offsets; "
        "smooth union; transforms; Newton; bbox and zero loops.",
    )
    fig.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.1)
    save(fig, __file__, "distance_and_csg")


if __name__ == "__main__":
    check_primitives()
    check_csg()
    check_offsets()
    smooth = check_smooth_union()
    check_transforms()
    hist = check_newton()
    check_bbox_and_loops()
    figure(smooth, hist)
