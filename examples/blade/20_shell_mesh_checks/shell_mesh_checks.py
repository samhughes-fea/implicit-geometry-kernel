"""The shell mesh: span stations, chordwise sizing, mesh lines, orientation and rib fills.

Draws Figure R9.1 of the R manual (docs/r/R9_shell_mesh.md) and checks:

1. target edge and span bays, (R9.1)-(R9.2);
2. a web that leaves the section is found by bisection, (R9.3);
3. the sagitta edge keeps the chord within epsilon of a circle, (R9.4);
4. windowed curvature: true on a circle, a kink spread over the window, (R9.5);
5. grading limits the growth of the edge, (R9.6);
6. edge measure and node counts, (R9.7)-(R9.8);
7. a partition mesh: shared web feet, held lines, span grading and the hard invariants, (R9.9)-(R9.11), (R9.14);
8. orientation: skin flux equals twice the enclosed volume, webs face the TE, ribs +z, (R9.13);
9. rib fill: seed spacing, cell areas and Euler's formula, (R9.12).

    python examples/blade/20_shell_mesh_checks/shell_mesh_checks.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from _shared.plotting import REFERENCE, SERIES, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, TransverseMember  # noqa: E402
from implicit_geometry_kernel.blade.topology.sizing import EPS_OVER_C  # noqa: E402
from implicit_geometry_kernel.mesh.fill.delaunay import _edge_clearance, _unstructured_seeds  # noqa: E402
from implicit_geometry_kernel.mesh.shell import shell_mesh  # noqa: E402
from implicit_geometry_kernel.mesh.shell.invariants import check_shell_mesh  # noqa: E402
from implicit_geometry_kernel.mesh.shell.partition_mesh import (  # noqa: E402
    _GRADING,
    _KAPPA_WINDOW,
    _edge_measure,
    _held_u,
    _mean_loop_edge,
    _ts_between,
    _window_kappa,
)
from implicit_geometry_kernel.mesh.shell.stations import (  # noqa: E402
    TIP_H_OVER_C,
    _anchor_stations,
    _drop_stations,
    _n_bays,
    _target_edge,
    mesh_stations,
)

H = 0.08


def blade(**kw):
    af = Airfoil.naca4("2418")
    st = [Station(z=0.0, chord=1.5, airfoil=af, twist=8.0), Station(z=6.0, chord=0.8, airfoil=af, twist=0.0)]
    kw.setdefault("transverse", (TransverseMember("rib_1", 3.0, 0.01),))
    return BladeGeometry(
        st,
        kw.pop("layout", "2D-F-LC-0R"),
        0.006,
        spar_cap_thickness=0.012,
        shear_web_thickness=0.006,
        pivot_fraction=0.3,
        **kw,
    )


def ring_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


# --- 1. target edge and span bays ------------------------------------------------------------


def check_span(geo):
    s0 = _target_edge(geo, 0.0, None, TIP_H_OVER_C)
    s1 = _target_edge(geo, 0.0, 0.05, TIP_H_OVER_C)
    bays = [(d, _n_bays(d, 0.1)) for d in (0.07, 0.1, 0.13, 0.35, 1.0)]
    for d, n in bays:
        assert n == max(1, int(np.ceil(d / 0.1 - 1e-9)))
        assert n == 1 or 0.05 < d / n <= 0.1 + 1e-12
    zs = mesh_stations(geo, h=H)
    anchors = _anchor_stations(geo)
    have = all(any(abs(a - z) < 1e-9 for z in zs) for a in anchors)
    worst = max(np.diff(zs)) / H
    assert abs(s0 - 0.05 * 1.5) < 1e-12 and s1 == 0.05 and have and worst <= 1 + 1e-9
    print(
        f"1. Target edge (R9.1): default 0.05 c = {s0:.4f} m at the root; min with h = 0.05 gives {s1}. Bays (R9.2): "
        + ", ".join(f"D = {d}: n = {n}" for d, n in bays)
        + f" for s = 0.1. Stations at h = {H}: {len(zs)}, anchors {anchors} kept, largest bay {worst:.3f} h"
    )
    return zs


# --- 2. drop station -----------------------------------------------------------------------


def check_drop():
    af = Airfoil.naca4("2418")
    st = [Station(z=0.0, chord=2.0, airfoil=af), Station(z=8.0, chord=0.6, airfoil=af)]
    geo = BladeGeometry(
        st, "2D-F-LE-0R", 0.006, spar_cap_thickness=0.012, shear_web_thickness=0.006, pivot_fraction=0.3
    )
    drops = _drop_stations(geo)
    zd = drops[0]
    before = len(geo.section_at(zd - 1e-4).dropped_web_indices)
    after = len(geo.section_at(zd + 1e-4).dropped_web_indices)
    assert len(drops) == 1 and after == before + 1
    print(
        f"2. Drop station (R9.3): an LE-frozen web on a chord tapering 2.0 to 0.6 m leaves the section at "
        f"z = {zd:.6f} m; dropped webs {before} at z - 1e-4 and {after} at z + 1e-4 (bracket 1e-6 of the span)"
    )


# --- 3-6. chordwise sizing ------------------------------------------------------------------


def check_sagitta():
    rows = []
    for c, radius in ((1.0, 0.02), (2.0, 0.1), (3.0, 0.5)):
        eps = EPS_OVER_C * c
        ell = np.sqrt(8 * eps / (1 / radius))
        sag = radius * (1 - np.cos(ell / (2 * radius)))
        rows.append((c, radius, ell, sag / eps))
    assert all(abs(r[3] - 1) < 0.02 for r in rows)
    print(
        "3. Sagitta edge (R9.4): chord-to-arc height over epsilon = "
        + ", ".join(f"{r[3]:.5f} (c = {r[0]}, R = {r[1]}, edge {r[2] * 1e3:.1f} mm)" for r in rows)
    )


def check_kappa():
    radius, c = 0.3, 1.0
    a = np.linspace(0, np.pi, 400)
    arc = np.column_stack([radius * np.cos(a), radius * np.sin(a)])
    k = _window_kappa(arc, _KAPPA_WINDOW * c)
    mid = k[50:-50]
    e1 = float(np.max(np.abs(mid * radius - 1)))
    alpha = np.radians(20.0)
    s = np.linspace(-0.2, 0.2, 401)
    kink = np.column_stack([s, np.where(s > 0, s * np.tan(alpha), 0.0)])
    kk = _window_kappa(kink, _KAPPA_WINDOW * c)
    peak = float(kk.max())
    assert e1 < 1e-3 and abs(peak - alpha / (_KAPPA_WINDOW * c)) / (alpha / (_KAPPA_WINDOW * c)) < 0.05
    print(
        f"4. Windowed curvature (R9.5), window {_KAPPA_WINDOW} c: circle R = {radius}: kappa R - 1 within {e1:.1e}; "
        f"a 20 deg kink gives {peak:.2f} 1/m = angle / window ({alpha / (_KAPPA_WINDOW * c):.2f}), not 1 / segment"
    )


def loop_points(geo, z):
    g = geo.section_at(z).shell_graph()
    return np.asarray(g.loop_polygon(), float), g


def check_grading(geo):
    """The edge each loop partition gets, computed on that partition alone, as the mesher does."""
    graph = geo.section_at(0.0).shell_graph()
    c = 1.5
    s_all, raw_all, edge_all, ends = [], [], [], []
    offset, slope = 0.0, 0.0
    for _name, xy, rev in graph.closed_loop_chain():
        pts = xy[::-1] if rev else xy
        seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        s = np.concatenate([[0], np.cumsum(seg)])
        kap = _window_kappa(pts, _KAPPA_WINDOW * c)
        sag = np.where(kap > 1e-12, np.sqrt(8 * EPS_OVER_C * c / np.maximum(kap, 1e-12)), np.inf)
        raw = np.minimum(H, sag)
        edge = raw.copy()
        for i in range(1, len(edge)):
            edge[i] = min(edge[i], edge[i - 1] + _GRADING * seg[i - 1])
        for i in range(len(edge) - 2, -1, -1):
            edge[i] = min(edge[i], edge[i + 1] + _GRADING * seg[i])
        tt, w = _edge_measure(pts, H, c)
        density = np.diff(w) / np.maximum(seg, 1e-15)
        assert np.allclose(density, 0.5 * (1 / edge[1:] + 1 / edge[:-1]), rtol=1e-12)
        assert np.all(edge <= raw + 1e-15)
        slope = max(slope, float(np.max(np.abs(np.diff(edge)) / np.maximum(seg, 1e-15))))
        s_all.append(s + offset)
        raw_all.append(raw)
        edge_all.append(edge)
        ends.append(offset)
        offset += float(s[-1])
    s, raw, edge = (np.concatenate(v) for v in (s_all, raw_all, edge_all))
    assert slope <= _GRADING + 1e-9
    print(
        f"5. Grading (R9.6): on each root loop partition the graded edge never exceeds the raw edge and changes at "
        f"most {slope:.4f} per metre of arc (limit {_GRADING}); smallest edge {edge.min() * 1e3:.1f} mm, at the "
        f"nose; the sharp TE is a partition end, so a mesh line and not a curvature"
    )
    return s, raw, edge, ends


def check_measure():
    pts = np.column_stack([np.linspace(0, 1.0, 11), np.zeros(11)])
    tt, w = _edge_measure(pts, 0.15, None)
    n = int(np.ceil(w[-1] - 1e-9))
    ts = _ts_between(tt, w, 0.0, 1.0, n)
    steps = np.diff(ts)
    assert abs(w[-1] - 1 / 0.15) < 1e-12 and n == 7 and np.allclose(steps, 1 / 7)
    ref = {"a": 0.1, "b": 0.3, "c": 0.7}
    here = {"a": 0.2, "c": 0.6}
    held = _held_u(ref, here, "b")
    assert abs(held - (0.2 + (0.3 - 0.1) / (0.7 - 0.1) * (0.6 - 0.2))) < 1e-12
    print(
        f"6. Edge measure (R9.7)-(R9.8): a 1 m run at h = 0.15 has W = {w[-1]:.4f} edges, n = ceil(W) = {n}, steps "
        f"1/{n} exactly. Held line (R9.9): b at 0.3 between a (0.1) and c (0.7) is held at {held:.4f} = 0.2 + (1/3)(0.6 - 0.2)"
    )


# --- 7-9. a partition mesh ------------------------------------------------------------------


def check_mesh(geo):
    mesh = shell_mesh(geo, h=H)
    assert set(mesh.constructor) <= {"partition", "interior", "inherited"}
    problems = check_shell_mesh(mesh)
    quads = sum(1 for e in mesh.elements if len(e) == 4)
    tris = sum(1 for e in mesh.elements if len(e) == 3)
    skin_web_quads = all(len(e) == 4 for e, s in zip(mesh.elements, mesh.surface) if s in ("skin", "web"))
    # every web end node is also a skin node
    skin_nodes = {i for e, s in zip(mesh.elements, mesh.surface) if s == "skin" for i in e}
    web_nodes = {i for e, s in zip(mesh.elements, mesh.surface) if s == "web" for i in e}
    zs = np.asarray(mesh.stations)
    shared = 0
    for z in zs:
        row = [i for i in web_nodes if abs(mesh.nodes[i, 2] - z) < 1e-9]
        shared += sum(1 for i in row if i in skin_nodes)
    # span grading: each bay within 1.25 * min(mean edge) or the step floor h/2
    long_bays = []
    for za, zb in zip(zs[:-1], zs[1:]):
        lim = 1.25 * min(
            _mean_loop_edge(geo.section_at(float(z)).shell_graph(), H, float(geo.interpolate(float(z))[0]))
            for z in (za, zb)
        )
        if zb - za > lim * (1 + 1e-9) and zb - za > H:
            long_bays.append((za, zb))
    assert not problems and skin_web_quads and shared == 2 * 2 * len(zs) and not long_bays
    print(
        f"7. Partition mesh at h = {H}: {len(mesh.nodes)} nodes, {quads} quads (all skin and web), {tris} rib and tip "
        f"triangles, {len(zs)} stations; both feet of both webs are skin nodes at every station ({shared} shared); "
        f"no bay longer than 1.25 mean edges (R9.10); hard invariants (R9.14): {len(problems)} problems"
    )
    return mesh


def check_orientation(geo, mesh):
    nodes = mesh.nodes
    skin = [i for i, s in enumerate(mesh.surface) if s == "skin"]
    flux = 0.0
    for i in skin:
        p = nodes[list(mesh.elements[i])]
        area = 0.5 * np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)
        mean = p.mean(axis=0)
        flux += float(area[:2] @ mean[:2])
    zs = np.asarray(mesh.stations)
    areas = []
    for z in zs:
        ring_xy, _ = loop_points(geo, float(z))
        areas.append(abs(ring_area(ring_xy[:-1])))
    volume = float(np.trapezoid(areas, zs))
    te = {float(z): geo.oml_polyline(float(z))[0] for z in zs}
    web_ok = 0
    webs = [i for i, s in enumerate(mesh.surface) if s == "web"]
    for i in webs:
        p = nodes[list(mesh.elements[i])]
        area = 0.5 * np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)
        z = float(zs[np.argmin(np.abs(zs - p[:, 2].mean()))])
        web_ok += float(area[:2] @ (te[z] - p.mean(axis=0))[:2]) > 0
    rib_ok = all(
        0.5 * np.sum(np.cross(nodes[list(e)], np.roll(nodes[list(e)], -1, axis=0)), axis=0)[2] > 0
        for e, s in zip(mesh.elements, mesh.surface)
        if s in ("rib", "tip_closure")
    )
    ratio = flux / (2 * volume)
    assert abs(ratio - 1) < 5e-3 and web_ok == len(webs) and rib_ok
    print(
        f"8. Orientation (R9.13): skin flux sum A.(x - c) = {flux:.5f} m^3 against 2V = {2 * volume:.5f} "
        f"(ratio {ratio:.5f}); {web_ok} of {len(webs)} web quads face the TE; every rib and tip triangle along +z"
    )


def check_rib(geo, mesh):
    nodes = mesh.nodes
    ribs = [i for i, s in enumerate(mesh.surface) if s == "rib"]
    tri_area = sum(
        0.5 * np.linalg.norm(np.cross(nodes[e[1]] - nodes[e[0]], nodes[e[2]] - nodes[e[0]]))
        for e in (mesh.elements[i] for i in ribs)
    )
    ring, _ = loop_points(geo, 3.0)
    loop_area = abs(ring_area(ring[:-1]))
    vs = {k for i in ribs for k in mesh.elements[i]}
    uses: dict = {}
    for i in ribs:
        e = mesh.elements[i]
        for a, b in zip(e, e[1:] + e[:1]):
            uses.setdefault(tuple(sorted((a, b))), []).append((a, b))
    euler = len(vs) - len(uses) + len(ribs)
    nxt = {d[0][0]: d[0][1] for d in uses.values() if len(d) == 1}
    start = next(iter(nxt))
    order, k = [start], nxt[start]
    while k != start:
        order.append(k)
        k = nxt[k]
    cell_area = abs(ring_area(nodes[order][:, :2]))
    a = np.linspace(0, 2 * np.pi, 60, endpoint=False)
    disk = np.column_stack([0.5 * np.cos(a), 0.3 * np.sin(a)])
    seeds = _unstructured_seeds(disk, 0.06)
    d = np.linalg.norm(seeds[:, None] - seeds[None], axis=2) + np.eye(len(seeds)) * 9
    clear = _edge_clearance(disk, seeds).min()
    assert abs(tri_area / cell_area - 1) < 1e-9 and euler == 1 and d.min() >= 0.06 and clear >= 0.62 * 0.06
    print(
        f"9. Rib fill (R9.12): rib triangles cover their boundary ring ({len(order)} skin nodes) at z = 3 to "
        f"{abs(tri_area / cell_area - 1):.1e} ({cell_area / loop_area:.4f} of the mid ring: the ring is inscribed); "
        f"V - E + F = {euler} (a disk); seeds on an ellipse at r = 0.06: closest pair {d.min():.4f} >= r, "
        f"clearance {clear:.4f} >= 0.62 r"
    )
    return seeds, disk


# --- figure ---------------------------------------------------------------------------------


def figure(geo, mesh, grading, seeds, disk):
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), gridspec_kw={"wspace": 0.25, "hspace": 0.38})
    ax = axes[0, 0]
    s, raw, edge, ends = grading
    ax.semilogy(s, np.full_like(s, H), "--", color=REFERENCE, lw=0.8, label="$h$")
    ax.semilogy(s, raw, color=SERIES[1], lw=1, label="$\\min(h, \\ell_\\mathrm{sag})$ (R9.4)")
    ax.semilogy(s, edge, color=SERIES[0], lw=1.4, label="graded edge (R9.6)")
    for e in ends[1:]:
        ax.axvline(e, color="0.8", lw=0.5, zorder=0)
    i = int(np.argmin(edge))
    ax.annotate("LE", (s[i], edge[i]), xytext=(6, -2), textcoords="offset points", fontsize=7)
    ax.set_xlabel("arc length round the root skin loop (m); grey: partition ends")
    ax.set_ylabel("edge (m)")
    ax.set_title("Chordwise edge on the root section", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="lower left")

    ax = axes[0, 1]
    z0 = float(mesh.stations[0])
    ring, graph = loop_points(geo, z0)
    ax.plot(ring[:, 0], ring[:, 1], color=REFERENCE, lw=0.6)
    for w in graph.webs:
        ax.plot(w.points[:, 0], w.points[:, 1], color=REFERENCE, lw=0.6)
    at = [i for i in range(len(mesh.nodes)) if abs(mesh.nodes[i, 2] - z0) < 1e-9]
    skin_nodes = {i for e, srf in zip(mesh.elements, mesh.surface) if srf == "skin" for i in e}
    web_nodes = {i for e, srf in zip(mesh.elements, mesh.surface) if srf == "web" for i in e}
    p = mesh.nodes[[i for i in at if i in skin_nodes and i not in web_nodes]]
    q = mesh.nodes[[i for i in at if i in web_nodes and i not in skin_nodes]]
    f = mesh.nodes[[i for i in at if i in web_nodes and i in skin_nodes]]
    ax.plot(p[:, 0], p[:, 1], ".", color=SERIES[0], ms=3, label="skin nodes")
    ax.plot(q[:, 0], q[:, 1], ".", color=SERIES[2], ms=3, label="web nodes")
    ax.plot(f[:, 0], f[:, 1], "o", mfc="none", color=SERIES[1], ms=6, label="web feet: shared")
    ax.set_aspect("equal")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$y$ (m)")
    ax.set_title("Root station of the partition mesh", fontsize=9)
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.4))

    ax = axes[1, 0]
    nodes = mesh.nodes
    for e, srf, part in zip(mesh.elements, mesh.surface, mesh.partition):
        if srf != "skin" or "suction" not in part:
            continue
        pts = nodes[list(e) + [e[0]]]
        if np.mean(pts[:, 2]) > 3.5:
            continue
        ax.plot(pts[:, 2], pts[:, 0], color=SERIES[0], lw=0.3)
    for z in mesh.stations:
        if z <= 3.5:
            ax.axvline(z, color="0.85", lw=0.4, zorder=0)
    ax.axvline(3.0, color=SERIES[1], lw=1.0, label="rib (anchor)")
    ax.set_xlabel("span $z$ (m)")
    ax.set_ylabel("$x$ (m)")
    ax.set_title("Suction-side skin quads, root to mid-span (planform)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="lower left")

    ax = axes[1, 1]
    from implicit_geometry_kernel.mesh.fill import triangulate_ring

    ring3 = np.column_stack([disk, np.zeros(len(disk))])
    tris = triangulate_ring(ring3, size=0.06, pattern="unstructured", split_boundary=False)
    for t in np.asarray(tris):
        ax.fill(t[:, 0], t[:, 1], facecolor="none", edgecolor=SERIES[0], lw=0.4)
    ax.plot(np.r_[disk[:, 0], disk[0, 0]], np.r_[disk[:, 1], disk[0, 1]], color=REFERENCE, lw=0.8)
    ax.plot(seeds[:, 0], seeds[:, 1], ".", color=SERIES[1], ms=3, label="seeds before smoothing")
    ax.set_aspect("equal")
    ax.set_title("Rib fill: seeds, Delaunay, smoothing (R9.12)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.05))
    ax.axis("off")
    note(
        fig,
        "Checks: target edge and bays; drop station; sagitta edge; windowed curvature; grading; edge measure; "
        "partition mesh and invariants; orientation; rib fill.",
    )
    fig.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.08)
    save(fig, __file__, "shell_mesh_checks")


if __name__ == "__main__":
    geo = blade()
    check_span(geo)
    check_drop()
    check_sagitta()
    check_kappa()
    grading = check_grading(geo)
    check_measure()
    mesh = check_mesh(geo)
    check_orientation(geo, mesh)
    seeds, disk = check_rib(geo, mesh)
    figure(geo, mesh, grading, seeds, disk)
