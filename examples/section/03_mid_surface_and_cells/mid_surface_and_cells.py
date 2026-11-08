"""The hand-off to meshing and beam models: wall midlines, the merged mid-surface and the cells.

A ``SectionCut`` holds three views of the same section:

``walls``          one midline per wall with its thickness, kind (skin, cap, web), side
                   and label; skins are split into bays at the webs, and each cap
                   records which skin vertices it is bonded to.
``cells``          the closed cells, each a polygon on the wall midlines with its area
                   and the labels of the walls around it.
``mid_surface()``  one graph of straight segments that a shell mesher or a thin-walled
                   beam model can use directly: skins and webs joined at shared nodes,
                   each cap a partition on the **skin mid** (``layers`` = (cap, skin),
                   inside first). Eccentricity is a laminate OFFSET downstream.

The script cuts a twisted NACA 2418 with box spar caps (2D-F) and checks the graph:

- Euler's formula for a connected plane graph: cells = segments - nodes + 1;
- every web joins two skin nodes and borders two cells; every skin bay borders one;
- the laminate normal z x (end - start) points out of the section on every skin
  segment and towards the trailing edge on every web segment;
- offset regions (cap footprints) lie on the skin midline, not t_cap/2 inside;
- each cell's area equals the shoelace area of its polygon.

    python examples/section/03_mid_surface_and_cells/mid_surface_and_cells.py
"""

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from shapely.geometry import LineString, Point, Polygon  # noqa: E402
from shapely.geometry.polygon import orient  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

from _shared.plotting import SERIES, WALL_COLOURS, draw_cut, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402

CHORD, TWIST, PITCH = 2.0, 10.0, 0.3
SKIN, CAP, WEB = 0.008, 0.03, 0.015


def cut_section():
    airfoil = Airfoil.naca4("2418")
    stations = [Station(0.0, CHORD, airfoil, twist=TWIST, pitch_axis=PITCH),
                Station(1.0, CHORD, airfoil, twist=TWIST, pitch_axis=PITCH)]
    return BladeGeometry(stations, "2D-F-LC-0R", SKIN, spar_cap_thickness=CAP, shear_web_thickness=WEB, pivot_fraction=0.3).section_at(0.5)


def shoelace(p):
    return 0.5 * abs(np.sum(p[:-1, 0] * p[1:, 1] - p[1:, 0] * p[:-1, 1]))


def check(cut, ms):
    P, segs = ms.points, ms.segments
    print("Walls")
    print(f"{'label':<36} {'kind':<5} {'side':<9} {'t [mm]':>7} {'length [m]':>10}  bonded to")
    for w in cut.walls:
        print(f"{w.label:<36} {w.kind:<5} {str(w.side):<9} {w.thickness * 1e3:7.1f} {w.length:10.4f}  "
              f"{', '.join(w.bonds) or '-'}")

    # Euler: a connected plane graph with V nodes and E edges has E - V + 1 bounded faces
    used = sorted({i for s in segs for i in (s.start, s.end)})
    euler = len(segs) - len(used) + 1
    print(f"\nMid-surface: {len(used)} nodes, {len(segs)} segments -> {euler} cells by Euler; the cut has {len(cut.cells)}")
    assert euler == len(cut.cells)

    degree = Counter(i for s in segs for i in (s.start, s.end))
    junctions = [i for i, d in degree.items() if d == 3]
    n_webs = len(cut.walls_of("web"))
    print(f"nodes of degree 3 (web-skin junctions): {len(junctions)} = 2 x {n_webs} webs")
    assert len(junctions) == 2 * n_webs and set(degree.values()) <= {2, 3}

    borders = Counter(label for c in cut.cells for label in c.wall_labels)
    for w in cut.walls_of("web"):
        assert borders[w.label] == 2, w.label
    for w in cut.walls_of("skin"):
        assert borders[w.label] == 1, w.label
    print("every web borders two cells and every skin bay one")

    # normals: z x d = (-dy, dx). Skins: compare with the outward normal of the nearest edge of
    # the region inside the skin midline (the union of the cells, oriented counter-clockwise).
    outline = orient(unary_union([Polygon(c.points[:, :2]) for c in cut.cells]), 1.0)
    ring = np.asarray(outline.exterior.coords)
    ra, rb = ring[:-1], ring[1:]
    rd = rb - ra
    outward = np.column_stack([rd[:, 1], -rd[:, 0]]) / np.linalg.norm(rd, axis=1)[:, None]
    chord_dir = np.array([np.cos(np.radians(TWIST)), -np.sin(np.radians(TWIST))])
    bad = 0
    for s in segs:
        d = P[s.end] - P[s.start]
        n = np.array([-d[1], d[0]]) / np.linalg.norm(d)
        if s.kind == "skin":
            m = 0.5 * (P[s.start] + P[s.end])
            u = np.clip(np.einsum("ij,ij->i", m - ra, rd) / np.einsum("ij,ij->i", rd, rd), 0.0, 1.0)
            k = int(np.argmin(np.linalg.norm(ra + u[:, None] * rd - m, axis=1)))
            bad += n @ outward[k] <= 0.0
        else:
            bad += n @ chord_dir <= 0.0
    print(f"segments whose normal points the wrong way: {bad} of {len(segs)}")
    assert bad == 0

    graph = cut.shell_graph()
    regions = graph.offset_regions()
    print(f"offset regions: {len(regions)}")
    for region in regions:
        print(f"  {region.name} eccentricity_m={region.eccentricity_m:.4f} "
              f"layers={region.layers}")
        assert region.requires_eccentricity
        assert region.eccentricity_m == 0.5 * CAP
    assert len(regions) == 2
    worst = 0.0
    skins = [w for w in cut.walls if w.kind == "skin"]
    for run in graph.loop:
        if not run.requires_eccentricity:
            continue
        host = next(w for w in skins if w.label in run.layers)
        skin = LineString(host.xy)
        for p in run.points[:, :2]:
            worst = max(worst, float(skin.distance(Point(p))))
    print(f"cap partitions: largest departure from the skin midline {worst * 1e6:.1f} µm")
    assert worst < 1e-9

    print(f"\n{'cell':<5} {'area [m^2]':>11} {'shoelace':>11}  walls")
    for c in cut.cells:
        a = shoelace(c.points[:, :2])
        print(f"{c.id:<5} {c.area:11.6f} {a:11.6f}  {', '.join(c.wall_labels)}")
        assert abs(a / c.area - 1) < 1e-9


def figure(cut, ms):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1.7, 1.0]})
    P = ms.points
    ax = axes[0]
    for c, colour in zip(cut.cells, SERIES[3:] + SERIES):
        ax.fill(c.points[:, 0], c.points[:, 1], color=colour, alpha=0.15, lw=0)
        cx, cy = Polygon(c.points[:, :2]).centroid.xy
        ax.text(cx[0], cy[0], f"cell {c.id}\n{c.area:.3f} m²", ha="center", va="center", fontsize=8)
    for s in ms.segments:
        colour = WALL_COLOURS["cap" if len(s.layers) == 2 else s.kind]
        ax.plot(P[[s.start, s.end], 0], P[[s.start, s.end], 1], color=colour, lw=1.4)
    every = ms.segments[::6]
    mids = np.array([0.5 * (P[s.start] + P[s.end]) for s in every])
    d = np.array([P[s.end] - P[s.start] for s in every])
    n = np.column_stack([-d[:, 1], d[:, 0]]) / np.linalg.norm(d, axis=1)[:, None]
    ax.quiver(mids[:, 0], mids[:, 1], n[:, 0], n[:, 1], color="k", scale=40, width=0.0018)
    ax.set_aspect("equal")
    ax.set_title("Skin-mid graph and cells; cap footprints coloured; arrows z × (end − start)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

    ax = axes[1]
    web = cut.walls_of("web")[0]
    top = web.xy[np.argmax(web.xy[:, 1])]
    draw_cut(ax, cut, alpha=0.3)
    for s in ms.segments:
        ax.plot(P[[s.start, s.end], 0], P[[s.start, s.end], 1], "k-", lw=1.0)
    ax.plot(P[:, 0], P[:, 1], "k.", ms=3)
    ax.set_xlim(top[0] - 0.12, top[0] + 0.12)
    ax.set_ylim(top[1] - 0.09, top[1] + 0.03)
    ax.set_title("Front web meets the suction cap: solids (coloured bands)\n"
                 "and the conformal skin-mid graph (black); OFFSET downstream")
    ax.set_xlabel("x [m]")
    note(fig, f"NACA 2418, 2 m chord, {TWIST:g}° twist, 2D-F; skin {SKIN * 1e3:g} mm, caps {CAP * 1e3:g} mm, "
              f"webs {WEB * 1e3:g} mm. Euler's formula, normals, offset regions and cell areas checked.")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, __file__, "mid_surface_and_cells")


if __name__ == "__main__":
    cut = cut_section()
    ms = cut.mid_surface()
    check(cut, ms)
    figure(cut, ms)
