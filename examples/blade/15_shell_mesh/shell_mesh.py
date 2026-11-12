"""The kernel's shell mesh: all-quad skin and webs, polygon ribs and tip closure.

``mesh.shell_mesh`` builds one conforming sheet: mapped quad strips on the skin
and web mid-surfaces, and a polygon mesh (every cell at least five nodes) on
each rib face and the tip closure, whose perimeter nodes are the skin and web
nodes of its station. Each element carries its partition and shell data.

Figures: a made-up two-web blade with four ribs, sized by element edge ``h``,
and (if the reference-blades catalogue is cloned) NREL 5 MW.

    python examples/blade/15_shell_mesh/shell_mesh.py
"""

from __future__ import annotations

import runpy
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from matplotlib.colors import to_rgba  # noqa: E402

from _shared.plotting import NAVY, RED, TEAL, add_faces_3d, note, plt, save  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402
from implicit_geometry_kernel.mesh import ShellMesh, shell_mesh  # noqa: E402

HERE = Path(__file__).resolve().parent
REF = HERE.parents[0] / "08_reference_blades" / "reference_blades.py"
GREY = "#9ca3af"
EDGE = "#1f2937"


def _colour(m: ShellMesh, i: int) -> str:
    if m.surface[i] == "web":
        return TEAL
    if m.surface[i] in ("rib", "tip_closure"):
        return GREY
    return RED if m.eccentricity[i] > 0.0 else NAVY


def _faces(m: ShellMesh, keep) -> tuple[list[np.ndarray], list[str]]:
    faces, cols = [], []
    for i, e in enumerate(m.elements):
        if keep(i):
            faces.append(m.nodes[list(e)][:, [2, 0, 1]])
            cols.append(_colour(m, i))
    return faces, cols


def _draw(ax, faces, cols, *, elev, azim, lw=0.15, aspect=None):
    rgba = [to_rgba(c, 0.9) for c in cols]
    add_faces_3d(ax, faces, facecolors=rgba, edgecolors=to_rgba(EDGE, 0.6), linewidths=lw, shade=True)
    pts = np.concatenate(faces)
    box = np.ptp(pts, axis=0)
    box[box == 0.0] = 1.0
    ax.set_xlim(pts[:, 0].min(), pts[:, 0].max())
    ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
    ax.set_zlim(pts[:, 2].min(), pts[:, 2].max())
    ax.set_box_aspect(box if aspect is None else aspect)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_zlabel("y [m]")
    ax.tick_params(labelsize=7)
    ax.view_init(elev=elev, azim=azim)


def min_angles(m: ShellMesh) -> np.ndarray:
    """Smallest corner angle of every element (degrees)."""
    out = []
    for e in m.elements:
        p = m.nodes[list(e)]
        a = np.roll(p, 1, axis=0) - p
        b = np.roll(p, -1, axis=0) - p
        c = (a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        out.append(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))).min())
    return np.array(out)


def report(label: str, m: ShellMesh) -> str:
    use = Counter(m.edge_use().values())
    z0 = m.stations[0]
    free = [(a, b) for (a, b), n in m.edge_use().items() if n == 1]
    root = sum(abs(m.nodes[a, 2] - z0) < 1e-9 and abs(m.nodes[b, 2] - z0) < 1e-9 for a, b in free)
    ang = min_angles(m)
    quads = list(m.quads)
    text = (f"{label}: {len(m.nodes)} nodes, {len(m.quads)} quads, {len(m.polygons)} polygons "
            f"(fewest nodes {min(len(m.elements[i]) for i in m.polygons)}); edge use {dict(sorted(use.items()))}; "
            f"open edges {len(free)} ({root} on the root); quads with a corner under 20°: "
            f"{int((ang[quads] < 20).sum())}")
    print(text)
    return text


def figure(m: ShellMesh, stem: str, title: str, text: str, window: tuple[float, float],
           azim: float = -70) -> None:
    fig = plt.figure(figsize=(12, 5.2))
    ax = fig.add_subplot(121, projection="3d")
    faces, cols = _faces(m, lambda i: m.surface[i] != "web")
    _draw(ax, faces, cols, elev=22, azim=-58, lw=0.06, aspect=(4.0, 1.2, 0.6))
    ax.set_title(f"{title}: skin (navy), caps (red), ribs and tip (grey)", fontsize=9)
    ax = fig.add_subplot(122, projection="3d")
    lo, hi = window
    zc = {i: float(np.mean(m.nodes[list(e)][:, 2])) for i, e in enumerate(m.elements)}
    yc = {i: float(np.mean(m.nodes[list(e)][:, 1])) for i, e in enumerate(m.elements)}

    def keep(i):
        # suction skin left out to show the webs and the rib / tip polygons
        return lo <= zc[i] <= hi and not (m.surface[i] == "skin" and yc[i] > 0.0)

    faces, cols = _faces(m, keep)
    _draw(ax, faces, cols, elev=30, azim=azim, lw=0.3, aspect=(1.6, 1.4, 0.7))
    ax.set_title(f"z {lo:.1f}–{hi:.1f} m without the suction skin: webs (teal), polygons (grey)", fontsize=9)
    note(fig, text)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, __file__, stem)


def angle_figure(meshes: list[tuple[str, ShellMesh]]) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    bins = np.linspace(0, 90, 46)
    for (label, m), c in zip(meshes, (NAVY, TEAL)):
        ang = min_angles(m)[list(m.quads)]
        ax.hist(ang, bins=bins, histtype="step", color=c, label=label, density=True)
    ax.set_xlabel("smallest corner angle of a quad [°]")
    ax.set_ylabel("density (log)")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    note(fig, "Skin and web quads. The few under 20° are transition quads in event bands (bug log BG-034); "
              "span spacing near h (K8) removes them.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, __file__, "quad_angles")


def main() -> int:
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=8.0),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
                        spar_cap_thickness=0.008, shear_web_thickness=0.004)
    m = shell_mesh(geo, h=0.08)
    rib = geo.transverse[1].z
    figure(m, "parametric", "2D-F-LC-4R, h = 0.08 m", report("2D-F-LC-4R", m), (rib - 0.6, rib + 0.6))
    meshes = [("2D-F-LC-4R, h = 0.08 m", m)]
    ns = runpy.run_path(str(REF))
    root = ns["catalogue_root"]()
    if root is None:
        print("No reference-blades catalogue; skip NREL 5 MW (see example 08).")
    else:
        _name, nrel, _ = ns["load_windio"](ns["catalogue_yaml"](root, "nrel5mw/nrel5mw.yaml"))
        mr = shell_mesh(nrel, 41)
        L = float(nrel.stations[-1].z)
        figure(mr, "nrel5mw", "NREL 5 MW, 41 stations", report("NREL 5 MW", mr), (L - 4.0, L), azim=-40)
        meshes.append(("NREL 5 MW, 41 stations", mr))
    angle_figure(meshes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
