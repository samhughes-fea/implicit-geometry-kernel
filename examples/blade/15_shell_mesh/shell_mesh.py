"""The kernel's shell mesh: all-quad skin and webs, triangle ribs and tip closure.

``mesh.shell_mesh`` (the package function) builds one conforming sheet: mapped quad strips on the skin
and web mid-surfaces, and an unstructured Delaunay triangle mesh on each rib
face and the tip closure, whose perimeter nodes are the skin and web nodes of
its station. Each element carries its partition and shell data.

With the reference-blades catalogue this writes two blades, NREL 5 MW and the
DOE RM1 tidal rotor, each at the 5% chord default. A folder holds the span and
tip pictures, the check tables, and the VTK set.

    python examples/blade/15_shell_mesh/shell_mesh.py
"""

from __future__ import annotations

import csv
import importlib.util
import logging
import runpy
import sys
import textwrap
import time
import warnings
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from implicit_geometry_kernel.plotting.figure_io import resolve_save_dpi  # noqa: E402
from implicit_geometry_kernel.plotting.palette import panel_colour  # noqa: E402
from _shared.plotting import (  # noqa: E402
    MESH_EDGE, TEXT_SECONDARY, mesh_edge_rgba, note, plt,
)
from implicit_geometry_kernel.mesh import ShellMesh, shell_mesh, write_mesh_tables, write_vtk  # noqa: E402
from implicit_geometry_kernel.profiling import LOGGER_NAME, STAGE, collecting  # noqa: E402
from implicit_geometry_kernel.mesh.fill import triangulate_ring  # noqa: E402

HERE = Path(__file__).resolve().parent
REF = HERE.parents[0] / "08_reference_blades" / "reference_blades.py"


def _colour(m: ShellMesh, i: int) -> str:
    """Same retro role colour as example 09: stable for a strip or rib-cell name."""
    return panel_colour(m.name[i])


_PARTITION_TRIS: dict[int, tuple[dict[int, list[np.ndarray]], set[int]]] = {}


def _partition_edge_size(m: ShellMesh, cells: list) -> dict[tuple[float, str], float]:
    """One spacing per rib or tip: the median length of that station's partition edges."""
    buckets: dict[tuple[float, str], list[np.ndarray]] = {}
    for z, surface, _name, _indices, outers, _holes in cells:
        for loop in outers:
            p = np.asarray(m.nodes[list(loop)], float)
            buckets.setdefault((z, surface), []).append(np.linalg.norm(np.roll(p, -1, 0) - p, axis=1))
    sizes = {}
    for key, parts in buckets.items():
        d = np.concatenate(parts)
        d = d[d > 1e-12]
        if len(d):
            sizes[key] = float(np.median(d))
    return sizes


def _partition_triangles(m: ShellMesh) -> tuple[dict[int, list[np.ndarray]], set[int]]:
    """Unstructured Delaunay fill of every rib and tip partition.

    The boundary vertices are the skin and web nodes of that station. The
    interior is filled at the median length of those edges. Sibling elements of
    the same cell are in the returned skip set so the cell is drawn once.
    """
    hit = _PARTITION_TRIS.get(id(m))
    if hit is not None:
        return hit
    emit: dict[int, list[np.ndarray]] = {}
    skip: set[int] = set()
    cells, _rows = _polygon_faces(m, None)
    sizes = _partition_edge_size(m, cells)
    for (z, surface, _name, indices, outers, _holes) in cells:
        tris: list[np.ndarray] = []
        h = sizes.get((z, surface))
        for loop in outers:
            tris.extend(triangulate_ring(
                m.nodes[list(loop)], size=h, pattern="unstructured", split_boundary=False,
            ))
        if not tris or not indices:
            continue
        emit[indices[0]] = tris
        skip.update(indices[1:])
    packed = (emit, skip)
    _PARTITION_TRIS[id(m)] = packed
    return packed


def _plot_xyz(xyz: np.ndarray) -> np.ndarray:
    """Stored ``(x, y, z)`` nodes as plot ``(z, x, y)``."""
    return np.asarray(xyz, float)[:, [2, 0, 1]]


def _stored_faces(m: ShellMesh, keep) -> tuple[np.ndarray | list[np.ndarray], list[str]]:
    """Stored skin and web quads, and the stored rib and tip triangles."""
    faces: list[np.ndarray] = []
    cols: list[str] = []
    nodes = m.nodes
    for i, elem in enumerate(m.elements):
        if not keep(i) or len(elem) < 3:
            continue
        faces.append(_plot_xyz(nodes[list(elem)]))
        cols.append(_colour(m, i))
    if faces and all(len(face) == 4 for face in faces):
        return np.stack(faces, axis=0), cols
    return faces, cols


# Each element edge is this many pixels at DPI, and the shorter side of a panel
# is at least MIN_SHORT_PX. The page is the equal-metre view, so those pixels
# fall on the quads and the ribs rather than on empty margin.
EDGE_PX = 24
MIN_SHORT_PX = 480
DPI = 300
_MARGIN_L = 0.12
_MARGIN_R = 0.08
_TITLE_IN = 0.36
_BAR_IN = 0.42
_GAP = 0.32


def _screen_basis(elev: float, azim: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Right, up, and towards the viewer, in plot metres.

    The eye matches mplot3d (elevation off the horizontal, azimuth in the plot
    plane). The vectors stay in metres, so the page is not letterboxed into a square.
    """
    elev_rad = np.deg2rad(float(elev))
    azim_rad = np.deg2rad(float(azim))
    towards = np.array([
        np.cos(elev_rad) * np.cos(azim_rad),
        np.cos(elev_rad) * np.sin(azim_rad),
        np.sin(elev_rad),
    ])
    towards = towards / np.linalg.norm(towards)
    right = np.cross(np.array([0.0, 0.0, 1.0]), towards)
    right = right / np.linalg.norm(right)
    up = np.cross(towards, right)
    return right, up, towards


def _as_faces(faces) -> list[np.ndarray]:
    if isinstance(faces, np.ndarray) and faces.ndim == 3:
        return [np.asarray(face, float) for face in faces]
    return [np.asarray(face, float) for face in faces]


def _project_view(faces, cols, right, up, towards) -> tuple[list[np.ndarray], list, np.ndarray]:
    """Screen polygons, far to near. A metre on the screen is a metre on the blade."""
    sheets = _as_faces(faces)
    if not sheets:
        return [], [], np.zeros((0, 2))
    depth = np.fromiter((float(face.mean(axis=0) @ towards) for face in sheets), float, count=len(sheets))
    order = np.argsort(depth)
    polygons: list[np.ndarray] = []
    colours: list = []
    for index in order:
        face = sheets[int(index)]
        polygons.append(np.column_stack((face @ right, face @ up)))
        colours.append(cols[int(index)])
    return polygons, colours, np.concatenate(polygons)


def _align_span(polygons: list[np.ndarray], cloud: np.ndarray, right, up):
    """Put the span along the page so the frame is the blade, not its diagonal box."""
    if len(cloud) == 0:
        return polygons, cloud
    centre = cloud.mean(axis=0)
    delta = cloud - centre
    _values, vectors = np.linalg.eigh(delta.T @ delta)
    axis = vectors[:, -1]
    span_on_screen = np.array([right[0], up[0]])
    if float(axis @ span_on_screen) < 0.0:
        axis = -axis
    perp = np.array([-axis[1], axis[0]])

    def turn(xy: np.ndarray) -> np.ndarray:
        shifted = xy - centre
        return np.column_stack((shifted @ axis, shifted @ perp))

    return [turn(poly) for poly in polygons], turn(cloud)


def _px_per_m(m: ShellMesh) -> float:
    """Pixels per metre so the smallest recipe edge is ``EDGE_PX`` pixels."""
    edge = min(m.recipe.edge) if m.recipe is not None and m.recipe.edge else 0.1
    return EDGE_PX / max(float(edge), 1e-3)


def _limits(cloud: np.ndarray) -> tuple[float, float, float, float]:
    span = np.ptp(cloud, axis=0)
    pad = np.maximum(0.03 * span, 1e-6)
    return (
        float(cloud[:, 0].min() - pad[0]), float(cloud[:, 0].max() + pad[0]),
        float(cloud[:, 1].min() - pad[1]), float(cloud[:, 1].max() + pad[1]),
    )


def _draw_ortho(ax, polygons, cols, limits) -> None:
    """One equal-metre panel. The axes rectangle already has the data's aspect."""
    if polygons:
        ax.add_collection(PolyCollection(
            polygons,
            facecolors=[to_rgba(c, 0.92) for c in cols],
            edgecolors=mesh_edge_rgba(len(polygons)),
            linewidths=0.4,
        ))
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _bar_length(limits) -> float:
    width = limits[1] - limits[0]
    length = 10.0 if width > 25.0 else 1.0
    if length > 0.45 * width:
        length = float(10 ** np.floor(np.log10(max(0.4 * width, 1e-6))))
    return length


def _scale_bar(fig, ax, length_m: float) -> None:
    pos = ax.get_position()
    width = ax.get_xlim()[1] - ax.get_xlim()[0]
    frac = (length_m / width) * pos.width
    y = pos.y0 - 0.16 / fig.get_figheight()
    fig.add_artist(Line2D(
        [pos.x0, pos.x0 + frac], [y, y],
        transform=fig.transFigure, color="black", lw=1.15, solid_capstyle="butt",
    ))
    fig.text(
        pos.x0 + 0.5 * frac, y - 0.01, f"{length_m:g} m",
        ha="center", va="top", fontsize=8, color="black", transform=fig.transFigure,
    )


def _note_block(width_in: float, text: str) -> tuple[str, float]:
    """Wrap a caption to the figure width. Returns the text and its height in inches."""
    cols = max(64, int(width_in * 13.5))
    wrapped = "\n".join(textwrap.wrap(text, width=cols))
    lines = wrapped.count("\n") + 1
    return wrapped, lines * (7.5 * 1.35) / 72 + 0.08


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
    text = (f"{label}: {len(m.nodes)} nodes, {len(m.quads)} quads, {len(m.triangles)} triangles; "
            f"edge use {dict(sorted(use.items()))}; "
            f"open edges {len(free)} ({root} on the root); quads with a corner under 20°: "
            f"{int((ang[quads] < 20).sum())}; skew max {float(m.element_quality['skew'].max()):.3f}; "
            f"volume ratio max {float(m.element_quality['volume_ratio'].max()):.2f}")
    print(text)
    return text


def _save(fig, name: str, *, folder: str | None = None, publish: bool = False) -> None:
    """Write ``figures/<name>.png``, or ``figures/<folder>/<name>.png`` when ``folder`` is set.

    Iteration saves use the package default DPI and a fixed pad. ``publish=True``
    writes 300 dpi with a tight box.
    """
    dpi = 300 if publish else resolve_save_dpi(presentation=False)
    out = HERE / "figures" if folder is None else HERE / "figures" / folder
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.png"
    if publish:
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.04)
    else:
        fig.savefig(path, dpi=dpi, pad_inches=0.15)
    plt.close(fig)
    print(f"wrote {path.relative_to(HERE.parents[1])}")


def _write_views(panels, caption: str, px_per_m: float, name: str, *,
                 folder: str | None, publish: bool, stack: str = "row") -> None:
    """One file. Each panel's rectangle is the projected mesh at ``px_per_m``.

    ``stack="column"`` places the panels one under the other, sharing the page
    width, so a leading-edge strip and a trailing-edge strip share the span.
    """
    inches = [
        ((limits[1] - limits[0]) * px_per_m / DPI, (limits[3] - limits[2]) * px_per_m / DPI)
        for _heading, _polys, _cols, limits in panels
    ]
    wrapped, cap_h = _note_block(
        _MARGIN_L + (max(w for w, _h in inches) if stack == "column" else sum(w for w, _h in inches)) + _MARGIN_R,
        caption,
    )
    if stack == "column":
        width = max(w for w, _h in inches)
        fig_w = _MARGIN_L + width + _MARGIN_R
        fig_h = cap_h + 0.06 + _BAR_IN + sum(h for _w, h in inches) + _TITLE_IN * len(panels) + 0.08 * (len(panels) - 1)
        fig = plt.figure(figsize=(fig_w, fig_h))
        y_in = fig_h - _TITLE_IN
        for index, ((heading, polygons, cols, limits), (_width, height)) in enumerate(zip(panels, inches)):
            y_in -= height
            ax = fig.add_axes([_MARGIN_L / fig_w, y_in / fig_h, width / fig_w, height / fig_h])
            _draw_ortho(ax, polygons, cols, limits)
            ax.set_title(heading, fontsize=9, pad=2)
            if index == len(panels) - 1:
                _scale_bar(fig, ax, _bar_length(limits))
            print(f"  {heading[:48]}: {width * DPI:.0f} x {height * DPI:.0f} px, {px_per_m:.0f} px/m")
            y_in -= _TITLE_IN + 0.08
    else:
        plot_h = max(height for _width, height in inches)
        fig_w = _MARGIN_L + sum(width for width, _height in inches) + _GAP * (len(panels) - 1) + _MARGIN_R
        fig_h = cap_h + 0.06 + _BAR_IN + plot_h + _TITLE_IN
        fig = plt.figure(figsize=(fig_w, fig_h))
        x_in = _MARGIN_L
        axes_bottom = cap_h + 0.06 + _BAR_IN
        for (heading, polygons, cols, limits), (width, height) in zip(panels, inches):
            ax = fig.add_axes([
                x_in / fig_w, axes_bottom / fig_h, width / fig_w, height / fig_h,
            ])
            _draw_ortho(ax, polygons, cols, limits)
            ax.set_title(heading, fontsize=9, pad=2)
            _scale_bar(fig, ax, _bar_length(limits))
            print(f"  {heading[:48]}: {width * DPI:.0f} x {height * DPI:.0f} px, {px_per_m:.0f} px/m")
            x_in += width + _GAP
    fig.text(0.012, 0.04 / fig.get_figheight(), wrapped, fontsize=7.5, color=TEXT_SECONDARY,
             ha="left", va="bottom", linespacing=1.25)
    _save(fig, name, folder=folder, publish=publish)


def figure(m: ShellMesh, stem: str, title: str, text: str, window: tuple[float, float],
           azim: float = -70, *, folder: str | None = None,
           split: bool = False, publish: bool = False) -> None:
    """Two equal-metre views. ``split`` writes each view to its own file.

    The camera is the usual elevation and azimuth. The span is then laid along
    the page and the figure is sized so one pixel is the same length in both
    directions. A square 3D axes would spend the page on the empty box.
    """
    full_faces, full_cols = _stored_faces(m, lambda i: m.surface[i] != "web")
    lo, hi = window
    zc = {i: float(np.mean(m.nodes[list(e)][:, 2])) for i, e in enumerate(m.elements)}

    def keep(i):
        # Suction skin is left out so the webs and the rib / tip cells stay visible.
        # The side is the partition name: a cambered trailing edge sits above the chord
        # on the pressure surface, and a tip centroid can sit a float past the station.
        suction = m.surface[i] == "skin" and "suction" in m.partition[i]
        return lo - 1e-6 <= zc[i] <= hi + 1e-6 and not suction

    close_faces, close_cols = _stored_faces(m, keep)
    caption = (text + " Skin and webs are stored quads. Ribs and the tip closure are "
               "the stored triangles. Equal metres, span along the page; the bar is the scale.")
    # The span camera sits on the leading-edge side. The opposite azimuth is the
    # trailing-edge side of the same skin, still from above.
    views = (
        (f"{title}: leading edge", full_faces, full_cols, 22.0, -58.0),
        (f"{title}: trailing edge", full_faces, full_cols, 22.0, -58.0 + 180.0),
        (f"z {lo:.1f}–{hi:.1f} m without the suction skin", close_faces, close_cols, 30.0, azim),
    )
    px_per_m = _px_per_m(m)
    projected = []
    for heading, faces, cols, elev, view_azim in views:
        right, up, towards = _screen_basis(elev, view_azim)
        polygons, colours, cloud = _project_view(faces, cols, right, up, towards)
        polygons, cloud = _align_span(polygons, cloud, right, up)
        projected.append((heading, polygons, colours, cloud))
    # Root at the same page position in both lengthwise views, then one shared
    # span range so the tip meets the tip.
    _le_heading, le_polygons, le_colours, le_cloud = projected[0]
    te_heading, te_polygons, te_colours, te_cloud = projected[1]
    # The opposite camera puts thickness toward the bottom of the page. Flip that
    # axis so the same edge stays up, then put the root on the same page position.
    te_polygons = [np.column_stack((poly[:, 0], -poly[:, 1])) for poly in te_polygons]
    te_cloud = te_cloud.copy()
    te_cloud[:, 1] *= -1.0
    shift = float(le_cloud[:, 0].min() - te_cloud[:, 0].min())
    te_polygons = [poly + (shift, 0.0) for poly in te_polygons]
    te_cloud[:, 0] += shift
    le_own = _limits(le_cloud)
    te_own = _limits(te_cloud)
    x0 = min(le_own[0], te_own[0])
    x1 = max(le_own[1], te_own[1])
    le_limits = (x0, x1, le_own[2], le_own[3])
    te_limits = (x0, x1, te_own[2], te_own[3])
    tip_heading, tip_polygons, tip_colours, tip_cloud = projected[2]
    tip_limits = _limits(tip_cloud)
    limits_for_scale = (le_limits, te_limits, tip_limits)
    for limits in limits_for_scale:
        short = min(limits[1] - limits[0], limits[3] - limits[2])
        px_per_m = max(px_per_m, MIN_SHORT_PX / max(short, 1e-6))
    if split:
        span_caption = caption + " Leading edge above, trailing edge below, root to the left."
        _write_views(
            [(_le_heading, le_polygons, le_colours, le_limits),
             (te_heading, te_polygons, te_colours, te_limits)],
            span_caption, px_per_m, f"{stem}_span", folder=folder, publish=publish, stack="column",
        )
        _write_views(
            [(tip_heading, tip_polygons, tip_colours, tip_limits)],
            caption, px_per_m, f"{stem}_tip", folder=folder, publish=publish,
        )
        return
    _write_views(
        [(f"{title}: each partition in the retro palette", le_polygons, le_colours, le_own),
         (tip_heading, tip_polygons, tip_colours, tip_limits)],
        caption, px_per_m, stem, folder=folder, publish=publish,
    )


def angle_figure(meshes: list[tuple[str, ShellMesh]], *, folder: str | None = None,
                 publish: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    bins = np.linspace(0, 90, 46)
    for label, m in meshes:
        ang = min_angles(m)[list(m.quads)]
        ax.hist(ang, bins=bins, histtype="step", color=panel_colour(label), label=label, density=True)
    ax.set_xlabel("smallest corner angle of a quad [°]")
    ax.set_ylabel("density (log)")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    note(fig, "Skin and web quads. Rib and tip partitions stay triangles and are not in this histogram. "
              "The few under 20° are transition quads in event bands (bug log BG-034); "
              "span spacing near h (K8) removes them.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    _save(fig, "quad_angles", folder=folder, publish=publish)


def _polygon_area(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _boundary_loops(m: ShellMesh, indices: list[int]) -> list[list[int]]:
    """Closed chains of edges that belong to one element of ``indices``."""
    count: dict[tuple[int, int], int] = {}
    direction: dict[tuple[int, int], tuple[int, int]] = {}
    for i in indices:
        elem = m.elements[i]
        for a, b in zip(elem, elem[1:] + elem[:1]):
            key = (a, b) if a < b else (b, a)
            count[key] = count.get(key, 0) + 1
            direction[key] = (a, b)
    succ: dict[int, list[int]] = {}
    for key, n in count.items():
        if n == 1:
            a, b = direction[key]
            succ.setdefault(a, []).append(b)
    pending = {(a, b) for a, nxt in succ.items() for b in nxt}
    loops: list[list[int]] = []
    while pending:
        start, b = pending.pop()
        loop = [start]
        while b != start:
            loop.append(b)
            choices = [c for c in succ.get(b, ()) if (b, c) in pending]
            if not choices:
                loop = []
                break
            c = choices[0]
            pending.remove((b, c))
            b = c
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def _section_xy(m: ShellMesh, ids) -> np.ndarray:
    return np.asarray(m.nodes[list(ids)][:, :2], float)


def _polygon_faces(m: ShellMesh, window: tuple[float, float] | None):
    """Rib and tip cells: element indices, outer loops, and any interior loops."""
    buckets: dict[tuple[float, str, str], list[int]] = {}
    for i, elem in enumerate(m.elements):
        surface = m.surface[i]
        if surface not in ("rib", "tip_closure"):
            continue
        z = round(float(np.mean(m.nodes[list(elem)][:, 2])), 6)
        if window is not None and not (window[0] - 1e-6 <= z <= window[1] + 1e-6):
            continue
        buckets.setdefault((z, surface, m.name[i]), []).append(i)
    cells = []
    rows: dict[tuple[float, str], list[int]] = {}
    for (z, surface, name), indices in sorted(buckets.items()):
        loops = _boundary_loops(m, indices)
        if not loops:
            continue
        areas = [_polygon_area(_section_xy(m, loop)) for loop in loops]
        outer_i = int(np.argmax(np.abs(areas)))
        sign = np.sign(areas[outer_i]) or 1.0
        outers = [loop for loop, area in zip(loops, areas) if area * sign > 1e-14]
        holes = [loop for loop, area in zip(loops, areas) if area * sign < -1e-14]
        cells.append((z, surface, name, indices, outers, holes))
        rows.setdefault((z, surface), []).append(len(cells) - 1)
    return cells, rows


def _cover(faces, ring_area: float) -> float:
    if ring_area <= 1e-18:
        return 1.0
    total = 0.0
    for face in faces:
        pts = np.asarray(face, float)
        total += abs(_polygon_area(pts[:, :2]))
    return total / ring_area


def _skin_node_ids(m: ShellMesh) -> set[int]:
    """Nodes of the skin midsurface. A rib meets the shell on this set."""
    ids: set[int] = set()
    for i, elem in enumerate(m.elements):
        if m.surface[i] == "skin":
            ids.update(int(n) for n in elem)
    return ids


def rib_compare(m: ShellMesh, stem: str, title: str, window: tuple[float, float] | None = None,
                *, folder: str | None = None, publish: bool = False) -> None:
    """Section view of each rib and tip partition beside a Gmsh fill of the same ring.

    The kernel panel is the unstructured Delaunay fill of that partition. Gmsh
    triangulates the same outer boundary (Frontal-Delaunay). White dots are the
    nodes the partition shares with the skin midsurface; the rib edge triangles
    are built from those nodes.
    """
    cells, rows = _polygon_faces(m, window)
    emit, _skip = _partition_triangles(m)
    skin_ids = _skin_node_ids(m)
    if not rows:
        print(f"{title}: no rib or tip faces to compare")
        return
    if importlib.util.find_spec("gmsh") is None:
        print("gmsh not installed; skip rib comparison")
        return
    from implicit_geometry_kernel.mesh.gmsh.plane import triangulate_rings

    rings: list[np.ndarray] = []
    owner: list[int] = []
    for i, (_z, _surface, _name, _indices, outers, _holes) in enumerate(cells):
        for loop in outers:
            rings.append(m.nodes[list(loop)])
            owner.append(i)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        filled = triangulate_rings(rings)
    tris: list[list[np.ndarray]] = [[] for _ in cells]
    for i, faces in enumerate(filled):
        tris[owner[i]].extend(faces)

    keys = sorted(rows)
    fig, axes = plt.subplots(len(keys), 2, figsize=(11.0, 1.72 * len(keys)), squeeze=False)
    kernel_cover: list[float] = []
    gmsh_cover: list[float] = []
    n_holes = 0
    for row, key in enumerate(keys):
        z, surface = key
        ids = rows[key]
        kind = "tip closure" if surface == "tip_closure" else "rib"
        pts = [_section_xy(m, m.elements[j]) for i in ids for j in cells[i][3]]
        cloud = np.concatenate([p for p in pts if len(p)])
        pad = 0.03 * max(float(np.ptp(cloud, axis=0).max()), 1e-6)
        limits = (cloud[:, 0].min() - pad, cloud[:, 0].max() + pad,
                  cloud[:, 1].min() - pad, cloud[:, 1].max() + pad)
        for col, src in enumerate(("kernel", "gmsh")):
            ax = axes[row, col]
            faces, colours = [], []
            for i in ids:
                _z, _surface, _name, indices, outers, holes = cells[i]
                colour = to_rgba(panel_colour(_name), 0.92)
                if src == "kernel":
                    owned = emit.get(indices[0], [])
                    batch = [np.asarray(t, float)[:, :2] for t in owned] or [
                        _section_xy(m, m.elements[j]) for j in indices]
                    ring_area = sum(abs(_polygon_area(_section_xy(m, loop))) for loop in outers)
                    kernel_cover.append(_cover(batch, ring_area))
                    n_holes += len(holes)
                else:
                    batch = [np.asarray(t, float)[:, :2] for t in tris[i]]
                    ring_area = sum(abs(_polygon_area(_section_xy(m, loop))) for loop in outers)
                    gmsh_cover.append(_cover(batch, ring_area))
                faces.extend(batch)
                colours.extend([colour] * len(batch))
            if faces:
                ax.add_collection(PolyCollection(
                    faces, facecolors=colours, edgecolors=mesh_edge_rgba(len(faces)), linewidths=0.35))
            for i in ids:
                for loop in cells[i][4]:
                    xy = _section_xy(m, list(loop) + [loop[0]])
                    ax.plot(xy[:, 0], xy[:, 1], color="black", lw=0.7, solid_capstyle="round")
                    shared = [n for n in loop if int(n) in skin_ids]
                    if shared:
                        nodes = _section_xy(m, shared)
                        ax.scatter(nodes[:, 0], nodes[:, 1], s=8, c="white",
                                   edgecolors=MESH_EDGE, linewidths=0.45, zorder=5)
                if src == "kernel":
                    for loop in cells[i][5]:
                        xy = _section_xy(m, list(loop) + [loop[0]])
                        ax.plot(xy[:, 0], xy[:, 1], color=MESH_EDGE, lw=1.1)
            ax.set_xlim(limits[0], limits[1])
            ax.set_ylim(limits[2], limits[3])
            ax.set_aspect("equal", adjustable="box")
            ax.tick_params(labelsize=7)
            ax.set_xlabel("x [m]")
            if col == 0:
                ax.set_ylabel("y [m]")
            side = "unstructured Delaunay" if src == "kernel" else "Gmsh triangles"
            ax.set_title(f"{kind}, z = {z:.2f} m — {side}", fontsize=9)
    def band(lo: float, hi: float) -> str:
        return f"{lo:.4f}" if hi - lo < 5e-5 else f"{lo:.4f}–{hi:.4f}"

    k0, k1 = min(kernel_cover), max(kernel_cover)
    g0, g1 = min(gmsh_cover), max(gmsh_cover)
    text = (f"{title}: unstructured Delaunay and a Gmsh Frontal-Delaunay fill of the same outer boundary. "
            f"White dots are the nodes shared with the skin midsurface; the rib edge triangles are built from them. "
            f"Covered area / boundary area {band(k0, k1)} (kernel), {band(g0, g1)} (Gmsh). "
            f"Interior boundary loops: {n_holes}.")
    print(text)
    note(fig, text)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    _save(fig, stem, folder=folder, publish=publish)


def _write_timings(directory: Path, reports) -> None:
    """One row per mesh step, beside the check tables."""
    path = directory / "timings.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["stage", "seconds", "cpu_seconds", "depth", "detail"],
        )
        writer.writeheader()
        for report in reports:
            detail = " ".join(f"{key}={value}" for key, value in report.measures.items())
            writer.writerow({
                "stage": report.stage,
                "seconds": f"{report.seconds:.6f}",
                "cpu_seconds": f"{report.cpu_seconds:.6f}",
                "depth": report.depth,
                "detail": detail,
            })
    print(f"wrote {path.relative_to(HERE.parents[1])}")


class _Echo(logging.Handler):
    """Print each mesh step as it finishes, so a long blade is not silent."""

    def emit(self, record: logging.LogRecord) -> None:
        name = getattr(record, "stage_name", None)
        seconds = getattr(record, "stage_seconds", None)
        if not isinstance(name, str) or not str(name).startswith("shell_"):
            return
        if isinstance(seconds, bool) or not isinstance(seconds, (int, float)):
            return
        depth = getattr(record, "stage_depth", 0)
        indent = "  " * (int(depth) if isinstance(depth, int) else 0)
        print(f"{indent}{name} {float(seconds):.1f}s", flush=True)


def _publish(mesh: ShellMesh, reports, stem: str, label: str, folder: str, length: float) -> None:
    """Span and tip pictures, the check tables, the step times, and the VTK set."""
    out = HERE / "figures" / folder
    write_mesh_tables(mesh, out)
    _write_timings(out, reports)
    write_vtk(mesh, out / f"{stem}.vtk", label)
    figure(
        mesh, stem, label, report(label, mesh), (length - 4.0, length),
        azim=-40, folder=folder, split=True, publish=True,
    )


def main() -> int:
    ns = runpy.run_path(str(REF))
    root = ns["catalogue_root"]()
    if root is None:
        print("No reference-blades catalogue; skip NREL 5 MW and RM1 (see example 08).")
        return 0
    blades = (
        ("nrel5mw/nrel5mw.yaml", "nrel_5mw", "nrel5mw", "NREL 5 MW, tip h/c = 5%"),
        ("tidal/rm1/RM1_Bir2011_web-ud_core.windio.yaml", "rm1", "rm1", "DOE RM1 tidal, tip h/c = 5%"),
    )
    for rel, folder, stem, label in blades:
        _name, geo, _z = ns["load_windio"](ns["catalogue_yaml"](root, rel))
        log = logging.getLogger(LOGGER_NAME)
        echo = _Echo(level=STAGE)
        log.addHandler(echo)
        t0 = time.perf_counter()
        try:
            with collecting() as handler:
                mesh = shell_mesh(geo)
            reports = handler.take()
        finally:
            log.removeHandler(echo)
        shell_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        _publish(mesh, reports, stem, label, folder, float(geo.stations[-1].z))
        print(f"{stem} shell_mesh {shell_s:.3f}s  figure {time.perf_counter() - t0:.3f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
