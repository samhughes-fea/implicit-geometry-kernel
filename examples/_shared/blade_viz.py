"""Shared figure tree, suction filter, loft live-web fracs, and planform sampling."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

from _shared.plotting import (
    NAVY, RED, SERIES, TEAL, draw_cap_solids, draw_shell, note, panel_colour, plt,
)
from implicit_geometry_kernel.blade.topology.network import wall_meta_from_cut
from implicit_geometry_kernel.blade.topology import (
    ShellTopologyTemplate, _chord_fracs, _cut_tc, _template_from_cut, role_atom,
    sample_region_polylines,
)
from implicit_geometry_kernel.contract import InfeasibleSectionGeometry
from implicit_geometry_kernel.integration.quadtree import build_shared_quadtree, leaf_bounds
from implicit_geometry_kernel.integration.quadtree.quadrature import (
    label_boundary_segments,
    label_fill_polygons,
)
from implicit_geometry_kernel.plotting.palette import region_color
from implicit_geometry_kernel.plotting.plot_cell_wall_loops import LOOP_STYLES
from implicit_geometry_kernel.section.boundary.cores import emit_constructive_polygons
from implicit_geometry_kernel.section.boundary.labels import (
    is_core_label,
    is_residual_label,
    is_wall_solid_label,
)
from implicit_geometry_kernel.section.boundary.polygon_ops import oriented_rings
from implicit_geometry_kernel.section.properties.backend import quadtree_spec_for_airfoil
from implicit_geometry_kernel.section.properties.section_evaluator import SectionFieldEvaluator

CROSS_SECTIONS = "cross_sections"
CSG_MESH = "csg_mesh"
REGIONS = "regions"
REGION_CUTAWAYS = "region_cutaways"
PLANFORM = "planform"
FULL_BLADE = "full_blade"
VIZ_DIRS = (CROSS_SECTIONS, CSG_MESH, REGIONS, REGION_CUTAWAYS, PLANFORM, FULL_BLADE)
N_SECTION_GRID = 16
_OCC_EXTERIOR, _OCC_CAVITY = "#f4f4f2", "#c5d4e0"
_CSG_QT_DEPTH = 10
_CSG_QT_MAX_CELLS = 150_000
_CSG_FIG_W = 15.6
_CSG_FIG_H_PAD = 2.2
_CSG_FIG_H_MIN = 4.0
_CSG_N_COL = 3
_CORE_ALPHA = 0.22
_WALL_ALPHA = 0.5
_BAY_LOOP_DASH = (0, (8.0, 3.5))
_UNION_KEY = "_union_walls"
VIEW_LE = (26.0, -62.0)
VIEW_TE = (26.0, 118.0)
_KEEP_ROOT = frozenset({
    "mids.obj", "mids.vtk", "cutaway.obj", "cutaway.vtk",
    "samples.npz", "oml.msh", "mesh_closeup.png",
})


def is_rib_role(name: str) -> bool:
    key = str(name).lower()
    return "rib" in key or key.startswith("tip")


def is_suction_role(name: str) -> bool:
    """True for suction-side loop strips; false for ribs, tip plates, webs, pressure."""
    key = str(name).lower()
    if "rib" in key or key.startswith("tip"):
        return False
    atom = role_atom(name)
    if len(atom) >= 2 and str(atom[1]) in ("su", "suction"):
        return True
    return key.endswith("_su") or "_su_" in key or "suction" in key


def live_web_fracs(cut) -> list[float]:
    """Twist-aware web chord fractions, matching ``templates_along``."""
    live = []
    for run in cut.shell_graph().webs:
        pts = np.asarray(run.points, float)[:, :2]
        live.append(float(np.clip(_chord_fracs(cut, pts[:1])[0], 0.0, 1.0)))
    return live


def span_section_zs(geo, n: int = N_SECTION_GRID) -> np.ndarray:
    """Equally spaced z from root to tip (inclusive)."""
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    return np.linspace(z0, z1, int(n))


def _section_xy_bounds(cut) -> tuple[float, float, float, float]:
    xs, ys = [cut.offset[0]], [cut.offset[1]]
    for wall in cut.walls:
        xy = np.asarray(wall.xy, float)
        if len(xy):
            xs.append(xy[:, 0])
            ys.append(xy[:, 1])
    x = np.concatenate([np.atleast_1d(v) for v in xs])
    y = np.concatenate([np.atleast_1d(v) for v in ys])
    return float(x.min()), float(x.max()), float(y.min()), float(y.max())


def _root_anchored_limits(cuts, pad_frac: float = 0.04) -> tuple[tuple[float, float], tuple[float, float]]:
    """Shared blade-frame window: root box, expanded if a later station sticks out."""
    boxes = [_section_xy_bounds(c) for c in cuts if c is not None]
    if not boxes:
        return (-1.0, 1.0), (-1.0, 1.0)
    x0, x1, y0, y1 = boxes[0]
    for xa, xb, ya, yb in boxes[1:]:
        x0, x1 = min(x0, xa), max(x1, xb)
        y0, y1 = min(y0, ya), max(y1, yb)
    pad = pad_frac * max(x1 - x0, y1 - y0, 1e-6)
    return (x0 - pad, x1 + pad), (y0 - pad, y1 + pad)


def _station_limits(cut, pad_frac: float = 0.04) -> tuple[tuple[float, float], tuple[float, float]]:
    x0, x1, y0, y1 = _section_xy_bounds(cut)
    pad = pad_frac * max(x1 - x0, y1 - y0, 1e-6)
    return (x0 - pad, x1 + pad), (y0 - pad, y1 + pad)


def _paint_section(ax, cut, z, xlim, ylim, *, show_axes: bool = False,
                   window: str | None = None) -> None:
    if cut is None:
        ax.set_title(f"z = {z:.1f} m\ninfeasible", fontsize=8, pad=6)
    else:
        draw_shell(ax, cut)
        draw_cap_solids(ax, cut)
        ax.plot(cut.offset[0], cut.offset[1], "k+", ms=5)
        extra = f"\n{window}" if window else ""
        ax.set_title(f"z = {z:.1f} m\nchord {cut.chord:.2f} m{extra}", fontsize=8, pad=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    if show_axes:
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        ax.tick_params(labelsize=7)
    else:
        ax.axis("off")


def _sample_cuts(geo):
    zs = span_section_zs(geo)
    cuts: list = []
    errors: list = []
    for z in zs:
        try:
            cuts.append(geo.section_at(float(z)))
            errors.append(None)
        except InfeasibleSectionGeometry as exc:
            cuts.append(None)
            errors.append(exc)
    return zs, cuts, errors


def figure_section_grid(stem: str, title: str, geo, *, figures: Path, repo: Path,
                        caption: str) -> None:
    """4×4 ``sections.png``, dual-scale ``01.png``–``16.png``, and 2×3 ``csg_mesh`` singles."""
    zs, cuts, errors = _sample_cuts(geo)
    xlim, ylim = _root_anchored_limits(cuts)
    fig, axes = plt.subplots(4, 4, figsize=(12.0, 13.2))
    fig.subplots_adjust(left=0.04, right=0.99, top=0.93, bottom=0.05, wspace=0.18, hspace=0.42)
    print(f"{'station':>12} {'z [m]':>8} {'chord':>8} {'twist':>8} {'walls':>6} {'cells':>6}  layout")
    for i, (ax, z, cut, exc) in enumerate(zip(axes.ravel(), zs, cuts, errors)):
        label = f"{i + 1}/{len(zs)}"
        _paint_section(ax, cut, z, xlim, ylim)
        if cut is None:
            print(f"{label:>12} {z:8.2f}  skipped ({exc})")
            continue
        dropped = f", dropped {len(cut.dropped_webs)}" if cut.dropped_webs else ""
        print(f"{label:>12} {z:8.2f} {cut.chord:8.2f} {cut.twist:8.2f} {len(cut.walls):6d} "
              f"{len(cut.cells):6d}  {cut.layout}{dropped}")
        if i == len(zs) - 1:
            kinds = sorted({w.kind for w in cut.walls})
            thick = {k: max((w.thickness for w in cut.walls if w.kind == k), default=0.0) for k in kinds}
            print(f"{'':>12} tip walls {kinds}  t={thick}  dropped={cut.dropped_webs}")
    fig.suptitle(f"{title} — 4×4 span cuts (skin / cap / web)", fontsize=10)
    note(fig, caption)
    save_fig(fig, figures, stem, CROSS_SECTIONS, "sections", repo=repo)
    for i, (z, cut) in enumerate(zip(zs, cuts), start=1):
        one, (ax_root, ax_st) = plt.subplots(1, 2, figsize=(8.4, 4.2))
        st_xlim, st_ylim = (xlim, ylim) if cut is None else _station_limits(cut)
        _paint_section(ax_root, cut, z, xlim, ylim, show_axes=True, window="root window")
        _paint_section(ax_st, cut, z, st_xlim, st_ylim, show_axes=True, window="station window")
        one.tight_layout()
        save_fig(one, figures, stem, CROSS_SECTIONS, f"{i:02d}", repo=repo)
    figure_csg_mesh(stem, geo, zs, cuts, figures=figures, repo=repo)


def _csg_mesh_figsize(xlim, ylim) -> tuple[float, float]:
    """2×3 page: fixed width, height from two equal-aspect station rows."""
    dx = max(float(xlim[1] - xlim[0]), 1e-9)
    dy = max(float(ylim[1] - ylim[0]), 1e-9)
    col_w = _CSG_FIG_W / float(_CSG_N_COL)
    return (_CSG_FIG_W, max(_CSG_FIG_H_MIN, 2.0 * col_w * (dy / dx) + _CSG_FIG_H_PAD))


def _finish_csg_axes(ax, xlim, ylim) -> None:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.tick_params(labelsize=7)


def _cut_oml_xy(geo, cut, z) -> np.ndarray:
    """OML polygon in the station plane, else wall xy, for ``quadtree_spec_for_airfoil``."""
    poly = np.asarray(geo.oml_polyline(float(z)), float)
    if poly.ndim == 2 and poly.shape[0] >= 3:
        return poly[:, :2]
    chunks = [np.asarray(w.xy, float) for w in cut.walls if len(np.asarray(w.xy, float))]
    if chunks:
        return np.vstack(chunks)
    return np.asarray([[cut.offset[0], cut.offset[1]]], dtype=float)


def _wall_tube_radius(cut) -> float:
    """Half-thickness of the thinnest positive wall (``PolylineBandSDF`` radius)."""
    ts = [float(w.thickness) for w in cut.walls if float(w.thickness) > 0.0]
    return 0.5 * min(ts) if ts else 1e-6


def _csg_quadtree_spec(geo, cut, z):
    """Airfoil domain with thin-band refine so wall tubes are visible at coarse cells."""
    chord = float(cut.chord)
    radius = _wall_tube_radius(cut)
    return quadtree_spec_for_airfoil(
        _cut_oml_xy(geo, cut, z),
        chord=chord,
        max_depth=_CSG_QT_DEPTH,
        max_cells=_CSG_QT_MAX_CELLS,
        min_cell_size=min(0.0015 * chord, 0.5 * radius),
        refine_band=radius,
        thin_band=True,
    )


def _leaf_rects(bounds: np.ndarray) -> np.ndarray:
    """Leaf corners ``(n, 4, 2)`` in QuadCell order."""
    x0, x1, y0, y1 = bounds.T
    return np.stack(
        [
            np.column_stack([x0, y0]),
            np.column_stack([x1, y0]),
            np.column_stack([x1, y1]),
            np.column_stack([x0, y1]),
        ],
        axis=1,
    )


def _sample_wall_meta(geo, cut, z, xlim, ylim):
    """Adaptive ``wall_meta`` on a thin-band quadtree (depth ≤ 10, min cell ≤ t/4)."""
    del xlim, ylim
    field, _table = wall_meta_from_cut(cut)
    spec = _csg_quadtree_spec(geo, cut, z)
    zc = float(z)

    def batch(xy: np.ndarray) -> dict[str, np.ndarray]:
        sample = field.evaluate_plane(np.asarray(xy, float).reshape(-1, 2), zc)
        return {"wall": np.asarray(sample.phi, float).reshape(-1)}

    mesh = build_shared_quadtree(spec, {"wall": lambda pts: batch(pts)["wall"]}, multi_sampler=batch)
    return leaf_bounds(mesh), mesh


def _map_naca_chord_xy(cut, xy: np.ndarray) -> np.ndarray:
    """Chord-frame NACA ``(..., 2)`` (host metres) into blade xy."""
    arr = np.asarray(xy, float)
    if arr.size == 0:
        return np.zeros((0, 2))
    k = float(getattr(cut, "scale", 1.0) or 1.0)
    flat = arr.reshape(-1, 2) / k
    local_x = flat[:, 0] - float(cut.pitch_axis) * float(cut.chord)
    blade = cut.to_blade(local_x, flat[:, 1])
    return np.asarray(blade, float).reshape(*arr.shape[:-1], 2)


def _map_host_xy(cut, xy: np.ndarray) -> np.ndarray:
    """Twisted-about-LE constructive rings into blade xy."""
    arr = np.asarray(xy, float)
    if arr.size == 0:
        return np.zeros((0, 2))
    t = float(getattr(getattr(cut, "host", None), "twist_rad", 0.0) or 0.0)
    flat = arr.reshape(-1, 2)
    if abs(t) > 1e-15:
        ct, st = float(np.cos(t)), float(np.sin(t))
        flat = np.column_stack([ct * flat[:, 0] + st * flat[:, 1], -st * flat[:, 0] + ct * flat[:, 1]])
    return _map_naca_chord_xy(cut, flat.reshape(arr.shape))


def _path_from_rings(rings: list[np.ndarray]) -> MplPath | None:
    verts_list: list[np.ndarray] = []
    codes_list: list[np.ndarray] = []
    for ring in rings:
        c = np.asarray(ring, float)
        if len(c) < 3:
            continue
        verts_list.append(c)
        codes = np.full(len(c), MplPath.LINETO, dtype=MplPath.code_type)
        codes[0] = MplPath.MOVETO
        codes[-1] = MplPath.CLOSEPOLY
        codes_list.append(codes)
    if not verts_list:
        return None
    return MplPath(np.vstack(verts_list), np.concatenate(codes_list))


def _constructive_draw_priority(label: str) -> int:
    if is_core_label(label):
        return 0
    if "skin" in label:
        return 1
    if "spar" in label:
        return 2
    if "shear_web" in label or label.startswith("web"):
        return 3
    return 4


def constructive_polygons_blade_xy(cut) -> dict[str, list[np.ndarray]]:
    """``emit_constructive_polygons(cut.host)`` rings mapped into blade xy."""
    host = getattr(cut, "host", None)
    if host is None:
        return {}
    out: dict[str, list[np.ndarray]] = {}
    for label, geom in emit_constructive_polygons(host).items():
        key = str(label)
        if is_residual_label(key):
            continue
        rings: list[np.ndarray] = []
        parts = getattr(geom, "geoms", None) or (geom,)
        for g in parts:
            if getattr(g, "is_empty", False):
                continue
            for ring in oriented_rings(g):
                xy = _map_host_xy(cut, np.asarray(ring, float))
                if len(xy) >= 3:
                    rings.append(xy)
        if rings:
            out[key] = rings
    return out


def _paint_wall_meta(ax, z, xlim, ylim, bounds, mesh) -> None:
    if bounds.size:
        ax.add_collection(
            PolyCollection(
                _leaf_rects(bounds), facecolors="none", edgecolors="0.75",
                linewidths=0.2, zorder=1,
            )
        )
    fills = label_fill_polygons(mesh, "wall")
    if fills:
        ax.add_collection(
            PolyCollection(fills, facecolors=NAVY, edgecolors="none", alpha=0.85, zorder=2)
        )
    segs = label_boundary_segments(mesh, "wall")
    if len(segs):
        ax.add_collection(LineCollection(segs, colors="k", linewidths=0.6, zorder=3))
    ax.set_title(f"z = {z:.1f} m\nwall occupancy (φ ≤ 0)", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def _paint_kind_solid(ax, geo, cut, z, xlim, ylim) -> None:
    """Same constructive fill as ``plot_section_hybrid``, in blade xy."""
    del geo
    ax.set_facecolor(_OCC_EXTERIOR)
    host = getattr(cut, "host", None)
    if host is not None:
        polys = emit_constructive_polygons(host)
        ordered = sorted(polys.items(), key=lambda kv: _constructive_draw_priority(str(kv[0])))
        for label, geom in ordered:
            lab = str(label)
            if is_residual_label(lab):
                continue
            color = region_color(lab)
            fill_a = _CORE_ALPHA if is_core_label(lab) else _WALL_ALPHA
            zord = 1 if is_core_label(lab) else 2
            parts = getattr(geom, "geoms", None) or (geom,)
            for g in parts:
                if getattr(g, "is_empty", False):
                    continue
                rings = [_map_host_xy(cut, np.asarray(r, float)) for r in oriented_rings(g)]
                path = _path_from_rings(rings)
                if path is None:
                    continue
                ax.add_patch(
                    PathPatch(
                        path, facecolor=color, edgecolor="none", alpha=fill_a,
                        rasterized=True, zorder=zord,
                    )
                )
    ax.set_title(f"z = {z:.1f} m\nconstructive polygons", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def _paint_bays(ax, cut, z, xlim, ylim) -> None:
    if not cut.cells:
        ax.set_title(f"z = {z:.1f} m\nno bays", fontsize=8, pad=6)
    else:
        for i, cell in enumerate(cut.cells):
            xy = np.asarray(cell.points, float)[:, :2]
            if len(xy) < 3:
                continue
            edge, fill = LOOP_STYLES[i % len(LOOP_STYLES)]
            ax.fill(
                xy[:, 0], xy[:, 1], facecolor=fill, edgecolor="none",
                alpha=0.55, zorder=3,
            )
            loop = xy
            if not np.allclose(loop[0], loop[-1]):
                loop = np.vstack([loop, loop[0]])
            ax.plot(
                loop[:, 0], loop[:, 1], color=edge, lw=1.2, ls=_BAY_LOOP_DASH,
                solid_capstyle="round", zorder=4,
            )
        ax.set_title(f"z = {z:.1f} m\ncell bays", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def _blade_xy_to_host_naca(cut, xy: np.ndarray) -> np.ndarray:
    """Inverse of ``_map_naca_chord_xy``: blade xy to host NACA metres."""
    arr = np.asarray(xy, float)
    if arr.size == 0:
        return np.zeros((0, 2))
    flat = arr.reshape(-1, 2)
    ox, oy = float(cut.offset[0]), float(cut.offset[1])
    t = np.radians(float(cut.twist))
    c, s = float(np.cos(t)), float(np.sin(t))
    x = flat[:, 0] - ox
    y = flat[:, 1] - oy
    local_x = c * x - s * y
    local_y = s * x + c * y
    naca = np.column_stack(
        [local_x + float(cut.pitch_axis) * float(cut.chord), local_y]
    )
    k = float(getattr(cut, "scale", 1.0) or 1.0)
    return (k * naca).reshape(*arr.shape[:-1], 2)


def _property_quadtree_spec(geo, cut, z):
    """Same OML box as occupancy; cavity-band refine (not thin-band tubes)."""
    return quadtree_spec_for_airfoil(
        _cut_oml_xy(geo, cut, z),
        chord=float(cut.chord),
        max_depth=_CSG_QT_DEPTH,
        thin_band=False,
    )


def _chord_frame_host(host):
    """Unrotated components so the property tree lives in the airfoil bbox."""
    comps = dict(getattr(host, "components_unrotated", None) or {})
    comps.update(dict(getattr(host, "spar_cap_components_unrotated", None) or {}))
    if not comps:
        return host

    class _ChordHost:
        twist_rad = 0.0
        airfoil = host.airfoil
        components = comps
        labels = list(comps)
        underlying_section = None

        def __getitem__(self, key):
            return comps[key]

        def __iter__(self):
            return iter(comps)

    wrapped = _ChordHost()
    wrapped.underlying_section = wrapped
    return wrapped


def _sample_property_qt(geo, cut, z):
    host = getattr(cut, "host", None)
    if host is None:
        raise ValueError("cut.host is missing")
    spec = _property_quadtree_spec(geo, cut, z)
    evaluator = SectionFieldEvaluator.from_section(_chord_frame_host(host))
    labels = list(evaluator.labels)
    if not labels:
        raise ValueError("cut.host has no property labels")
    inner = evaluator.make_multi_sampler(labels)
    samplers = evaluator.make_samplers(labels)

    def batch(xy: np.ndarray) -> dict[str, np.ndarray]:
        return inner(_blade_xy_to_host_naca(cut, xy))

    wrapped = {
        key: (lambda pts, fn=fn: fn(_blade_xy_to_host_naca(cut, pts)))
        for key, fn in samplers.items()
    }
    return build_shared_quadtree(spec, wrapped, multi_sampler=batch)


def _attach_union_phi(mesh) -> bool:
    phis = []
    for key, arr in (getattr(mesh, "phi_by_label", None) or {}).items():
        if is_wall_solid_label(str(key)):
            phis.append(np.asarray(arr, float))
    if not phis:
        return False
    stacked = np.stack(phis, axis=0)
    mesh.phi_by_label[_UNION_KEY] = np.min(stacked, axis=0)
    return True


def _paint_qt_leaves(ax, mesh, *, fills, color, alpha: float, zorder: int) -> None:
    bounds = leaf_bounds(mesh)
    if bounds.size:
        ax.add_collection(
            PolyCollection(
                _leaf_rects(bounds),
                facecolors="none", edgecolors="0.75", linewidths=0.2, zorder=1,
            )
        )
    if fills:
        ax.add_collection(
            PolyCollection(fills, facecolors=color, edgecolors="none", alpha=alpha, zorder=zorder)
        )


def _paint_property_union(ax, cut, z, xlim, ylim, mesh) -> None:
    del cut
    ax.set_facecolor(_OCC_EXTERIOR)
    if _attach_union_phi(mesh):
        _paint_qt_leaves(
            ax, mesh,
            fills=label_fill_polygons(mesh, _UNION_KEY),
            color=NAVY, alpha=0.85, zorder=2,
        )
        segs = label_boundary_segments(mesh, _UNION_KEY)
        if len(segs):
            ax.add_collection(LineCollection(segs, colors="k", linewidths=0.6, zorder=3))
    ax.set_title(f"z = {z:.1f} m\nproperty quadtree (union)", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def _paint_property_by_label(ax, cut, z, xlim, ylim, mesh) -> None:
    del cut
    ax.set_facecolor(_OCC_EXTERIOR)
    bounds = leaf_bounds(mesh)
    if bounds.size:
        ax.add_collection(
            PolyCollection(
                _leaf_rects(bounds),
                facecolors="none", edgecolors="0.75", linewidths=0.2, zorder=1,
            )
        )
    labels = [str(k) for k in (getattr(mesh, "phi_by_label", None) or {}) if k != _UNION_KEY]
    labels.sort(key=_constructive_draw_priority)
    for lab in labels:
        if is_residual_label(lab):
            continue
        fills = label_fill_polygons(mesh, lab)
        if not fills:
            continue
        color = region_color(lab)
        alpha = _CORE_ALPHA if is_core_label(lab) else 0.75
        zord = 2 if is_core_label(lab) else 3
        ax.add_collection(
            PolyCollection(fills, facecolors=color, edgecolors="none", alpha=alpha, zorder=zord)
        )
    ax.set_title(f"z = {z:.1f} m\nproperty occupancy (by label)", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def _conformal_samples(geo, cut, z, te_blunt: bool | None = None) -> dict:
    tmpl = _template_from_cut(geo.region_graph, cut, float(z), 2, te_blunt=te_blunt)
    return sample_region_polylines(
        cut, cut.shell_graph(), tmpl, geo.region_graph, size_field=True,
    )


def _te_xy_from_samples(samples: dict) -> np.ndarray:
    pts = []
    for name in ("te_su", "te_pr"):
        xy = np.asarray(samples.get(name, []), float)
        if len(xy):
            pts.append(xy[-1, :2])
    if pts:
        return np.vstack(pts)
    chunks = [np.asarray(p, float)[:, :2] for p in samples.values() if len(np.asarray(p))]
    if not chunks:
        return np.zeros((1, 2))
    allxy = np.vstack(chunks)
    return allxy[int(np.argmax(allxy[:, 0]))].reshape(1, 2)


def _te_mesh_limits(cut, samples: dict) -> tuple[tuple[float, float], tuple[float, float]]:
    """Square window around TE keys, pad ``0.08 c`` (y at least ``0.04 c``)."""
    c = max(float(getattr(cut, "chord", 0.0) or 0.0), 1e-9)
    xy = _te_xy_from_samples(samples)
    x0, x1 = float(xy[:, 0].min()), float(xy[:, 0].max())
    y0, y1 = float(xy[:, 1].min()), float(xy[:, 1].max())
    pad = 0.08 * c
    y_pad = max(pad, 0.04 * c)
    x0, x1 = x0 - pad, x1 + pad
    y0, y1 = y0 - y_pad, y1 + y_pad
    side = max(x1 - x0, y1 - y0, 1e-9)
    cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    h = 0.5 * side
    return (cx - h, cx + h), (cy - h, cy + h)


def _paint_mesh_nodes(ax, geo, cut, z, xlim, ylim, samples=None, *,
                      title: str | None = None, marker_scale: float = 1.0) -> None:
    if samples is None:
        samples = _conformal_samples(geo, cut, z)
    s_mid = 10.0 * float(marker_scale)
    s_end = 22.0 * float(marker_scale)
    for name, pts in samples.items():
        xy = np.asarray(pts, float)
        if len(xy) < 1:
            continue
        colour = panel_colour(name)
        ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=0.9, zorder=2)
        ax.scatter(xy[1:-1, 0], xy[1:-1, 1], s=s_mid, c="k", zorder=3, linewidths=0)
        ax.scatter(xy[[0, -1], 0], xy[[0, -1], 1], s=s_end, c="k", zorder=4, linewidths=0)
    ax.set_title(title or f"z = {z:.1f} m\nmid-surface nodes (h)", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def _blank_csg_axes(axes, z: float) -> None:
    for ax in axes:
        ax.set_title(f"z = {z:.1f} m\ninfeasible", fontsize=8, pad=6)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        ax.tick_params(labelsize=7)


def figure_csg_mesh(stem: str, geo, zs, cuts, *, figures: Path, repo: Path) -> None:
    """``csg_mesh/01.png``–``16.png`` plus ``te_01.png``–``te_16.png`` TE zooms."""
    for i, (z, cut) in enumerate(zip(zs, cuts), start=1):
        samples = None
        if cut is None:
            size = (_CSG_FIG_W, 10.0)
        else:
            size = _csg_mesh_figsize(*_station_limits(cut))
            samples = _conformal_samples(geo, cut, z)
        fig, axes = plt.subplots(2, 3, figsize=size)
        ax_csg, ax_phi, ax_prop_u = axes[0]
        ax_bay, ax_mesh, ax_prop_l = axes[1]
        flat = [ax_csg, ax_phi, ax_prop_u, ax_bay, ax_mesh, ax_prop_l]
        if cut is None:
            _blank_csg_axes(flat, z)
        else:
            xlim, ylim = _station_limits(cut)
            _paint_kind_solid(ax_csg, geo, cut, z, xlim, ylim)
            try:
                bounds, occ = _sample_wall_meta(geo, cut, z, xlim, ylim)
            except ValueError:
                _blank_csg_axes((ax_phi,), z)
            else:
                _paint_wall_meta(ax_phi, z, xlim, ylim, bounds, occ)
            try:
                prop_mesh = _sample_property_qt(geo, cut, z)
            except (ValueError, RuntimeError):
                _blank_csg_axes((ax_prop_u, ax_prop_l), z)
            else:
                _paint_property_union(ax_prop_u, cut, z, xlim, ylim, prop_mesh)
                _paint_property_by_label(ax_prop_l, cut, z, xlim, ylim, prop_mesh)
            _paint_bays(ax_bay, cut, z, xlim, ylim)
            _paint_mesh_nodes(ax_mesh, geo, cut, z, xlim, ylim, samples)
        fig.subplots_adjust(left=0.05, right=0.99, top=0.94, bottom=0.06, wspace=0.22, hspace=0.32)
        save_fig(fig, figures, stem, CSG_MESH, f"{i:02d}", repo=repo)
        if cut is None or not samples:
            continue
        te_lim = _te_mesh_limits(cut, samples)
        fig_te, ax_te = plt.subplots(1, 1, figsize=(5.6, 5.6))
        _paint_mesh_nodes(
            ax_te, geo, cut, z, *te_lim, samples,
            title=f"z = {z:.1f} m\nmid-surface nodes (h) TE",
            marker_scale=1.6,
        )
        fig_te.subplots_adjust(left=0.14, right=0.97, top=0.90, bottom=0.12)
        save_fig(fig_te, figures, stem, CSG_MESH, f"te_{i:02d}", repo=repo)


def landmarks(geo) -> list[tuple[str, float, str]]:
    """Root, cylinder–airfoil transition, max chord, 70% span, tip."""
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = np.linspace(z0, z1, 241)
    chords = np.array([geo.interpolate(float(z))[0] for z in zs])
    tcs = []
    for z in zs:
        *_, (_x, upper, lower) = geo.interpolate(float(z))
        tcs.append(float(np.max(upper - lower)))
    tcs = np.asarray(tcs)
    z_maxc = float(zs[int(np.argmax(chords))])
    trans = z0 + 0.12 * (z1 - z0)
    thick = tcs >= 0.88
    if thick.any():
        start = int(np.argmax(thick))
        later = np.where((~thick) & (np.arange(len(zs)) > start))[0]
        if later.size:
            trans = float(zs[int(later[0])])
    return [
        ("root", z0, "root"),
        ("transition", trans, "transition"),
        ("max_chord", z_maxc, "max chord"),
        ("span70", z0 + 0.70 * (z1 - z0), "70% span"),
        ("tip", z1, "tip"),
    ]


def save_fig(fig, root: Path, stem: str, category: str, name: str, *,
             dpi: int = 150, repo: Path | None = None) -> Path:
    out = Path(root) / stem / category
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.png"
    try:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    except ValueError:
        fig.savefig(path, dpi=dpi)
    plt.close(fig)
    if repo is not None:
        print(f"wrote {path.relative_to(repo)}")
    else:
        print(f"wrote {path}")
    return path


def tidy_legacy_pngs(blade_dir: Path) -> None:
    """Remove flat PNGs left at the stem root after the subdirectory move."""
    if not blade_dir.is_dir():
        return
    for path in blade_dir.iterdir():
        if path.is_file() and path.suffix.lower() == ".png" and path.name not in _KEEP_ROOT:
            path.unlink()


def _midline_centroid_x(cut) -> float | None:
    a_sum = sx = 0.0
    for wall in cut.walls:
        p, t = wall.xy, float(wall.thickness)
        if len(p) < 2 or t <= 0.0:
            continue
        a, b = p[:-1], p[1:]
        length = np.linalg.norm(b - a, axis=1)
        a_sum += t * float(length.sum())
        sx += t * float(np.sum(length * 0.5 * (a[:, 0] + b[:, 0])))
    if a_sum <= 1e-18:
        return None
    return sx / a_sum


def _cell_centroid_x(cut) -> float | None:
    if not cut.cells:
        return None
    area = 0.0
    sx = 0.0
    for cell in cut.cells:
        a = abs(float(cell.area))
        if a <= 1e-18:
            continue
        xy = np.asarray(cell.points, float)[:, :2]
        area += a
        sx += a * float(np.mean(xy[:, 0]))
    if area <= 1e-18:
        return None
    return sx / area


def sample_planform(geo, n: int = 80) -> dict:
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = np.linspace(z0, z1, n)
    le, te, pitch, centroid, shear = [], [], [], [], []
    webs: dict[str, list[tuple[float, float]]] = {}
    for zi in zs:
        ring = geo.oml_polyline(float(zi))
        le.append(ring[np.argmin(ring[:, 0])])
        te.append(ring[np.argmax(ring[:, 0])])
        _, _, _, offset, _ = geo.interpolate(float(zi))
        pitch.append((offset[0], offset[1], float(zi)))
        try:
            cut = geo.section_at(float(zi))
        except InfeasibleSectionGeometry:
            continue
        cx = _midline_centroid_x(cut)
        if cx is not None:
            centroid.append((cx, float(zi)))
        sc = _cell_centroid_x(cut)
        if sc is not None:
            shear.append((sc, float(zi)))
        for run in cut.shell_graph().webs:
            xy = np.asarray(run.points, float)[:, :2]
            webs.setdefault(str(run.name), []).append((float(np.mean(xy[:, 0])), float(zi)))
    ribs = []
    for member in getattr(geo, "transverse", ()):
        z = float(member.z)
        ring = geo.oml_polyline(z)
        x_le = float(ring[np.argmin(ring[:, 0]), 0])
        x_te = float(ring[np.argmax(ring[:, 0]), 0])
        ribs.append((z, x_le, x_te, str(member.name)))
    return {
        "z": zs,
        "le": np.asarray(le),
        "te": np.asarray(te),
        "pitch": np.asarray(pitch),
        "centroid": np.asarray(centroid) if centroid else np.zeros((0, 2)),
        "shear": np.asarray(shear) if shear else np.zeros((0, 2)),
        "webs": {k: np.asarray(v) for k, v in webs.items()},
        "ribs": ribs,
    }


def _landmark_ticks(ax, marks, y_top: float) -> None:
    for _key, zi, label in marks:
        ax.axvline(zi, color="#9ca3af", lw=0.6)
        ax.text(zi, y_top, f" {label}", fontsize=7, color="#6b7280", rotation=90, va="top")


def _finish_planform(ax, fig, title: str, marks, te_x) -> None:
    te_x = np.asarray(te_x, float)
    te_x = te_x[np.isfinite(te_x)]
    y_top = float(np.max(te_x)) if te_x.size else 1.0
    _landmark_ticks(ax, marks, y_top)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_title(title)
    n = len(ax.get_legend_handles_labels()[0])
    ncol = min(max(n, 1), 4)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.28),
        ncol=ncol,
        fontsize=8,
        frameon=False,
        borderaxespad=0.0,
    )
    try:
        fig.tight_layout(rect=(0, 0.12, 1, 1))
    except ValueError:
        pass


def figure_planforms(stem: str, title: str, geo, marks, *, figures: Path, repo: Path,
                    n: int = 80) -> None:
    data = sample_planform(geo, n=n)
    figure_planform_axes(stem, title, marks, data, figures=figures, repo=repo)
    figure_planform_structure(stem, title, marks, data, figures=figures, repo=repo)


def figure_planform_axes(stem: str, title: str, marks, data: dict, *, figures: Path,
                         repo: Path) -> None:
    le, te, pitch = data["le"], data["te"], data["pitch"]
    fig, ax = plt.subplots(figsize=(10.5, 3.9))
    ax.plot(le[:, 2], le[:, 0], color=NAVY, lw=1.2, label="LE")
    ax.plot(te[:, 2], te[:, 0], color=TEAL, lw=1.2, label="TE")
    ax.fill_between(le[:, 2], le[:, 0], te[:, 0], color=NAVY, alpha=0.08, lw=0)
    ax.plot(pitch[:, 2], pitch[:, 0], color=RED, lw=1.0, ls="--", label="pitch axis")
    if len(data["centroid"]):
        ax.plot(data["centroid"][:, 1], data["centroid"][:, 0], color="#c45c26", lw=1.0,
                label="centroid")
    if len(data["shear"]):
        ax.plot(data["shear"][:, 1], data["shear"][:, 0], color="#6e4a7a", lw=1.0, ls=":",
                label="shear centre")
    _finish_planform(ax, fig, f"{title} — planform axes (blade-frame x–z)", marks, te[:, 0])
    save_fig(fig, figures, stem, PLANFORM, "axes", repo=repo)


def figure_planform_structure(stem: str, title: str, marks, data: dict, *, figures: Path,
                              repo: Path) -> None:
    le, te, pitch = data["le"], data["te"], data["pitch"]
    fig, ax = plt.subplots(figsize=(10.5, 3.9))
    ax.plot(le[:, 2], le[:, 0], color=NAVY, lw=1.2, label="LE")
    ax.plot(te[:, 2], te[:, 0], color=TEAL, lw=1.2, label="TE")
    ax.fill_between(le[:, 2], le[:, 0], te[:, 0], color=NAVY, alpha=0.08, lw=0)
    ax.plot(pitch[:, 2], pitch[:, 0], color=RED, lw=0.8, ls="--", alpha=0.55, label="pitch axis")
    for i, (name, pts) in enumerate(sorted(data["webs"].items())):
        if len(pts) < 2:
            continue
        ax.plot(pts[:, 1], pts[:, 0], color=SERIES[i % len(SERIES)], lw=1.15, label=name)
    rib_label = False
    for z, x_le, x_te, _name in data.get("ribs", ()):
        ax.plot([z, z], [x_le, x_te], color="#111827", lw=1.05,
                label="ribs" if not rib_label else None)
        rib_label = True
    _finish_planform(ax, fig, f"{title} — planform structure (webs, ribs)", marks, te[:, 0])
    save_fig(fig, figures, stem, PLANFORM, "structure", repo=repo)
