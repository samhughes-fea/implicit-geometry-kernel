"""Catalogue blade concepts on a made-up planform (not a published reference blade).

A ``2D-F`` matrix (``LC``/``LE`` × ``0R``/``4R``/``10R``) plus ``PA-10R``,
three-web ``LC``/``LE``/``PA``, and ``4D-F-PA-10R`` on one synthetic PCHIP loft.
Wall strips use the authored ``WallRegionGraph`` and ``loft_mapped_strips``.

    python examples/blade/09_blade_concepts/blade_concepts.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from _shared.blade_viz import (  # noqa: E402
    FULL_BLADE,
    REGION_CUTAWAYS,
    REGIONS,
    VIEW_LE,
    VIEW_TE,
    figure_planforms,
    figure_section_grid,
    is_suction_role,
    landmarks,
    live_web_fracs,
    save_fig as save_viz,
    tidy_legacy_pngs,
)
from _shared.debug_mesh import dump_debug_mesh  # noqa: E402
from _shared.plotting import note, panel_colour, plt  # noqa: E402
from implicit_geometry_kernel import (  # noqa: E402
    Airfoil,
    BladeConcept,
    BladeGeometry,
    Station,
    loft_mapped_strips,
)
from implicit_geometry_kernel.blade.topology import (  # noqa: E402
    ShellTopologyTemplate,
    _cut_tc,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
    section_cut_near,
)
KEYS = (
    "2D-F-LC-0R", "2D-F-LC-4R", "2D-F-LC-10R",
    "2D-F-LE-0R", "2D-F-LE-4R", "2D-F-LE-10R",
    "2D-F-PA-10R",
    "3D-F-LC-10R", "3D-F-LE-10R", "3D-F-PA-10R",
    "4D-F-PA-10R",
)
N_SPAN = 41
N_CHORD = 24
CLOSEUP_STRIPS = 8
SKIN, CAP, WEB_T = 0.004, 0.008, 0.004
FIGURES = Path(__file__).resolve().parent / "figures"
SCRIPT = Path(__file__).resolve()


def stations() -> list[Station]:
    root, foil = Airfoil.ellipse(1.0), Airfoil.naca4("0024")
    return [
        Station(0.0, 1.80, root, twist=12.0, pitch_axis=0.50),
        Station(1.0, 1.88, Airfoil.blend(root, foil, 0.15), twist=11.5, pitch_axis=0.47),
        Station(2.5, 2.05, Airfoil.blend(root, foil, 0.50), twist=10.0, pitch_axis=0.40),
        Station(4.0, 2.28, Airfoil.blend(root, foil, 0.85), twist=8.8, pitch_axis=0.33),
        Station(5.0, 2.40, foil, twist=8.0, pitch_axis=0.30),
        Station(6.5, 2.20, Airfoil.naca4("0022"), twist=6.5, pitch_axis=0.30),
        Station(8.5, 1.85, Airfoil.naca4("0020"), twist=4.8, pitch_axis=0.30),
        Station(11.0, 1.50, Airfoil.naca4("0018"), twist=2.8, pitch_axis=0.30),
        Station(14.0, 1.10, Airfoil.naca4("0015"), twist=1.0, pitch_axis=0.30),
        Station(17.0, 0.90, Airfoil.naca4("0013"), twist=-0.6, pitch_axis=0.30),
        Station(20.0, 0.75, Airfoil.naca4("0012"), twist=-2.0, pitch_axis=0.30),
    ]


def build(key: str) -> BladeGeometry:
    st = stations()
    return BladeGeometry.from_concept(
        BladeConcept(key),
        st,
        SKIN,
        spar_cap_thickness=CAP,
        shear_web_thickness=WEB_T,
        span_kind="pchip",
        pivot_fraction=([s.z for s in st], [s.pitch_axis for s in st]),
    )


def _marks(geo: BladeGeometry) -> list[tuple[str, float, str]]:
    return landmarks(geo)


def _station_zs(geo: BladeGeometry) -> list[float]:
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = list(np.linspace(z0, z1, N_SPAN))
    for z in geo.forced_stations():
        if z0 - 1e-9 <= float(z) <= z1 + 1e-9:
            zs.append(float(z))
    return sorted({round(float(z), 9) for z in zs})


def _wall_stations(geo: BladeGeometry) -> list[tuple[float, object, dict]]:
    out = []
    for z in _station_zs(geo):
        cut = section_cut_near(geo, float(z))
        if cut is None:
            continue
        sg = cut.shell_graph()
        tmpl = ShellTopologyTemplate.from_region_graph(
            geo.region_graph, float(z), n_chord=N_CHORD,
            live_web_fracs=live_web_fracs(cut) or None, t_c=_cut_tc(cut),
            dropped_web_indices=tuple(cut.dropped_web_indices),
        )
        samples = sample_region_polylines(cut, sg, tmpl, geo.region_graph)
        if samples:
            out.append((float(z), tmpl, samples))
    return out


def _wall_faces(stations: list, geo: BladeGeometry | None = None, *,
                skip_suction: bool = False) -> tuple[list[np.ndarray], list[str]]:
    faces, roles = [], []
    for (_z0, t0, a), (_z1, t1, b) in zip(stations[:-1], stations[1:]):
        for name, face, _collapsing in loft_mapped_strips(a, b, t0, t1):
            if skip_suction and is_suction_role(name):
                continue
            faces.append(face)
            roles.append(name)
    if geo is not None and stations:
        z0, z1 = float(stations[0][0]), float(stations[-1][0])
        seen = set()
        for z, tmpl, samples in stations:
            if not (z0 - 1e-9 <= float(z) <= z1 + 1e-9):
                continue
            for member in geo.transverse:
                if abs(float(member.z) - float(z)) > 1e-8:
                    continue
                key = (str(member.name), round(float(z), 9))
                if key in seen:
                    continue
                seen.add(key)
                for name, ring in rib_cell_polygons_from_samples(tmpl, samples, member):
                    faces.append(ring)
                    roles.append(name)
    return faces, roles


def _span_xyz(points: np.ndarray) -> np.ndarray:
    p = np.asarray(points, float)
    return np.column_stack([p[..., 2], p[..., 0], p[..., 1]])


def _finite_span_faces(
    faces: list[np.ndarray], kinds: list[str]
) -> tuple[list[np.ndarray], list[str]]:
    shown, keep = [], []
    for face, kind in zip(faces, kinds):
        xyz = _span_xyz(face)
        if xyz.size == 0 or xyz.ndim != 2 or xyz.shape[0] < 3 or xyz.shape[1] != 3:
            continue
        if not np.all(np.isfinite(xyz)) or float(np.max(np.abs(xyz))) > 1.0e6:
            continue
        shown.append(xyz)
        keep.append(kind)
    return shown, keep


def _face_rgba(kinds: list[str], alpha: float = 0.88) -> np.ndarray:
    return np.array([to_rgba(panel_colour(k), alpha) for k in kinds])


def _tight(fig, **kwargs) -> None:
    try:
        fig.tight_layout(**kwargs)
    except ValueError:
        pass


def _legend_kinds(fig, kinds: list[str], *, ncol: int | None = None, fontsize: int = 7) -> None:
    seen = []
    for k in kinds:
        if k not in seen:
            seen.append(k)
    if not seen:
        return
    if ncol is None:
        ncol = min(6, max(1, len(seen)))
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=panel_colour(k)) for k in seen],
               seen, loc="lower center", ncol=ncol, fontsize=fontsize)


def _finish_3d(ax, pts: np.ndarray, *, elev: float, azim: float, zoom: float = 1.0) -> None:
    pts = np.asarray(pts, float)
    if pts.size == 0 or not np.all(np.isfinite(pts)):
        lo, hi = np.array([-1.0, -1.0, -1.0]), np.array([1.0, 1.0, 1.0])
    else:
        lo = pts.min(axis=0)
        hi = pts.max(axis=0)
    box = hi - lo
    pad = np.where(box == 0.0, 1.0, 0.0)
    lo = lo - 0.5 * pad
    hi = hi + 0.5 * pad
    box = hi - lo
    ax.set_xlim(float(lo[0]), float(hi[0]))
    ax.set_ylim(float(lo[1]), float(hi[1]))
    ax.set_zlim(float(lo[2]), float(hi[2]))
    ax.set_box_aspect(box, zoom=float(zoom))
    ax.set_xlabel("z [m]", labelpad=4)
    ax.set_ylabel("x [m]", labelpad=4)
    ax.set_zlabel("y [m]", labelpad=4)
    ax.tick_params(labelsize=7, pad=1)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
    ax.view_init(elev=elev, azim=azim)


def _draw_walls(ax, faces: list[np.ndarray], kinds: list[str], *,
                elev: float = VIEW_LE[0], azim: float = VIEW_LE[1],
                zoom: float = 1.0) -> None:
    shown, kinds = _finite_span_faces(faces, kinds)
    if not shown:
        _finish_3d(ax, np.zeros((0, 3), float), elev=elev, azim=azim, zoom=zoom)
        return
    rgba = _face_rgba(kinds)
    try:
        ax.add_collection3d(Poly3DCollection(
            shown, facecolors=rgba, shade=True,
            edgecolors=np.broadcast_to(to_rgba("#9ca3af", 0.45), (len(shown), 4)).copy(),
            linewidths=0.12,
        ))
    except ValueError:
        pass
    _finish_3d(ax, np.concatenate(shown), elev=elev, azim=azim, zoom=zoom)


def _strip_window(stations: list, zc: float) -> list:
    if len(stations) < 2:
        return stations
    zs = np.array([s[0] for s in stations])
    i = int(np.argmin(np.abs(zs - zc)))
    i0 = max(0, i - CLOSEUP_STRIPS // 2)
    i1 = min(len(stations) - 1, i0 + CLOSEUP_STRIPS)
    i0 = max(0, i1 - CLOSEUP_STRIPS)
    if i0 >= i1:
        i0, i1 = max(0, i - 1), min(len(stations) - 1, i + 1)
    return stations[i0: i1 + 1]


def save_fig(fig, key: str, name: str, *, category: str) -> Path:
    return save_viz(fig, FIGURES, key, category, name, repo=SCRIPT.parents[2])


def figure_sections(key: str, geo: BladeGeometry) -> None:
    figure_section_grid(
        key, key, geo, figures=FIGURES, repo=SCRIPT.parents[2],
        caption="16 equally spaced stations, root to tip, one root-anchored blade-frame window. "
                "Singles are root window | station window. Manufacturing gauges (skin / cap / web).",
    )


def figure_closeups(key: str, geo: BladeGeometry, walls: list, marks: list) -> None:
    fig = plt.figure(figsize=(16.5, 3.8))
    panel_kinds: list[str] = []
    for i, (_k, zc, label) in enumerate(marks, start=1):
        ax = fig.add_subplot(1, len(marks), i, projection="3d")
        crop = _strip_window(walls, zc)
        faces, kinds = _wall_faces(crop, geo=geo)
        z0, z1 = crop[0][0], crop[-1][0]
        _draw_walls(ax, faces, kinds, elev=VIEW_LE[0], azim=VIEW_LE[1])
        panel_kinds.extend(kinds)
        ax.set_title(f"{label}\nz={z0:.1f}–{z1:.1f} m", fontsize=8)
    fig.suptitle(f"{key} — mapped loft close-ups (true scale)", fontsize=10)
    _legend_kinds(fig, panel_kinds)
    note(fig, "Authored WallRegionGraph; web ends are shared skin-loop keys. tip_closure is a station face.")
    _tight(fig, rect=(0, 0.10, 1, 0.92))
    save_fig(fig, key, "closeups", category=FULL_BLADE)
    te_strip = plt.figure(figsize=(16.5, 3.8))
    te_kinds: list[str] = []
    for i, (_k, zc, label) in enumerate(marks, start=1):
        ax = te_strip.add_subplot(1, len(marks), i, projection="3d")
        crop = _strip_window(walls, zc)
        faces, kinds = _wall_faces(crop, geo=geo)
        z0, z1 = crop[0][0], crop[-1][0]
        _draw_walls(ax, faces, kinds, elev=VIEW_TE[0], azim=VIEW_TE[1])
        te_kinds.extend(kinds)
        ax.set_title(f"{label}\nz={z0:.1f}–{z1:.1f} m", fontsize=8)
    te_strip.suptitle(f"{key} — mapped loft close-ups TE (true scale)", fontsize=10)
    _legend_kinds(te_strip, te_kinds)
    note(te_strip, "Trailing-edge view of the same five landmark strips.")
    _tight(te_strip, rect=(0, 0.10, 1, 0.92))
    save_fig(te_strip, key, "closeups_te", category=FULL_BLADE)
    for name, zc, label in marks:
        crop = _strip_window(walls, zc)
        faces, kinds = _wall_faces(crop, geo=geo)
        z0, z1 = crop[0][0], crop[-1][0]
        one = plt.figure(figsize=(8.5, 6.1))
        ax = one.add_subplot(111, projection="3d")
        _draw_walls(ax, faces, kinds, elev=VIEW_LE[0], azim=VIEW_LE[1])
        ax.set_title(f"{key} — {label}, true scale, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(one, kinds)
        note(one, f"{max(len(crop) - 1, 0)} mapped loft strips (region roles + ribs).")
        _tight(one, rect=(0, 0.10, 1, 0.96))
        save_fig(one, key, name, category=REGIONS)
        te = plt.figure(figsize=(8.5, 6.1))
        ax_te = te.add_subplot(111, projection="3d")
        _draw_walls(ax_te, faces, kinds, elev=VIEW_TE[0], azim=VIEW_TE[1])
        ax_te.set_title(f"{key} — {label} TE, true scale, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(te, kinds)
        note(te, f"{max(len(crop) - 1, 0)} mapped loft strips, trailing-edge view.")
        _tight(te, rect=(0, 0.10, 1, 0.96))
        save_fig(te, key, f"{name}_te", category=REGIONS)


def figure_cutaways(key: str, geo: BladeGeometry, walls: list, marks: list) -> None:
    for name, zc, label in marks:
        crop = _strip_window(walls, zc)
        faces, kinds = _wall_faces(crop, geo=geo, skip_suction=True)
        if not faces:
            continue
        z0, z1 = crop[0][0], crop[-1][0]
        fig = plt.figure(figsize=(8.5, 6.1))
        ax = fig.add_subplot(111, projection="3d")
        _draw_walls(ax, faces, kinds, elev=VIEW_LE[0], azim=VIEW_LE[1])
        ax.set_title(f"{key} — {label}, suction skin removed, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(fig, kinds)
        note(fig, "Cutaway of mapped region strips (pressure, caps, webs, ribs).")
        _tight(fig, rect=(0, 0.10, 1, 0.96))
        save_fig(fig, key, name, category=REGION_CUTAWAYS)
        te = plt.figure(figsize=(8.5, 6.1))
        ax_te = te.add_subplot(111, projection="3d")
        _draw_walls(ax_te, faces, kinds, elev=VIEW_TE[0], azim=VIEW_TE[1])
        ax_te.set_title(f"{key} — {label} TE, suction skin removed, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(te, kinds)
        note(te, "Cutaway trailing-edge view.")
        _tight(te, rect=(0, 0.10, 1, 0.96))
        save_fig(te, key, f"{name}_te", category=REGION_CUTAWAYS)


def _elevation_polys(faces: list[np.ndarray], kinds: list[str], *, depth_sign: float = 1.0):
    kept = []
    for face, kind in zip(faces, kinds):
        arr = np.asarray(face, float)
        if arr.ndim != 2 or arr.shape[1] < 3 or not np.all(np.isfinite(arr)):
            continue
        if float(np.max(np.abs(arr))) > 1.0e6:
            continue
        kept.append((arr, kind))
    faces = [item[0] for item in kept]
    kinds = [item[1] for item in kept]
    polys = [np.column_stack([np.asarray(f, float)[:, 2], np.asarray(f, float)[:, 1]]) for f in faces]
    depth = depth_sign * np.array([float(np.mean(np.asarray(f, float)[:, 0])) for f in faces])
    all_uv = np.vstack(polys) if polys else np.zeros((0, 2))
    span = float(np.ptp(all_uv[:, 0])) if len(all_uv) else 1.0
    cross = float(np.ptp(all_uv[:, 1])) if len(all_uv) else 1.0
    factor = 1.0 if cross < 1e-9 else max(1.0, 0.16 * span / cross)
    mid = float(np.median(all_uv[:, 1])) if len(all_uv) else 0.0
    uv = [p.copy() for p in polys]
    for p in uv:
        p[:, 1] = mid + (p[:, 1] - mid) * factor
    return uv, np.argsort(depth), factor, span, kinds


def _save_elevation(key: str, name: str, faces: list[np.ndarray], kinds: list[str], *,
                    title: str, caption: str, depth_sign: float = 1.0,
                    legend: bool = True) -> None:
    if not faces:
        return
    uv, order, factor, span, kinds = _elevation_polys(faces, kinds, depth_sign=depth_sign)
    if not uv:
        return
    stretched = np.vstack(uv)
    if stretched.size == 0 or not np.all(np.isfinite(stretched)):
        return
    fig_h = 13.0 * float(np.ptp(stretched[:, 1])) / max(span, 1e-6) + (1.35 if legend else 0.55)
    if not np.isfinite(fig_h):
        return
    fig, ax = plt.subplots(figsize=(13.0, max(3.2, fig_h)))
    ax.add_collection(PolyCollection(
        [uv[i] for i in order], facecolors=_face_rgba(kinds, 0.92)[order],
        edgecolors="none", rasterized=True,
    ))
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontsize=10)
    if legend:
        _legend_kinds(fig, kinds)
        _tight(fig, rect=(0, 0.08, 1, 1))
    else:
        _tight(fig, rect=(0, 0.04, 1, 1))
    note(fig, caption)
    save_fig(fig, key, name, category=FULL_BLADE)


def _camera_basis(elev: float, azim: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Screen-right, screen-up, toward-camera axes (matplotlib elev/azim, z vertical)."""
    e, a = np.deg2rad(elev), np.deg2rad(azim)
    toward = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)], float)
    toward /= np.linalg.norm(toward)
    right = np.cross(np.array([0.0, 0.0, 1.0]), toward)
    right /= np.linalg.norm(right)
    up = np.cross(toward, right)
    return right, up, toward


def _save_corner_cutaway(key: str, name: str, faces: list[np.ndarray], kinds: list[str], *,
                         title: str, caption: str, elev: float, azim: float) -> None:
    """Orthographic corner view, cropped to the blade (no 3-D letterboxing)."""
    if not faces:
        return
    right, up, toward = _camera_basis(elev, azim)
    shown, kinds = _finite_span_faces(faces, kinds)
    if not shown:
        return
    raw, depth = [], np.empty(len(shown), float)
    for i, pts in enumerate(shown):
        raw.append(np.column_stack([pts @ right, pts @ up]))
        depth[i] = float(np.mean(pts @ toward))
    z_span = np.array([float(np.mean(p[:, 0])) for p in shown])
    cents = np.array([np.mean(p, axis=0) for p in raw])
    if not np.all(np.isfinite(cents)):
        return
    root = cents[z_span <= z_span.min() + 1e-9].mean(axis=0)
    tip = cents[z_span >= z_span.max() - 1e-9].mean(axis=0)
    axis = tip - root
    ang = float(np.arctan2(axis[1], axis[0]))
    c, s = np.cos(-ang), np.sin(-ang)
    rot = np.array([[c, -s], [s, c]], float)
    uv = [p @ rot.T for p in raw]
    order = np.argsort(depth)
    stacked = np.vstack(uv)
    if stacked.size == 0 or not np.all(np.isfinite(stacked)):
        return
    span_u = float(np.ptp(stacked[:, 0]))
    span_v = float(np.ptp(stacked[:, 1]))
    pad_u, pad_v = 0.015 * span_u, 0.06 * span_v
    fig_w = 14.5
    data_h = fig_w * (span_v + 2.0 * pad_v) / max(span_u + 2.0 * pad_u, 1e-6)
    fig, ax = plt.subplots(figsize=(fig_w, data_h + 0.72))
    ax.add_collection(PolyCollection(
        [uv[i] for i in order], facecolors=_face_rgba(kinds, 0.92)[order],
        edgecolors=to_rgba("#6b7280", 0.28),
        linewidths=0.08, rasterized=True,
    ))
    ax.set_aspect("equal")
    ax.set_xlim(float(stacked[:, 0].min()) - pad_u, float(stacked[:, 0].max()) + pad_u)
    ax.set_ylim(float(stacked[:, 1].min()) - pad_v, float(stacked[:, 1].max()) + pad_v)
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontsize=10)
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.09, top=0.88)
    note(fig, caption)
    save_fig(fig, key, name, category=FULL_BLADE)


def figure_full_cutaway(key: str, faces: list[np.ndarray], kinds: list[str]) -> None:
    """Full-span blade, suction skin omitted, looking down from a corner."""
    note_txt = "suction skin omitted; pressure, webs, ribs; corner view looking down"
    _save_corner_cutaway(
        key, "cutaway", faces, kinds,
        title=f"{key} — full blade, suction skin removed",
        caption=note_txt, elev=38.0, azim=-52.0,
    )
    _save_corner_cutaway(
        key, "cutaway_te", faces, kinds,
        title=f"{key} — full blade TE, suction skin removed",
        caption=note_txt + "; trailing-edge corner",
        elev=38.0, azim=128.0,
    )


def figure_elevation(key: str, title: str, faces: list[np.ndarray], kinds: list[str]) -> None:
    _, _, factor, _, _ = _elevation_polys(faces, kinds)
    note_txt = "true scale" if factor < 1.05 else f"flapwise drawn at {factor:.0f}×"
    _save_elevation(
        key, "elevation", faces, kinds,
        title=f"{title} — elevation ({note_txt})",
        caption="Mapped loft by authored region. Catalogue concept on a made-up chord/twist; not a reference blade.",
        legend=True,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--smoke", action="store_true",
                        help="one catalogue key (CI). Default is the full KEYS gallery.")
    parser.add_argument("--keys", nargs="+", metavar="TOKEN",
                        help="subset of KEYS, e.g. 2D-F-PA-10R 3D-F-LC-10R")
    args = parser.parse_args(argv)
    keys = list(KEYS)
    if args.smoke:
        keys = ["2D-F-PA-10R"]
    elif args.keys:
        unknown = [k for k in args.keys if k not in KEYS]
        if unknown:
            parser.error(f"unknown key(s) {unknown}; choose from {list(KEYS)}")
        keys = list(args.keys)
    for key in keys:
        geo = build(key)
        marks = _marks(geo)
        print(f"\n{key}  regions={len(geo.region_graph.regions)}  "
              f"webs={len(geo.region_graph.webs)}")
        tidy_legacy_pngs(FIGURES / key)
        figure_sections(key, geo)
        figure_planforms(key, key, geo, marks, figures=FIGURES, repo=SCRIPT.parents[2], n=41)
        walls = _wall_stations(geo)
        if not walls:
            print(f"{key}: no feasible loft stations")
            continue
        figure_closeups(key, geo, walls, marks)
        figure_cutaways(key, geo, walls, marks)
        faces, kinds = _wall_faces(walls, geo=geo)
        dump_debug_mesh(FIGURES / key, walls, faces, kinds, label=key, repo=SCRIPT.parents[2], geo=geo)
        figure_elevation(key, key, faces, kinds)
        cut_faces, cut_kinds = _wall_faces(walls, geo=geo, skip_suction=True)
        dump_debug_mesh(
            FIGURES / key, walls, cut_faces, cut_kinds,
            label=f"{key} cutaway", repo=SCRIPT.parents[2], geo=geo,
            stem="cutaway", write_samples=False,
        )
        figure_full_cutaway(key, cut_faces, cut_kinds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
