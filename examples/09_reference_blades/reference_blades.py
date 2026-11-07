"""Cut published reference blades (windIO) with the kernel; loft the OML in 3-D.

The private catalogue https://github.com/samhughes-fea/reference-blades is cloned
locally as ``reference_blades/`` (gitignored; do not commit a submodule gitlink).

    git clone https://github.com/samhughes-fea/reference-blades.git reference_blades
    python examples/09_reference_blades/reference_blades.py

Optional Gmsh OML probe (IEA 15 MW only)::

    pip install -e ".[plots,reference,mesh]"
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from _shared.plotting import NAVY, RED, TEAL, WALL_COLOURS, draw_cut, note, plt  # noqa: E402
from implicit_geometry_kernel import Airfoil, BladeGeometry, Station, shared_loft_names  # noqa: E402
from implicit_geometry_kernel.blade.windio import (  # noqa: E402
    cap_chord_windows_from_windio_structure,
    web_positions_from_windio_structure,
)
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule  # noqa: E402
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout  # noqa: E402
from implicit_geometry_kernel.contract import InfeasibleSectionGeometry  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
CASES = (
    ("nrel5mw", "nrel5mw/nrel5mw.yaml", "NREL 5 MW (SNL61.5 lineage)"),
    ("iea15", "iea15/IEA-15-240-RWT.yaml", "IEA Wind 15 MW"),
    ("iea22", "iea22/IEA-22-280-RWT.yaml", "IEA Wind 22 MW"),
    ("iea3p4", "iea3.4_iea10/IEA-3.4-130-RWT.yaml", "IEA Wind 3.4 MW"),
    ("iea10", "iea3.4_iea10/IEA-10-198-RWT.yaml", "IEA Wind 10 MW"),
)
LOFT_STATIONS = 161
RING_POINTS = 96
WALL_POINTS = 48
CLOSEUP_HALF_CHORDS = 0.4
CLOSEUP_STRIPS = 8
GMSH_STEM = "iea15"
FIGURES = Path(__file__).resolve().parent / "figures"


def catalogue_root() -> Path | None:
    env = os.environ.get("REFERENCE_BLADES", "").strip()
    candidates = [Path(env)] if env else []
    candidates.extend((REPO / "reference_blades", REPO.parent / "reference-blades"))
    for path in candidates:
        if (path / "nrel5mw" / "nrel5mw.yaml").is_file():
            return path
    return None


def _curve(block: dict) -> tuple[np.ndarray, np.ndarray]:
    return np.asarray(block["grid"], float), np.asarray(block["values"], float)


def _interp(block: dict, s: float) -> float:
    grid, values = _curve(block)
    return float(np.interp(s, grid, values))


def _sum_tables(tables: list[dict]) -> dict | None:
    if not tables:
        return None
    grid = np.unique(np.concatenate([_curve(t)[0] for t in tables]))
    values = np.zeros_like(grid)
    for table in tables:
        values += np.interp(grid, *_curve(table))
    return {"grid": grid, "values": values}


def _is_core(name: str) -> bool:
    text = name.lower()
    return any(w in text for w in ("foam", "filler", "balsa", "sandwich")) or text.endswith("_core")


def _is_skin_layer(name: str) -> bool:
    key = name.lower()
    if _is_core(key) or key.startswith("dp"):
        return False
    return "shell" in key or "skin" in key or "triax" in key


def _layer_tables(layers: list, *, needles: tuple[str, ...] | None = None,
                  require_web: bool | None = None, predicate=None) -> dict | None:
    matches = []
    for layer in layers:
        name = str(layer.get("name", ""))
        if _is_core(name):
            continue
        tagged = layer.get("web") not in (None, "")
        if require_web is True and not tagged:
            continue
        if require_web is False and tagged:
            continue
        key = name.lower()
        if predicate is not None:
            if not predicate(name):
                continue
        elif needles is None or not any(n in key for n in needles):
            continue
        matches.append(layer["thickness"])
    return _sum_tables(matches)


def _web_gauge(layers: list) -> dict | None:
    by_web: dict[str, list] = {}
    for layer in layers:
        web = layer.get("web")
        if not web or _is_core(str(layer.get("name", ""))):
            continue
        by_web.setdefault(str(web), []).append(layer["thickness"])
    summed = [_sum_tables(tables) for tables in by_web.values() if tables]
    summed = [t for t in summed if t is not None]
    if not summed:
        return None
    grid = np.unique(np.concatenate([_curve(t)[0] for t in summed]))
    stacked = np.column_stack([np.interp(grid, *_curve(t)) for t in summed])
    return {"grid": grid, "values": stacked.max(axis=1)}


def _z_table(z_of_s: np.ndarray, block: dict) -> tuple[np.ndarray, np.ndarray]:
    s, values = _curve(block)
    z = np.interp(s, z_of_s[:, 0], z_of_s[:, 1])
    order = np.argsort(z)
    z, values = z[order], np.asarray(values, float)[order]
    _, idx = np.unique(np.round(z, 12), return_index=True)
    return z[idx], values[idx]


def _pitch_table(z_of_s: np.ndarray, outer: dict) -> tuple[np.ndarray, np.ndarray]:
    s_c, chord = _curve(outer["chord"])
    s_p, pitch_m = _curve(outer["section_offset_y"])
    s = np.unique(np.concatenate([s_c, s_p]))
    c = np.interp(s, s_c, chord)
    p = np.interp(s, s_p, pitch_m)
    frac = np.clip(p / np.maximum(c, 1e-12), 0.05, 0.95)
    z = np.interp(s, z_of_s[:, 0], z_of_s[:, 1])
    order = np.argsort(z)
    z, frac = z[order], frac[order]
    _, idx = np.unique(np.round(z, 12), return_index=True)
    return z[idx], frac[idx]


def _offset_table(z_of_s: np.ndarray, axis: dict) -> tuple[np.ndarray, np.ndarray]:
    s_x, prebend = _curve(axis["x"])
    s_y, sweep = _curve(axis["y"])
    s = np.unique(np.concatenate([z_of_s[:, 0], s_x, s_y]))
    z = np.interp(s, z_of_s[:, 0], z_of_s[:, 1])
    xy = np.column_stack([np.interp(s, s_y, sweep), np.interp(s, s_x, prebend)])
    order = np.argsort(z)
    z, xy = z[order], xy[order]
    _, idx = np.unique(np.round(z, 12), return_index=True)
    return z[idx], xy[idx]


def _thickness_at_z(table: dict | None, z_of_s, floor: float):
    def gauge(z: float) -> float:
        if table is None:
            return floor
        s = float(np.interp(z, z_of_s[:, 1], z_of_s[:, 0]))
        return max(floor, _interp(table, s))

    return gauge


def _airfoil(entry: dict) -> Airfoil:
    coords = entry["coordinates"]
    points = np.column_stack([np.asarray(coords["x"], float), np.asarray(coords["y"], float)])
    return Airfoil(points, str(entry["name"]))


def _layout_from_structure(structure: dict) -> SectionLayout:
    webs = web_positions_from_windio_structure(structure) or ()
    windows = cap_chord_windows_from_windio_structure(structure)
    n = len(webs)
    if n >= 2:
        cap = CapPlacementRule.continuous_box()
        fam = "D"
    elif n == 1:
        cap = CapPlacementRule.per_web_pads()
        fam = "C"
    elif windows:
        cap = CapPlacementRule.from_chord_windows(windows)
        fam = "B"
    else:
        cap = CapPlacementRule.none()
        fam = "A"
    return SectionLayout(
        web_chord_fracs=tuple(float(f) for f in webs),
        orientation="flapwise",
        cap_placement=cap,
        structural_family=fam,
    )


def load_windio(path: Path) -> tuple[str, BladeGeometry, np.ndarray]:
    import yaml

    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    blade = data["components"]["blade"]
    outer = blade["outer_shape"]
    ref_z_g, ref_z_v = _curve(blade["reference_axis"]["z"])
    z_of_s = np.column_stack([ref_z_g, ref_z_v])
    airfoils = {a["name"]: _airfoil(a) for a in data["airfoils"]}
    stations = []
    seen_z = set()
    for row in outer["airfoils"]:
        s = float(row["spanwise_position"])
        z = float(np.interp(s, ref_z_g, ref_z_v))
        key = round(z, 6)
        if key in seen_z:
            continue
        seen_z.add(key)
        chord = _interp(outer["chord"], s)
        twist = _interp(outer["twist"], s)
        pitch_m = _interp(outer["section_offset_y"], s)
        pitch = float(np.clip(pitch_m / chord, 0.05, 0.95))
        prebend = _interp(blade["reference_axis"]["x"], s)
        sweep = _interp(blade["reference_axis"]["y"], s)
        foil = airfoils[str(row["name"])]
        stations.append(Station(z, chord, foil, twist=twist, pitch_axis=pitch, offset=(sweep, prebend)))
    if len(stations) < 2:
        raise ValueError(f"{path} has fewer than two airfoil stations")

    structure = blade.get("structure", {})
    layers = list(structure.get("layers", []))
    webs = list(structure.get("webs", []))
    skin_t = _layer_tables(layers, require_web=False, predicate=_is_skin_layer)
    cap_t = _layer_tables(layers, needles=("spar_cap",), require_web=False)
    web_t = _web_gauge(layers)
    layout = _layout_from_structure(structure)
    web_positions = web_positions_from_windio_structure(structure)
    cap_windows = cap_chord_windows_from_windio_structure(structure)
    geo = BladeGeometry(
        stations,
        layout,
        _thickness_at_z(skin_t, z_of_s, 0.0025),
        spar_cap_thickness=_thickness_at_z(cap_t, z_of_s, 0.005),
        shear_web_thickness=_thickness_at_z(web_t, z_of_s, 0.004),
        web_positions=web_positions,
        cap_windows=(cap_windows if not web_positions else None),
        chord=_z_table(z_of_s, outer["chord"]),
        twist=_z_table(z_of_s, outer["twist"]),
        pitch_axis=_pitch_table(z_of_s, outer),
        offset=_offset_table(z_of_s, blade["reference_axis"]),
    )
    return str(data.get("name", path.stem)), geo, z_of_s


def _blade_dir(stem: str) -> Path:
    out = FIGURES / stem
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_fig(fig, stem: str, name: str, *, dpi: int = 150) -> Path:
    path = _blade_dir(stem) / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path.relative_to(REPO)}")
    return path


def _resample_closed(ring: np.ndarray, n: int = RING_POINTS) -> np.ndarray:
    pts = np.asarray(ring, float)
    if np.linalg.norm(pts[0] - pts[-1]) > 1e-12:
        pts = np.vstack([pts, pts[0]])
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    u = np.concatenate([[0.0], np.cumsum(seg)])
    if u[-1] == 0.0:
        return pts[:n]
    u /= u[-1]
    s = np.linspace(0.0, 1.0, n, endpoint=False)
    return np.column_stack([np.interp(s, u, pts[:, i]) for i in range(pts.shape[1])])


def _oml_rings(geo: BladeGeometry, z_of_s: np.ndarray, n_span: int = LOFT_STATIONS) -> list[np.ndarray]:
    z = np.linspace(float(z_of_s[0, 1]), float(z_of_s[-1, 1]), n_span)
    return [_resample_closed(geo.oml_polyline(zi)) for zi in z]


def _landmarks(geo: BladeGeometry) -> list[tuple[str, float, str]]:
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


def figure_sections(stem: str, title: str, geo: BladeGeometry,
                    marks: list[tuple[str, float, str]]) -> None:
    fig, axes = plt.subplots(1, len(marks), figsize=(15.5, 3.4))
    print(f"\n{title}")
    print(f"{'station':>12} {'z [m]':>8} {'chord':>8} {'twist':>8} {'walls':>6} {'cells':>6}  layout")
    for ax, (_key, z, label) in zip(np.atleast_1d(axes), marks):
        try:
            cut = geo.section_at(z)
        except InfeasibleSectionGeometry as exc:
            ax.set_title(f"{label}\ninfeasible", fontsize=8)
            ax.axis("off")
            print(f"{label:>12} {z:8.2f}  skipped ({exc})")
            continue
        draw_cut(ax, cut)
        ax.plot(0.0, 0.0, "k+", ms=5)
        ax.set_title(f"{label}  z={z:.1f} m\nchord {cut.chord:.2f} m, {cut.twist:.1f}°", fontsize=8)
        ax.axis("off")
        dropped = f", dropped {len(cut.dropped_webs)}" if cut.dropped_webs else ""
        print(f"{label:>12} {z:8.2f} {cut.chord:8.2f} {cut.twist:8.2f} {len(cut.walls):6d} "
              f"{len(cut.cells):6d}  {cut.layout}{dropped}")
    fig.suptitle(f"{title} — {geo.layout} wall midlines", fontsize=10)
    note(fig, "Outer shape, gauges, and web chord fractions from the windIO file; layout is a catalogue key.")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save_fig(fig, stem, "sections")


def figure_planform(stem: str, title: str, geo: BladeGeometry, z_of_s: np.ndarray,
                    marks: list[tuple[str, float, str]]) -> None:
    z = np.linspace(float(z_of_s[0, 1]), float(z_of_s[-1, 1]), 80)
    le, te, pitch = [], [], []
    for zi in z:
        ring = geo.oml_polyline(zi)
        le.append(ring[np.argmin(ring[:, 0])])
        te.append(ring[np.argmax(ring[:, 0])])
        _, _, _, offset, _ = geo.interpolate(zi)
        pitch.append((offset[0], offset[1], zi))
    le, te, pitch = map(np.asarray, (le, te, pitch))
    fig, ax = plt.subplots(figsize=(10.5, 3.2))
    ax.plot(le[:, 2], le[:, 0], color=NAVY, lw=1.2, label="LE")
    ax.plot(te[:, 2], te[:, 0], color=TEAL, lw=1.2, label="TE")
    ax.plot(pitch[:, 2], pitch[:, 0], color=RED, lw=1.0, ls="--", label="pitch axis")
    ax.fill_between(z, le[:, 0], te[:, 0], color=NAVY, alpha=0.08, lw=0)
    y_top = float(np.max(te[:, 0]))
    for _key, zi, label in marks:
        ax.axvline(zi, color="#9ca3af", lw=0.6)
        ax.text(zi, y_top, f" {label}", fontsize=7, color="#6b7280", rotation=90, va="top")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_title(f"{title} — planform (blade-frame x–z)")
    ax.legend(loc="upper right", fontsize=8)
    note(fig, "LE/TE from oml_polyline; ticks mark root, transition, max chord, 70% span, tip.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save_fig(fig, stem, "planform")


def _span_xyz(points: np.ndarray) -> np.ndarray:
    """Plot coordinates (z, x, y) so span is the figure's long axis; lengths unchanged."""
    p = np.asarray(points, float)
    return np.column_stack([p[..., 2], p[..., 0], p[..., 1]])


def _resample_open(pts: np.ndarray, n: int) -> np.ndarray:
    pts = np.asarray(pts, float)
    if len(pts) == 0:
        return pts
    if len(pts) == 1:
        return np.repeat(pts, n, axis=0)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    u = np.concatenate([[0.0], np.cumsum(seg)])
    if u[-1] == 0.0:
        return np.repeat(pts[:1], n, axis=0)
    u /= u[-1]
    s = np.linspace(0.0, 1.0, n)
    return np.column_stack([np.interp(s, u, pts[:, i]) for i in range(pts.shape[1])])


def _wall_stations(geo: BladeGeometry, z_of_s: np.ndarray, n_span: int = LOFT_STATIONS
                   ) -> list[tuple[float, dict[str, tuple[np.ndarray, str, str | None]]]]:
    zs = np.linspace(float(z_of_s[0, 1]), float(z_of_s[-1, 1]), n_span)
    out = []
    for z in zs:
        try:
            cut = geo.section_at(float(z))
        except InfeasibleSectionGeometry:
            continue
        walls = {}
        for run in cut.shell_graph().runs():
            walls[run.name] = (_resample_open(run.points, WALL_POINTS), run.kind, run.name)
        if walls:
            out.append((float(z), walls))
    return out


def _kinds_from_wall_meta(geo: BladeGeometry, faces: list[np.ndarray], kinds: list[str]) -> list[str]:
    """Recolour loft quads by querying ``wall_meta`` at centroids (kind from FeatureTable)."""
    if not faces:
        return kinds
    field = geo.wall_meta()
    cents = np.mean(np.stack(faces), axis=1)
    meta = field.evaluate(cents).meta
    out = list(kinds)
    for i, fid in enumerate(meta):
        if 0 <= int(fid) < len(field.table):
            out[i] = field.table[int(fid)].kind
    return out


def _wall_faces(stations: list, *, skip_suction: bool = False) -> tuple[list[np.ndarray], list[str]]:
    faces, kinds = [], []
    for (_z0, a), (_z1, b) in zip(stations[:-1], stations[1:]):
        for label in shared_loft_names(a.keys(), b.keys()):
            pa, kind, _name = a[label]
            pb = b[label][0]
            if skip_suction and "suction" in label:
                continue
            for i in range(len(pa) - 1):
                faces.append(np.stack([pa[i], pa[i + 1], pb[i + 1], pb[i]]))
                kinds.append(kind)
    return faces, kinds


def _face_rgba(kinds: list[str], alpha: float = 0.88) -> np.ndarray:
    return np.array([to_rgba(WALL_COLOURS.get(k, NAVY), alpha) for k in kinds])


def _legend_kinds(fig, kinds: list[str]) -> None:
    seen = [k for k in ("skin", "cap", "web") if k in set(kinds)]
    if not seen:
        return
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=WALL_COLOURS[k]) for k in seen],
               seen, loc="lower center", ncol=len(seen), fontsize=8)


def _oml_quads(rings: list[np.ndarray]) -> list[np.ndarray]:
    faces = []
    n = rings[0].shape[0]
    for a, b in zip(rings[:-1], rings[1:]):
        for i in range(n):
            j = (i + 1) % n
            faces.append(np.stack([a[i], a[j], b[j], b[i]]))
    return faces


def _crop_rings(geo: BladeGeometry, rings: list[np.ndarray], zc: float) -> tuple[float, list[np.ndarray]]:
    zc = float(np.clip(zc, float(rings[0][0, 2]), float(rings[-1][0, 2])))
    half = CLOSEUP_HALF_CHORDS * max(geo.interpolate(zc)[0], 1e-6)
    crop = [r for r in rings if abs(float(r[0, 2]) - zc) <= half]
    if len(crop) < 2:
        i = int(np.argmin([abs(float(r[0, 2]) - zc) for r in rings]))
        crop = rings[max(0, i - 1): min(len(rings), i + 2)]
    return half, crop


def _strip_window(rings: list[np.ndarray], zc: float, n_strips: int = CLOSEUP_STRIPS) -> list[np.ndarray]:
    """A few loft intervals around ``zc`` so quad edges stay visible (tidalshell recipe)."""
    zs = np.array([float(r[0, 2]) for r in rings])
    i = int(np.argmin(np.abs(zs - zc)))
    i0 = max(0, i - n_strips // 2)
    i1 = min(len(rings) - 1, i0 + n_strips)
    i0 = max(0, i1 - n_strips)
    if i0 >= i1:
        i0, i1 = max(0, i - 1), min(len(rings) - 1, i + 1)
    return rings[i0: i1 + 1]


def _strip_wall_stations(stations: list, zc: float, n_strips: int = CLOSEUP_STRIPS) -> list:
    if len(stations) < 2:
        return stations
    zs = np.array([s[0] for s in stations])
    i = int(np.argmin(np.abs(zs - zc)))
    i0 = max(0, i - n_strips // 2)
    i1 = min(len(stations) - 1, i0 + n_strips)
    i0 = max(0, i1 - n_strips)
    if i0 >= i1:
        i0, i1 = max(0, i - 1), min(len(stations) - 1, i + 1)
    return stations[i0: i1 + 1]


def _finish_3d(ax, pts: np.ndarray, *, elev: float, azim: float) -> None:
    box = np.ptp(pts, axis=0)
    box[box == 0.0] = 1.0
    ax.set_box_aspect(box)
    ax.set_xlabel("z [m]", labelpad=6)
    ax.set_ylabel("x [m]", labelpad=6)
    ax.set_zlabel("y [m]", labelpad=6)
    ax.tick_params(labelsize=7, pad=1)
    ax.view_init(elev=elev, azim=azim)


def _draw_oml(ax, rings: list[np.ndarray], geo: BladeGeometry, *, elev: float, azim: float,
              mesh_edges: bool = False, every_ring: bool = False) -> None:
    shown = [_span_xyz(r) for r in rings]
    faces = _oml_quads(shown)
    rgba = np.broadcast_to(to_rgba(NAVY, 0.88), (len(faces), 4)).copy()
    kw = dict(facecolors=rgba, shade=True)
    if mesh_edges:
        kw["edgecolors"] = np.broadcast_to(to_rgba("#9ca3af", 0.55), (len(faces), 4)).copy()
        kw["linewidths"] = 0.15
    else:
        kw["linewidths"] = 0.0
    ax.add_collection3d(Poly3DCollection(faces, **kw))
    n = len(shown)
    picks = range(n) if every_ring else sorted({0, n // 4, n // 2, (3 * n) // 4, n - 1})
    for i in picks:
        closed = np.vstack([shown[i], shown[i][0]])
        ax.plot(closed[:, 0], closed[:, 1], closed[:, 2], color=RED, lw=0.7)
    pa = []
    for ring in rings:
        zi = float(ring[0, 2])
        _, _, _, offset, _ = geo.interpolate(zi)
        pa.append((zi, offset[0], offset[1]))
    pa = np.asarray(pa)
    ax.plot(pa[:, 0], pa[:, 1], pa[:, 2], color=TEAL, lw=1.2)
    _finish_3d(ax, np.concatenate(shown), elev=elev, azim=azim)


def _draw_walls(ax, faces: list[np.ndarray], kinds: list[str], *, elev: float, azim: float,
                mesh_edges: bool = False) -> None:
    shown = [_span_xyz(f) for f in faces]
    rgba = _face_rgba(kinds)
    kw: dict = dict(facecolors=rgba, shade=True)
    if mesh_edges:
        kw["edgecolors"] = np.broadcast_to(to_rgba("#9ca3af", 0.45), (len(shown), 4)).copy()
        kw["linewidths"] = 0.12
    else:
        kw["linewidths"] = 0.0
    ax.add_collection3d(Poly3DCollection(shown, **kw))
    _finish_3d(ax, np.concatenate(shown), elev=elev, azim=azim)


def _exaggerate(v: np.ndarray, span: float) -> tuple[np.ndarray, float]:
    cross = float(np.ptp(v))
    factor = 1.0 if cross < 1e-9 else max(1.0, 0.16 * span / cross)
    mid = float(np.median(v))
    return mid + (v - mid) * factor, factor


def _project_faces(faces: list[np.ndarray], mode: str) -> tuple[np.ndarray, np.ndarray, str]:
    corners = np.stack(faces)
    x, y, z = corners[:, :, 0], corners[:, :, 1], corners[:, :, 2]
    if mode == "elevation":
        return np.stack([z, y], axis=-1), x.mean(axis=1), "flapwise"
    if mode == "plan":
        return np.stack([z, x], axis=-1), -y.mean(axis=1), "edgewise"
    return np.stack([z, 0.45 * x + 0.90 * y], axis=-1), (-0.90 * x + 0.45 * y).mean(axis=1), "cross-section"


def figure_projected(stem: str, title: str, faces: list[np.ndarray], kinds: list[str]) -> None:
    fc = _face_rgba(kinds, 0.92)
    for mode in ("elevation", "plan", "iso"):
        uv, depth, word = _project_faces(faces, mode)
        span = float(np.ptp(uv[:, :, 0]))
        uv = uv.copy()
        uv[:, :, 1], factor = _exaggerate(uv[:, :, 1], span)
        order = np.argsort(depth)
        fig_w = 13.0
        fig_h = fig_w * float(np.ptp(uv[:, :, 1])) / max(span, 1e-6) + 1.35
        fig, ax = plt.subplots(figsize=(fig_w, max(3.2, fig_h)))
        ax.add_collection(PolyCollection(uv[order], facecolors=fc[order], edgecolors="none", rasterized=True))
        ax.set_aspect("equal")
        ax.autoscale_view()
        ax.set_axis_off()
        note_txt = "true scale" if factor < 1.05 else f"{word} drawn at {factor:.0f}×"
        ax.set_title(f"{title} — {mode} ({note_txt})", loc="left", fontsize=10)
        _legend_kinds(fig, kinds)
        note(fig, "Wall strips from section_at, coloured by kind. loft3d.png is the navy OML default.")
        fig.tight_layout(rect=(0, 0.08, 1, 1))
        save_fig(fig, stem, mode, dpi=200)


def figure_loft3d(stem: str, title: str, geo: BladeGeometry, rings: list[np.ndarray]) -> None:
    fig = plt.figure(figsize=(12.0, 4.2))
    ax = fig.add_subplot(111, projection="3d")
    _draw_oml(ax, rings, geo, elev=14, azim=-72)
    ax.set_title("full span, true metres")
    fig.suptitle(f"{title} — OML loft (navy default)", fontsize=10)
    note(fig, "Single OML from oml_polyline. Coloured wall strips are elevation/plan/iso and close-ups.")
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    save_fig(fig, stem, "loft3d", dpi=200)


def figure_closeups(stem: str, title: str, wall_stations: list,
                    marks: list[tuple[str, float, str]]) -> None:
    fig = plt.figure(figsize=(16.5, 3.8))
    panel_kinds: list[str] = []
    for i, (_key, zc, label) in enumerate(marks, start=1):
        ax = fig.add_subplot(1, len(marks), i, projection="3d")
        crop = _strip_wall_stations(wall_stations, zc)
        faces, kinds = _wall_faces(crop)
        z0, z1 = crop[0][0], crop[-1][0]
        if faces:
            _draw_walls(ax, faces, kinds, elev=8, azim=0, mesh_edges=True)
            panel_kinds.extend(kinds)
        ax.set_title(f"{label}\nz={z0:.1f}–{z1:.1f} m", fontsize=8)
    fig.suptitle(f"{title} — true-scale close-ups", fontsize=10)
    _legend_kinds(fig, panel_kinds)
    note(fig, "Root, transition, max chord, 70% span, tip. Wall strips coloured by kind; lengths true m.")
    fig.tight_layout(rect=(0, 0.10, 1, 0.92))
    save_fig(fig, stem, "closeups", dpi=200)
    for key, zc, label in marks:
        crop = _strip_wall_stations(wall_stations, zc)
        faces, kinds = _wall_faces(crop)
        z0, z1 = crop[0][0], crop[-1][0]
        one = plt.figure(figsize=(8.5, 6.1))
        ax = one.add_subplot(111, projection="3d")
        if faces:
            _draw_walls(ax, faces, kinds, elev=26, azim=-62, mesh_edges=True)
        ax.set_title(f"{title} — {label}, true scale, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(one, kinds)
        note(one, f"{max(len(crop) - 1, 0)} loft strips of section_at walls (skin / cap / web).")
        one.tight_layout(rect=(0, 0.10, 1, 0.96))
        save_fig(one, stem, key, dpi=200)


def figure_cutaway(stem: str, title: str, wall_stations: list, zc: float) -> None:
    crop = _strip_wall_stations(wall_stations, zc)
    faces, kinds = _wall_faces(crop, skip_suction=True)
    if not faces:
        return
    z0, z1 = crop[0][0], crop[-1][0]
    fig = plt.figure(figsize=(8.5, 6.1))
    ax = fig.add_subplot(111, projection="3d")
    _draw_walls(ax, faces, kinds, elev=26, azim=-62, mesh_edges=True)
    ax.set_title(f"{title} — suction skin removed, z = {z0:.2f}–{z1:.2f} m")
    _legend_kinds(fig, kinds)
    note(fig, "Cutaway of section_at wall strips (pressure skin, caps, webs). Not a tidalshell mesh.")
    fig.tight_layout(rect=(0, 0.10, 1, 0.96))
    save_fig(fig, stem, "cutaway", dpi=200)


def _gmsh_faces(gmsh, zmid: float, half: float) -> list[np.ndarray]:
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = np.asarray(coords, float).reshape(-1, 3)
    lookup = {int(t): xyz[i] for i, t in enumerate(tags)}
    n_per = {2: 3, 3: 4}
    faces = []
    types, _, node_tags = gmsh.model.mesh.getElements(2)
    for etype, nodes in zip(types, node_tags):
        k = n_per.get(int(etype))
        if k is None:
            continue
        for row in np.asarray(nodes, int).reshape(-1, k):
            pts = np.array([lookup[int(n)] for n in row])
            if abs(float(pts[:, 2].mean()) - zmid) <= half:
                faces.append(_span_xyz(pts))
    return faces


def figure_mesh_closeup(stem: str, faces: list[np.ndarray], half: float, n_nodes: int, n_elem: int) -> None:
    if not faces:
        print("Gmsh close-up: no triangles in the z-window")
        return
    fig = plt.figure(figsize=(7.2, 5.4))
    ax = fig.add_subplot(111, projection="3d")
    rgba = np.broadcast_to(to_rgba(NAVY, 0.82), (len(faces), 4)).copy()
    edges = np.broadcast_to(to_rgba("#d1d5db", 0.7), (len(faces), 4)).copy()
    ax.add_collection3d(Poly3DCollection(faces, facecolors=rgba, edgecolors=edges, linewidths=0.2, shade=True))
    _finish_3d(ax, np.concatenate(faces), elev=8, azim=0)
    ax.set_title(f"IEA 15 MW Gmsh OML at max chord |Δz| ≤ {half:.2f} m")
    note(fig, f"ThruSections probe ({n_nodes} nodes, {n_elem} faces globally) — not shell-fem wall strips.")
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    save_fig(fig, stem, "mesh_closeup", dpi=200)


def gmsh_oml_probe(rings: list[np.ndarray], stem: str, max_chord: float, geo: BladeGeometry) -> None:
    if importlib.util.find_spec("gmsh") is None:
        print("gmsh not installed; skip OML mesh probe")
        return
    import gmsh

    marks = _landmarks(geo)
    z_maxc = next(z for k, z, _ in marks if k == "max_chord")
    half, _ = _crop_rings(geo, rings, z_maxc)
    zmid = z_maxc
    dz = abs(float(rings[-1][0, 2] - rings[0][0, 2])) / max(len(rings) - 1, 1)
    size = min(0.02 * max(max_chord, 1.0), dz)
    out = _blade_dir(stem) / "oml.msh"
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("oml")
        occ = gmsh.model.occ
        loops = []
        for ring in rings:
            pts = ring if np.linalg.norm(ring[0] - ring[-1]) > 1e-9 * max_chord else ring[:-1]
            tags = [occ.addPoint(float(p[0]), float(p[1]), float(p[2]), size) for p in pts]
            lines = [occ.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))]
            loops.append(occ.addCurveLoop(lines))
        for a, b in zip(loops[:-1], loops[1:]):
            occ.addThruSections([a, b], makeSolid=False, makeRuled=True)
        occ.synchronize()
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", size)
        gmsh.model.mesh.generate(2)
        n_nodes = gmsh.model.mesh.getNodes()[0].size
        n_elem = sum(len(t) for t in gmsh.model.mesh.getElements(2)[1])
        gmsh.write(str(out))
        print(f"Gmsh OML probe: {n_nodes} nodes, {n_elem} surface elements -> {out.relative_to(REPO)}")
        print("  (OML shell only — not the conformal wall-strip mesh in tidalshell / shell-fem)")
        figure_mesh_closeup(stem, _gmsh_faces(gmsh, zmid, half), half, n_nodes, n_elem)
    finally:
        gmsh.finalize()


def main() -> int:
    root = catalogue_root()
    if root is None:
        print("No reference-blades catalogue found. Clone it (gitignored, not a committed submodule):")
        print("  git clone https://github.com/samhughes-fea/reference-blades.git reference_blades")
        print("or set REFERENCE_BLADES to that checkout.")
        return 0
    print(f"catalogue: {root}")
    iea15_rings = None
    iea15_geo = None
    iea15_chord = 1.0
    for stem, rel, title in CASES:
        path = root / rel
        if not path.is_file():
            print(f"missing {path}")
            continue
        name, geo, z_of_s = load_windio(path)
        heading = f"{title} — {name}"
        marks = _landmarks(geo)
        figure_sections(stem, heading, geo, marks)
        figure_planform(stem, heading, geo, z_of_s, marks)
        rings = _oml_rings(geo, z_of_s)
        walls = _wall_stations(geo, z_of_s)
        faces, kinds = _wall_faces(walls)
        if faces:
            kinds = _kinds_from_wall_meta(geo, faces, kinds)
            figure_projected(stem, heading, faces, kinds)
            figure_closeups(stem, heading, walls, marks)
            z_maxc = next(z for k, z, _ in marks if k == "max_chord")
            figure_cutaway(stem, heading, walls, z_maxc)
        else:
            print(f"{stem}: no feasible section_at stations; span views stay OML-only")
        figure_loft3d(stem, heading, geo, rings)
        if stem == GMSH_STEM:
            iea15_rings = rings
            iea15_geo = geo
            zs = np.linspace(geo.stations[0].z, geo.stations[-1].z, 21)
            iea15_chord = max(geo.interpolate(float(zi))[0] for zi in zs)
    if iea15_rings is not None and iea15_geo is not None:
        gmsh_oml_probe(iea15_rings, GMSH_STEM, iea15_chord, iea15_geo)
    else:
        print("IEA 15 MW not loaded; skip Gmsh probe.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
