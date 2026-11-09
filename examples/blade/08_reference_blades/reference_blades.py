"""Cut published reference blades (windIO) with the kernel; loft the OML in 3-D.

The private catalogue https://github.com/samhughes-fea/reference-blades is cloned
locally as ``reference_blades/`` (gitignored; do not commit a submodule gitlink).

    git clone https://github.com/samhughes-fea/reference-blades.git reference_blades
    python examples/blade/08_reference_blades/reference_blades.py

Optional Gmsh wall-strip mesh (IEA 15 MW only)::

    pip install -e ".[plots,reference,mesh]"
"""

from __future__ import annotations

import importlib.util
import os
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
from _shared.debug_mesh import as_quad, dump_debug_mesh  # noqa: E402
from _shared.plotting import (  # noqa: E402
    MESH_LW_3D,
    MESH_LW_FLAT,
    NAVY,
    RED,
    TEAL,
    mesh_edge_rgba,
    mesh_linewidths,
    note,
    panel_colour,
    plt,
)
from implicit_geometry_kernel import (  # noqa: E402
    Airfoil,
    BladeGeometry,
    Station,
    WallRegionGraph,
    loft_mapped_strips,
)
from implicit_geometry_kernel.blade.topology import (  # noqa: E402
    ShellTopologyTemplate,
    _cut_tc,
    rib_cell_polygons_from_samples,
    sample_region_polylines,
    section_cut_near,
)
from implicit_geometry_kernel.mesh.polygon import fill_rib_cells  # noqa: E402
from implicit_geometry_kernel.blade.windio import (  # noqa: E402
    cap_window_table_from_windio_structure,
    transverse_members_from_windio,
    web_positions_fn_from_windio_structure,
    web_positions_from_windio_structure,
)
from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule  # noqa: E402
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout  # noqa: E402

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
    kernel_root = Path(__file__).resolve().parents[3]
    candidates.extend(
        (
            kernel_root / "reference_blades",
            kernel_root.parent / "reference-blades",
            REPO / "reference_blades",
            REPO.parent / "reference-blades",
        )
    )
    for path in candidates:
        if not path:
            continue
        for marker in (
            path / "wind" / "nrel5mw" / "nrel5mw.yaml",
            path / "nrel5mw" / "nrel5mw.yaml",
        ):
            if marker.is_file():
                return path
    return None


def catalogue_yaml(root: Path, rel: str) -> Path:
    """Resolve a blade YAML; catalogue folders are ``wind/<blade>/``."""
    wind = root / "wind" / rel
    if wind.is_file():
        return wind
    return root / rel


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
    windows = cap_window_table_from_windio_structure(structure)
    n = len(webs)
    if windows is not None:
        z, vals = windows
        pair = ((float(vals[0, 0]), float(vals[0, 1])),)
        cap = CapPlacementRule.from_chord_windows(pair)
        fam = "B" if n == 0 else "D"
    elif n >= 2:
        cap = CapPlacementRule.continuous_box()
        fam = "D"
    elif n == 1:
        cap = CapPlacementRule.per_web_pads()
        fam = "C"
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
    web_positions = web_positions_fn_from_windio_structure(structure, z_of_s)
    cap_windows = cap_window_table_from_windio_structure(structure, z_of_s)
    geo = BladeGeometry(
        stations,
        layout,
        _thickness_at_z(skin_t, z_of_s, 0.0025),
        spar_cap_thickness=_thickness_at_z(cap_t, z_of_s, 0.005),
        shear_web_thickness=_thickness_at_z(web_t, z_of_s, 0.004),
        web_positions=web_positions,
        cap_windows=cap_windows,
        chord=_z_table(z_of_s, outer["chord"]),
        twist=_z_table(z_of_s, outer["twist"]),
        pitch_axis=_pitch_table(z_of_s, outer),
        offset=_offset_table(z_of_s, blade["reference_axis"]),
        region_graph=WallRegionGraph.from_windio_structure(structure, z_of_s),
        transverse=transverse_members_from_windio(structure, z_of_s),
    )
    return str(data.get("name", path.stem)), geo, z_of_s


def _blade_dir(stem: str) -> Path:
    out = FIGURES / stem
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_fig(fig, stem: str, name: str, *, dpi: int = 150, category: str | None = None) -> Path:
    if category:
        return save_viz(fig, FIGURES, stem, category, name, dpi=dpi, repo=REPO)
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
    return landmarks(geo)


def figure_sections(stem: str, title: str, geo: BladeGeometry) -> None:
    print(f"\n{title}")
    figure_section_grid(
        stem, title, geo, figures=FIGURES, repo=REPO,
        caption="16 equally spaced stations, root to tip, one root-anchored blade-frame window. "
                "Singles are root window | station window. Manufacturing gauges from the windIO file.",
    )


def _span_xyz(points: np.ndarray) -> np.ndarray:
    """Plot coordinates (z, x, y) so span is the figure's long axis; lengths unchanged."""
    p = np.asarray(points, float)
    return np.column_stack([p[..., 2], p[..., 0], p[..., 1]])


def _station_zs(geo: BladeGeometry, z_of_s: np.ndarray, n_span: int) -> list[float]:
    z0, z1 = float(z_of_s[0, 1]), float(z_of_s[-1, 1])
    zs = list(np.linspace(z0, z1, n_span))
    for z in geo.forced_stations():
        if z0 - 1e-9 <= float(z) <= z1 + 1e-9:
            zs.append(float(z))
    return sorted({round(float(z), 9) for z in zs})


def _wall_stations(geo: BladeGeometry, z_of_s: np.ndarray, n_span: int = LOFT_STATIONS
                   ) -> list[tuple[float, object, dict[str, np.ndarray]]]:
    out = []
    for z in _station_zs(geo, z_of_s, n_span):
        cut = section_cut_near(geo, float(z))
        if cut is None:
            continue
        sg = cut.shell_graph()
        live = live_web_fracs(cut)
        tmpl = ShellTopologyTemplate.from_region_graph(
            geo.region_graph, float(getattr(cut, "z", z)), n_chord=WALL_POINTS,
            live_web_fracs=live or None,
            t_c=_cut_tc(cut),
            dropped_web_indices=tuple(cut.dropped_web_indices),
        )
        samples = sample_region_polylines(cut, sg, tmpl, geo.region_graph)
        if samples:
            out.append((float(getattr(cut, "z", z)), tmpl, samples))
    return out


def _is_suction_role(name: str) -> bool:
    return is_suction_role(name)


_WALL_CACHE: dict[tuple, tuple[list[np.ndarray], list[str]]] = {}


def _wall_faces(stations: list, *, skip_suction: bool = False, geo: BladeGeometry | None = None
                ) -> tuple[list[np.ndarray], list[str]]:
    key = (
        id(geo) if geo is not None else 0,
        skip_suction,
        tuple(round(float(s[0]), 9) for s in stations),
    )
    hit = _WALL_CACHE.get(key)
    if hit is not None:
        return hit
    faces, roles = [], []
    for (_z0, t0, a), (_z1, t1, b) in zip(stations[:-1], stations[1:]):
        for name, face, _collapsing in loft_mapped_strips(a, b, t0, t1):
            if skip_suction and _is_suction_role(name):
                continue
            faces.append(face)
            roles.append(name)
    if geo is not None and stations:
        z0, z1 = float(stations[0][0]), float(stations[-1][0])
        seen = set()
        cells: list[tuple[str, np.ndarray]] = []
        for member in geo.transverse:
            zm = float(member.z)
            if not (z0 - 1e-9 <= zm <= z1 + 1e-9):
                continue
            _z, tmpl, samples = min(stations, key=lambda s: abs(float(s[0]) - zm))
            name = str(member.name)
            if name in seen:
                continue
            seen.add(name)
            cells.extend(rib_cell_polygons_from_samples(tmpl, samples, member))
        for name, face in fill_rib_cells(cells):
            faces.append(face)
            roles.append(name)
    _WALL_CACHE[key] = (faces, roles)
    return faces, roles


def _face_rgba(kinds: list[str], alpha: float = 0.88) -> np.ndarray:
    return np.array([to_rgba(panel_colour(k), alpha) for k in kinds])


def _legend_kinds(fig, kinds: list[str]) -> None:
    seen = []
    for k in kinds:
        if k not in seen:
            seen.append(k)
    if not seen:
        return
    ncol = min(6, max(1, len(seen)))
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=panel_colour(k)) for k in seen],
               seen, loc="lower center", ncol=ncol, fontsize=7)


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
        kw["edgecolors"] = mesh_edge_rgba(len(faces))
        kw["linewidths"] = MESH_LW_3D
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
        kw["edgecolors"] = mesh_edge_rgba(len(shown))
        kw["linewidths"] = mesh_linewidths(kinds, base=MESH_LW_3D)
    else:
        kw["linewidths"] = 0.0
    ax.add_collection3d(Poly3DCollection(shown, **kw))
    _finish_3d(ax, np.concatenate(shown), elev=elev, azim=azim)


def _exaggerate(v: np.ndarray, span: float) -> tuple[np.ndarray, float]:
    cross = float(np.ptp(v))
    factor = 1.0 if cross < 1e-9 else max(1.0, 0.16 * span / cross)
    mid = float(np.median(v))
    return mid + (v - mid) * factor, factor


def _as_quad(face: np.ndarray) -> np.ndarray:
    return as_quad(face)


def _dump_debug_mesh(stem: str, walls: list, faces: list[np.ndarray], kinds: list[str],
                     geo: BladeGeometry | None = None) -> None:
    dump_debug_mesh(_blade_dir(stem), walls, faces, kinds, label=stem, repo=REPO, geo=geo)


def _project_faces(faces: list[np.ndarray], mode: str) -> tuple[list[np.ndarray], np.ndarray, str]:
    uvs, depths = [], []
    word = "flapwise"
    for f in faces:
        pts = np.asarray(f, float)
        x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
        if mode == "elevation":
            uv = np.column_stack([z, y])
            d = float(x.mean())
            word = "flapwise"
        elif mode == "plan":
            uv = np.column_stack([z, x])
            d = float(-y.mean())
            word = "edgewise"
        else:
            uv = np.column_stack([z, 0.45 * x + 0.90 * y])
            d = float((-0.90 * x + 0.45 * y).mean())
            word = "cross-section"
        uvs.append(uv)
        depths.append(d)
    return uvs, np.asarray(depths, float), word


def figure_projected(stem: str, title: str, faces: list[np.ndarray], kinds: list[str]) -> None:
    fc = _face_rgba(kinds, 0.92)
    for mode in ("elevation", "plan", "iso"):
        uv, depth, word = _project_faces(faces, mode)
        stacked = np.vstack(uv) if uv else np.zeros((0, 2))
        span = float(np.ptp(stacked[:, 0])) if len(stacked) else 1.0
        stretched = [p.copy() for p in uv]
        y_all = stacked[:, 1] if len(stacked) else np.zeros(1)
        y_all, factor = _exaggerate(y_all, span)
        mid = float(np.median(stacked[:, 1])) if len(stacked) else 0.0
        for p in stretched:
            p[:, 1] = mid + (p[:, 1] - mid) * factor
        order = np.argsort(depth)
        fig_w = 13.0
        fig_h = fig_w * float(np.ptp(np.vstack(stretched)[:, 1])) / max(span, 1e-6) + 1.35
        fig, ax = plt.subplots(figsize=(fig_w, max(3.2, fig_h)))
        ax.add_collection(PolyCollection(
            [stretched[i] for i in order], facecolors=fc[order],
            edgecolors=mesh_edge_rgba(len(order)),
            linewidths=[mesh_linewidths(kinds, base=MESH_LW_FLAT)[i] for i in order],
            rasterized=True,
        ))
        ax.set_aspect("equal")
        ax.autoscale_view()
        ax.set_axis_off()
        note_txt = "true scale" if factor < 1.05 else f"{word} drawn at {factor:.0f}×"
        ax.set_title(f"{title} — {mode} ({note_txt})", loc="left", fontsize=10)
        _legend_kinds(fig, kinds)
        note(fig, "Mapped loft by authored region (RETRO_PANELS).")
        fig.tight_layout(rect=(0, 0.08, 1, 1))
        save_fig(fig, stem, mode, dpi=200, category=FULL_BLADE)


def figure_closeups(stem: str, title: str, wall_stations: list,
                    marks: list[tuple[str, float, str]], geo: BladeGeometry) -> None:
    fig = plt.figure(figsize=(16.5, 3.8))
    panel_kinds: list[str] = []
    for i, (_key, zc, label) in enumerate(marks, start=1):
        ax = fig.add_subplot(1, len(marks), i, projection="3d")
        crop = _strip_wall_stations(wall_stations, zc)
        faces, kinds = _wall_faces(crop, geo=geo)
        z0, z1 = crop[0][0], crop[-1][0]
        if faces:
            _draw_walls(ax, faces, kinds, elev=VIEW_LE[0], azim=VIEW_LE[1], mesh_edges=True)
            panel_kinds.extend(kinds)
        ax.set_title(f"{label}\nz={z0:.1f}–{z1:.1f} m", fontsize=8)
    fig.suptitle(f"{title} — true-scale close-ups", fontsize=10)
    _legend_kinds(fig, panel_kinds)
    note(fig, "Root, transition, max chord, 70% span, tip. Colour by region role; tip_closure at tip.")
    fig.tight_layout(rect=(0, 0.10, 1, 0.92))
    save_fig(fig, stem, "closeups", dpi=200, category=FULL_BLADE)
    te_strip = plt.figure(figsize=(16.5, 3.8))
    te_kinds: list[str] = []
    for i, (_key, zc, label) in enumerate(marks, start=1):
        ax = te_strip.add_subplot(1, len(marks), i, projection="3d")
        crop = _strip_wall_stations(wall_stations, zc)
        faces, kinds = _wall_faces(crop, geo=geo)
        z0, z1 = crop[0][0], crop[-1][0]
        if faces:
            _draw_walls(ax, faces, kinds, elev=VIEW_TE[0], azim=VIEW_TE[1], mesh_edges=True)
            te_kinds.extend(kinds)
        ax.set_title(f"{label}\nz={z0:.1f}–{z1:.1f} m", fontsize=8)
    te_strip.suptitle(f"{title} — true-scale close-ups TE", fontsize=10)
    _legend_kinds(te_strip, te_kinds)
    note(te_strip, "Trailing-edge view of the same five landmark strips.")
    te_strip.tight_layout(rect=(0, 0.10, 1, 0.92))
    save_fig(te_strip, stem, "closeups_te", dpi=200, category=FULL_BLADE)
    for key, zc, label in marks:
        crop = _strip_wall_stations(wall_stations, zc)
        faces, kinds = _wall_faces(crop, geo=geo)
        z0, z1 = crop[0][0], crop[-1][0]
        one = plt.figure(figsize=(8.5, 6.1))
        ax = one.add_subplot(111, projection="3d")
        if faces:
            _draw_walls(ax, faces, kinds, elev=VIEW_LE[0], azim=VIEW_LE[1], mesh_edges=True)
        ax.set_title(f"{title} — {label}, true scale, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(one, kinds)
        note(one, f"{max(len(crop) - 1, 0)} mapped loft strips (region roles + ribs).")
        one.tight_layout(rect=(0, 0.10, 1, 0.96))
        save_fig(one, stem, key, dpi=200, category=REGIONS)
        te = plt.figure(figsize=(8.5, 6.1))
        ax_te = te.add_subplot(111, projection="3d")
        if faces:
            _draw_walls(ax_te, faces, kinds, elev=VIEW_TE[0], azim=VIEW_TE[1], mesh_edges=True)
        ax_te.set_title(f"{title} — {label} TE, true scale, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(te, kinds)
        note(te, f"{max(len(crop) - 1, 0)} mapped loft strips, trailing-edge view.")
        te.tight_layout(rect=(0, 0.10, 1, 0.96))
        save_fig(te, stem, f"{key}_te", dpi=200, category=REGIONS)


def figure_cutaways(stem: str, title: str, wall_stations: list,
                    marks: list[tuple[str, float, str]], geo: BladeGeometry) -> None:
    for key, zc, label in marks:
        crop = _strip_wall_stations(wall_stations, zc)
        faces, kinds = _wall_faces(crop, skip_suction=True, geo=geo)
        if not faces:
            continue
        z0, z1 = crop[0][0], crop[-1][0]
        fig = plt.figure(figsize=(8.5, 6.1))
        ax = fig.add_subplot(111, projection="3d")
        _draw_walls(ax, faces, kinds, elev=VIEW_LE[0], azim=VIEW_LE[1], mesh_edges=True)
        ax.set_title(f"{title} — {label}, suction skin removed, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(fig, kinds)
        note(fig, "Suction skin omitted; camera on the suction side looking onto inner pressure skin.")
        fig.tight_layout(rect=(0, 0.10, 1, 0.96))
        save_fig(fig, stem, key, dpi=200, category=REGION_CUTAWAYS)
        te = plt.figure(figsize=(8.5, 6.1))
        ax_te = te.add_subplot(111, projection="3d")
        _draw_walls(ax_te, faces, kinds, elev=VIEW_TE[0], azim=VIEW_TE[1], mesh_edges=True)
        ax_te.set_title(f"{title} — {label} TE, suction skin removed, z = {z0:.2f}–{z1:.2f} m")
        _legend_kinds(te, kinds)
        note(te, "Same omitted-suction mesh; trailing-edge corner, still looking from the suction side.")
        te.tight_layout(rect=(0, 0.10, 1, 0.96))
        save_fig(te, stem, f"{key}_te", dpi=200, category=REGION_CUTAWAYS)


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
    ax.set_title(f"IEA 15 MW Gmsh wall strips at max chord |Δz| ≤ {half:.2f} m")
    note(fig, f"Transfinite wall-strip mesh ({n_nodes} nodes, {n_elem} faces globally).")
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    save_fig(fig, stem, "mesh_closeup", dpi=200)


def gmsh_shell_probe(geo: BladeGeometry, stem: str) -> None:
    if importlib.util.find_spec("gmsh") is None:
        print("gmsh not installed; skip wall-strip mesh probe")
        return
    from implicit_geometry_kernel.mesh.gmsh import write_shell_msh

    marks = _landmarks(geo)
    z_maxc = next(z for k, z, _ in marks if k == "max_chord")
    z0, z1 = float(geo.stations[0].z), float(geo.stations[-1].z)
    zs = np.linspace(z0, z1, 9)
    out = _blade_dir(stem) / "shell.msh"
    stats = write_shell_msh(geo, out, zs, n_chord=6)
    print(f"Gmsh wall-strip mesh: {stats['n_nodes']} nodes, {stats['n_elements']} faces "
          f"-> {out.relative_to(REPO)}")
    half = 0.15 * max(float(geo.interpolate(z_maxc)[0]), 1e-3)
    faces = [_span_xyz(f) for f in stats.get("faces", [])
             if abs(float(np.mean(np.asarray(f)[:, 2])) - z_maxc) <= half]
    figure_mesh_closeup(stem, faces, half, stats["n_nodes"], stats["n_elements"])


def main() -> int:
    root = catalogue_root()
    if root is None:
        print("No reference-blades catalogue found. Clone it (gitignored, not a committed submodule):")
        print("  git clone https://github.com/samhughes-fea/reference-blades.git reference_blades")
        print("or set REFERENCE_BLADES to that checkout.")
        return 0
    print(f"catalogue: {root}")
    iea15_geo = None
    for stem, rel, title in CASES:
        path = catalogue_yaml(root, rel)
        if not path.is_file():
            print(f"missing {path}")
            continue
        name, geo, z_of_s = load_windio(path)
        heading = f"{title} — {name}"
        marks = _landmarks(geo)
        tidy_legacy_pngs(_blade_dir(stem))
        figure_sections(stem, heading, geo)
        figure_planforms(stem, heading, geo, marks, figures=FIGURES, repo=REPO)
        walls = _wall_stations(geo, z_of_s)
        faces, kinds = _wall_faces(walls, geo=geo)
        if walls:
            _dump_debug_mesh(stem, walls, faces, kinds, geo)
        if faces:
            figure_projected(stem, heading, faces, kinds)
            figure_closeups(stem, heading, walls, marks, geo)
            figure_cutaways(stem, heading, walls, marks, geo)
        else:
            print(f"{stem}: no feasible section_at stations; span views stay OML-only")
        if stem == GMSH_STEM:
            iea15_geo = geo
    if iea15_geo is not None:
        gmsh_shell_probe(iea15_geo, GMSH_STEM)
    else:
        print("IEA 15 MW not loaded; skip Gmsh probe.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
