"""Shared figure tree, suction filter, loft live-web fracs, and planform sampling."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from _shared.plotting import (
    NAVY, RED, SERIES, TEAL, draw_cap_solids, draw_shell, note, panel_colour, plt,
)
from implicit_geometry_kernel.blade.network import wall_meta_from_cut
from implicit_geometry_kernel.blade.topology import (
    ShellTopologyTemplate, _chord_fracs, _cut_tc, role_atom, sample_region_polylines,
)
from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

CROSS_SECTIONS = "cross_sections"
CSG_MESH = "csg_mesh"
REGIONS = "regions"
REGION_CUTAWAYS = "region_cutaways"
PLANFORM = "planform"
FULL_BLADE = "full_blade"
VIZ_DIRS = (CROSS_SECTIONS, CSG_MESH, REGIONS, REGION_CUTAWAYS, PLANFORM, FULL_BLADE)
N_SECTION_GRID = 16
CSG_MESH_N_CHORD = 8
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
    """4×4 ``sections.png``, dual-scale ``01.png``–``16.png``, and ``csg_mesh`` singles."""
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


def _finish_csg_axes(ax, xlim, ylim) -> None:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.tick_params(labelsize=7)


def _paint_wall_meta(ax, fig, cut, z, xlim, ylim) -> None:
    field, _table = wall_meta_from_cut(cut)
    xs = np.linspace(xlim[0], xlim[1], 80)
    ys = np.linspace(ylim[0], ylim[1], 50)
    X, Y = np.meshgrid(xs, ys)
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, float(z))])
    phi = np.asarray(field.evaluate(pts).phi, float).reshape(X.shape)
    cf = ax.contourf(X, Y, phi, levels=16, cmap="RdBu_r")
    ax.contour(X, Y, phi, levels=[0.0], colors="k", linewidths=0.7)
    fig.colorbar(cf, ax=ax, shrink=0.72, pad=0.02, label="φ / m")
    ax.set_title(f"z = {z:.1f} m\nwall_meta φ", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def _paint_mesh_nodes(ax, geo, cut, z, xlim, ylim) -> None:
    tmpl = ShellTopologyTemplate.from_region_graph(
        geo.region_graph, float(z), n_chord=CSG_MESH_N_CHORD,
        live_web_fracs=live_web_fracs(cut) or None, t_c=_cut_tc(cut),
        dropped_web_indices=tuple(cut.dropped_web_indices),
    )
    samples = sample_region_polylines(cut, cut.shell_graph(), tmpl, geo.region_graph)
    for name, pts in samples.items():
        xy = np.asarray(pts, float)
        if len(xy) < 1:
            continue
        colour = panel_colour(name)
        ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=0.9, zorder=2)
        ax.scatter(xy[1:-1, 0], xy[1:-1, 1], s=10, c="k", zorder=3, linewidths=0)
        ax.scatter(xy[[0, -1], 0], xy[[0, -1], 1], s=22, c="k", zorder=4, linewidths=0)
    ax.set_title(f"z = {z:.1f} m\nconformal nodes (n_chord={CSG_MESH_N_CHORD})", fontsize=8, pad=6)
    _finish_csg_axes(ax, xlim, ylim)


def figure_csg_mesh(stem: str, geo, zs, cuts, *, figures: Path, repo: Path) -> None:
    """``csg_mesh/01.png``–``16.png``: wall_meta φ vs authored station mesh nodes."""
    for i, (z, cut) in enumerate(zip(zs, cuts), start=1):
        fig, (ax_csg, ax_mesh) = plt.subplots(1, 2, figsize=(8.4, 4.2))
        if cut is None:
            ax_csg.set_title(f"z = {z:.1f} m\ninfeasible", fontsize=8, pad=6)
            ax_mesh.set_title(f"z = {z:.1f} m\ninfeasible", fontsize=8, pad=6)
            for ax in (ax_csg, ax_mesh):
                ax.set_aspect("equal", adjustable="box")
                ax.set_xlabel("x / m")
                ax.set_ylabel("y / m")
                ax.tick_params(labelsize=7)
        else:
            xlim, ylim = _station_limits(cut)
            try:
                _paint_wall_meta(ax_csg, fig, cut, z, xlim, ylim)
            except ValueError:
                ax_csg.set_title(f"z = {z:.1f} m\ninfeasible", fontsize=8, pad=6)
                _finish_csg_axes(ax_csg, xlim, ylim)
            _paint_mesh_nodes(ax_mesh, geo, cut, z, xlim, ylim)
        fig.tight_layout()
        save_fig(fig, figures, stem, CSG_MESH, f"{i:02d}", repo=repo)


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


def _finish_planform(ax, fig, title: str, caption: str, marks, te_x) -> None:
    te_x = np.asarray(te_x, float)
    te_x = te_x[np.isfinite(te_x)]
    y_top = float(np.max(te_x)) if te_x.size else 1.0
    _landmark_ticks(ax, marks, y_top)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    note(fig, caption)
    try:
        fig.tight_layout(rect=(0, 0.06, 1, 1))
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
    fig, ax = plt.subplots(figsize=(10.5, 3.2))
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
    _finish_planform(
        ax, fig, f"{title} — planform axes (blade-frame x–z)",
        "LE/TE from oml_polyline. Pitch from interpolate offset; centroid from wall midlines; "
        "shear centre is the cell-area centroid where cells exist.",
        marks, te[:, 0],
    )
    save_fig(fig, figures, stem, PLANFORM, "axes", repo=repo)


def figure_planform_structure(stem: str, title: str, marks, data: dict, *, figures: Path,
                              repo: Path) -> None:
    le, te, pitch = data["le"], data["te"], data["pitch"]
    fig, ax = plt.subplots(figsize=(10.5, 3.2))
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
    _finish_planform(
        ax, fig, f"{title} — planform structure (webs, ribs)",
        "Web loci from conformal shell_graph() midlines at each sample z; dropped webs stop. "
        "Ribs are chordwise lines at each TransverseMember z (LE to TE).",
        marks, te[:, 0],
    )
    save_fig(fig, figures, stem, PLANFORM, "structure", repo=repo)
