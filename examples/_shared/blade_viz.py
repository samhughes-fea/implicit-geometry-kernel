"""Shared figure tree, suction filter, loft live-web fracs, and planform sampling."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from _shared.plotting import NAVY, RED, SERIES, TEAL, note, plt
from implicit_geometry_kernel.blade.topology import _chord_fracs, role_atom
from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

CROSS_SECTIONS = "cross_sections"
REGIONS = "regions"
REGION_CUTAWAYS = "region_cutaways"
PLANFORM = "planform"
FULL_BLADE = "full_blade"
VIZ_DIRS = (CROSS_SECTIONS, REGIONS, REGION_CUTAWAYS, PLANFORM, FULL_BLADE)
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
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
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
    y_top = float(np.max(te_x))
    _landmark_ticks(ax, marks, y_top)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    note(fig, caption)
    fig.tight_layout(rect=(0, 0.06, 1, 1))


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
