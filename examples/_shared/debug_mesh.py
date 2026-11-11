"""OBJ / VTK / NPZ dumps of a mapped loft (same arrays as the PNG strips)."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def as_quad(face: np.ndarray) -> np.ndarray:
    pts = np.asarray(face, float)
    if len(pts) == 3:
        return np.vstack([pts, pts[-1]])
    return pts[:4]


def seam_stats(walls: list) -> tuple[float, float]:
    """Max consecutive-loop-endpoint gap and max web-end-to-skin gap (SHARED = 0)."""
    max_loop, max_web = 0.0, 0.0
    for _z, tmpl, samples in walls:
        for side in ("suction", "pressure"):
            chain = [b for b in tmpl.loop if b.side == side and not b.wrap]
            for a, b in zip(chain, chain[1:]):
                if a.name not in samples or b.name not in samples:
                    continue
                gap = float(np.linalg.norm(samples[a.name][-1] - samples[b.name][0]))
                max_loop = max(max_loop, gap)
        skin_bits = [samples[b.name] for b in tmpl.loop if b.name in samples]
        if not skin_bits:
            continue
        skin = np.vstack(skin_bits)
        for web in tmpl.webs:
            if web.name not in samples:
                continue
            pts = np.asarray(samples[web.name], float)
            for end in (pts[0], pts[-1]):
                gap = float(np.min(np.linalg.norm(skin - end, axis=1)))
                max_web = max(max_web, gap)
    return max_loop, max_web


def _te_xyz(samples: dict, oml_te=None) -> np.ndarray | None:
    """TE sample: nearest OML ridge, else max chord-frac from LE (not blade-x)."""
    pts = []
    wrap = samples.get("te")
    if wrap is not None and len(wrap):
        pts.append(np.asarray(wrap, float))
    for key in ("te_pr", "te_su"):
        if key in samples and len(samples[key]):
            arr = np.asarray(samples[key], float)
            pts.append(arr)
    if not pts:
        return None
    stacked = np.vstack(pts)
    if oml_te is not None:
        o = np.asarray(oml_te, float).reshape(-1)[:2]
        return stacked[int(np.argmin(np.linalg.norm(stacked[:, :2] - o, axis=1)))]
    le = None
    if "le_su" in samples and len(samples["le_su"]):
        le = np.asarray(samples["le_su"], float)[-1]
    elif "le_pr" in samples and len(samples["le_pr"]):
        le = np.asarray(samples["le_pr"], float)[-1]
    if le is not None:
        d = np.linalg.norm(stacked[:, :2] - np.asarray(le, float)[:2], axis=1)
        return stacked[int(np.argmax(d))]
    return stacked[int(np.argmax(stacked[:, 0]))]


def span_metrics(walls: list, faces: np.ndarray, roles: np.ndarray, geo=None) -> dict:
    """Spanwise loft checks: TE chatter, wrap size, coverage holes, root web count."""
    if not walls:
        return {
            "max_te_chatter": 0.0, "max_te_oml": 0.0, "max_wrap_frac": 0.0,
            "n_coverage_holes": 0, "n_web_clusters": 0,
            "z_tip": 0.0, "web_z_max": {},
        }
    tes, chords, wraps, te_oml, oml_tes = [], [], [], [], []
    for zi, _tmpl, samples in walls:
        oml_te_pt = None
        if geo is not None:
            from implicit_geometry_kernel.blade.topology import _oml_te_xyz, section_cut_near

            cut = section_cut_near(geo, float(zi))
            ridge = _oml_te_xyz(cut) if cut is not None else None
            if ridge is None:
                oml = np.asarray(geo.oml_polyline(float(zi)), float)
                ridge = oml[int(np.argmax(oml[:, 0]))]
            oml_te_pt = np.asarray(ridge, float)[:2].copy()
            oml_tes.append(oml_te_pt)
        te = _te_xyz(samples, oml_te=oml_te_pt)
        le = None
        if "le_su" in samples:
            le = np.asarray(samples["le_su"], float)[-1]
        elif "le_pr" in samples:
            le = np.asarray(samples["le_pr"], float)[-1]
        chord = 1.0
        if te is not None and le is not None:
            chord = max(float(np.linalg.norm(te[:2] - le[:2])), 1e-9)
        chords.append(chord)
        tes.append(te[:2] if te is not None else np.zeros(2))
        wrap = samples.get("te")
        if wrap is None or len(wrap) < 2:
            wraps.append(0.0)
        else:
            wpts = np.asarray(wrap, float)
            slen = float(np.sum(np.linalg.norm(np.diff(wpts, axis=0), axis=1)))
            wraps.append(slen / chord)
        if geo is not None:
            if te is not None:
                te_oml.append(float(np.linalg.norm(te[:2] - oml_tes[-1])) / chord)
            else:
                te_oml.append(0.0)
        else:
            oml_tes.append(tes[-1].copy())
            te_oml.append(0.0)
    tes = np.asarray(tes)
    oml_tes = np.asarray(oml_tes)
    chords_a = np.asarray(chords)
    chatter = 0.0
    if len(tes) > 1:
        jumps = np.linalg.norm(np.diff(tes, axis=0) - np.diff(oml_tes, axis=0), axis=1)
        denom = np.maximum(0.5 * (chords_a[:-1] + chords_a[1:]), 1e-9)
        chatter = float(np.max(jumps / denom))
    fz = faces[:, 0, 2] if len(faces) else np.zeros(0)
    holes = 0
    for i in range(len(walls) - 1):
        n0, n1 = set(walls[i][2]), set(walls[i + 1][2])
        common = [nm for nm in (n0 & n1) if not str(nm).startswith("web") and str(nm) != "te"]
        z0, z1 = float(walls[i][0]), float(walls[i + 1][0])
        mask = (fz >= z0 - 1e-9) & (fz <= z1 + 1e-9) if len(fz) else np.zeros(0, bool)
        present = set(roles[mask].tolist()) if np.any(mask) else set()
        if any(nm not in present for nm in common):
            holes += 1
    webs0 = [n for n in walls[0][2] if str(n).startswith("web") or str(n).endswith("_web")
             or str(n) in ("fore_web", "rear_web", "web_te")]
    if not webs0:
        webs0 = [n for n in walls[0][2] if "web" in str(n).lower()]
    web_z_max: dict[str, float] = {}
    for zi, _tmpl, samples in walls:
        for name in samples:
            if "web" in str(name).lower():
                key = str(name)
                web_z_max[key] = max(web_z_max.get(key, float(zi)), float(zi))
    return {
        "max_te_chatter": chatter,
        "max_te_oml": float(max(te_oml) if te_oml else 0.0),
        "max_wrap_frac": float(max(wraps) if wraps else 0.0),
        "n_coverage_holes": int(holes),
        "n_web_clusters": int(len(set(webs0))),
        "z_tip": float(walls[-1][0]),
        "web_z_max": web_z_max,
    }


def _mesh_faces(faces: list[np.ndarray], kinds: list[str]) -> tuple[list[np.ndarray], list[str]]:
    mesh, roles = [], []
    for face, kind in zip(faces, kinds):
        pts = np.asarray(face, float)
        if len(pts) < 3:
            continue
        mesh.append(pts)
        roles.append(str(kind))
    return mesh, roles


def write_loft_obj(path: Path, faces: list[np.ndarray], kinds: list[str]) -> Path:
    """OBJ groups: skins/webs as quads, ribs/tip as n-gons (same arrays as the PNGs)."""
    mesh, roles = _mesh_faces(faces, kinds)
    groups: dict[str, list[np.ndarray]] = {}
    for face, kind in zip(mesh, roles):
        groups.setdefault(kind, []).append(face)
    lines: list[str] = []
    base = 0
    for name, flist in groups.items():
        lines.append(f"g {name}")
        for face in flist:
            for p in face:
                lines.append(f"v {float(p[0]):g} {float(p[1]):g} {float(p[2]):g}")
            n = len(face)
            lines.append("f " + " ".join(str(base + i + 1) for i in range(n)))
            base += n
    path = Path(path)
    path.write_text(("\n".join(lines) + "\n") if lines else "", encoding="utf-8")
    return path


def _role_rgb(role: str) -> tuple[float, float, float]:
    """Same RETRO_PANELS hash as the PNG loft (0–1 RGB for VTK COLOR_SCALARS)."""
    from matplotlib.colors import to_rgb

    from _shared.plotting import panel_colour

    return tuple(float(c) for c in to_rgb(panel_colour(role)))


def write_loft_vtk(path: Path, faces: list[np.ndarray], kinds: list[str]) -> Path:
    """Legacy ASCII VTK POLYDATA: PNG panel RGB plus a ``role`` name field."""
    mesh, roles = _mesh_faces(faces, kinds)
    n_pts = int(sum(len(f) for f in mesh))
    n_cells = len(mesh)
    conn = int(sum(len(f) + 1 for f in mesh))
    lines = [
        "# vtk DataFile Version 3.0",
        "mapped loft",
        "ASCII",
        "DATASET POLYDATA",
        f"POINTS {n_pts} float",
    ]
    for face in mesh:
        for p in face:
            lines.append(f"{float(p[0]):.9g} {float(p[1]):.9g} {float(p[2]):.9g}")
    lines.append(f"POLYGONS {n_cells} {conn}")
    base = 0
    for face in mesh:
        n = len(face)
        lines.append(f"{n} " + " ".join(str(base + i) for i in range(n)))
        base += n
    if n_cells:
        lines.append(f"CELL_DATA {n_cells}")
        lines.append("COLOR_SCALARS rgb 3")
        for kind in roles:
            r, g, b = _role_rgb(kind)
            lines.append(f"{r:.6f} {g:.6f} {b:.6f}")
        lines.append("FIELD attributes 1")
        lines.append(f"role 1 {n_cells} string")
        for kind in roles:
            lines.append(kind.replace(" ", "_"))
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _loft_quads(faces: list[np.ndarray], kinds: list[str]) -> tuple[np.ndarray, np.ndarray]:
    loft_faces, loft_kinds = [], []
    for face, kind in zip(faces, kinds):
        pts = np.asarray(face, float)
        key = str(kind)
        if len(pts) > 4 or "rib" in key.lower() or key.startswith("tip"):
            continue
        loft_faces.append(pts)
        loft_kinds.append(kind)
    if loft_faces:
        return np.stack([as_quad(f) for f in loft_faces]), np.asarray(loft_kinds, dtype="U64")
    return np.zeros((0, 4, 3), float), np.zeros((0,), dtype="U64")


def _rib_ngons(walls: list, geo) -> tuple[list[np.ndarray], list[str]]:
    """Authored rib n-gons and the tip closure at loft stations (OBJ/VTK; figures may ear-clip)."""
    from implicit_geometry_kernel.blade.topology import rib_cell_polygons_from_samples

    out_f, out_k = [], []
    if geo is None or not walls:
        return out_f, out_k
    z0, z1 = float(walls[0][0]), float(walls[-1][0])
    seen: set[str] = set()
    faces = [(float(m.z), m, str(m.name)) for m in geo.transverse]
    faces.append((float(geo.tip_closure_z), str(geo.tip_closure.name), str(geo.tip_closure.name)))
    for zm, member, key in faces:
        if not (z0 - 1e-9 <= zm <= z1 + 1e-9):
            continue
        _z, tmpl, samples = min(walls, key=lambda s, z=zm: abs(float(s[0]) - z))
        if key in seen:
            continue
        seen.add(key)
        for name, ring in rib_cell_polygons_from_samples(tmpl, samples, member):
            out_f.append(ring)
            out_k.append(name)
    return out_f, out_k


def _dump_faces(walls: list, faces: list[np.ndarray], kinds: list[str], geo=None
                ) -> tuple[list[np.ndarray], list[str]]:
    loft_f, loft_k, rib_pass_f, rib_pass_k = [], [], [], []
    for face, kind in zip(faces, kinds):
        key = str(kind)
        if "rib" in key.lower() or key.startswith("tip"):
            rib_pass_f.append(face)
            rib_pass_k.append(kind)
        else:
            loft_f.append(face)
            loft_k.append(kind)
    if geo is not None:
        rib_f, rib_k = _rib_ngons(walls, geo)
        return loft_f + rib_f, loft_k + rib_k
    return loft_f + rib_pass_f, loft_k + rib_pass_k


def dump_debug_mesh(out: Path, walls: list, faces: list[np.ndarray], kinds: list[str],
                    *, label: str, repo: Path, geo=None, stem: str = "mids",
                    write_samples: bool = True) -> None:
    """Write ``{stem}.obj`` / ``{stem}.vtk`` and optionally ``samples.npz``."""
    if not walls:
        return
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    mesh, roles = _mesh_faces(*_dump_faces(walls, faces, kinds, geo))
    obj_path = write_loft_obj(out / f"{stem}.obj", mesh, roles)
    vtk_path = write_loft_vtk(out / f"{stem}.vtk", mesh, roles)
    groups = {k for k in roles}
    print(f"wrote {obj_path.relative_to(repo)}")
    print(f"wrote {vtk_path.relative_to(repo)}")
    if not write_samples:
        print(f"{label}: {len(mesh)} {stem} faces, {len(groups)} groups")
        return
    z = np.array([float(w[0]) for w in walls], float)
    names = np.empty(len(walls), dtype=object)
    xyz = np.empty(len(walls), dtype=object)
    for i, (_zi, _t, samples) in enumerate(walls):
        names[i] = tuple(samples)
        xyz[i] = {k: np.asarray(v, float) for k, v in samples.items()}
    n_chord = int(walls[0][1].n_chord)
    face_arr, quad_roles = _loft_quads(faces, kinds)
    polygons = np.empty(len(mesh), dtype=object)
    for i, face in enumerate(mesh):
        polygons[i] = np.asarray(face, float)
    max_loop, max_web = seam_stats(walls)
    metrics = span_metrics(walls, face_arr, quad_roles, geo=geo)
    npz_path = out / "samples.npz"
    np.savez_compressed(
        npz_path,
        z=z,
        n_chord=np.int32(n_chord),
        names=names,
        xyz=xyz,
        faces=face_arr,
        roles=quad_roles,
        polygons=polygons,
        polygon_roles=np.asarray(roles, dtype="U64"),
        max_loop_seam=np.float64(max_loop),
        max_web_skin=np.float64(max_web),
        max_te_chatter=np.float64(metrics["max_te_chatter"]),
        max_te_oml=np.float64(metrics["max_te_oml"]),
        max_wrap_frac=np.float64(metrics["max_wrap_frac"]),
        n_coverage_holes=np.int32(metrics["n_coverage_holes"]),
        n_web_clusters=np.int32(metrics["n_web_clusters"]),
        z_tip=np.float64(metrics["z_tip"]),
        web_z_max=np.asarray(metrics["web_z_max"], dtype=object),
    )
    print(f"{label}: {len(walls)} stations, {len(groups)} groups, "
          f"{len(mesh)} polygons, max_loop_seam={max_loop:.3g} max_web_skin={max_web:.3g} "
          f"te_chatter={metrics['max_te_chatter']:.3g} te_oml={metrics['max_te_oml']:.3g} "
          f"wrap_frac={metrics['max_wrap_frac']:.3g} holes={metrics['n_coverage_holes']}")
    print(f"wrote {npz_path.relative_to(repo)}")
