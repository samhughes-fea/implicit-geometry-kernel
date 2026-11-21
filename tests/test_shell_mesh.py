"""The kernel's shell mesh: all-quad skin and webs, triangle ribs and tip, one conforming sheet."""

from collections import Counter

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, PlyDrop, Station, Stringer, TransverseMember
from implicit_geometry_kernel.mesh.quality import mesh_check
from implicit_geometry_kernel.mesh.shell import TIP_H_OVER_C, mesh_stations, shell_mesh


def _geo(key="2D-F-LC-4R"):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, key, 0.004, pivot_fraction=0.3,
                         spar_cap_thickness=0.008, shear_web_thickness=0.004)


_REFINED = {"h": 0.1, "region_h": {"le_su": 0.05, "le_pr": 0.05}, "refine": [(0.0, 2.0, 0.05)]}


@pytest.fixture(scope="module", params=[{"n_span": 11}, {"h": 0.08}, _REFINED], ids=["n_span", "h", "refined"])
def mesh(request):
    return _geo(), shell_mesh(_geo(), **request.param)


def test_elements_are_quads_and_triangles(mesh):
    _geo_, m = mesh
    for e, surface in zip(m.elements, m.surface):
        assert len(set(e)) == len(e)
        if surface in ("skin", "web"):
            assert len(e) == 4
        else:
            assert len(e) == 3
    assert {"skin", "web", "rib", "tip_closure"} == set(m.surface)


def test_one_conforming_sheet_open_only_at_the_root(mesh):
    _geo_, m = mesh
    use = m.edge_use()
    free = [e for e, n in use.items() if n == 1]
    assert free
    for a, b in free:
        assert m.nodes[a, 2] == 0.0 and m.nodes[b, 2] == 0.0


def test_rib_and_tip_perimeters_are_skin_and_web_nodes(mesh):
    _geo_, m = mesh
    shell_nodes = {k for e, s in zip(m.elements, m.surface) if s in ("skin", "web") for k in e}
    for kind in ("rib", "tip_closure"):
        edges: Counter = Counter()
        for e, s in zip(m.elements, m.surface):
            if s == kind:
                for a, b in zip(e, e[1:] + e[:1]):
                    edges[(min(a, b), max(a, b))] += 1
        perimeter = {k for (a, b), n in edges.items() if n == 1 for k in (a, b)}
        assert perimeter and perimeter <= shell_nodes, kind


def test_mesh_quality_is_reported_per_partition(mesh):
    _geo_, m = mesh
    assert m.quality
    names = {row.partition for row in m.quality}
    assert names == set(m.partition)
    kinds = {}
    for part, surface in zip(m.partition, m.surface):
        kinds.setdefault(part, set()).add(surface)
    for row in m.quality:
        assert row.n > 0
        assert row.aspect_min >= 1.0
        assert row.aspect_max >= row.aspect_min
        assert row.jacobian_max <= 1.0 + 1e-12
        if kinds[row.partition] <= {"skin", "web"}:
            assert row.inverted == 0
            assert row.jacobian_min > 0.0


def test_elements_carry_their_partition_data(mesh):
    geo, m = mesh
    assert np.all(m.thickness > 0.0)
    loop = {r.name for z in m.stations for r in geo.section_at(z).shell_graph().loop}
    for surface, part, ecc, layers, lt in zip(m.surface, m.partition, m.eccentricity, m.layers, m.layer_thicknesses):
        assert len(layers) == len(lt)
        if surface == "skin":
            assert part in loop
            assert (ecc > 0.0) == part.startswith("spar_cap")
    ribs = {r.name: r.thickness for r in geo.transverse}
    for surface, part, t in zip(m.surface, m.partition, m.thickness):
        if surface == "rib":
            assert t == ribs[part]
        if surface == "tip_closure":
            assert t == pytest.approx(geo.tip_closure_thickness())


def test_stations_include_ribs_and_tip():
    geo = _geo()
    zs = mesh_stations(geo, 5)
    for r in geo.transverse:
        assert any(abs(z - r.z) < 1e-9 for z in zs)
    assert zs[0] == 0.0 and zs[-1] == pytest.approx(geo.tip_closure_z)
    with pytest.raises(ValueError):
        mesh_stations(geo)


def _catalogue():
    import os
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    env = os.environ.get("REFERENCE_BLADES_DIR")
    for path in ([Path(env)] if env else []) + [root / "reference_blades", root.parent / "reference-blades"]:
        if (path / "wind" / "nrel5mw" / "nrel5mw.yaml").is_file() or (path / "nrel5mw" / "nrel5mw.yaml").is_file():
            return path
    return None


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
@pytest.mark.parametrize("rel", ["nrel5mw/nrel5mw.yaml", "iea22/IEA-22-280-RWT.yaml"])
def test_reference_blade_mesh_is_all_quad_and_conforming(rel):
    import runpy
    from pathlib import Path

    ns = runpy.run_path(str(Path(__file__).resolve().parents[1] / "examples" / "blade" / "08_reference_blades"
                            / "reference_blades.py"))
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_catalogue(), rel))
    m = shell_mesh(geo, 21)
    assert all(len(e) == 4 for e, s in zip(m.elements, m.surface) if s in ("skin", "web"))
    assert all(len(e) == 3 for e, s in zip(m.elements, m.surface) if s in ("rib", "tip_closure"))
    z_root = float(geo.stations[0].z)
    web = {k for e, s in zip(m.elements, m.surface) if s == "web" for k in e}
    for (a, b), n in m.edge_use().items():
        if n != 1:
            continue
        # open only along the root, or along a web's start or end line
        on_root = abs(m.nodes[a, 2] - z_root) < 1e-9 and abs(m.nodes[b, 2] - z_root) < 1e-9
        assert on_root or (a in web and b in web)


@pytest.mark.skipif(_catalogue() is None, reason="reference-blades catalogue not cloned")
def test_nrel_te_cap_seam_shares_endpoints_on_probe_stations():
    """BG-049: cap_su_1–te_su (and pressure twin) meet; no off-root opens in the band."""
    import runpy
    from pathlib import Path

    from implicit_geometry_kernel.mesh.shell import _freeze_sagitta, _sample_stations
    from implicit_geometry_kernel.mesh.tables import off_root_open_edges

    ns = runpy.run_path(str(Path(__file__).resolve().parents[1] / "examples" / "blade" / "08_reference_blades"
                            / "reference_blades.py"))
    _name, geo, _ = ns["load_windio"](ns["catalogue_yaml"](_catalogue(), "nrel5mw/nrel5mw.yaml"))
    zs = list(mesh_stations(geo, h_over_c=TIP_H_OVER_C))
    zs, _ = _freeze_sagitta(geo, zs, None, TIP_H_OVER_C)
    zs = [float(z) for z in zs if 27.5 <= float(z) <= 38.5]
    _cuts, tmpls, smps = _sample_stations(geo, zs, 8, None, TIP_H_OVER_C, None)
    for tmpl, smp in zip(tmpls, smps):
        for side in ("suction", "pressure"):
            prev = None
            for branch in tmpl.loop:
                if branch.wrap or branch.side != side or branch.name not in smp:
                    continue
                if prev is not None:
                    gap = float(np.linalg.norm(smp[prev][-1] - smp[branch.name][0]))
                    assert gap <= 1e-12, (prev, branch.name, gap)
                prev = branch.name
    # Match the sample sizing: bare zs= turns the tip-h/c default off.
    mesh = shell_mesh(geo, zs=zs, h_over_c=TIP_H_OVER_C)
    # Probe free ends leave station opens at zs[0]/zs[-1]; the TE seam must not.
    groups = off_root_open_edges(mesh)
    assert not any(row["kind"] == "spanwise" for row in groups)
    for row in groups:
        assert row["kind"] == "station"
        assert min(abs(row["z0"] - zs[0]), abs(row["z0"] - zs[-1])) <= 1e-6


def test_vtk_keeps_triangles(tmp_path, mesh):
    from implicit_geometry_kernel.mesh import write_vtk

    _geo_, m = mesh
    path = tmp_path / "shell.vtk"
    parts = write_vtk(m, path)
    text = path.read_text().splitlines()
    i = text.index(f"CELL_TYPES {len(m.elements)}")
    types = [int(t) for t in text[i + 1:i + 1 + len(m.elements)]]
    assert types.count(9) == len(m.quads) and types.count(5) == len(m.triangles)
    assert types.count(7) == len(m.polygons)
    assert set(parts) == set(m.partition)
    assert f"POINTS {len(m.nodes)} double" in text


def test_gmsh_path_is_deprecated():
    pytest.importorskip("gmsh")
    from implicit_geometry_kernel.mesh.gmsh import triangulate_ring

    with pytest.deprecated_call():
        triangulate_ring(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]))


def _corner_angles(face):
    out = []
    for i in range(len(face)):
        a, b = face[i - 1] - face[i], face[(i + 1) % len(face)] - face[i]
        out.append(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1.0, 1.0))))
    return out


def test_sized_quads_have_no_slivers():
    m = shell_mesh(_geo(), h=0.08)
    angles = [a for i in m.quads for a in _corner_angles(m.nodes[list(m.elements[i])])]
    # A tapering station edge is a little off a right angle. 20° is the sliver line.
    assert min(angles) > 20.0 and max(angles) < 165.0


def test_leading_edge_and_span_refinement():
    geo = _geo()
    m = shell_mesh(geo, **_REFINED)
    zs = np.array(m.stations)
    assert np.all(np.diff(zs[zs <= 2.0 + 1e-9]) <= 0.05 + 1e-9)
    assert np.diff(zs[zs >= 2.0]).max() > 0.05
    def root_edges(prefix):
        out = []
        for i in m.quads:
            if m.name[i].startswith(prefix):
                e = m.nodes[list(m.elements[i])]
                for k in range(4):
                    if e[k, 2] == 0.0 and e[k - 1, 2] == 0.0:
                        out.append(float(np.linalg.norm(e[k] - e[k - 1])))
        return out

    le, cap = root_edges("le_"), root_edges("cap_")
    assert le and max(le) <= 0.05 * (1 + 1e-9) and max(cap) > 0.05
    with pytest.raises(ValueError):
        shell_mesh(geo, n_span=5, region_h={"le_su": 0.05})
    with pytest.raises(ValueError):
        mesh_stations(geo, 5, refine=[(0.0, 1.0, 0.0)])


def _area_vector(p):
    return 0.5 * np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)


def test_normals_out_of_the_skin_to_the_te_on_webs_and_along_z_on_ribs(mesh):
    _geo_, m = mesh
    skin = [i for i, s in enumerate(m.surface) if s == "skin"]
    runs: dict = {}
    for i in skin:
        e = m.elements[i]
        for a, b in zip(e, e[1:] + e[:1]):
            runs.setdefault((min(a, b), max(a, b)), []).append((a, b))
    for key, uses in runs.items():
        if len(uses) == 2:
            assert uses[0] != uses[1], key  # neighbours run a shared edge in opposite directions
    flux = 0.0
    for i in skin:
        p = m.nodes[list(m.elements[i])]
        flux += float(_area_vector(p)[:2] @ (p.mean(axis=0)[:2] - np.array([0.0, 0.0])))
    assert flux > 0.0
    for e, s in zip(m.elements, m.surface):
        a = _area_vector(m.nodes[list(e)])
        if s == "web":
            assert a[0] > 0.0  # untwisted blade: the TE is +x
        elif s in ("rib", "tip_closure"):
            assert a[2] > 0.0


def _plain(scale: float = 1.0, length: float = 4.0, chord0: float = 2.0, chord1: float = 1.0,
           layout: str = "0A-LC-0R", **kwargs):
    stations = [
        Station(0.0, chord0 * scale, Airfoil.naca4("0012"), pitch_axis=0.5),
        Station(length * scale, chord1 * scale, Airfoil.naca4("0012"), pitch_axis=0.5),
    ]
    return BladeGeometry(stations, layout, 0.004, pivot_fraction=0.5, **kwargs)


def _skin_edges(m, surface: str = "skin"):
    """Spanwise edges (dz) and edges that lie in a station."""
    span, chord = [], []
    for i in m.quads:
        if m.surface[i] != surface:
            continue
        e = m.nodes[list(m.elements[i])]
        for k in range(4):
            d = e[k] - e[k - 1]
            (span if abs(d[2]) > 1e-8 else chord).append(float(np.linalg.norm(d)))
    return span, chord


def test_h_over_c_scales_with_the_blade():
    fine, coarse = _plain(1.0), _plain(2.0)
    a = shell_mesh(fine, h_over_c=0.3)
    b = shell_mesh(coarse, h_over_c=0.3)
    assert len(a.stations) == len(b.stations)
    assert np.asarray(b.stations) == pytest.approx(np.asarray(a.stations) * 2.0, abs=1e-6)

    def ratios(m, geo):
        zs = [s.z for s in geo.stations]
        cs = [s.chord for s in geo.stations]
        out = []
        for i in m.quads:
            if m.surface[i] != "skin":
                continue
            e = m.nodes[list(m.elements[i])]
            for k in range(4):
                p, q = e[k], e[k - 1]
                if abs(p[2] - q[2]) > 1e-8:
                    continue
                chord = float(np.interp(p[2], zs, cs))
                out.append(float(np.linalg.norm(p - q)) / chord)
        return np.asarray(out)

    ra, rb = ratios(a, fine), ratios(b, coarse)
    assert ra.size and rb.size
    assert ra.max() <= 0.3 * (1 + 1e-6) and rb.max() <= 0.3 * (1 + 1e-6)
    assert np.median(ra) == pytest.approx(np.median(rb), rel=0.05)


def test_h_bounds_the_edge_and_the_sharp_te_stays_a_shared_corner():
    h = 0.4
    geo = _plain(length=6.0, chord1=1.5)
    m = shell_mesh(geo, h=h)
    span, chord = _skin_edges(m)
    assert span and chord
    # The station step is h. On a taper the 3D edge is a little longer than that step.
    assert max(span) <= h * 1.02 and min(span) >= 0.5 * h * (1 - 1e-6)
    assert max(chord) <= h * 1.02
    capped = shell_mesh(geo, h=0.25, h_over_c=0.5)
    span_c, chord_c = _skin_edges(capped)
    assert max(span_c) <= 0.25 * 1.02 and max(chord_c) <= 0.25 * 1.02

    skin = [i for i in m.quads if m.surface[i] == "skin"]
    on_root = {k for i in skin for k in m.elements[i] if abs(m.nodes[k, 2]) < 1e-8}
    te = max(on_root, key=lambda k: float(m.nodes[k, 0]))
    sides = []
    for i in skin:
        elem = m.elements[i]
        if te not in elem:
            continue
        others = m.nodes[[k for k in elem if k != te]]
        y = others[:, 1]
        assert np.all(y >= -1e-5) or np.all(y <= 1e-5)
        same = [k for k in elem if k != te and abs(m.nodes[k, 2] - m.nodes[te, 2]) < 1e-8]
        assert same
        length = min(float(np.linalg.norm(m.nodes[k] - m.nodes[te])) for k in same)
        assert 0.25 * h <= length <= 1.5 * h
        sides.append(float(np.mean(y)))
    assert any(s > 1e-6 for s in sides) and any(s < -1e-6 for s in sides)


def test_a_rib_does_not_leave_a_sliver_against_a_regular_station():
    geo = _plain(length=10.0, chord1=1.0, transverse=(TransverseMember("rib", 2.003, 0.012),))
    zs = np.asarray(mesh_stations(geo, h=0.5))
    i = int(np.argmin(np.abs(zs - 2.003)))
    assert zs[i] == pytest.approx(2.003, abs=1e-6)
    assert zs[i] - zs[i - 1] >= 0.25 - 1e-6
    assert zs[i + 1] - zs[i] >= 0.25 - 1e-6


def test_a_dropped_web_ends_on_its_station_and_the_tip_uses_the_webs_that_remain():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, pivot_fraction=0.3,
                        spar_cap_thickness=0.008, shear_web_thickness=0.004,
                        web_offsets=(0.4, 1.2))
    m = shell_mesh(geo, h=1.2)
    lost = geo.section_at(10.0).dropped_web_indices
    assert lost
    idx = int(lost[0])
    lo, hi = 0.0, 10.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if idx in geo.section_at(mid).dropped_web_indices:
            hi = mid
        else:
            lo = mid
    assert any(abs(z - hi) < 1e-3 for z in m.stations)
    reach: dict[str, float] = {}
    for elem, surface, name in zip(m.elements, m.surface, m.name):
        if surface != "web":
            continue
        reach[name] = max(reach.get(name, 0.0), float(m.nodes[list(elem), 2].max()))
    assert any(z <= hi + 1e-2 for z in reach.values())
    assert any(z > hi + 0.2 for z in reach.values())
    tip_cells = {m.name[i] for i, s in enumerate(m.surface) if s == "tip_closure"}
    assert len(tip_cells) == 2


def test_a_short_cap_keeps_short_edges():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0018"), pitch_axis=0.3),
        Station(4.0, 1.6, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(stations, "1B-F-LC-0R", 0.004, pivot_fraction=0.3,
                        spar_cap_thickness=0.008, shear_web_thickness=0.004,
                        cap_windows=((0.46, 0.50),))
    m = shell_mesh(geo, h=0.5)
    edges = []
    for i in m.quads:
        if not str(m.name[i]).startswith("cap"):
            continue
        e = m.nodes[list(m.elements[i])]
        for k in range(4):
            if abs(e[k, 2] - e[k - 1, 2]) < 1e-8:
                edges.append(float(np.linalg.norm(e[k] - e[k - 1])))
    assert edges and max(edges) < 0.25


def test_flapwise_web_feet_are_skin_nodes():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("2412"), pitch_axis=0.3, twist=12.0),
        Station(4.0, 1.6, Airfoil.naca4("2412"), pitch_axis=0.3, twist=4.0),
    ]
    geo = BladeGeometry(stations, "2D-F-LC-0R", 0.004, pivot_fraction=0.3,
                        spar_cap_thickness=0.008, shear_web_thickness=0.004)
    m = shell_mesh(geo, h=0.8)
    skin = {k for e, s in zip(m.elements, m.surface) if s == "skin" for k in e}
    by_z: dict[float, set[int]] = {}
    for e, s in zip(m.elements, m.surface):
        if s != "web":
            continue
        for k in e:
            if k in skin:
                by_z.setdefault(round(float(m.nodes[k, 2]), 6), set()).add(k)
    assert by_z
    assert all(len(ids) >= 2 for ids in by_z.values())


def test_ply_drop_and_stringer_add_stations_and_no_new_skin_regions():
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.5),
        Station(8.0, 1.0, Airfoil.naca4("0012"), pitch_axis=0.5),
    ]
    drop = PlyDrop("pad", 0.2, 0.45, 3.15, 3.40, 0.002)
    stringer = Stringer("stiff", 0.35, 0.01, 0.02, extent=(2.15, 5.4))
    bare = BladeGeometry(stations, "0A-LC-0R", 0.004, pivot_fraction=0.5)
    full = BladeGeometry(stations, "0A-LC-0R", 0.004, pivot_fraction=0.5,
                         ply_drops=(drop,), stringers=(stringer,))
    zs = mesh_stations(full, h=1.0)
    for z in (2.15, 3.15, 3.40, 5.4):
        assert any(abs(s - z) < 1e-6 for s in zs)
    m = shell_mesh(full, h=1.0)
    base = shell_mesh(bare, h=1.0)
    assert "stringer" not in m.surface
    assert all("stiff" not in n and "pad" not in n for n in m.name)
    skin = {n for n, s in zip(m.name, m.surface) if s == "skin"}
    assert skin == {n for n, s in zip(base.name, base.surface) if s == "skin"}


def _chordwise_over_chord(m, geo):
    zs = [s.z for s in geo.stations]
    cs = [s.chord for s in geo.stations]
    out = []
    for i in m.quads:
        if m.surface[i] != "skin":
            continue
        e = m.nodes[list(m.elements[i])]
        for k in range(4):
            p, q = e[k], e[k - 1]
            if abs(p[2] - q[2]) > 1e-8:
                continue
            chord = float(np.interp(p[2], zs, cs))
            out.append((float(p[2]), float(np.linalg.norm(p - q)) / chord))
    return out


def test_the_default_edge_is_five_percent_of_chord_and_every_face_uses_the_sheet():
    geo = _geo()
    m = shell_mesh(geo)
    assert m.recipe is not None
    assert m.recipe.tip_h_over_c == TIP_H_OVER_C
    assert m.recipe.h is None and m.recipe.h_over_c is None and m.recipe.n_span is None
    assert len(m.recipe.edge) == len(m.stations)
    tip_c = float(np.interp(m.stations[-1], [s.z for s in geo.stations], [s.chord for s in geo.stations]))
    assert m.recipe.edge[-1] == pytest.approx(TIP_H_OVER_C * tip_c)
    ratios = _chordwise_over_chord(m, geo)
    assert ratios
    assert max(r for _z, r in ratios) <= TIP_H_OVER_C * 1.02
    tip = [r for z, r in ratios if abs(z - m.stations[-1]) < 1e-6]
    assert tip and np.median(tip) == pytest.approx(TIP_H_OVER_C, rel=0.25)

    sheet = {k for e, s in zip(m.elements, m.surface) if s in ("skin", "web") for k in e}
    assert sheet and max(sheet) < m.n_source
    interior = {k for e, s in zip(m.elements, m.surface) if s in ("rib", "tip_closure") for k in e} - sheet
    assert all(k >= m.n_source for k in interior)
    check = mesh_check(m)
    assert check.boundary_in_source
    assert check.open_edges > 0
    assert check.quad_angle_min == pytest.approx(float(np.min(m.element_quality["angle_min"][list(m.quads)])))
    for kind in ("rib", "tip_closure"):
        assert any(s == kind for s in m.surface)

    q = m.element_quality
    assert len(q) == len(m.elements)
    assert np.all(np.isfinite(q["aspect"]))
    assert np.all((q["skew"] >= 0.0) & (q["skew"] <= 1.0))
    assert np.allclose(q["shell_volume"], q["area"] * m.thickness)
    assert np.all(q["volume_ratio"] >= 1.0 - 1e-9)
    i = m.quads[0]
    face = m.nodes[list(m.elements[i])]
    assert q["angle_min"][i] == pytest.approx(min(_corner_angles(face)), abs=1e-6)
    for row in m.quality:
        idx = [k for k, part in enumerate(m.partition) if part == row.partition]
        assert row.aspect_min == pytest.approx(float(np.min(q["aspect"][idx])))
        assert row.skew_max == pytest.approx(float(np.max(q["skew"][idx])))


def test_an_explicit_size_or_count_leaves_the_chord_default_off(tmp_path):
    geo = _plain(length=10.0, chord1=1.0)
    sized = shell_mesh(geo, h=1.5)
    assert sized.recipe is not None and sized.recipe.tip_h_over_c is None and sized.recipe.h == 1.5
    steps = np.diff(sized.stations)
    assert steps.size and max(steps) <= 1.5 * 1.02 and min(steps) >= 0.75 - 1e-6

    counted = shell_mesh(_geo(), 11)
    assert counted.recipe is not None and counted.recipe.tip_h_over_c is None and counted.recipe.edge == ()
    assert counted.stations == mesh_stations(_geo(), 11)
    assert mesh_check(counted).boundary_in_source
    from implicit_geometry_kernel.mesh import write_vtk
    path = tmp_path / "shell.vtk"
    write_vtk(sized, path)
    text = path.read_text()
    assert "SCALARS aspect double 1" in text and "SCALARS shell_volume double 1" in text
    assert "SCALARS skew double 1" in text and "SCALARS volume_ratio double 1" in text
