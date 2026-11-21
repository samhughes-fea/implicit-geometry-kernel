"""Shell mesh layers: station staircase, on-section quads, shared web feet, rib vertices."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.blade.topology.map_loft import _step_quads
from implicit_geometry_kernel.mesh.shell import shell_mesh
from implicit_geometry_kernel.mesh.tables import mesh_failures, write_mesh_tables
from implicit_geometry_kernel.mesh.export import write_vtk


def _concept():
    """Example 15's self-similar two-web blade."""
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=8.0),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )


def _taper():
    """Constant edge, chord falling by five, so a region's arc leaves the hold band."""
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3, twist=8.0),
        Station(2.0, 0.4, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    return BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )


def _angles(face: np.ndarray) -> list[float]:
    out = []
    for i in range(len(face)):
        a, b = face[i - 1] - face[i], face[(i + 1) % len(face)] - face[i]
        out.append(float(np.degrees(np.arccos(np.clip(
            a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1.0, 1.0)))))
    return out


def _site_halves(faces) -> tuple[list, list]:
    """Step-patch faces split into the two sites by their centroid along the strip."""
    sites = [f for f in faces if f.constructor == "step_patch"]
    xs = [float(np.mean(f.face[:, 0])) for f in sites]
    mid = 0.5 * (min(xs) + max(xs))
    left, right = [], []
    for face, x in zip(sites, xs):
        angles = tuple(round(a, 6) for a in sorted(_angles(face.face)))
        (left if x < mid - 1e-9 else right).append(angles)
    return left, right


def test_flat_two_edge_patch_keeps_its_corners():
    coarse = np.column_stack([np.linspace(0, 1, 5), np.zeros(5), np.zeros(5)])
    fine = np.column_stack([np.linspace(0, 1, 7), np.zeros(7), np.full(7, 1.0)])
    stats: dict = {}
    faces = _step_quads("panel", coarse, fine, stats=stats)
    assert faces
    assert (stats["site0"], stats["site1"]) == (1, 2)
    assert {f.constructor for f in faces} <= {"step_patch", "midspan_split"}
    assert sum(f.constructor == "step_patch" for f in faces) == 8
    angles = [a for f in faces for a in _angles(f.face)]
    assert min(angles) >= 20.0 - 1e-6
    assert max(angles) <= 160.0 + 1e-6
    assert all("centroid" not in f.origins for f in faces)
    left, right = _site_halves(faces)
    assert left and right and sorted(left) == sorted(right)
    rev_stats: dict = {}
    rev = _step_quads("panel", fine, coarse, stats=rev_stats)
    assert (rev_stats["site0"], rev_stats["site1"]) == (1, 2)
    rev_angles = [a for f in rev for a in _angles(f.face)]
    assert min(rev_angles) >= 20.0 - 1e-6
    assert sorted(_site_halves(rev)[0]) == sorted(_site_halves(rev)[1])


def test_two_coarse_edges_are_both_sites():
    coarse = np.column_stack([np.linspace(0, 1, 3), np.zeros(3), np.zeros(3)])
    fine = np.column_stack([np.linspace(0, 1, 5), np.zeros(5), np.ones(5)])
    stats: dict = {}
    faces = _step_quads("te", coarse, fine, stats=stats)
    assert (stats["site0"], stats["site1"]) == (0, 1)
    assert faces and {f.constructor for f in faces} == {"step_patch"}
    assert all("centroid" not in f.origins for f in faces)
    angles = [a for f in faces for a in _angles(f.face)]
    assert min(angles) >= 20.0 - 1e-6
    left, right = _site_halves(faces)
    assert left and right and sorted(left) == sorted(right)


def test_self_similar_skin_stays_ruled_when_the_nose_fits():
    geo = _concept()
    zs = np.linspace(0.0, 10.0, 5)
    mesh = shell_mesh(geo, h_over_c=0.05, zs=zs)
    skin = [s for s in mesh.strips if not s.region.startswith("web") and s.constructor != "omitted"]
    assert skin
    assert {s.constructor for s in skin} == {"ruled"}
    assert all(s.n0 == s.n1 for s in skin)
    assert "centroid" not in mesh.node_origin


@pytest.fixture(scope="module")
def taper():
    return shell_mesh(_taper(), h=0.12)


def test_a_halving_arc_keeps_one_count(taper):
    """A region that exists on every station of the taper takes one block count."""
    shared = [s for s in taper.strips if s.n0 and s.n1]
    assert shared
    assert all(s.constructor != "step_patch" for s in shared)
    by: dict[str, set[int]] = {}
    for decision in taper.counts:
        by.setdefault(decision.region, set()).add(decision.count)
    assert by and all(len(counts) == 1 for counts in by.values())
    for strip in shared:
        if strip.n0 == strip.n1:
            assert strip.constructor in ("ruled", "midspan_split")
    assert float(np.max(taper.section_distance)) < 1e-8
    for i in taper.quads:
        if float(taper.element_quality["angle_min"][i]) < 20.0:
            assert taper.constructor[i] in ("fan", "transition")
    assert mesh_failures(taper) == []


def test_web_foot_matches_the_skin_and_a_rib_uses_source_nodes(taper):
    skin = {k for e, s in zip(taper.elements, taper.surface) if s == "skin" for k in e}
    web = {k for e, s in zip(taper.elements, taper.surface) if s == "web" for k in e}
    feet = [k for k in skin & web if taper.node_origin[k] in ("station", "span_mid")]
    assert feet
    for node in feet:
        skin_bays = {taper.bay[i] for i, e in enumerate(taper.elements)
                     if node in e and taper.surface[i] == "skin"}
        web_bays = {taper.bay[i] for i, e in enumerate(taper.elements)
                    if node in e and taper.surface[i] == "web"}
        assert skin_bays & web_bays
    for bay in {s.bay for s in taper.strips if s.trigger}:
        for strip in taper.strips:
            if strip.bay != bay or strip.constructor == "omitted" or not (strip.n0 and strip.n1):
                continue
            if strip.n0 == strip.n1:
                assert strip.constructor == "midspan_split"
            else:
                assert strip.constructor == "fan" and min(strip.n0, strip.n1) == 1
    assert any(s == "rib" for s in taper.surface)
    for element, surface, constructor in zip(taper.elements, taper.surface, taper.constructor):
        if surface not in ("rib", "tip_closure"):
            continue
        for node in element:
            if node < taper.n_source:
                assert taper.node_origin[node] in ("station", "span_mid")
            else:
                assert constructor == "interior"
                assert taper.node_origin[node] == "interior"
    assert all(len(e) == 4 for e, s in zip(taper.elements, taper.surface) if s in ("skin", "web"))
    assert all(len(e) == 3 for e, s in zip(taper.elements, taper.surface) if s in ("rib", "tip_closure"))


def test_mesh_tables_and_vtk_carry_the_record(tmp_path, taper):
    write_mesh_tables(taper, tmp_path)
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "quads under 20 deg: 0" in summary
    assert "open edges off the root: 0" in summary
    assert "t-junctions: 0" in summary
    assert "inverted elements: 0" in summary
    for name in ("nodes.csv", "elements.csv", "edges.csv", "counts.csv", "stations.csv",
                 "strips.csv", "corners.csv", "loops.csv", "failures.csv", "health_counts.csv"):
        assert (tmp_path / name).exists()
    header = (tmp_path / "elements.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "constructor" in header and "count0" in header
    strips = (tmp_path / "strips.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "chord_ratio" in strips and "site0" in strips and "site1" in strips
    corners = (tmp_path / "corners.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "element_id" in corners and "node_id" in corners and "angle_deg" in corners
    assert "yaw_deg" in corners and "pitch_deg" in corners and "roll_deg" in corners
    edge_header = (tmp_path / "edges.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "dihedral_deg" in edge_header and "dx" in edge_header
    assert "chain_yaw_deg" in edge_header and "chain_pitch_deg" in edge_header
    assert "diagonal_min" in header
    assert "shallowest quad corner:" in summary
    assert "health hits:" in summary
    write_vtk(taper, tmp_path / "shell.vtk")
    assert (tmp_path / "shell.legend.txt").stat().st_size > 0
    assert "health_code" in (tmp_path / "shell.legend.txt").read_text(encoding="utf-8")
    edges = (tmp_path / "shell.edges.vtk").read_text(encoding="utf-8")
    assert f"CELL_TYPES {len(taper.edges)}" in edges
    assert "SCALARS meet int 1" in edges
    assert "SCALARS health_code int 1" in edges
    face = (tmp_path / "shell.vtk").read_text(encoding="utf-8")
    assert f"CELL_TYPES {len(taper.elements)}" in face
    assert "SCALARS health_code int 1" in face
    corner_vtk = (tmp_path / "shell.corners.vtk").read_text(encoding="utf-8")
    assert f"CELL_TYPES {len(taper.corners)}" in corner_vtk
    assert "SCALARS angle_deg double 1" in corner_vtk
    assert "SCALARS element_id int 1" in corner_vtk


def test_a_block_takes_the_finer_station_count():
    """Both stations of a region take the larger edge count, whatever the bay length."""
    from implicit_geometry_kernel.blade.topology.region_counts import region_counts
    from implicit_geometry_kernel.blade.topology.template import LoopBranch, ShellTopologyTemplate, WebBranch

    def station(z):
        loop = (LoopBranch("panel", ("TE",), ("anchor", "le"), 3, "suction"),)
        web = (WebBranch("web_0", 0, ("anchor", "w", "suction"), ("anchor", "w", "pressure"), 3),)
        return ShellTopologyTemplate(z=z, signature=(), loop=loop, webs=web, vertices=())

    def samples(z, panel):
        def line(length):
            return np.column_stack([np.linspace(0.0, length, 2), np.zeros(2), np.full(2, z)])
        return {"panel": line(panel), "web_0": line(0.4)}

    def panel_at(zs, chords):
        log = []
        region_counts(
            [station(zs[0]), station(zs[1])],
            [samples(zs[0], 1.0), samples(zs[1], 0.5)],
            0.2, zs=zs, chords=chords, log=log,
        )
        return [row for row in log if row.region == "panel"]

    mild = panel_at((0.0, 0.5), (1.0, 0.95))
    assert mild[0].reason == "length"
    assert mild[1].reason == "block" and mild[1].count == mild[0].count
    short = panel_at((0.0, 0.02), (1.0, 0.95))
    assert short[1].reason == "block" and short[1].count == short[0].count
    jumped = panel_at((0.0, 0.35), (1.0, 0.05))
    assert jumped[1].reason == "block" and jumped[1].count == jumped[0].count


def test_an_oversized_edge_takes_the_block_count():
    """The station that needs more edges sets the count for the whole block."""
    from implicit_geometry_kernel.blade.topology.region_counts import region_counts
    from implicit_geometry_kernel.blade.topology.template import LoopBranch, ShellTopologyTemplate, WebBranch

    def station(z):
        loop = (
            LoopBranch("panel", ("TE",), ("anchor", "le"), 3, "suction"),
            LoopBranch("other", ("TE",), ("anchor", "le"), 3, "pressure"),
        )
        web = (WebBranch("web_0", 0, ("anchor", "w", "suction"), ("anchor", "w", "pressure"), 3),)
        return ShellTopologyTemplate(z=z, signature=(), loop=loop, webs=web, vertices=())

    def line(length, z):
        return np.column_stack([np.linspace(0.0, length, 2), np.zeros(2), np.full(2, z)])

    log = []
    region_counts(
        [station(0.0), station(0.05)],
        [
            {"panel": line(0.4, 0.0), "other": line(0.4, 0.0), "web_0": line(0.4, 0.0)},
            {"panel": line(1.2, 0.05), "other": line(0.4, 0.05), "web_0": line(0.4, 0.05)},
        ],
        0.2, zs=(0.0, 0.05), chords=(1.0, 0.95), log=log,
    )
    panel = [row for row in log if row.region == "panel"]
    other = [row for row in log if row.region == "other"]
    assert panel[0].reason == "block" and panel[1].reason == "length"
    assert panel[0].count == panel[1].count
    assert other[1].count == other[0].count


def test_a_closing_run_keeps_the_fan_budget():
    """A short closing arc limits the shared count to the edges its fan can carry."""
    from implicit_geometry_kernel.blade.topology.region_counts import _block_counts

    asked = [{"panel": 6, "web_0": 4}, {"panel": 6, "web_0": 4}, {"web_0": 4}]
    lengths = [{"panel": 1.0, "web_0": 0.4}, {"panel": 0.05, "web_0": 0.4}, {"web_0": 0.4}]
    out = _block_counts(asked, lengths, (0.0, 1.0, 1.07), {"panel"})
    assert [row["panel"] for row in out[:2]] == [2, 2]
    assert [row["web_0"] for row in out] == [4, 4, 4]


def test_the_size_ceiling_wins_over_the_fan_budget():
    """The fan is refused when its count would put a station's edges over ``h``."""
    from implicit_geometry_kernel.blade.topology.region_counts import region_counts
    from implicit_geometry_kernel.blade.topology.template import LoopBranch, ShellTopologyTemplate, WebBranch

    def station(z, names):
        loop = tuple(LoopBranch(name, ("TE",), ("anchor", "le"), 3, "suction") for name in names)
        web = (WebBranch("web_0", 0, ("anchor", "w", "suction"), ("anchor", "w", "pressure"), 3),)
        return ShellTopologyTemplate(z=z, signature=(), loop=loop, webs=web, vertices=())

    def line(length, z):
        return np.column_stack([np.linspace(0.0, length, 2), np.zeros(2), np.full(2, z)])

    log = []
    region_counts(
        [station(0.0, ("panel",)), station(1.0, ("panel",)), station(1.07, ())],
        [
            {"panel": line(1.0, 0.0), "web_0": line(0.4, 0.0)},
            {"panel": line(0.05, 1.0), "web_0": line(0.4, 1.0)},
            {"web_0": line(0.4, 1.07)},
        ],
        0.2, zs=(0.0, 1.0, 1.07), chords=(1.0, 0.95, 0.9), log=log,
    )
    panel = [row for row in log if row.region == "panel"]
    # A fan of 2 would give 0.5 m edges at the root station, over h = 0.2 m.
    assert [row.count for row in panel] == [6, 6]
    assert all(row.edge <= 0.2 + 1e-12 for row in panel)
    assert panel[1].reason == "block"


def _scaled(z1: float, chord: float, h: float):
    """Two stations, explicit, so the bay is the whole span and the mesh stays small."""
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0012"), pitch_axis=0.3),
        Station(z1, chord, Airfoil.naca4("0012"), pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )
    return shell_mesh(geo, h=h, zs=(0.0, z1))


def test_a_steep_shell_bay_stays_ruled():
    mesh = _scaled(0.08, 0.6, 0.4)
    shared = [s for s in mesh.strips if s.n0 and s.n1 and s.n0 == s.n1]
    assert shared
    for strip in shared:
        assert strip.constructor == "ruled"
        assert strip.site0 < 0 and strip.site1 < 0
    assert not any(s.constructor == "step_patch" for s in mesh.strips)


def test_a_thin_trailing_face_is_the_ridge():
    """A face under ``2 * h * tan(10°)`` is the shared ridge. The corners have not met."""
    from implicit_geometry_kernel.blade.topology.map_walk import (
        _skin_mid_loop,
        te_blunt_eps,
        te_face_open,
        te_opening,
        te_wrap_for_mesh,
    )

    foil = Airfoil.naca4("0012").with_te_thickness(0.004)
    stations = [
        Station(0.0, 2.0, foil, pitch_axis=0.3),
        Station(0.4, 0.6, foil, pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )
    root = geo.section_at(0.0)
    root_pts = _skin_mid_loop(root.shell_graph(), 0.0)
    opening = te_opening(root, root_pts)
    assert opening < te_blunt_eps(root)
    assert te_face_open(root, root_pts)
    assert not te_wrap_for_mesh(root, root_pts, 0.5)
    mesh = shell_mesh(geo, h=0.5, zs=(0.0, 0.4))
    assert [s for s in mesh.strips if s.region == "te"] == []


def test_a_closing_fan_stays_on_the_triangle():
    """The apex fan does not send its samples to a section that no longer has the face."""
    from implicit_geometry_kernel.blade.topology.map_loft import _transition_quads

    pa = np.array([[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.4, 0.0, 0.0]])
    pb = np.array([[0.2, 0.0, 1.0]])

    def snap(_name, _xyz, _at_end):
        return np.array([9.0, 9.0, 9.0])

    faces = _transition_quads("te", pa, pb, on_section=snap)
    cloud = np.vstack([face.face for face in faces])
    assert len(faces) > 0
    assert float(np.max(np.abs(cloud[:, 1]))) < 1e-9
    assert float(np.max(cloud[:, 0])) <= 0.4 + 1e-9
    assert not np.any(np.all(np.abs(cloud - 9.0) < 1e-9, axis=1))


def test_a_closing_flatback_leaves_a_legal_trailing_face():
    """A flatback that closes keeps a legal trailing face, and the fan uses the neighbour seam.

    The last open station sits just above ``2 * h * tan(10°)`` and the next
    station is sharp, so the close is a fan of one bay. Its spanwise corners
    are nodes of the neighbouring strips.
    """
    from implicit_geometry_kernel.mesh.tables import trailing_face_failures

    foil = Airfoil.naca4("0012")
    h = 0.45
    stations = [
        Station(0.0, 0.8, foil.with_te_thickness(0.55), pitch_axis=0.3),
        Station(0.5, 6.0, foil.with_te_thickness(0.08), pitch_axis=0.3),
        Station(0.95, 6.0, foil.with_te_thickness(0.030), pitch_axis=0.3),
        Station(1.40, 6.0, foil, pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )
    mesh = shell_mesh(geo, h=h, zs=tuple(s.z for s in stations))
    te = [s for s in mesh.strips if s.region == "te"]
    assert any(s.constructor == "fan" for s in te)
    assert trailing_face_failures(mesh) == []
    assert mesh.junctions == ()
    assert all(lp.is_root for lp in mesh.loops)
    fan_nodes = {n for i, e in enumerate(mesh.elements) if mesh.constructor[i] == "fan" for n in e}
    other = {n for i, e in enumerate(mesh.elements) if mesh.constructor[i] != "fan" for n in e}
    shared = fan_nodes & other
    assert shared
    assert all(mesh.node_origin[n] != "centroid" for n in shared)


def test_a_long_step_shares_the_taper():
    """A long flat region takes one block count and the bay stays ruled.

    The root edge would be several times the bay at the coarse count. The block
    uses the finer station, so both ends match and the strip is ruled quads.
    """
    from implicit_geometry_kernel.mesh.tables import transition_band_failures

    foil = Airfoil.naca4("0018").with_te_thickness(0.06)
    stations = [
        Station(0.0, 2.0, foil, pitch_axis=0.3),
        Station(0.05, 2.16, foil, pitch_axis=0.3),
    ]
    geo = BladeGeometry(
        stations, "2D-F-LC-4R", 0.004, pivot_fraction=0.3,
        spar_cap_thickness=0.008, shear_web_thickness=0.004,
    )
    mesh = shell_mesh(geo, h=0.20, zs=(0.0, 0.05))
    webs = [c for c in mesh.counts if c.region == "web_1"]
    assert len(webs) == 2 and webs[0].count == webs[1].count
    strips = [s for s in mesh.strips if s.region == "web_1"]
    assert len(strips) == 1 and strips[0].constructor == "ruled"
    assert strips[0].n0 == strips[0].n1
    assert not any(s.constructor == "step_patch" for s in mesh.strips)
    assert transition_band_failures(mesh) == []


def test_a_mild_shell_bay_is_ruled():
    mesh = _scaled(0.4, 0.6, 0.5)
    shared = [s for s in mesh.strips if s.n0 and s.n1]
    assert shared
    assert not any(s.constructor == "step_patch" for s in shared)
    for strip in shared:
        assert strip.n0 == strip.n1
        assert strip.constructor == "ruled"
    assert all(float(mesh.element_quality["angle_min"][i]) >= 20.0 - 1e-6 for i in mesh.quads)


def test_corner_rows_name_the_node_they_sit_on(taper):
    nodes = np.asarray(taper.nodes, float)
    by: dict[int, list] = {}
    for corner in taper.corners:
        element_id = int(corner["element_id"])
        element = taper.elements[element_id]
        assert element[int(corner["corner"])] == int(corner["node_id"])
        assert np.allclose((corner["x"], corner["y"], corner["z"]), nodes[int(corner["node_id"])])
        by.setdefault(element_id, []).append(corner)
    for i, element in enumerate(taper.elements):
        rows = by[i]
        assert len(rows) == len(element)
        assert min(float(row["angle_deg"]) for row in rows) == pytest.approx(
            float(taper.element_quality["angle_min"][i]), abs=1e-6)
        if len(element) == 4 and float(taper.element_quality["angle_min"][i]) < 20.0:
            assert len([row for row in rows if float(row["angle_deg"]) < 20.0]) == 1


def test_three_rotations_separate_a_flat_corner_from_a_fold():
    """Yaw is the in-plane turn, pitch leaves the element plane, roll is the fold about the edge."""
    from types import SimpleNamespace

    from implicit_geometry_kernel.mesh.tables import build_corners, build_edges

    flat = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
    ])
    sheet = SimpleNamespace(
        nodes=flat, elements=((0, 1, 2, 3),), surface=("skin",), name=("panel",),
        constructor=("ruled",), bay=(0,), node_origin=(),
    )
    corner = next(c for c in build_corners(sheet) if int(c["node_id"]) == 0)
    assert float(corner["angle_deg"]) == pytest.approx(90.0, abs=1e-6)
    assert float(corner["yaw_deg"]) == pytest.approx(90.0, abs=1e-6)
    assert float(corner["pitch_deg"]) == pytest.approx(0.0, abs=1e-6)
    assert float(corner["roll_deg"]) == pytest.approx(0.0, abs=1e-6)
    assert np.allclose((corner["dx_in"], corner["dy_in"], corner["dz_in"]), (-1.0, 0.0, 0.0))
    assert np.allclose((corner["dx_out"], corner["dy_out"], corner["dz_out"]), (0.0, 1.0, 0.0))

    bent = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.86602540378, 0.0, -0.5],
        [1.86602540378, 1.0, -0.5],
    ])
    pair = SimpleNamespace(
        nodes=bent, elements=((0, 1, 2, 3), (2, 1, 4, 5)), surface=("skin", "skin"),
        name=("te", "te"), constructor=("ruled", "ruled"), bay=(0, 0), node_origin=(),
    )
    edges = {(int(e["a"]), int(e["b"])): e for e in build_edges(pair)}
    hinge = edges[(1, 2)]
    assert float(hinge["dihedral_deg"]) == pytest.approx(30.0, abs=0.05)
    assert all(float(e["dihedral_deg"]) == pytest.approx(0.0, abs=1e-6) for key, e in edges.items() if key != (1, 2))
    leaving = next(c for c in build_corners(pair) if int(c["element_id"]) == 0 and int(c["node_id"]) == 1)
    assert float(leaving["yaw_deg"]) == pytest.approx(90.0, abs=0.05)
    assert float(leaving["roll_deg"]) == pytest.approx(30.0, abs=0.05)
    assert abs(float(leaving["pitch_deg"])) < 1.0

    warped = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]])
    twist = SimpleNamespace(
        nodes=warped, elements=((0, 1, 2, 3),), surface=("skin",), name=("panel",),
        constructor=("ruled",), bay=(0,), node_origin=(),
    )
    lifted = next(c for c in build_corners(twist) if int(c["node_id"]) == 1)
    assert float(lifted["yaw_deg"]) == pytest.approx(90.0, abs=0.05)
    assert float(lifted["pitch_deg"]) == pytest.approx(np.degrees(np.arctan(0.5)), abs=0.05)
    assert float(lifted["roll_deg"]) == pytest.approx(np.degrees(np.arctan(0.5)), abs=0.05)

    joined = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0], [2.0, 1.0, 0.0],
        [1.0, 0.0, 1.0], [1.0, 1.0, 1.0],
    ])
    triple = SimpleNamespace(
        nodes=joined,
        elements=((0, 1, 2, 3), (1, 4, 5, 2), (1, 2, 7, 6)),
        surface=("skin", "skin", "web"),
        name=("a", "b", "web"),
        constructor=("ruled", "ruled", "ruled"),
        bay=(0, 0, 0),
        node_origin=(),
    )
    shared = next(e for e in build_edges(triple) if {int(e["a"]), int(e["b"])} == {1, 2})
    assert int(shared["use"]) == 3
    assert float(shared["dihedral_deg"]) == pytest.approx(0.0, abs=1e-6)
    skin_roll = next(c for c in build_corners(triple) if int(c["element_id"]) == 0 and int(c["node_id"]) == 1)
    web_roll = next(c for c in build_corners(triple) if int(c["element_id"]) == 2 and int(c["node_id"]) == 1)
    assert float(skin_roll["roll_deg"]) == pytest.approx(0.0, abs=1e-6)
    assert float(web_roll["roll_deg"]) == pytest.approx(-90.0, abs=0.05)


def test_a_forced_shallow_corner_points_at_its_node():
    from types import SimpleNamespace

    from implicit_geometry_kernel.mesh.quality import element_quality
    from implicit_geometry_kernel.mesh.tables import build_corners

    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.99, 0.001, 0.0],
        [0.0, 1.0, 0.0],
    ])
    elements = ((0, 1, 2, 3),)
    quality = element_quality(nodes, elements, [0.004])
    corners = build_corners(SimpleNamespace(
        nodes=nodes, elements=elements, name=("te",), surface=("skin",),
        constructor=("step_patch",), bay=(4,),
    ))
    sharp = [c for c in corners if float(c["angle_deg"]) < 20.0]
    assert len(sharp) == 1
    corner = sharp[0]
    assert int(corner["node_id"]) == 1
    assert int(corner["corner"]) == 1
    assert np.allclose((corner["x"], corner["y"], corner["z"]), nodes[int(corner["node_id"])])
    assert float(quality["angle_min"][0]) == pytest.approx(float(corner["angle_deg"]))
    groups = mesh_failures(SimpleNamespace(
        nodes=nodes, elements=elements, element_quality=quality,
        section_distance=np.zeros(len(nodes)), edges=(), constructor=("step_patch",),
        bay=(4,), z0=(0.0,), z1=(0.05,), corners=corners,
    ))
    assert len(groups) == 1
    assert groups[0]["corner"] == int(corner["corner"])
    assert groups[0]["constructor"] == "step_patch"
