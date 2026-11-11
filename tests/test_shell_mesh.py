"""The kernel's shell mesh: all-quad skin and webs, polygon ribs and tip, one conforming sheet."""

from collections import Counter

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.mesh.shell_mesh import mesh_stations, shell_mesh


def _geo(key="2D-F-LC-4R"):
    stations = [
        Station(0.0, 2.0, Airfoil.naca4("0024"), pitch_axis=0.3),
        Station(10.0, 1.0, Airfoil.naca4("0018"), pitch_axis=0.3),
    ]
    return BladeGeometry(stations, key, 0.004, pivot_fraction=0.3,
                         spar_cap_thickness=0.008, shear_web_thickness=0.004)


@pytest.fixture(scope="module", params=[{"n_span": 11}, {"h": 0.08}], ids=["n_span", "h"])
def mesh(request):
    return _geo(), shell_mesh(_geo(), **request.param)


def test_elements_are_quads_and_polygons(mesh):
    _geo_, m = mesh
    for e, surface in zip(m.elements, m.surface):
        assert len(set(e)) == len(e)
        if surface in ("skin", "web"):
            assert len(e) == 4
        else:
            assert len(e) >= 5
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
    assert all(len(e) >= 5 for e, s in zip(m.elements, m.surface) if s in ("rib", "tip_closure"))
    z_root = float(geo.stations[0].z)
    web = {k for e, s in zip(m.elements, m.surface) if s == "web" for k in e}
    for (a, b), n in m.edge_use().items():
        if n != 1:
            continue
        # open only along the root, or along a web's start or end line
        on_root = abs(m.nodes[a, 2] - z_root) < 1e-9 and abs(m.nodes[b, 2] - z_root) < 1e-9
        assert on_root or (a in web and b in web)


def test_vtk_keeps_polygons(tmp_path, mesh):
    from implicit_geometry_kernel.mesh import write_vtk

    _geo_, m = mesh
    path = tmp_path / "shell.vtk"
    parts = write_vtk(m, path)
    text = path.read_text().splitlines()
    i = text.index(f"CELL_TYPES {len(m.elements)}")
    types = [int(t) for t in text[i + 1:i + 1 + len(m.elements)]]
    assert types.count(9) == len(m.quads) and types.count(7) == len(m.polygons)
    assert set(parts) == set(m.partition)
    assert f"POINTS {len(m.nodes)} double" in text


def test_gmsh_path_is_deprecated():
    pytest.importorskip("gmsh")
    from implicit_geometry_kernel.mesh.gmsh import triangulate_ring

    with pytest.deprecated_call():
        triangulate_ring(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]))
