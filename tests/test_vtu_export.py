"""Binary XML VTK export of the shell mesh."""

import numpy as np
import pytest

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station
from implicit_geometry_kernel.mesh import shell_mesh, write_vtk, write_vtu
from implicit_geometry_kernel.mesh.export import read_vtu_arrays


@pytest.fixture(scope="module")
def mesh():
    foil = Airfoil.naca4("0018")
    stations = [Station(0.0, 1.0, foil, pitch_axis=0.3), Station(4.0, 0.7, foil, pitch_axis=0.3)]
    geo = BladeGeometry(stations, "2D-F-LC-1R", 0.004, spar_cap_thickness=0.008,
                        shear_web_thickness=0.004, pivot_fraction=0.3)
    return shell_mesh(geo, h=0.25)


def test_round_trip(mesh, tmp_path):
    path = tmp_path / "blade.vtu"
    partitions = write_vtu(mesh, path)
    data = read_vtu_arrays(path)
    np.testing.assert_array_equal(data["Points"], np.asarray(mesh.nodes, float))
    sizes = np.array([len(e) for e in mesh.elements])
    np.testing.assert_array_equal(data["offsets"], np.cumsum(sizes))
    np.testing.assert_array_equal(data["connectivity"], [i for e in mesh.elements for i in e])
    assert set(data["types"].tolist()) <= {5, 9, 7}
    np.testing.assert_array_equal(data["thickness"], mesh.thickness)
    np.testing.assert_array_equal(data["eccentricity"], mesh.eccentricity)
    assert [partitions[k] for k in data["partition"]] == list(mesh.partition)
    meta = data["_meta"]
    assert meta["provenance"]["digest"] == mesh.digest()
    assert meta["partition"] == list(partitions)


def test_smaller_than_ascii(mesh, tmp_path):
    write_vtu(mesh, tmp_path / "a.vtu")
    write_vtk(mesh, tmp_path / "a.vtk")
    assert (tmp_path / "a.vtu").stat().st_size < (tmp_path / "a.vtk").stat().st_size


def test_reads_in_vtk_when_available(mesh, tmp_path):
    vtk = pytest.importorskip("vtk")
    path = tmp_path / "b.vtu"
    write_vtu(mesh, path)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    assert grid.GetNumberOfPoints() == len(mesh.nodes)
    assert grid.GetNumberOfCells() == len(mesh.elements)


def test_benchmark_script_runs(tmp_path):
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "bench_shell_mesh", Path(__file__).resolve().parents[1] / "scripts" / "bench_shell_mesh.py")
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    rows = bench.run([0.4])
    assert {r["mesher"] for r in rows} == {"ruled", "partition"}
    assert all(r["elements"] > 0 and len(r["digest"]) == 64 for r in rows)
