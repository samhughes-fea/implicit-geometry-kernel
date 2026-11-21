"""The 0.9.2 flat ``mesh.*`` module paths still import, with a warning, until 0.9.5."""

from __future__ import annotations

import importlib

import pytest

import implicit_geometry_kernel.mesh as mesh

FLAT = {
    "shell_mesh": ("ShellMesh", "implicit_geometry_kernel.mesh.shell.shell_mesh"),
    "vtk": ("write_vtk", "implicit_geometry_kernel.mesh.export.vtk"),
    "polygon": ("polygon_mesh_ring", "implicit_geometry_kernel.mesh.fill.polygon"),
    "polygon_voronoi": ("fill_rib_cells", "implicit_geometry_kernel.mesh.fill.polygon_voronoi"),
    "delaunay": ("triangulate_ring", "implicit_geometry_kernel.mesh.fill.delaunay"),
}


@pytest.mark.parametrize("old", sorted(FLAT))
def test_a_flat_path_warns_and_forwards(old):
    name, new = FLAT[old]
    legacy = importlib.import_module(f"implicit_geometry_kernel.mesh.{old}")
    with pytest.warns(DeprecationWarning, match="removed in 0.9.5"):
        value = getattr(legacy, name)
    assert value is getattr(importlib.import_module(new), name)


def test_importing_the_flat_module_keeps_the_shell_mesh_function():
    importlib.import_module("implicit_geometry_kernel.mesh.shell_mesh")
    assert callable(mesh.shell_mesh)
    assert mesh.shell_mesh is importlib.import_module("implicit_geometry_kernel.mesh.shell").shell_mesh


def test_an_unknown_name_is_an_attribute_error():
    legacy = importlib.import_module("implicit_geometry_kernel.mesh.vtk")
    with pytest.raises(AttributeError):
        legacy.not_a_name  # noqa: B018
