"""Conforming blade shell mesh (quads on skin/webs, triangles on ribs/tip)."""

from implicit_geometry_kernel.mesh.shell.shell_mesh import (
    TIP_H_OVER_C,
    MeshRecipe,
    ShellMesh,
    SURFACES,
    mesh_stations,
    orient_shell_elements,
    shell_mesh,
)

from implicit_geometry_kernel.mesh.shell.partition_mesh import PartitionMeshError, partition_mesh

__all__ = [
    "PartitionMeshError",
    "partition_mesh",
    "TIP_H_OVER_C",
    "MeshRecipe",
    "ShellMesh",
    "SURFACES",
    "mesh_stations",
    "orient_shell_elements",
    "shell_mesh",
]


def __getattr__(name: str):
    """Expose ``shell_mesh`` module helpers (including private) for tests."""
    import importlib

    mod = importlib.import_module("implicit_geometry_kernel.mesh.shell.shell_mesh")
    try:
        return getattr(mod, name)
    except AttributeError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
