"""Planar fills: polygon Voronoi and Steiner Delaunay."""

from implicit_geometry_kernel.mesh.fill.delaunay import triangulate_ring
from implicit_geometry_kernel.mesh.fill.polygon import (
    fill_rib_cells, mesh_rib_cells, mesh_rib_cells_many, point_in_ring,
    points_in_ring, polygon_mesh_ring,
)

__all__ = [
    "fill_rib_cells", "mesh_rib_cells", "mesh_rib_cells_many", "point_in_ring",
    "points_in_ring", "polygon_mesh_ring", "triangulate_ring",
]


def __getattr__(name: str):
    """Expose fill helpers (including private) for compatibility imports."""
    import importlib
    for module in (
        "implicit_geometry_kernel.mesh.fill.polygon",
        "implicit_geometry_kernel.mesh.fill.delaunay",
        "implicit_geometry_kernel.mesh.fill.polygon_voronoi",
    ):
        mod = importlib.import_module(module)
        if hasattr(mod, name):
            return getattr(mod, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
