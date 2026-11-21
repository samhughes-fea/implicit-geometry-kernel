"""Mesh tables: edges, corners, seams, failure adapters, CSV dump."""

from implicit_geometry_kernel.mesh.tables.tables import (
    CONSTRUCTORS, CORNER_DTYPE, EDGE_DTYPE, KINDS, MEETS, ORIGINS,
    Corner, MeshEdge, OpenLoop, Seam, StationDecision, StripDecision,
    assemble, build_corners, build_edges, mesh_failures, near_misses,
    off_root_open_edges, t_junctions, trailing_face_failures,
    transition_band_failures, valence_per_node, write_mesh_tables,
)

__all__ = [
    "CONSTRUCTORS", "CORNER_DTYPE", "EDGE_DTYPE", "KINDS", "MEETS", "ORIGINS",
    "Corner", "MeshEdge", "OpenLoop", "Seam", "StationDecision", "StripDecision",
    "assemble", "build_corners", "build_edges", "mesh_failures", "near_misses",
    "off_root_open_edges", "t_junctions", "trailing_face_failures",
    "transition_band_failures", "valence_per_node", "write_mesh_tables",
]


def __getattr__(name: str):
    """Expose tables helpers (including private) for compatibility imports."""
    import implicit_geometry_kernel.mesh.tables.tables as tables
    try:
        return getattr(tables, name)
    except AttributeError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
