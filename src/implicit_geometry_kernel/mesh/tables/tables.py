"""Mesh tables on a shell mesh: edges, open loops, the failure query, and the CSV dump.

Quality stays in :data:`QUALITY_DTYPE`. These tables are the mesh itself: node ids,
coordinates, connectivity, and the decisions that built each face.
"""
from __future__ import annotations


import numpy as np

from implicit_geometry_kernel.mesh.tables.geometry import (
    build_corners,
    build_edges,
    diagonal_mins,
)
from implicit_geometry_kernel.mesh.tables.checks import (
    _chain_open,
    _ruled_valence_defects,
    near_misses,
    t_junctions,
    valence_per_node,
)
from implicit_geometry_kernel.mesh.tables.records import (  # noqa: F401 (re-exported helpers)
    CONSTRUCTORS,
    CORNER_DTYPE,
    Corner,
    EDGE_DTYPE,
    KINDS,
    MEETS,
    MeshEdge,
    ORIGINS,
    OpenLoop,
    Seam,
    StationDecision,
    StripDecision,
    _edge_users,
)
from implicit_geometry_kernel.mesh.tables.geometry import (  # noqa: F401 (re-exported helpers)
    _by_size,
    _chain_turns,
    _corner_angles_block,
    _diagonal_min,
    _dihedral,
    _dihedral_block,
    _dihedral_pair,
    _element_normals,
    _face_normal,
    _kind,
    _meet,
    _newell_units,
    _other_face,
    _trajectory_angles,
    _trajectory_block,
    _unit,
    _unit_rows,
)
from implicit_geometry_kernel.mesh.tables.checks import (  # noqa: F401 (re-exported helpers)
    _cell_groups,
    _segment_hits,
    mesh_failures,
    trailing_face_failures,
    transition_band_failures,
)
from implicit_geometry_kernel.mesh.tables.dump import (  # noqa: F401 (re-exported helpers)
    _csv,
    _csv_columns,
    _shallowest_quad_corner,
    _summary,
    off_root_open_edges,
    write_mesh_tables,
)


"""One undirected edge. ``e0``..``e2`` are the element ids, ``-1`` when unused."""

"""One element corner. Columns match :class:`Corner`."""


def assemble(mesh, z_root: float) -> dict:
    """Edges, loops, corners, valence, near-misses, and T-junctions for a mesh whose elements are already set."""
    edges = build_edges(mesh)
    nodes = np.asarray(mesh.nodes, float)
    valence = valence_per_node(mesh, edges)
    return {
        "edges": edges,
        "loops": _chain_open(edges, nodes, float(z_root)),
        "corners": build_corners(mesh, edges),
        "valence": valence,
        "near_miss": near_misses(nodes),
        "junctions": t_junctions(nodes, edges),
        "valence_defects": _ruled_valence_defects(mesh, edges, valence),
        "diagonal_min": diagonal_mins(mesh),
    }


__all__ = [
    "CONSTRUCTORS", "CORNER_DTYPE", "EDGE_DTYPE", "KINDS", "MEETS", "ORIGINS",
    "Corner", "MeshEdge", "OpenLoop", "Seam", "StationDecision", "StripDecision",
    "assemble", "build_corners", "build_edges", "mesh_failures", "near_misses",
    "off_root_open_edges", "t_junctions", "trailing_face_failures",
    "transition_band_failures", "valence_per_node", "write_mesh_tables",
]
