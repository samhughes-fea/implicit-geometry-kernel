"""implicit_geometry_kernel — geometry of blade and wing cross-sections.

The kernel describes the geometry of blade and wing sections and lofts it along the
span. It holds nothing about materials, plies, stiffness, loads or solvers; the
laminate and section models build those on top of it.

Public API
----------

========== ==================================================================
Airfoils   Airfoil (from coordinates, ``Airfoil.naca4``, ``Airfoil.ellipse``)
Catalogue  SYSTEM_TYPE_KEYS and resolve_system_type (WX-Y-Z-NR tokens)
Concept    BladeConcept (catalogue or parametric layout plus spanwise extras)
Station    Station (span position, chord, twist, pitch axis, offset),
           SkinRegion, ShearWeb
Loft       BladeGeometry (layout lofted between stations; optional span
           tables), interpolate(z), section_at(z) -> SectionCut,
           oml_polyline, oml_sdf. Scalars interpolate linearly; airfoils
           stay 2-D linear blends.
Section    implicit_geometry_kernel.section: one-station SectionSpec /
           SectionGeometry
Feature    FieldSample.meta ids, GeometrySolver, wall_meta / wall_mid_sdf,
           shell_graph / ShellGraph, LockedLoftPlan / resample_locked,
           WallRegionGraph / loft_mapped_strips.
========== ==================================================================

Rare topology helpers (``TopologyDiff``, ``ShellTopologyTemplate``, …) live under
``implicit_geometry_kernel.blade`` / ``blade.topology``. Importing them from the
package root still works in 0.9.x with a ``DeprecationWarning``.

Subpackages
-----------
blade         airfoils, stations, lofting and cut sections
section       one 2-D station: layout, walls, midlines, junctions, cells, properties
sdf           signed-distance modelling, AirfoilSDF, GeometrySolver
frames        geometric coordinate frames
integration   quadtree and shoelace integrals
contract      ShellMidlineStrip, CellWallLoop (public DTOs)
extract       in-plane zero loops and feature curves (not a mesh)
catalogue     WX-Y-Z-NR tokens
mesh          conforming shell mesh (``shell_mesh``); Gmsh comparison stays on ``mesh.gmsh``

Dependencies: numpy and numba; system/conda GEOS (``libgeos_c``) for polygon
offset and overlay; scipy (``[splines]``) for spline-resampled airfoils.
Matplotlib (``[plots]``) is ``implicit_geometry_kernel.plotting``, not imported
by geometry modules. Offset and overlay call GEOS through
``implicit_geometry_kernel.geos``. Public products are numpy.
"""

from __future__ import annotations

import warnings

from implicit_geometry_kernel.blade import (
    Airfoil,
    Anchor,
    BladeConcept,
    BladeGeometry,
    GeometricIntegrals,
    geometric_integrals,
    ShearWeb,
    LockedLoftPlan,
    SectionSample,
    SectionCut,
    ShellGraph,
    SkinRegion,
    Station,
    StationNetwork,
    Stringer,
    SpanwiseScheme,
    TipClosure,
    TransverseMember,
    WallRegionGraph,
    WebMember,
    PlyDrop,
    loft_mapped_strips,
    read_planform_csv,
    write_windio,
    resample_locked,
    rib_face_polygons,
    section_at,
    shell_graph,
    shell_loft_zones,
    station_network,
    tip_closure_polygons,
    transverse_members_from_windio,
    wall_meta_from_cut,
)
from implicit_geometry_kernel.sdf import (
    FeatureTable,
    FieldSample,
    GeometrySolver,
    LoftedPolygonSDF,
    NO_FEATURE,
    PriorityUnionNode,
    TagNode,
)
from implicit_geometry_kernel.catalogue import (
    SYSTEM_TYPE_KEYS,
    SYSTEM_TYPE_LAYOUT_KEYS,
    resolve_system_type,
)

__version__ = "0.9.4"

# Topology helpers demoted from the barrel in 0.9.x (K5). Still importable from
# ``blade`` / ``blade.topology``; root import warns once per name this process.
_DEPRECATED_ROOT = {
    "CellPolygon": "implicit_geometry_kernel.blade",
    "MidSurface": "implicit_geometry_kernel.blade",
    "PartitionEdge": "implicit_geometry_kernel.blade",
    "Region": "implicit_geometry_kernel.blade",
    "SectionWall": "implicit_geometry_kernel.blade",
    "ShellRun": "implicit_geometry_kernel.blade",
    "ShellTopologyTemplate": "implicit_geometry_kernel.blade",
    "TopologyDiff": "implicit_geometry_kernel.blade",
    "TopologyEvent": "implicit_geometry_kernel.blade",
    "partition_edge_curve": "implicit_geometry_kernel.blade",
    "refine_section_samples": "implicit_geometry_kernel.blade",
    "shared_loft_names": "implicit_geometry_kernel.blade",
    "signatures_compatible": "implicit_geometry_kernel.blade",
    "shell_strips_from_graph": "implicit_geometry_kernel.blade",
    "station_network_from_meta": "implicit_geometry_kernel.blade",
    "topology_diff": "implicit_geometry_kernel.blade",
    "zone_region_graphs": "implicit_geometry_kernel.blade",
}

__all__ = [
    "Airfoil", "Anchor", "BladeConcept", "BladeGeometry", "FeatureTable", "FieldSample",
    "GeometrySolver", "LoftedPolygonSDF", "LockedLoftPlan", "NO_FEATURE", "PriorityUnionNode",
    "SYSTEM_TYPE_KEYS", "SYSTEM_TYPE_LAYOUT_KEYS",
    "GeometricIntegrals", "geometric_integrals",
    "SectionCut", "ShearWeb", "ShellGraph", "SkinRegion",
    "Station", "StationNetwork", "Stringer", "TagNode",
    "SpanwiseScheme", "TipClosure", "TransverseMember", "WallRegionGraph", "WebMember", "PlyDrop",
    "loft_mapped_strips", "read_planform_csv", "write_windio", "resolve_system_type",
    "SectionSample",
    "resample_locked", "rib_face_polygons",
    "section_at",
    "shell_graph", "shell_loft_zones",
    "station_network", "tip_closure_polygons",
    "transverse_members_from_windio", "wall_meta_from_cut",
]


def __getattr__(name: str):
    home = _DEPRECATED_ROOT.get(name)
    if home is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    warnings.warn(
        f"importing {name!r} from implicit_geometry_kernel is deprecated; "
        f"import it from {home} instead (removed from the top-level barrel in 0.9.x, K5)",
        DeprecationWarning,
        stacklevel=2,
    )
    import importlib

    module = importlib.import_module(home)
    value = getattr(module, name)
    globals()[name] = value
    return value
