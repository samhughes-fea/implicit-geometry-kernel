"""implicit_geometry_kernel — geometry of blade and wing cross-sections.

The kernel describes the geometry of blade and wing sections and lofts it along the
span. It holds nothing about materials, plies, stiffness, loads or solvers; the
laminate and section models build those on top of it.

Public API
----------
Airfoils   Airfoil (from coordinates, ``Airfoil.naca4``, ``Airfoil.ellipse``)
Catalogue  SYSTEM_TYPE_KEYS and resolve_system_type (WX-Y-Z-NR tokens)
Concept    BladeConcept (catalogue or parametric layout plus spanwise extras)
Station    Station (span position, chord, twist, pitch axis, offset),
           SkinRegion, ShearWeb
Loft       BladeGeometry (layout lofted between stations; optional span tables),
           interpolate(z), section_at(z) -> SectionCut, oml_polyline, oml_sdf.
           Scalars interpolate linearly; airfoils stay 2-D linear blends.
Section    implicit_geometry_kernel.section — one-station SectionSpec / SectionGeometry
Feature    FieldSample.meta ids, GeometrySolver, wall_meta / wall_mid_sdf,
           shell_graph / partition_edges / station_network,
           shell_loft_zones / LockedLoftPlan / resample_locked,
           WallRegionGraph / loft_mapped_strips.

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

Dependencies: numpy and numba; system/conda GEOS (``libgeos_c``) for polygon
offset and overlay; scipy (``[splines]``) for spline-resampled airfoils.
Matplotlib (``[plots]``) is ``implicit_geometry_kernel.plotting``, not imported
by geometry modules. Offset and overlay call GEOS through
``implicit_geometry_kernel.geos``. Public products are numpy.
"""

from implicit_geometry_kernel.blade import (
    Airfoil,
    Anchor,
    BladeConcept,
    BladeGeometry,
    CellPolygon,
    MidSurface,
    SectionCut,
    SectionWall,
    GeometricIntegrals,
    geometric_integrals,
    ShearWeb,
    LockedLoftPlan,
    SectionSample,
    PartitionEdge,
    Region,
    ShellGraph,
    ShellRun,
    ShellTopologyTemplate,
    TopologyDiff,
    TopologyEvent,
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
    partition_edge_curve,
    resample_locked,
    rib_face_polygons,
    section_at,
    shared_loft_names,
    signatures_compatible,
    shell_graph,
    shell_loft_zones,
    shell_strips_from_graph,
    topology_diff,
    transverse_members_from_windio,
    zone_region_graphs,
    station_network,
    station_network_from_meta,
    tip_closure_polygons,
    wall_meta_from_cut,
    refine_section_samples,
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

__version__ = "0.9.2"

__all__ = [
    "Airfoil", "Anchor", "BladeConcept", "BladeGeometry", "CellPolygon", "FeatureTable", "FieldSample", "GeometrySolver",
    "LoftedPolygonSDF", "LockedLoftPlan", "MidSurface", "NO_FEATURE", "PriorityUnionNode", "SYSTEM_TYPE_KEYS",
    "SYSTEM_TYPE_LAYOUT_KEYS",
    "GeometricIntegrals", "geometric_integrals",
    "PartitionEdge", "Region", "SectionCut", "SectionWall", "ShearWeb", "ShellGraph", "ShellRun",
    "ShellTopologyTemplate", "SkinRegion",
    "Station", "StationNetwork", "Stringer", "TagNode", "TopologyDiff", "TopologyEvent",
    "SpanwiseScheme", "TipClosure", "TransverseMember", "WallRegionGraph", "WebMember", "PlyDrop",
    "loft_mapped_strips", "partition_edge_curve", "read_planform_csv", "refine_section_samples", "write_windio", "resolve_system_type",
    "SectionSample",
    "resample_locked", "rib_face_polygons",
    "section_at",
    "shared_loft_names", "signatures_compatible", "shell_graph", "shell_loft_zones",
    "shell_strips_from_graph",     "station_network", "station_network_from_meta", "tip_closure_polygons", "topology_diff",
    "transverse_members_from_windio", "wall_meta_from_cut", "zone_region_graphs",
]
