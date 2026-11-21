"""Authored wall-region graph and mapped loft: named regions keep their identity along the span."""

from implicit_geometry_kernel.blade.topology.fields import CutWallMetaSDF, StackThickness, span_thickness
from implicit_geometry_kernel.blade.topology.graph import (
    Anchor,
    Region,
    WallRegionGraph,
    WebMember,
)
from implicit_geometry_kernel.blade.topology.signatures import (
    role_atom,
    signatures_compatible,
)
from implicit_geometry_kernel.blade.topology.template import (
    LoopBranch,
    ShellTopologyTemplate,
    WebBranch,
    ZonedShellTopology,
    zone_region_graphs,
)
from implicit_geometry_kernel.blade.topology.map import (
    place_station_keys,
    sample_region_polylines,
    station_key_points,
)
from implicit_geometry_kernel.blade.topology.sizing import (
    resample_open_h,
    size_field_h,
    station_live_thickness,
)
from implicit_geometry_kernel.blade.topology.map_loft import loft_mapped_strips
from implicit_geometry_kernel.blade.topology.map_walk import (
    hold_te_blunt,
    te_blunt_eps,
    te_blunt_flags_for_cuts,
    te_blunt_of_cut,
    te_is_blunt,
    te_opening,
    te_ridge_point,
)
from implicit_geometry_kernel.blade.topology.ribs import (
    ear_clip_ring,
    rib_cell_polygons_from_samples,
    rib_cell_triangles_from_samples,
    rib_face_polygons,
    rib_face_triangles,
    section_cut_near,
    templates_along,
    tip_closure_polygons,
)

from implicit_geometry_kernel.blade.topology.network import (
    StationNetwork,
    WallBranch,
    station_network,
    station_network_from_meta,
    wall_meta_from_cut,
)

__all__ = [
    "Anchor",
    "CutWallMetaSDF",
    "StationNetwork",
    "StackThickness",
    "WallBranch",
    "LoopBranch",
    "Region",
    "ShellTopologyTemplate",
    "WallRegionGraph",
    "WebBranch",
    "WebMember",
    "ZonedShellTopology",
    "hold_te_blunt",
    "loft_mapped_strips",
    "ear_clip_ring",
    "rib_cell_polygons_from_samples",
    "rib_cell_triangles_from_samples",
    "rib_face_polygons",
    "rib_face_triangles",
    "tip_closure_polygons",
    "place_station_keys",
    "station_key_points",
    "resample_open_h",
    "role_atom",
    "sample_region_polylines",
    "size_field_h",
    "span_thickness",
    "station_live_thickness",
    "station_network",
    "station_network_from_meta",
    "section_cut_near",
    "signatures_compatible",
    "wall_meta_from_cut",
    "te_blunt_eps",
    "te_blunt_flags_for_cuts",
    "te_blunt_of_cut",
    "te_is_blunt",
    "te_opening",
    "te_ridge_point",
    "templates_along",
    "zone_region_graphs",
]
