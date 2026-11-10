"""Authored wall-region graph and mapped loft (tidalshell identity, no mesher)."""

from implicit_geometry_kernel.blade.topology.fields import CutWallMetaSDF, StackThickness, span_thickness
from implicit_geometry_kernel.blade.topology.graph import (
    Anchor,
    Region,
    WallRegionGraph,
    WebMember,
    _cluster_web_members,
    _cut_tc,
    _regularise_anchors,
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
    _key_name,
    place_station_keys,
    sample_region_polylines,
)
from implicit_geometry_kernel.blade.topology.sizing import (
    resample_open_h,
    size_field_h,
    station_live_thickness,
)
from implicit_geometry_kernel.blade.topology.map_loft import (
    _quad_bowtie,
    loft_mapped_strips,
)
from implicit_geometry_kernel.blade.topology.map_walk import (
    _chord_fracs,
    _oml_te_xyz,
    _te_wrap_arc,
    hold_te_blunt,
    te_blunt_eps,
    te_blunt_flags_for_cuts,
    te_blunt_of_cut,
    te_is_blunt,
    te_opening,
)
from implicit_geometry_kernel.blade.topology.ribs import (
    ear_clip_ring,
    rib_cell_polygons_from_samples,
    rib_cell_triangles_from_samples,
    rib_face_polygons,
    rib_face_triangles,
    section_cut_near,
    templates_along,
)

# Private helpers used by BladeGeometry.write_mid_surfaces
from implicit_geometry_kernel.blade.topology.ribs import _template_from_cut  # noqa: F401
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
    "_chord_fracs",
    "_cluster_web_members",
    "_cut_tc",
    "_key_name",
    "_oml_te_xyz",
    "_quad_bowtie",
    "_regularise_anchors",
    "_te_wrap_arc",
    "_template_from_cut",
    "hold_te_blunt",
    "loft_mapped_strips",
    "ear_clip_ring",
    "rib_cell_polygons_from_samples",
    "rib_cell_triangles_from_samples",
    "rib_face_polygons",
    "rib_face_triangles",
    "place_station_keys",
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
    "templates_along",
    "zone_region_graphs",
]
