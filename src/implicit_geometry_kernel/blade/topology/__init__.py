"""Authored wall-region graph and mapped loft (tidalshell identity, no mesher)."""

from implicit_geometry_kernel.blade.topology.graph import (
    Anchor,
    Region,
    WallRegionGraph,
    WebMember,
    _cluster_web_members,
    _cut_tc,
    _regularise_anchors,
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
    _chord_fracs,
    _key_name,
    _oml_te_xyz,
    _quad_bowtie,
    _te_wrap_arc,
    loft_mapped_strips,
    place_station_keys,
    sample_region_polylines,
)
from implicit_geometry_kernel.blade.topology.ribs import (
    rib_cell_polygons_from_samples,
    rib_face_polygons,
    rib_face_triangles,
    section_cut_near,
    templates_along,
)

# Private helpers used by BladeGeometry.write_mid_surfaces
from implicit_geometry_kernel.blade.topology.ribs import _template_from_cut  # noqa: F401

__all__ = [
    "Anchor",
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
    "loft_mapped_strips",
    "rib_cell_polygons_from_samples",
    "rib_face_polygons",
    "rib_face_triangles",
    "place_station_keys",
    "role_atom",
    "sample_region_polylines",
    "section_cut_near",
    "signatures_compatible",
    "templates_along",
    "zone_region_graphs",
]
