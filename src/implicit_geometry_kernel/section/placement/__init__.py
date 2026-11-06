"""Former path of the through-thickness anchoring modules in ``section.layout``; kept as an alias until 0.6.0."""

from implicit_geometry_kernel.section.layout.anchoring import (
    AnchoredWallGeometry,
    airfoil_inward_depth_for_midsurface,
    anchored_wall_geometry,
    anchored_wall_geometry_from_policy,
    default_construction_for_kind,
    midsurface_offset_from_reference_m,
    physical_face_offsets_from_reference_m,
)
from implicit_geometry_kernel.section.layout.collision import (
    PlacementGeometryError,
    assert_anchored_wall_feasible,
    find_placement_collisions,
)
from implicit_geometry_kernel.section.layout.handoff import (
    placement_stamp_for_kind,
    require_placement_stamp,
    stamp_strip_placement,
    stamp_strip_placement_fields,
)

__all__ = [
    "AnchoredWallGeometry",
    "PlacementGeometryError",
    "airfoil_inward_depth_for_midsurface",
    "anchored_wall_geometry",
    "anchored_wall_geometry_from_policy",
    "assert_anchored_wall_feasible",
    "default_construction_for_kind",
    "find_placement_collisions",
    "midsurface_offset_from_reference_m",
    "physical_face_offsets_from_reference_m",
    "placement_stamp_for_kind",
    "require_placement_stamp",
    "stamp_strip_placement",
    "stamp_strip_placement_fields",
]
