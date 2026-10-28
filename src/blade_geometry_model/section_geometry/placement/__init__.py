"""Anchored wall geometry: placement → physical faces and midsurface (WP3).

``section_layup`` owns *where* plies sit through the thickness
(``PlyPlacementPolicy``). This package owns the *geometry* realisation:
which physical face stays fixed, where the wall midsurface sits, and how far
each free face is from the reference.

Construction tokens are plain strings matching the layup enums so this package
never imports ``material_utilities`` (section_geometry G-row boundary).
"""

from __future__ import annotations

from blade_geometry_model.section_geometry.placement.anchoring import (
    AnchoredWallGeometry,
    airfoil_inward_depth_for_midsurface,
    anchored_wall_geometry,
    anchored_wall_geometry_from_policy,
    default_construction_for_kind,
    midsurface_offset_from_reference_m,
    physical_face_offsets_from_reference_m,
)
from blade_geometry_model.section_geometry.placement.collision import (
    PlacementGeometryError,
    assert_anchored_wall_feasible,
    find_placement_collisions,
)
from blade_geometry_model.section_geometry.placement.handoff import (
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
