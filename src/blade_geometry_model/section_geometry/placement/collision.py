"""Placement collision checks — raise, do not repair (WP3 D3.3)."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from blade_structure_utilities.geometry_utilities.section_geometry.placement.anchoring import (
    AnchoredWallGeometry,
)


class PlacementGeometryError(ValueError):
    """A placement realisation is geometrically invalid."""


@dataclass(frozen=True)
class PlacementCollision:
    """One detected collision / inversion under a placement change."""

    code: str
    message: str
    wall_kind: str | None = None


def assert_anchored_wall_feasible(wall: AnchoredWallGeometry) -> None:
    """Refuse non-positive thickness or inverted free-face offsets."""
    label = wall.placement_policy_id or f"{wall.reference_face}/{wall.build_direction}"
    if not (wall.thickness_m > 0.0):
        raise PlacementGeometryError(
            f"placement {label!r}: thickness_m must be positive, got {wall.thickness_m!r}"
        )
    free = float(wall.free_face_offset_m)
    if wall.second_free_face_offset_m is None:
        if free <= 0.0:
            raise PlacementGeometryError(
                f"placement {label!r}: free face offset must be "
                f"positive for a one-sided wall, got {free!r}"
            )
        return
    second = float(wall.second_free_face_offset_m)
    if free * second >= 0.0:
        raise PlacementGeometryError(
            f"placement {label!r}: paired free faces must sit on "
            f"opposite sides of the midsurface, got {free!r} and {second!r}"
        )
    if abs(free) + abs(second) <= 0.0:
        raise PlacementGeometryError(
            f"placement {label!r}: web faces collapsed to zero width"
        )


def find_placement_collisions(
    walls: Sequence[AnchoredWallGeometry],
    *,
    skin_cavity_clearance_m: float | None = None,
    cap_into_cavity_m: float | None = None,
    cavity_depth_m: float | None = None,
) -> tuple[PlacementCollision, ...]:
    """Detect cap/skin intrusion, negative cavity, and wall inversions."""
    findings: list[PlacementCollision] = []
    for wall in walls:
        try:
            assert_anchored_wall_feasible(wall)
        except PlacementGeometryError as exc:
            findings.append(
                PlacementCollision(
                    code="wall_inversion",
                    message=str(exc),
                    wall_kind=wall.reference_face,
                )
            )
    if (
        skin_cavity_clearance_m is not None
        and cavity_depth_m is not None
        and float(skin_cavity_clearance_m) > float(cavity_depth_m) + 1e-15
    ):
        findings.append(
            PlacementCollision(
                code="negative_cavity",
                message=(
                    f"skin cavity clearance {skin_cavity_clearance_m!r} exceeds "
                    f"available cavity depth {cavity_depth_m!r}"
                ),
                wall_kind="skin",
            )
        )
    if (
        cap_into_cavity_m is not None
        and cavity_depth_m is not None
        and float(cap_into_cavity_m) > float(cavity_depth_m) + 1e-15
    ):
        findings.append(
            PlacementCollision(
                code="cap_skin_separation",
                message=(
                    f"cap cavity build {cap_into_cavity_m!r} exceeds available "
                    f"cavity depth {cavity_depth_m!r}"
                ),
                wall_kind="cap",
            )
        )
    return tuple(findings)


def raise_on_placement_collisions(
    walls: Iterable[AnchoredWallGeometry],
    **kwargs: float | None,
) -> None:
    """Raise :class:`PlacementGeometryError` when any collision is found."""
    hits = find_placement_collisions(tuple(walls), **kwargs)
    if hits:
        detail = "; ".join(f"{h.code}: {h.message}" for h in hits)
        raise PlacementGeometryError(detail)


__all__ = [
    "PlacementCollision",
    "PlacementGeometryError",
    "assert_anchored_wall_feasible",
    "find_placement_collisions",
    "raise_on_placement_collisions",
]
