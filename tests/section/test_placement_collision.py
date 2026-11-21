"""Placement collision checks raise and report; they never repair."""

from __future__ import annotations

import pytest

from implicit_geometry_kernel.section.layout.anchoring import AnchoredWallGeometry
from implicit_geometry_kernel.section.layout.collision import (
    PlacementGeometryError,
    assert_anchored_wall_feasible,
    find_placement_collisions,
    raise_on_placement_collisions,
)


def _wall(t: float = 0.004, free: float = 0.004, second: float | None = None) -> AnchoredWallGeometry:
    return AnchoredWallGeometry(
        reference_face="oml",
        build_direction="inward",
        thickness_m=t,
        reference_offset_m=0.0,
        free_face_offset_m=free,
        midsurface_offset_m=0.5 * t,
        second_free_face_offset_m=second,
    )


def test_feasible_one_sided_and_paired_walls_pass() -> None:
    assert_anchored_wall_feasible(_wall())
    assert_anchored_wall_feasible(_wall(free=0.002, second=-0.002))


@pytest.mark.parametrize(
    ("wall", "match"),
    [
        (_wall(t=0.0), "thickness_m must be positive"),
        (_wall(free=-0.001), "one-sided wall"),
        (_wall(free=0.002, second=0.002), "opposite sides"),
    ],
)
def test_infeasible_walls_raise(wall: AnchoredWallGeometry, match: str) -> None:
    with pytest.raises(PlacementGeometryError, match=match):
        assert_anchored_wall_feasible(wall)


def test_findings_cover_inversion_negative_cavity_and_cap_intrusion() -> None:
    hits = find_placement_collisions(
        [_wall(), _wall(t=-1.0)],
        skin_cavity_clearance_m=0.02,
        cap_into_cavity_m=0.03,
        cavity_depth_m=0.01,
    )
    assert [h.code for h in hits] == ["wall_inversion", "negative_cavity", "cap_skin_separation"]
    assert [h.wall_kind for h in hits] == ["oml", "skin", "cap"]


def test_no_findings_when_everything_fits() -> None:
    assert find_placement_collisions([_wall()], skin_cavity_clearance_m=0.005, cavity_depth_m=0.01) == ()
    raise_on_placement_collisions([_wall()], cavity_depth_m=0.01)


def test_raise_joins_every_finding() -> None:
    with pytest.raises(PlacementGeometryError, match="negative_cavity: .*; cap_skin_separation: "):
        raise_on_placement_collisions(
            [_wall()], skin_cavity_clearance_m=0.02, cap_into_cavity_m=0.03, cavity_depth_m=0.01
        )
