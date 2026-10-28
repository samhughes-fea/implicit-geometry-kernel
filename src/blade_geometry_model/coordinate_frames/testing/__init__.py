"""Test helpers for coordinate_frames (not a public engineering API)."""

from .invariance import (
    FRAME_B,
    assert_planar_pierce_rotates,
    assert_scalar_invariant_under_rotation,
    rotate_midlines_in_b,
    rotate_points_in_b,
)

__all__ = [
    "FRAME_B",
    "assert_planar_pierce_rotates",
    "assert_scalar_invariant_under_rotation",
    "rotate_midlines_in_b",
    "rotate_points_in_b",
]
