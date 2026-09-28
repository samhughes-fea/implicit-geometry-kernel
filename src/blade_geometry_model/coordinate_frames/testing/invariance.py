"""Shared planar frame-invariance checks for B-consuming packages.

A quantity computed on midlines in canonical B must transform correctly when
those midlines are rotated in the section plane: scalars (areas, ``GJ``) are
invariant; planar pierces (centroids, shear centre) rotate with the points.

B-consuming packages adopt these helpers rather than inventing their own
rotation checks. Inputs and outputs are tagged ``FRAME_B``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

from blade_structure_utilities.geometry_utilities.coordinate_frames.frame_tags import (
    FRAME_B,
    require_frame_tag,
)

__all__ = [
    "FRAME_B",
    "assert_planar_pierce_rotates",
    "assert_scalar_invariant_under_rotation",
    "rotate_midlines_in_b",
    "rotate_points_in_b",
]


def rotate_points_in_b(
    points: np.ndarray,
    angle_rad: float,
    *,
    frame: str = FRAME_B,
) -> np.ndarray:
    """Rotate ``(N, 2)`` section-plane points by ``angle_rad`` about the origin in B."""
    require_frame_tag(frame)
    if frame != FRAME_B:
        raise ValueError(f"planar section rotation is defined in {FRAME_B!r}; got {frame!r}")
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] < 2:
        raise ValueError("points must be (N, 2) in the section plane")
    c, s = float(np.cos(angle_rad)), float(np.sin(angle_rad))
    x, y = pts[:, 0], pts[:, 1]
    out = np.column_stack([c * x - s * y, s * x + c * y])
    return out


def _copy_strip_with_points(strip: Any, points: np.ndarray) -> Any:
    if isinstance(strip, Mapping):
        out = dict(strip)
        for key in (
            "midline_twisted_pitch_axis_blade_frame",
            "xy",
            "points",
            "midline_xy",
            "coords",
        ):
            if key in out:
                out[key] = points
                return out
        out["xy"] = points
        return out
    # Frozen dataclass / ShellMidlineStrip: rebuild via replace when available.
    from dataclasses import is_dataclass, replace

    if is_dataclass(strip) and not isinstance(strip, type):
        for field_name in (
            "midline_twisted_pitch_axis_blade_frame",
            "xy",
            "points",
            "midline_xy",
        ):
            if hasattr(strip, field_name):
                return replace(strip, **{field_name: points})
    raise TypeError(f"cannot rotate strip of type {type(strip)!r}")


def rotate_midlines_in_b(
    midlines: Sequence[Any],
    angle_rad: float,
    *,
    frame: str = FRAME_B,
) -> list[Any]:
    """Return midlines with every polyline rotated in ``FRAME_B``."""
    require_frame_tag(frame)
    rotated: list[Any] = []
    for strip in midlines:
        pts = None
        for attr in (
            "midline_twisted_pitch_axis_blade_frame",
            "xy",
            "points",
            "midline_xy",
            "coords",
        ):
            val = getattr(strip, attr, None)
            if val is None and isinstance(strip, Mapping):
                val = strip.get(attr)
            if val is not None:
                arr = np.asarray(val, dtype=float)
                if arr.ndim == 2 and arr.shape[1] >= 2:
                    pts = arr[:, :2]
                    break
        if pts is None:
            raise TypeError(f"cannot extract points from {type(strip)!r}")
        rotated.append(_copy_strip_with_points(strip, rotate_points_in_b(pts, angle_rad, frame=frame)))
    return rotated


def assert_scalar_invariant_under_rotation(
    midlines: Sequence[Any],
    *,
    scalar: Callable[[Sequence[Any]], float],
    angles_rad: Sequence[float] = (0.3, -0.7, 1.1),
    rtol: float = 1e-9,
    atol: float = 0.0,
    frame: str = FRAME_B,
) -> None:
    """``scalar(midlines)`` must be unchanged under planar rotation in B."""
    require_frame_tag(frame)
    base = float(scalar(midlines))
    for angle in angles_rad:
        got = float(scalar(rotate_midlines_in_b(midlines, float(angle), frame=frame)))
        if not np.isclose(got, base, rtol=rtol, atol=atol):
            raise AssertionError(
                f"scalar not invariant under B rotation by {angle} rad: "
                f"base={base!r}, rotated={got!r} (frame={frame!r})"
            )


def assert_planar_pierce_rotates(
    midlines: Sequence[Any],
    *,
    pierce: Callable[[Sequence[Any]], tuple[float, float]],
    angles_rad: Sequence[float] = (0.3, -0.7, 1.1),
    rtol: float = 1e-9,
    atol: float = 1e-12,
    frame: str = FRAME_B,
) -> None:
    """``pierce(midlines)`` must rotate with the midlines in B."""
    require_frame_tag(frame)
    x0, y0 = pierce(midlines)
    for angle in angles_rad:
        c, s = float(np.cos(angle)), float(np.sin(angle))
        expected = (c * x0 - s * y0, s * x0 + c * y0)
        got = pierce(rotate_midlines_in_b(midlines, float(angle), frame=frame))
        if not (
            np.isclose(got[0], expected[0], rtol=rtol, atol=atol)
            and np.isclose(got[1], expected[1], rtol=rtol, atol=atol)
        ):
            raise AssertionError(
                f"pierce did not rotate under B rotation by {angle} rad: "
                f"expected={expected!r}, got={got!r} (frame={frame!r})"
            )
