"""2-D point transforms for S/B outlines and S-prime centroid shifts.

Used by Scale airfoil stacking and ContourStressRecovery consumers.
"""

from __future__ import annotations

from typing import Union

import numpy as np


def rotate_points_2d(
    points: np.ndarray,
    angle_deg: Union[float, np.ndarray],
) -> np.ndarray:
    """Rotate 2-D points counter-clockwise about the origin.

    Accepts (N, 2) with a scalar angle, or (n_eval, n_pts, 2) with
    (n_eval,) angles.
    """
    points = np.asarray(points, dtype=np.float64)
    angle_deg_arr = np.asarray(angle_deg, dtype=np.float64)
    angle_rad = np.deg2rad(angle_deg_arr)
    cos_t = np.cos(angle_rad)
    sin_t = np.sin(angle_rad)

    if points.ndim == 2:
        if cos_t.ndim != 0:
            raise ValueError(
                "rotate_points_2d with (N, 2) points requires a scalar angle_deg"
            )
        rot = np.array([[cos_t, -sin_t], [sin_t, cos_t]], dtype=np.float64)
        return points @ rot.T

    if points.ndim == 3:
        if cos_t.ndim == 0:
            cos_t = np.full(points.shape[0], cos_t)
            sin_t = np.full(points.shape[0], sin_t)
        if cos_t.shape != (points.shape[0],):
            raise ValueError(
                f"angle_deg length {cos_t.shape} must match points batch "
                f"dimension {points.shape[0]}"
            )
        x = points[..., 0]
        y = points[..., 1]
        return np.stack(
            [
                x * cos_t[:, None] - y * sin_t[:, None],
                x * sin_t[:, None] + y * cos_t[:, None],
            ],
            axis=-1,
        )

    raise ValueError(
        f"points must have shape (N, 2) or (n_eval, n_pts, 2); got {points.shape}"
    )


def transform_points_S_to_B(
    points: np.ndarray,
    twist_deg: Union[float, np.ndarray],
) -> np.ndarray:
    """Station (S) to Blade (B) via inverse twist rotation."""
    return rotate_points_2d(points, -np.asarray(twist_deg, dtype=np.float64))


def transform_points_B_to_S(
    points: np.ndarray,
    twist_deg: Union[float, np.ndarray],
) -> np.ndarray:
    """Blade (B) to Station (S) via forward twist rotation."""
    return rotate_points_2d(points, twist_deg)


def transform_points_S_to_Sprime(
    points: np.ndarray,
    centroid_S: np.ndarray,
) -> np.ndarray:
    """S to S-prime by translating the origin onto the section centroid."""
    points = np.asarray(points, dtype=np.float64)
    centroid_S = np.asarray(centroid_S, dtype=np.float64)
    return points - centroid_S


def transform_points_Sprime_to_S(
    points: np.ndarray,
    centroid_S: np.ndarray,
) -> np.ndarray:
    """S-prime to S (inverse of transform_points_S_to_Sprime)."""
    points = np.asarray(points, dtype=np.float64)
    centroid_S = np.asarray(centroid_S, dtype=np.float64)
    return points + centroid_S


def transform_points_aerofoil_to_S(
    points: np.ndarray,
    pivot_xy: np.ndarray,
) -> np.ndarray:
    """Chord-scaled aerofoil (x, y) to Station frame via axis swap + pivot."""
    points = np.asarray(points, dtype=np.float64)
    pivot_xy = np.asarray(pivot_xy, dtype=np.float64)
    swapped = np.empty_like(points)
    swapped[..., 0] = points[..., 1]
    swapped[..., 1] = points[..., 0]
    return swapped - np.array([pivot_xy[1], pivot_xy[0]], dtype=np.float64)


def rR_from_zL_scalars(
    zL: Union[float, np.ndarray],
    *,
    R: float,
    L: float,
) -> Union[float, np.ndarray]:
    """Convert B-frame z/L to R-frame r/R from scalars alone.

    Assumes the root sits at r = R - L. Prefer ScaledGeometry.rR_grid when available.
    """
    zL_arr = np.asarray(zL, dtype=np.float64)
    return (zL_arr * float(L) + float(R) - float(L)) / float(R)
