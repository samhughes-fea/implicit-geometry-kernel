"""Wall thickness between the S' and B frames."""

from typing import Union

import numpy as np


def transform_thickness_B_to_Sprime(
    thickness_B: Union[float, np.ndarray],
    twist_deg: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    """
    Transform a scalar section depth from B frame to S' frame.

    This is a low-order geometric interpretation helper used when comparing a
    depth defined along the B-frame flapwise axis with a depth measured along the
    S'-frame chord-normal axis.
    """
    thickness_B = np.asarray(thickness_B, dtype=np.float64)
    twist_deg = np.asarray(twist_deg, dtype=np.float64)
    if thickness_B.shape != twist_deg.shape and thickness_B.ndim != 0 and twist_deg.ndim != 0:
        raise ValueError(
            f"thickness_B and twist_deg must have same shape. "
            f"Got thickness_B.shape={thickness_B.shape}, twist_deg.shape={twist_deg.shape}"
        )
    twist_rad = np.deg2rad(twist_deg)
    cos_twist = np.cos(twist_rad)
    cos_twist_safe = np.where(np.abs(cos_twist) < 1e-10, 1.0, cos_twist)
    return thickness_B / cos_twist_safe

def transform_thickness_Sprime_to_B(
    thickness_Sprime: Union[float, np.ndarray],
    twist_deg: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    """
    Transform a scalar section depth from S' frame to B frame.

    Inverse of :func:`transform_thickness_B_to_Sprime`.
    """
    thickness_Sprime = np.asarray(thickness_Sprime, dtype=np.float64)
    twist_deg = np.asarray(twist_deg, dtype=np.float64)
    if thickness_Sprime.shape != twist_deg.shape and thickness_Sprime.ndim != 0 and twist_deg.ndim != 0:
        raise ValueError(
            f"thickness_Sprime and twist_deg must have same shape. "
            f"Got thickness_Sprime.shape={thickness_Sprime.shape}, twist_deg.shape={twist_deg.shape}"
        )
    twist_rad = np.deg2rad(twist_deg)
    return thickness_Sprime * np.cos(twist_rad)
