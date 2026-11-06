"""Force vectors between the S and B frames."""

from typing import Tuple, Union

import numpy as np
from implicit_geometry_kernel.frames.tensors.moments import transform_moments_B_to_S


def transform_forces_B_to_S(
    V_XB: Union[float, np.ndarray],
    V_YB: Union[float, np.ndarray] = 0.0,
    twist_deg: Union[float, np.ndarray] = 0.0,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Transform force/shear vectors from Blade frame (B) to Station frame (S).

    Applies the same 2D rotation as :func:`transform_moments_B_to_S`.
    The mathematics are identical for any 2D vector; this method provides
    correct parameter naming for force / shear-force quantities.

    Coordinate System Convention:
    - **Blade frame (B)**: Fixed frame, does NOT rotate with blade twist
      - X_B axis: Flapwise direction (perpendicular to rotor plane)
      - Y_B axis: Edgewise direction (parallel to rotor plane)
    - **Station frame (S)**: Rotates with blade twist at each station
      - X_S axis: Chord-normal direction (perpendicular to chord line)
      - Y_S axis: Chord-parallel direction (along chord line)
    - **Twist convention**: Positive twist_deg indicates CCW rotation from B to S

    Transformation Formulas (rotation matrix for angle θ = twist_deg):
    - V_XS = V_XB × cos(θ) - V_YB × sin(θ)
    - V_YS = V_XB × sin(θ) + V_YB × cos(θ)

    Parameters
    ----------
    V_XB : float or np.ndarray
        Force component along X_B (flapwise) [N]
    V_YB : float or np.ndarray, optional
        Force component along Y_B (edgewise) [N]. Default 0.0.
    twist_deg : float or np.ndarray, optional
        Twist angle [degrees]. Positive = CCW rotation from B to S.
        Can be scalar or array matching length of V_XB. Default 0.0.

    Returns
    -------
    tuple
        (V_XS, V_YS) - Force components in Station frame (S) [N]
        - V_XS: Force along X_S (chord-normal direction)
        - V_YS: Force along Y_S (chord-parallel direction)
        Returns arrays if inputs are arrays, scalars if inputs are scalars.
    """
    V_XS, V_YS = transform_moments_B_to_S(
        M_YB=V_XB, M_XB=V_YB, twist_deg=twist_deg,
    )
    return V_XS, V_YS

def transform_forces_S_to_B(
    V_XS: Union[float, np.ndarray],
    V_YS: Union[float, np.ndarray] = 0.0,
    twist_deg: Union[float, np.ndarray] = 0.0,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Transform force/shear vectors from Station frame (S) to Blade frame (B).

    This is the inverse of :func:`transform_forces_B_to_S`.
    Uses inverse rotation (-θ).

    Parameters
    ----------
    V_XS : float or np.ndarray
        Force component along X_S (chord-normal) [N]
    V_YS : float or np.ndarray, optional
        Force component along Y_S (chord-parallel) [N]. Default 0.0.
    twist_deg : float or np.ndarray, optional
        Twist angle [degrees]. Positive = CCW rotation from B to S.
        Default 0.0.

    Returns
    -------
    tuple
        (V_XB, V_YB) - Force components in Blade frame (B) [N]
    """
    return transform_forces_B_to_S(V_XS, V_YS, -twist_deg)
