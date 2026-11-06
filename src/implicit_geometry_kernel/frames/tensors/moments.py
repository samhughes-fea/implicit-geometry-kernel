"""Bending and torsion moment vectors between the S and B frames."""

from typing import Tuple, Union

import numpy as np


def transform_moments_B_to_S(
    M_YB: Union[float, np.ndarray],
    M_XB: Union[float, np.ndarray] = 0.0,
    twist_deg: Union[float, np.ndarray] = 0.0,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Transform bending moments from Blade frame (B) to Station frame (S).
    
    Coordinate System Convention:
    - **Blade frame (B)**: Fixed frame, does NOT rotate with blade twist
      - X_B axis: Flapwise direction (perpendicular to rotor plane)
      - Y_B axis: Edgewise direction (parallel to rotor plane)
      - Z_B axis: Spanwise direction (root to tip)
    - **Station frame (S)**: Rotates with blade twist at each station
      - X_S axis: Chord-normal direction (perpendicular to chord line)
      - Y_S axis: Chord-parallel direction (along chord line)
      - Z_S axis: Spanwise direction (same as Z_B)
    - **Twist convention**: Positive twist_deg indicates CCW rotation from B to S
    
    Physical Interpretation for Flapwise Loading:
    - Input: M_YB (flapwise moment about Y_B axis in B frame)
    - When blade is twisted, this B-frame moment projects onto two S-frame components:
      - M_YS = M_YB × cos(θ): Moment about Y_S (chord-parallel axis) → creates stress in X_S (chord-normal)
      - M_XS = M_YB × sin(θ): Moment about X_S (chord-normal axis) → creates stress in Y_S (chord-parallel)
    - For stress calculations, we use M_YS because fiber distances are measured in X_S direction
    
    Transformation Formulas (rotation matrix for angle θ = twist_deg):
    - M_YS = M_YB × cos(θ) - M_XB × sin(θ)
    - M_XS = M_YB × sin(θ) + M_XB × cos(θ)
    
    For flapwise-only loading (M_XB ≈ 0):
    - M_YS ≈ M_YB × cos(θ)  (primary component for stress calculation)
    - M_XS ≈ M_YB × sin(θ)  (secondary component, usually small for small twist)
    
    Parameters
    ----------
    M_YB : float or np.ndarray
        Bending moment about Y_B axis (flapwise, edgewise-perpendicular) in B frame [N·m]
    M_XB : float or np.ndarray, optional
        Bending moment about X_B axis (edgewise, flapwise-perpendicular) in B frame [N·m].
        Default 0.0 (flapwise loading only).
    twist_deg : float or np.ndarray, optional
        Twist angle [degrees]. Positive = CCW rotation from B to S frame.
        Can be scalar or array matching length of M_YB. Default 0.0.
    
    Returns
    -------
    tuple
        (M_YS, M_XS) - Bending moments in Station frame (S) [N·m]
        - M_YS: Moment about Y_S axis (chord-parallel) - used for stress calculation
        - M_XS: Moment about X_S axis (chord-normal) - typically small for flapwise loading
        Returns arrays if inputs are arrays, scalars if inputs are scalars.
    
    Examples
    --------
    Flapwise loading with 10° twist:
    >>> transform = CoordinateFrameTransforms()
    >>> M_YB = 1000.0  # N·m
    >>> M_YS, M_XS = transform.transform_moments_B_to_S(M_YB, M_XB=0.0, twist_deg=10.0)
    >>> # M_YS ≈ 985 N·m (cos(10°) ≈ 0.985)
    >>> # M_XS ≈ 174 N·m (sin(10°) ≈ 0.174)
    
    Notes
    -----
    - This method is the complement to transform_second_moments_S_to_B
    - For stress calculation: σ = M_YS × c_XS / I_YS (all in S frame)
    - Moments and fiber distances must be in the same frame for correct stress
    """
    # Convert to numpy arrays for vectorized operations
    M_YB = np.asarray(M_YB, dtype=np.float64)
    M_XB = np.asarray(M_XB, dtype=np.float64) if not isinstance(M_XB, (int, float)) or M_XB != 0.0 else np.zeros_like(M_YB)
    twist_deg = np.asarray(twist_deg, dtype=np.float64) if not isinstance(twist_deg, (int, float)) or twist_deg != 0.0 else np.zeros_like(M_YB)
    
    # Ensure all arrays have compatible shapes
    if M_XB.shape != M_YB.shape:
        if M_XB.size == 1:
            M_XB = np.full_like(M_YB, float(M_XB))
        else:
            raise ValueError(f"M_XB shape {M_XB.shape} incompatible with M_YB shape {M_YB.shape}")
    if twist_deg.shape != M_YB.shape:
        if twist_deg.size == 1:
            twist_deg = np.full_like(M_YB, float(np.asarray(twist_deg).reshape(-1)[0]))
        elif (
            twist_deg.ndim == 1
            and M_YB.ndim >= 1
            and twist_deg.shape[0] == M_YB.shape[0]
        ):
            # Per-station twist with trailing load axes (n_span,) -> (n_span, ...)
            twist_deg = np.broadcast_to(
                twist_deg.reshape(twist_deg.shape[0], *([1] * (M_YB.ndim - 1))),
                M_YB.shape,
            ).copy()
        elif (
            twist_deg.ndim == M_YB.ndim
            and twist_deg.shape[0] == M_YB.shape[0]
            and all(
                t == 1 or t == m
                for t, m in zip(twist_deg.shape[1:], M_YB.shape[1:])
            )
        ):
            twist_deg = np.broadcast_to(twist_deg, M_YB.shape).copy()
        else:
            raise ValueError(f"twist_deg shape {twist_deg.shape} incompatible with M_YB shape {M_YB.shape}")
    
    # Input validation: Check for finite values
    if not np.all(np.isfinite(M_YB)):
        if M_YB.ndim == 0:
            raise ValueError(f"Invalid M_YB value (must be finite): {float(M_YB)}")
        invalid_indices = np.where(~np.isfinite(M_YB))[0]
        raise ValueError(f"Invalid M_YB values (must be finite): {M_YB[invalid_indices]} at indices {invalid_indices}")
    
    if not np.all(np.isfinite(M_XB)):
        if M_XB.ndim == 0:
            raise ValueError(f"Invalid M_XB value (must be finite): {float(M_XB)}")
        invalid_indices = np.where(~np.isfinite(M_XB))[0]
        raise ValueError(f"Invalid M_XB values (must be finite): {M_XB[invalid_indices]} at indices {invalid_indices}")
    
    if not np.all(np.isfinite(twist_deg)):
        if twist_deg.ndim == 0:
            raise ValueError(f"Invalid twist_deg value (must be finite): {float(twist_deg)}")
        invalid_indices = np.where(~np.isfinite(twist_deg))[0]
        raise ValueError(f"Invalid twist_deg values (must be finite): {twist_deg[invalid_indices]} at indices {invalid_indices}")
    
    # Convert twist to radians
    twist_rad = np.deg2rad(twist_deg)
    
    # Compute trigonometric values
    cos_theta = np.cos(twist_rad)
    sin_theta = np.sin(twist_rad)
    
    # Apply rotation transformation: R(θ) @ [M_YB, M_XB]^T = [M_YS, M_XS]^T
    # Standard 2D rotation matrix for coordinate transformation
    M_YS = M_YB * cos_theta - M_XB * sin_theta
    M_XS = M_YB * sin_theta + M_XB * cos_theta
    
    # Return scalars if inputs were scalars
    if isinstance(M_YB, (int, float)) or M_YB.ndim == 0:
        return float(M_YS), float(M_XS)
    return M_YS, M_XS

def transform_moments_S_to_B(
    M_YS: Union[float, np.ndarray],
    M_XS: Union[float, np.ndarray] = 0.0,
    twist_deg: Union[float, np.ndarray] = 0.0,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Transform bending moments from Station frame (S) to Blade frame (B).
    
    This is the inverse of transform_moments_B_to_S. Uses inverse rotation (-θ).
    
    Transformation Formulas (inverse rotation for angle θ = -twist_deg):
    - M_YB = M_YS × cos(θ) + M_XS × sin(θ)
    - M_XB = -M_YS × sin(θ) + M_XS × cos(θ)
    
    Parameters
    ----------
    M_YS : float or np.ndarray
        Moment about Y_S axis (chord-parallel) in S frame [N·m]
    M_XS : float or np.ndarray, optional
        Moment about X_S axis (chord-normal) in S frame [N·m]. Default 0.0.
    twist_deg : float or np.ndarray, optional
        Twist angle [degrees]. Positive = CCW rotation from B to S.
        Default 0.0.
    
    Returns
    -------
    tuple
        (M_YB, M_XB) - Bending moments in Blade frame (B) [N·m]
    
    Notes
    -----
    This uses the inverse rotation matrix: R(-θ) = R(θ)^T
    For rotation matrices: R(-θ) = R(θ)^T (transpose = inverse)
    """
    # Use negative twist for inverse transformation
    # R(-θ) @ [M_YS, M_XS]^T = R(θ)^T @ [M_YS, M_XS]^T = [M_YB, M_XB]^T
    M_YB, M_XB = transform_moments_B_to_S(M_YS, M_XS, -twist_deg)
    return M_YB, M_XB

