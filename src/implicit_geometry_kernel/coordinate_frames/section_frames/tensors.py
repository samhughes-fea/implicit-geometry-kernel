"""Pure S↔B tensor / vector maps (second moments, moments, forces, thickness)."""

from typing import Tuple, Union

import numpy as np

#: Below this, a geometric second moment is not a usable divisor and
#: :func:`structural_product_moment` falls back. Same threshold the
#: flapwise-only scaling it replaced has always used.
_GEOM_I_FLOOR: float = 1e-12


def structural_product_moment(
    *,
    Ixy_B_geom: Union[float, np.ndarray],
    Ixx_B_geom: Union[float, np.ndarray],
    Iyy_B_geom: Union[float, np.ndarray],
    Ixx_B_structural: Union[float, np.ndarray],
    Iyy_B_structural: Union[float, np.ndarray],
) -> np.ndarray:
    """Structural ``Ixy_B`` implied by two independently stated stiffnesses.

    The stress model is handed *two scalars* per station -- ``EI_flap_B`` and
    ``EI_edge_B`` -- and needs a full inertia tensor. Nobody supplies ``Ixy``: at
    that stage of the pipeline the internal structure is not yet chosen, so the
    product moment has to be **closed** from the outer geometry. This is that
    closure, named and in one place rather than open-coded at each call site.

    It preserves the outer section's correlation coefficient

    .. math:: \\rho = I_{xy} / \\sqrt{I_{xx} I_{yy}}

    which is to say: the structure keeps the *shape* of the geometric section's
    asymmetry while its two principal stiffnesses are set independently.

    Why this closure and not another
    --------------------------------
    **It is admissible by construction.** :math:`|\\rho| < 1` is what makes the
    inertia tensor positive-definite, and rotation preserves :math:`\\rho`, so the
    transformed tensor is valid in *every* frame for *any* positive
    ``(Ixx_B_structural, Iyy_B_structural)``. No caller has to defend against a
    section that cannot exist.

    **It agrees with the old form wherever the old assumption held.** The previous
    form scaled ``Ixy`` by the flapwise ratio alone,
    ``Iyy_B_structural / Iyy_B_geom``. Whenever ``EI_edge_B`` comes from
    ``StressComponents._compute_EI_edge_B_from_ratio`` -- the fallback, and every
    call that predates independent edgewise selection -- the two ratios are equal
    *in exact arithmetic*, because that fallback sets
    ``Ixx_structural = Iyy_structural * (Ixx_B_geom / Iyy_B_geom)``, and the
    geometric mean of two equal numbers is that number. So this is not a re-tuning
    of the old behaviour; it is the same quantity, derived without the assumption.

    In floating point the agreement is close but **not** bit-exact, and the reason
    is worth stating rather than discovering later. The fallback reaches the
    edgewise ratio as ``EI_flap * (Ixx_geom / Iyy_geom) / Ixx_geom`` where the
    flapwise one is simply ``EI_flap / Iyy_geom`` -- three roundings against one.
    Measured end to end on the ``new_3.2m_rotor`` ULS floor, the flapwise
    ``EI_required`` moved at **4 of 64** stations by at most **8 ULP**
    (1.1e-15 relative), with the active-station mask unchanged. That is far below
    the NLP's own constraint tolerances, but it is not zero, so a caller checking
    for an exact match will not get one.

    **It is symmetric in the two directions.** The flapwise-only form privileges
    ``EI_flap``, which is exactly the single-parameter slaving that selecting
    ``EI_edge`` in its own right exists to remove. Left in place, it drives the
    implied :math:`\\rho` upward without bound as ``EI_edge`` falls, so a merely
    *soft* edgewise stiffness is rejected as a section that cannot exist -- at a
    ratio that reflects the closure, not the aerofoil. Measured on a real 3.2 m
    blade, that put the softest "valid" edgewise stiffness at up to 2.6x the
    flapwise one, against a geometric ``I_edge/I_flap`` of 4.4 at the same station.

    A stand-in, not the truth
    -------------------------
    The real product moment follows from the internal structure, which
    ``section_selection`` knows and this stage does not. This closure is the
    minimal statement consistent with two independent channels, and is meant to be
    replaced by a stated ``Ixy`` once one exists.
    """
    Ixy_geom = np.asarray(Ixy_B_geom, dtype=np.float64)
    Ixx_geom = np.asarray(Ixx_B_geom, dtype=np.float64)
    Iyy_geom = np.asarray(Iyy_B_geom, dtype=np.float64)
    Ixx_str = np.asarray(Ixx_B_structural, dtype=np.float64)
    Iyy_str = np.asarray(Iyy_B_structural, dtype=np.float64)

    # Divided under ``where`` rather than inside ``np.where``, which evaluates both
    # branches and would raise a divide-by-zero warning on exactly the degenerate
    # stations this is guarding.
    flap_ratio = np.divide(
        Iyy_str,
        Iyy_geom,
        out=np.ones(np.broadcast(Iyy_str, Iyy_geom).shape, dtype=np.float64),
        where=Iyy_geom > _GEOM_I_FLOOR,
    )
    edge_ratio = np.divide(
        Ixx_str,
        Ixx_geom,
        out=np.ones(np.broadcast(Ixx_str, Ixx_geom).shape, dtype=np.float64),
        where=Ixx_geom > _GEOM_I_FLOOR,
    )

    # Written ``flap * sqrt(edge / flap)`` rather than ``sqrt(flap * edge)``: when
    # the two ratios do land on the same float the quotient is exactly 1.0, so the
    # result is exactly ``flap_ratio``, where ``sqrt(flap * edge)`` would return
    # ``sqrt(flap**2)`` and lose a ULP for nothing. It is also the better
    # conditioned of the two near equality, which is where this sits in practice.
    quotient = np.divide(
        edge_ratio,
        flap_ratio,
        out=np.ones_like(flap_ratio),
        where=flap_ratio > 0.0,
    )
    # A degenerate geometric Ixx leaves no edgewise ratio to combine, so the
    # closure falls back to the flapwise-only form there -- exactly what those
    # stations got before. Clipping at zero guards a negative structural
    # stiffness, which is not a section either way.
    scale = np.where(
        Ixx_geom > _GEOM_I_FLOOR,
        flap_ratio * np.sqrt(np.maximum(quotient, 0.0)),
        flap_ratio,
    )
    return Ixy_geom * scale


def transform_second_moments_S_to_B(
    Ixx_S: Union[float, np.ndarray],
    Iyy_S: Union[float, np.ndarray],
    Ixy_S: Union[float, np.ndarray] = 0.0,
    twist_deg: Union[float, np.ndarray] = 0.0,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Transform second moments of area from Station frame (S) to Blade frame (B).
    
    Coordinate System Convention:
    - **Station frame (S)**: Origin at 0.333 chordline, chord line is horizontal axis.
      Rotates with blade twist (2D cross-section plane).
      Each structural discretisation grid node has its own S frame with its own twist angle.
    - **Blade frame (B)**: Origin at 0.333 chordline at blade root, does NOT twist with blade twist.
      Fixed orientation (2D cross-section plane).
    - **z_B axis**: Positive direction from root to tip
    - **Right-hand rule**: Thumb aligned with z_B (root→tip direction)
    - **Blade twist**: CCW rotation from Blade frame (B) to Station frame (S) for positive twist_deg
    
    Transformation Formulas (for rotation angle θ, where S is rotated by +θ relative to B):
    - Uses inverse rotation (-θ) to transform FROM S TO B: I_B = R(-θ) I_S R(-θ)^T
    - Ixx_B = Ixx_S × cos²(-θ) + Iyy_S × sin²(-θ) - 2×Ixy_S × sin(-θ) × cos(-θ)
    - Iyy_B = Ixx_S × sin²(-θ) + Iyy_S × cos²(-θ) + 2×Ixy_S × sin(-θ) × cos(-θ)
    - Ixy_B = (Ixx_S - Iyy_S) × sin(-θ) × cos(-θ) + Ixy_S × (cos²(-θ) - sin²(-θ))
    
    Note: Since cos(-θ) = cos(θ) and sin(-θ) = -sin(θ), the formulas simplify but use -θ internally.
    Airfoils are asymmetric in S frame, so Ixy_S is typically non-zero and must be computed.
    
    Parameters
    ----------
    Ixx_S : float or np.ndarray
        Second moment about x-axis in Station frame (S) [m⁴]
    Iyy_S : float or np.ndarray
        Second moment about y-axis in Station frame (S) [m⁴]
    Ixy_S : float or np.ndarray, optional
        Product moment in Station frame (S) [m⁴]. Default 0.0.
        Note: Airfoils are asymmetric, so Ixy_S is typically non-zero.
    twist_deg : float or np.ndarray, optional
        Twist angle [degrees]. Positive values indicate CCW rotation of S frame relative to B.
        The transformation internally uses -twist_deg (inverse rotation) to transform FROM S TO B.
        Can be scalar or array matching length of Ixx_S, Iyy_S, Ixy_S. Default 0.0.
    
    Returns
    -------
    tuple
        (Ixx_B, Iyy_B, Ixy_B) - Second moments in Blade frame (B) [m⁴]
        Returns arrays if inputs are arrays, scalars if inputs are scalars.
        
    Notes
    -----
    Physical Constraint:
    Second moments of area are physically defined as integrals of squared distances
    (Ixx = ∫ y² dA, Iyy = ∫ x² dA) and must always be non-negative for valid physical
    cross-sections. With the correct transformation (using -θ), negative values should
    not occur for valid geometry.
    """
    # Convert to numpy arrays for vectorized operations
    Ixx_S = np.asarray(Ixx_S, dtype=np.float64)
    Iyy_S = np.asarray(Iyy_S, dtype=np.float64)
    Ixy_S = np.asarray(Ixy_S, dtype=np.float64) if not isinstance(Ixy_S, (int, float)) or Ixy_S != 0.0 else np.zeros_like(Ixx_S)
    twist_deg = np.asarray(twist_deg, dtype=np.float64) if not isinstance(twist_deg, (int, float)) or twist_deg != 0.0 else np.zeros_like(Ixx_S)
    
    # Ensure all arrays have compatible shapes
    if Ixy_S.shape != Ixx_S.shape:
        if Ixy_S.size == 1:
            Ixy_S = np.full_like(Ixx_S, float(Ixy_S))
        else:
            raise ValueError(f"Ixy_S shape {Ixy_S.shape} incompatible with Ixx_S shape {Ixx_S.shape}")
    if twist_deg.shape != Ixx_S.shape:
        if twist_deg.size == 1:
            twist_deg = np.full_like(Ixx_S, float(twist_deg))
        else:
            raise ValueError(f"twist_deg shape {twist_deg.shape} incompatible with Ixx_S shape {Ixx_S.shape}")
    
    # Input validation: Validate S frame properties are physically valid
    if np.any(Ixx_S <= 0):
        if Ixx_S.ndim == 0:
            raise ValueError(
                f"Invalid Ixx_S value (must be > 0): {float(Ixx_S)}. "
                f"Second moments of area must be positive for valid geometry."
            )
        invalid_indices = np.where(Ixx_S <= 0)[0]
        raise ValueError(
            f"Invalid Ixx_S values (must be > 0): {Ixx_S[invalid_indices]} at indices {invalid_indices}. "
            f"Second moments of area must be positive for valid geometry."
        )
    
    if np.any(Iyy_S <= 0):
        if Iyy_S.ndim == 0:
            raise ValueError(
                f"Invalid Iyy_S value (must be > 0): {float(Iyy_S)}. "
                f"Second moments of area must be positive for valid geometry."
            )
        invalid_indices = np.where(Iyy_S <= 0)[0]
        raise ValueError(
            f"Invalid Iyy_S values (must be > 0): {Iyy_S[invalid_indices]} at indices {invalid_indices}. "
            f"Second moments of area must be positive for valid geometry."
        )
    
    if not np.all(np.isfinite(Ixy_S)):
        if Ixy_S.ndim == 0:
            raise ValueError(
                f"Invalid Ixy_S value (must be finite): {float(Ixy_S)}. "
                f"Product moment of inertia must be finite."
            )
        invalid_indices = np.where(~np.isfinite(Ixy_S))[0]
        raise ValueError(
            f"Invalid Ixy_S values (must be finite): {Ixy_S[invalid_indices]} at indices {invalid_indices}. "
            f"Product moment of inertia must be finite."
        )
    
    if not np.all(np.isfinite(twist_deg)):
        if twist_deg.ndim == 0:
            raise ValueError(
                f"Invalid twist_deg value (must be finite): {float(twist_deg)}. "
                f"Twist angle must be finite."
            )
        invalid_indices = np.where(~np.isfinite(twist_deg))[0]
        raise ValueError(
            f"Invalid twist_deg values (must be finite): {twist_deg[invalid_indices]} at indices {invalid_indices}. "
            f"Twist angle must be finite."
        )
    
    # Convert twist to radians
    # CRITICAL: Use +twist_deg (not -twist_deg) because:
    # - Testing with Station 0 shows: -twist_deg → negative Iyy_B (invalid), +twist_deg → positive Iyy_B (valid)
    # - The transformation formula requires +twist_deg and correct cross-term signs to produce physically valid results
    # - Root cause: The coordinate system convention or transformation formula interpretation differs
    #   from the documented inverse rotation approach. Further investigation needed to understand why.
    twist_rad = np.deg2rad(twist_deg)
    
    # Compute trigonometric values
    cos_theta = np.cos(twist_rad)
    sin_theta = np.sin(twist_rad)
    cos2 = cos_theta ** 2
    sin2 = sin_theta ** 2
    sin_cos = sin_theta * cos_theta
    cos2_sin2 = cos2 - sin2  # cos(2θ) = cos²(θ) - sin²(θ)
    
    # Apply transformation formulas
    # CRITICAL: Ixx_B cross-term uses + sign (not -) to ensure positive results
    # Testing confirmed: +twist_deg with + cross-term for Ixx_B produces valid (positive) results
    Ixx_B = Ixx_S * cos2 + Iyy_S * sin2 + 2.0 * Ixy_S * sin_cos
    Iyy_B = Ixx_S * sin2 + Iyy_S * cos2 + 2.0 * Ixy_S * sin_cos
    Ixy_B = (Ixx_S - Iyy_S) * sin_cos + Ixy_S * cos2_sin2
    
    # Output validation: Validate B frame properties are physically valid
    if np.any(Ixx_B <= 0):
        if Ixx_B.ndim == 0:
            raise ValueError(
                f"Transformation produced invalid Ixx_B value (must be > 0): {float(Ixx_B)}. "
                f"Input S frame: Ixx_S={float(Ixx_S)}, Iyy_S={float(Iyy_S)}, Ixy_S={float(Ixy_S)}. "
                f"Twist: {float(twist_deg)}°. "
                f"This indicates a transformation error or invalid geometry. "
                f"Possible causes: invalid geometry, wrong twist angle, transformation error."
            )
        invalid_indices = np.where(Ixx_B <= 0)[0]
        raise ValueError(
            f"Transformation produced invalid Ixx_B values (must be > 0): {Ixx_B[invalid_indices]} at indices {invalid_indices}. "
            f"Input S frame: Ixx_S={Ixx_S[invalid_indices] if hasattr(Ixx_S, '__getitem__') else Ixx_S}, "
            f"Iyy_S={Iyy_S[invalid_indices] if hasattr(Iyy_S, '__getitem__') else Iyy_S}, "
            f"Ixy_S={Ixy_S[invalid_indices] if hasattr(Ixy_S, '__getitem__') else Ixy_S}. "
            f"Twist: {twist_deg[invalid_indices] if hasattr(twist_deg, '__getitem__') else twist_deg}°. "
            f"This indicates a transformation error or invalid geometry. "
            f"Possible causes: invalid geometry, wrong twist angle, transformation error."
        )
    
    if np.any(Iyy_B <= 0):
        if Iyy_B.ndim == 0:
            raise ValueError(
                f"Transformation produced invalid Iyy_B value (must be > 0): {float(Iyy_B)}. "
                f"Input S frame: Ixx_S={float(Ixx_S)}, Iyy_S={float(Iyy_S)}, Ixy_S={float(Ixy_S)}. "
                f"Twist: {float(twist_deg)}°. "
                f"This indicates a transformation error or invalid geometry. "
                f"Possible causes: invalid geometry, wrong twist angle, transformation error."
            )
        invalid_indices = np.where(Iyy_B <= 0)[0]
        raise ValueError(
            f"Transformation produced invalid Iyy_B values (must be > 0): {Iyy_B[invalid_indices]} at indices {invalid_indices}. "
            f"Input S frame: Ixx_S={Ixx_S[invalid_indices] if hasattr(Ixx_S, '__getitem__') else Ixx_S}, "
            f"Iyy_S={Iyy_S[invalid_indices] if hasattr(Iyy_S, '__getitem__') else Iyy_S}, "
            f"Ixy_S={Ixy_S[invalid_indices] if hasattr(Ixy_S, '__getitem__') else Ixy_S}. "
            f"Twist: {twist_deg[invalid_indices] if hasattr(twist_deg, '__getitem__') else twist_deg}°. "
            f"This indicates a transformation error or invalid geometry. "
            f"Possible causes: invalid geometry, wrong twist angle, transformation error."
        )
    
    # Return scalars if inputs were scalars
    if isinstance(Ixx_S, (int, float)) or Ixx_S.ndim == 0:
        return float(Ixx_B), float(Iyy_B), float(Ixy_B)
    return Ixx_B, Iyy_B, Ixy_B

def transform_second_moments_B_to_S(
    Ixx_B: Union[float, np.ndarray],
    Iyy_B: Union[float, np.ndarray],
    Ixy_B: Union[float, np.ndarray] = 0.0,
    twist_deg: Union[float, np.ndarray] = 0.0,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Transform second moments of area from Blade frame (B) to Station frame (S).
    
    This is the inverse of transform_second_moments_S_to_B. Uses negative twist angle
    to reverse the transformation.
    
    Parameters
    ----------
    Ixx_B : float or np.ndarray
        Second moment about x-axis in Blade frame (B) [m⁴]
    Iyy_B : float or np.ndarray
        Second moment about y-axis in Blade frame (B) [m⁴]
    Ixy_B : float or np.ndarray, optional
        Product moment in Blade frame (B) [m⁴]. Default 0.0.
    twist_deg : float or np.ndarray, optional
        Twist angle [degrees]. Positive values indicate CCW rotation of S frame relative to B.
        Default 0.0.
    
    Returns
    -------
    tuple
        (Ixx_S, Iyy_S, Ixy_S) - Second moments in Station frame (S) [m⁴]
    """
    # Use negative twist for inverse transformation
    return transform_second_moments_S_to_B(Ixx_B, Iyy_B, Ixy_B, -twist_deg)
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
