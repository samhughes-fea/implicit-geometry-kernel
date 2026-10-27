"""
Verify consistency between coordinate rotation matrix and moment transformation.

Key question: If coordinates are rotated by R(θ), should moments transform by R(-θ) or R(θ)?
"""

import numpy as np


def rotation_matrix(theta_rad):
    """Standard 2D rotation matrix (CCW rotation)"""
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    return np.array([[cos_t, -sin_t], [sin_t, cos_t]])


def moment_tensor(Ixx, Iyy, Ixy):
    """Create moment of inertia tensor"""
    return np.array([[Ixx, Ixy], [Ixy, Iyy]])


def transform_moments_matrix(Ixx_S, Iyy_S, Ixy_S, theta_deg):
    """
    Transform moments using matrix approach: I_B = R(-θ) I_S R(-θ)^T
    
    If coordinates are rotated by R(θ) to create S frame, then to transform
    moments FROM S TO B, we use R(-θ) (inverse rotation).
    """
    theta_rad = np.deg2rad(theta_deg)
    R = rotation_matrix(-theta_rad)  # Inverse rotation
    I_S = moment_tensor(Ixx_S, Iyy_S, Ixy_S)
    I_B = R @ I_S @ R.T
    return I_B[0, 0], I_B[1, 1], I_B[0, 1]  # Ixx_B, Iyy_B, Ixy_B


def transform_moments_formula_positive(Ixx_S, Iyy_S, Ixy_S, theta_deg):
    """Current formula-based approach using +θ"""
    theta_rad = np.deg2rad(theta_deg)
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    cos2 = cos_t ** 2
    sin2 = sin_t ** 2
    sin_cos = sin_t * cos_t
    
    Ixx_B = Ixx_S * cos2 + Iyy_S * sin2 - 2.0 * Ixy_S * sin_cos
    Iyy_B = Ixx_S * sin2 + Iyy_S * cos2 + 2.0 * Ixy_S * sin_cos
    Ixy_B = (Ixx_S - Iyy_S) * sin_cos + Ixy_S * (cos2 - sin2)
    
    return Ixx_B, Iyy_B, Ixy_B


def transform_moments_formula_negative(Ixx_S, Iyy_S, Ixy_S, theta_deg):
    """Formula-based approach using -θ"""
    theta_rad = np.deg2rad(-theta_deg)
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    cos2 = cos_t ** 2
    sin2 = sin_t ** 2
    sin_cos = sin_t * cos_t
    
    Ixx_B = Ixx_S * cos2 + Iyy_S * sin2 - 2.0 * Ixy_S * sin_cos
    Iyy_B = Ixx_S * sin2 + Iyy_S * cos2 + 2.0 * Ixy_S * sin_cos
    Ixy_B = (Ixx_S - Iyy_S) * sin_cos + Ixy_S * (cos2 - sin2)
    
    return Ixx_B, Iyy_B, Ixy_B


print("=" * 80)
print("Verification: Matrix vs Formula Approach")
print("=" * 80)

# Test case with asymmetric shape
Ixx_S = 0.100
Iyy_S = 0.066
Ixy_S = 0.020
theta_deg = 30.0

print(f"\nInput (S frame): Ixx_S={Ixx_S:.6f}, Iyy_S={Iyy_S:.6f}, Ixy_S={Ixy_S:.6f}")
print(f"Twist angle: {theta_deg} deg")

# Matrix approach (using R(-θ))
Ixx_B_mat, Iyy_B_mat, Ixy_B_mat = transform_moments_matrix(Ixx_S, Iyy_S, Ixy_S, theta_deg)

# Formula with +θ
Ixx_B_pos, Iyy_B_pos, Ixy_B_pos = transform_moments_formula_positive(Ixx_S, Iyy_S, Ixy_S, theta_deg)

# Formula with -θ
Ixx_B_neg, Iyy_B_neg, Ixy_B_neg = transform_moments_formula_negative(Ixx_S, Iyy_S, Ixy_S, theta_deg)

print(f"\nMatrix approach (R(-theta)):")
print(f"  Ixx_B = {Ixx_B_mat:.6f}")
print(f"  Iyy_B = {Iyy_B_mat:.6f}")
print(f"  Ixy_B = {Ixy_B_mat:.6f}")

print(f"\nFormula with +theta:")
print(f"  Ixx_B = {Ixx_B_pos:.6f}")
print(f"  Iyy_B = {Iyy_B_pos:.6f}")
print(f"  Ixy_B = {Ixy_B_pos:.6f}")
print(f"  Matches matrix? {abs(Ixx_B_pos - Ixx_B_mat) < 1e-6 and abs(Iyy_B_pos - Iyy_B_mat) < 1e-6}")

print(f"\nFormula with -theta:")
print(f"  Ixx_B = {Ixx_B_neg:.6f}")
print(f"  Iyy_B = {Iyy_B_neg:.6f}")
print(f"  Ixy_B = {Ixy_B_neg:.6f}")
print(f"  Matches matrix? {abs(Ixx_B_neg - Ixx_B_mat) < 1e-6 and abs(Iyy_B_neg - Iyy_B_mat) < 1e-6}")

print(f"\n>>> Matrix approach uses R(-theta), so formula should use -theta for consistency")
print(f">>> If formula with -theta matches matrix, then current implementation (+theta) is WRONG")


