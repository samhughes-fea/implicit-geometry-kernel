"""
Test script to verify the sign convention for moment of inertia transformation.

This script tests whether the transformation should use +θ or -θ when transforming
from Station frame (S) to Blade frame (B).

Test case: Simple rectangle (symmetric, Ixy=0) at 45° twist
- If transformation is correct, Ixx_B and Iyy_B should swap at 45°
- Sum Ixx_B + Iyy_B should be invariant
- Both should remain non-negative
"""

import numpy as np

def transform_second_moments_S_to_B_positive_theta(Ixx_S, Iyy_S, Ixy_S, theta_deg):
    """Current implementation using +θ"""
    theta_rad = np.deg2rad(theta_deg)
    cos_theta = np.cos(theta_rad)
    sin_theta = np.sin(theta_rad)
    cos2 = cos_theta ** 2
    sin2 = sin_theta ** 2
    sin_cos = sin_theta * cos_theta
    
    Ixx_B = Ixx_S * cos2 + Iyy_S * sin2 - 2.0 * Ixy_S * sin_cos
    Iyy_B = Ixx_S * sin2 + Iyy_S * cos2 + 2.0 * Ixy_S * sin_cos
    Ixy_B = (Ixx_S - Iyy_S) * sin_cos + Ixy_S * (cos2 - sin2)
    
    return Ixx_B, Iyy_B, Ixy_B

def transform_second_moments_S_to_B_negative_theta(Ixx_S, Iyy_S, Ixy_S, theta_deg):
    """Alternative using -θ (inverse rotation)"""
    theta_rad = np.deg2rad(-theta_deg)  # Use negative angle
    cos_theta = np.cos(theta_rad)
    sin_theta = np.sin(theta_rad)
    cos2 = cos_theta ** 2
    sin2 = sin_theta ** 2
    sin_cos = sin_theta * cos_theta
    
    Ixx_B = Ixx_S * cos2 + Iyy_S * sin2 - 2.0 * Ixy_S * sin_cos
    Iyy_B = Ixx_S * sin2 + Iyy_S * cos2 + 2.0 * Ixy_S * sin_cos
    Ixy_B = (Ixx_S - Iyy_S) * sin_cos + Ixy_S * (cos2 - sin2)
    
    return Ixx_B, Iyy_B, Ixy_B

# Test case: Rectangle 2m × 1m (width × height)
# In untwisted frame (B): Ixx_B = (1/12) × 2 × 1³ = 0.1667 m⁴, Iyy_B = (1/12) × 1 × 2³ = 0.6667 m⁴
# For a rectangle, Ixy = 0 (symmetric)

# If we rotate coordinates by +45° to create S frame:
# In S frame at 45°: Ixx_S ≈ Iyy_S (should be equal for 45° rotation of rectangle)
# Actually, for a rectangle rotated 45°, Ixx_S = Iyy_S = (Ixx_B + Iyy_B) / 2 = 0.4167 m⁴

# Test with known values
print("=" * 80)
print("Test: Rectangle 2m x 1m rotated by 45 deg")
print("=" * 80)

# Moments in S frame (computed from 45° rotated rectangle)
# For a rectangle, when rotated 45°, Ixx_S = Iyy_S = (Ixx_B + Iyy_B) / 2
Ixx_S = 0.4167  # Approximate
Iyy_S = 0.4167  # Approximate
Ixy_S = 0.0  # Rectangle is symmetric

theta_deg = 45.0

print(f"\nInput (S frame): Ixx_S={Ixx_S:.6f}, Iyy_S={Iyy_S:.6f}, Ixy_S={Ixy_S:.6f}")
print(f"Twist angle: {theta_deg} deg")

# Test with +θ (current implementation)
Ixx_B_pos, Iyy_B_pos, Ixy_B_pos = transform_second_moments_S_to_B_positive_theta(
    Ixx_S, Iyy_S, Ixy_S, theta_deg
)

# Test with -θ (alternative)
Ixx_B_neg, Iyy_B_neg, Ixy_B_neg = transform_second_moments_S_to_B_negative_theta(
    Ixx_S, Iyy_S, Ixy_S, theta_deg
)

print(f"\nUsing +theta (current):")
print(f"  Ixx_B = {Ixx_B_pos:.6f}")
print(f"  Iyy_B = {Iyy_B_pos:.6f}")
print(f"  Ixy_B = {Ixy_B_pos:.6f}")
print(f"  Sum = {Ixx_B_pos + Iyy_B_pos:.6f}")
print(f"  Both non-negative? {Ixx_B_pos >= 0 and Iyy_B_pos >= 0}")

print(f"\nUsing -theta (alternative):")
print(f"  Ixx_B = {Ixx_B_neg:.6f}")
print(f"  Iyy_B = {Iyy_B_neg:.6f}")
print(f"  Ixy_B = {Ixy_B_neg:.6f}")
print(f"  Sum = {Ixx_B_neg + Iyy_B_neg:.6f}")
print(f"  Both non-negative? {Ixx_B_neg >= 0 and Iyy_B_neg >= 0}")

# Expected: For 45° rotation, Ixx_B should equal Iyy_B, and both should be positive
# The correct transformation should give Ixx_B ≈ Iyy_B ≈ 0.4167

print(f"\nExpected (for 45 deg rotation back): Ixx_B ~= Iyy_B ~= 0.4167")
print(f"  +theta gives: Ixx_B={Ixx_B_pos:.6f}, Iyy_B={Iyy_B_pos:.6f} (match? {abs(Ixx_B_pos - 0.4167) < 0.01 and abs(Iyy_B_pos - 0.4167) < 0.01})")
print(f"  -theta gives: Ixx_B={Ixx_B_neg:.6f}, Iyy_B={Iyy_B_neg:.6f} (match? {abs(Ixx_B_neg - 0.4167) < 0.01 and abs(Iyy_B_neg - 0.4167) < 0.01})")

# Test with 0° (should give identity)
print("\n" + "=" * 80)
print("Test: 0 deg twist (should be identity)")
print("=" * 80)

Ixx_S_0 = 0.1667
Iyy_S_0 = 0.6667
Ixy_S_0 = 0.0

Ixx_B_pos_0, Iyy_B_pos_0, _ = transform_second_moments_S_to_B_positive_theta(
    Ixx_S_0, Iyy_S_0, Ixy_S_0, 0.0
)
Ixx_B_neg_0, Iyy_B_neg_0, _ = transform_second_moments_S_to_B_negative_theta(
    Ixx_S_0, Iyy_S_0, Ixy_S_0, 0.0
)

print(f"Input: Ixx_S={Ixx_S_0:.6f}, Iyy_S={Iyy_S_0:.6f}")
print(f"+theta: Ixx_B={Ixx_B_pos_0:.6f}, Iyy_B={Iyy_B_pos_0:.6f}")
print(f"-theta: Ixx_B={Ixx_B_neg_0:.6f}, Iyy_B={Iyy_B_neg_0:.6f}")
print(f"Both should equal input at 0 deg: {abs(Ixx_B_pos_0 - Ixx_S_0) < 1e-6 and abs(Iyy_B_pos_0 - Iyy_S_0) < 1e-6}")


