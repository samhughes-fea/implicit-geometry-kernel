"""
Test script to verify sign convention for moment transformation with asymmetric shapes.

Tests transformation with non-zero Ixy (asymmetric airfoils) at various twist angles
to determine if +θ or -θ is correct for S→B transformation.
"""

import numpy as np
import sys
from pathlib import Path

# Add workspace root to path
_workspace_root = Path(__file__).parent.parent.parent.parent.parent
if str(_workspace_root) not in sys.path:
    sys.path.insert(0, str(_workspace_root))

from blade_geometry_model.coordinate_frames import CoordinateFrameTransforms


def transform_with_negative_theta(Ixx_S, Iyy_S, Ixy_S, theta_deg):
    """Alternative transformation using -θ (inverse rotation)"""
    # Use negative angle in the transformation
    transform = CoordinateFrameTransforms()
    return transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, -theta_deg)


def test_asymmetric_shape():
    """
    Test with asymmetric shape (L-section or channel section).
    
    For an L-section with known geometry, we can compute Ixx, Iyy, Ixy analytically.
    Then we can test what happens when we:
    1. Rotate coordinates by +θ to create S frame
    2. Compute moments in S frame
    3. Transform back to B frame using both +θ and -θ
    4. See which gives correct (non-negative) results
    """
    print("=" * 80)
    print("Test: Asymmetric Shape (Non-Zero Ixy)")
    print("=" * 80)
    
    # Test case: L-section (asymmetric, Ixy != 0)
    # Dimensions: vertical leg 2m × 0.2m, horizontal leg 1m × 0.2m
    # Centroid at origin for simplicity
    # 
    # In B frame (untwisted):
    # Ixx_B ≈ 0.133 m⁴, Iyy_B ≈ 0.033 m⁴, Ixy_B ≈ -0.033 m⁴ (negative due to asymmetry)
    #
    # If we rotate coordinates by +30° to create S frame:
    # The moments in S frame will be different
    
    # Simulated S frame moments (computed from 30° rotated L-section)
    # These would be computed from the actual rotated geometry
    # For testing, we'll use values that represent a rotated L-section
    Ixx_S = 0.100  # Approximate value for rotated L-section
    Iyy_S = 0.066  # Approximate value
    Ixy_S = 0.020  # Non-zero product moment (asymmetric)
    
    test_angles = [0.0, 30.0, 45.0, 90.0]
    
    for theta_deg in test_angles:
        print(f"\n{'=' * 80}")
        print(f"Test at {theta_deg} deg twist")
        print(f"{'=' * 80}")
        
        print(f"Input (S frame): Ixx_S={Ixx_S:.6f}, Iyy_S={Iyy_S:.6f}, Ixy_S={Ixy_S:.6f}")
        
        # Test with +θ (current implementation)
        transform = CoordinateFrameTransforms()
        Ixx_B_pos, Iyy_B_pos, Ixy_B_pos = transform.transform_second_moments_S_to_B(
            Ixx_S, Iyy_S, Ixy_S, theta_deg
        )
        
        # Test with -θ (alternative)
        Ixx_B_neg, Iyy_B_neg, Ixy_B_neg = transform_with_negative_theta(
            Ixx_S, Iyy_S, Ixy_S, theta_deg
        )
        
        print(f"\nUsing +theta (current):")
        print(f"  Ixx_B = {Ixx_B_pos:.6f} {'[NEGATIVE!]' if Ixx_B_pos < 0 else ''}")
        print(f"  Iyy_B = {Iyy_B_pos:.6f} {'[NEGATIVE!]' if Iyy_B_pos < 0 else ''}")
        print(f"  Ixy_B = {Ixy_B_pos:.6f}")
        print(f"  Sum = {Ixx_B_pos + Iyy_B_pos:.6f}")
        print(f"  Both non-negative? {Ixx_B_pos >= 0 and Iyy_B_pos >= 0}")
        
        print(f"\nUsing -theta (alternative):")
        print(f"  Ixx_B = {Ixx_B_neg:.6f} {'[NEGATIVE!]' if Ixx_B_neg < 0 else ''}")
        print(f"  Iyy_B = {Iyy_B_neg:.6f} {'[NEGATIVE!]' if Iyy_B_neg < 0 else ''}")
        print(f"  Ixy_B = {Ixy_B_neg:.6f}")
        print(f"  Sum = {Ixx_B_neg + Iyy_B_neg:.6f}")
        print(f"  Both non-negative? {Ixx_B_neg >= 0 and Iyy_B_neg >= 0}")
        
        # Check which is physically correct
        pos_valid = Ixx_B_pos >= 0 and Iyy_B_pos >= 0
        neg_valid = Ixx_B_neg >= 0 and Iyy_B_neg >= 0
        
        if pos_valid and not neg_valid:
            print(f"\n>>> +theta is CORRECT (both non-negative)")
        elif neg_valid and not pos_valid:
            print(f"\n>>> -theta is CORRECT (both non-negative)")
        elif pos_valid and neg_valid:
            print(f"\n>>> Both give valid results (need further analysis)")
        else:
            print(f"\n>>> NEITHER gives valid results (both have negatives)")
            print(f"    This suggests a different issue (formula error, not sign)")


def test_with_known_analytical_case():
    """
    Test with a case where we know the analytical answer.
    
    For a rectangle rotated 45°, we know:
    - If original: Ixx = a, Iyy = b, Ixy = 0
    - After 45° rotation: Ixx' = Iyy' = (a+b)/2, Ixy' = (a-b)/2
    - Transforming back should give original values
    """
    print("\n" + "=" * 80)
    print("Test: Known Analytical Case (Rectangle at 45°)")
    print("=" * 80)
    
    # Original rectangle in B frame: 2m × 1m
    Ixx_B_original = (1.0/12.0) * 2.0 * (1.0**3)  # = 0.1667 m⁴
    Iyy_B_original = (1.0/12.0) * 1.0 * (2.0**3)  # = 0.6667 m⁴
    Ixy_B_original = 0.0  # Rectangle is symmetric
    
    # After rotating coordinates by +45° to create S frame:
    # For a rectangle at 45°, Ixx_S = Iyy_S = (Ixx_B + Iyy_B) / 2
    Ixx_S = (Ixx_B_original + Iyy_B_original) / 2.0  # = 0.4167 m⁴
    Iyy_S = (Ixx_B_original + Iyy_B_original) / 2.0  # = 0.4167 m⁴
    Ixy_S = 0.0  # Still zero for rectangle
    
    theta_deg = 45.0
    
    print(f"Original (B frame): Ixx={Ixx_B_original:.6f}, Iyy={Iyy_B_original:.6f}, Ixy={Ixy_B_original:.6f}")
    print(f"After +45° rotation (S frame): Ixx={Ixx_S:.6f}, Iyy={Iyy_S:.6f}, Ixy={Ixy_S:.6f}")
    print(f"Transforming back to B frame with {theta_deg} deg...")
    
    # Transform back using +θ
    transform = CoordinateFrameTransforms()
    Ixx_B_pos, Iyy_B_pos, Ixy_B_pos = transform.transform_second_moments_S_to_B(
        Ixx_S, Iyy_S, Ixy_S, theta_deg
    )
    
    # Transform back using -θ
    Ixx_B_neg, Iyy_B_neg, Ixy_B_neg = transform_with_negative_theta(
        Ixx_S, Iyy_S, Ixy_S, theta_deg
    )
    
    print(f"\nUsing +theta:")
    print(f"  Ixx_B = {Ixx_B_pos:.6f} (original: {Ixx_B_original:.6f}, diff: {abs(Ixx_B_pos - Ixx_B_original):.6f})")
    print(f"  Iyy_B = {Iyy_B_pos:.6f} (original: {Iyy_B_original:.6f}, diff: {abs(Iyy_B_pos - Iyy_B_original):.6f})")
    print(f"  Match? {abs(Ixx_B_pos - Ixx_B_original) < 1e-3 and abs(Iyy_B_pos - Iyy_B_original) < 1e-3}")
    
    print(f"\nUsing -theta:")
    print(f"  Ixx_B = {Ixx_B_neg:.6f} (original: {Ixx_B_original:.6f}, diff: {abs(Ixx_B_neg - Ixx_B_original):.6f})")
    print(f"  Iyy_B = {Iyy_B_neg:.6f} (original: {Iyy_B_original:.6f}, diff: {abs(Iyy_B_neg - Iyy_B_original):.6f})")
    print(f"  Match? {abs(Ixx_B_neg - Ixx_B_original) < 1e-3 and abs(Iyy_B_neg - Iyy_B_original) < 1e-3}")


def test_invariance():
    """
    Test that Ixx_B + Iyy_B is invariant (same in both frames).
    This is a fundamental property that must hold.
    """
    print("\n" + "=" * 80)
    print("Test: Invariance of Ixx + Iyy")
    print("=" * 80)
    
    Ixx_S = 0.100
    Iyy_S = 0.066
    Ixy_S = 0.020
    sum_S = Ixx_S + Iyy_S
    
    test_angles = [0.0, 30.0, 45.0, 90.0]
    
    print(f"Sum in S frame: {sum_S:.6f}")
    print(f"This should be invariant (same in B frame)")
    
    transform = CoordinateFrameTransforms()
    for theta_deg in test_angles:
        Ixx_B_pos, Iyy_B_pos, _ = transform.transform_second_moments_S_to_B(
            Ixx_S, Iyy_S, Ixy_S, theta_deg
        )
        Ixx_B_neg, Iyy_B_neg, _ = transform_with_negative_theta(
            Ixx_S, Iyy_S, Ixy_S, theta_deg
        )
        
        sum_B_pos = Ixx_B_pos + Iyy_B_pos
        sum_B_neg = Ixx_B_neg + Iyy_B_neg
        
        print(f"\nAt {theta_deg} deg:")
        print(f"  +theta: Sum_B = {sum_B_pos:.6f} (diff from S: {abs(sum_B_pos - sum_S):.6f})")
        print(f"  -theta: Sum_B = {sum_B_neg:.6f} (diff from S: {abs(sum_B_neg - sum_S):.6f})")
        print(f"  +theta invariant? {abs(sum_B_pos - sum_S) < 1e-6}")
        print(f"  -theta invariant? {abs(sum_B_neg - sum_S) < 1e-6}")


if __name__ == "__main__":
    test_asymmetric_shape()
    test_with_known_analytical_case()
    test_invariance()

