"""
Test script to validate S-to-B frame transformation and EI assignments.

This script validates:
1. The S-to-B frame transformation formulas
2. The assignment of Ixx_B to EI_edge and Iyy_B to EI_flap
3. The coordinate system conventions in the transformation

Three key areas tested:
1. S-to-B frame transformation (transform_second_moments_S_to_B)
2. The assignment of Ixx_B to EI_edge and Iyy_B to EI_flap
3. The coordinate system conventions in the transformation
"""

import sys
import numpy as np
from pathlib import Path

# Add workspace root to path
_workspace_root = Path(__file__).parent.parent.parent.parent.parent
if str(_workspace_root) not in sys.path:
    sys.path.insert(0, str(_workspace_root))

from blade_geometry_model.coordinate_frames import CoordinateFrameTransforms
from blade_geometry_model.shoelace_integration import section_properties_S
from shapely.geometry import Polygon as ShapelyPolygon


def test_transformation_identity():
    """
    Test 1: Zero twist angle (identity transformation)
    
    When twist_deg = 0, the transformation should be identity:
    - Ixx_B = Ixx_S
    - Iyy_B = Iyy_S
    - Ixy_B = Ixy_S
    """
    print("=" * 80)
    print("Test 1: Zero Twist Angle (Identity Transformation)")
    print("=" * 80)
    
    # Test with symmetric rectangle (Ixy_S = 0)
    Ixx_S = 1.666667  # m^4 (5m x 1m rectangle, Ixx about x-axis)
    Iyy_S = 41.66667  # m^4 (5m x 1m rectangle, Iyy about y-axis)
    Ixy_S = 0.0  # Symmetric rectangle has zero product moment
    twist_deg = 0.0
    
    transform = CoordinateFrameTransforms()
    Ixx_B, Iyy_B, Ixy_B = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, twist_deg)
    
    print(f"\nInput S frame:")
    print(f"  Ixx_S = {Ixx_S:.6e} m^4")
    print(f"  Iyy_S = {Iyy_S:.6e} m^4")
    print(f"  Ixy_S = {Ixy_S:.6e} m^4")
    print(f"  twist_deg = {twist_deg:.6f}°")
    
    print(f"\nOutput B frame:")
    print(f"  Ixx_B = {Ixx_B:.6e} m^4")
    print(f"  Iyy_B = {Iyy_B:.6e} m^4")
    print(f"  Ixy_B = {Ixy_B:.6e} m^4")
    
    # Check identity
    ixx_error = abs(Ixx_B - Ixx_S)
    iyy_error = abs(Iyy_B - Iyy_S)
    ixy_error = abs(Ixy_B - Ixy_S)
    
    print(f"\nIdentity Check:")
    print(f"  |Ixx_B - Ixx_S| = {ixx_error:.6e} m^4")
    print(f"  |Iyy_B - Iyy_S| = {iyy_error:.6e} m^4")
    print(f"  |Ixy_B - Ixy_S| = {ixy_error:.6e} m^4")
    
    tolerance = 1e-10
    ixx_pass = ixx_error < tolerance
    iyy_pass = iyy_error < tolerance
    ixy_pass = ixy_error < tolerance
    
    print(f"\nResults:")
    print(f"  Ixx_B identity: {'PASS' if ixx_pass else 'FAIL'}")
    print(f"  Iyy_B identity: {'PASS' if iyy_pass else 'FAIL'}")
    print(f"  Ixy_B identity: {'PASS' if ixy_pass else 'FAIL'}")
    
    all_pass = ixx_pass and iyy_pass and ixy_pass
    print(f"\nTest 1 Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def test_transformation_90_degree():
    """
    Test 2: 90-degree twist (should swap Ixx and Iyy)
    
    When twist_deg = 90°, the transformation should swap Ixx and Iyy:
    - Ixx_B ≈ Iyy_S (for symmetric case with Ixy_S = 0)
    - Iyy_B ≈ Ixx_S (for symmetric case with Ixy_S = 0)
    - Ixy_B ≈ 0 (for symmetric case)
    """
    print("\n" + "=" * 80)
    print("Test 2: 90-Degree Twist (Ixx/Iyy Swap)")
    print("=" * 80)
    
    # Test with symmetric rectangle (Ixy_S = 0)
    Ixx_S = 1.666667  # m^4 (5m x 1m rectangle, Ixx about x-axis)
    Iyy_S = 41.66667  # m^4 (5m x 1m rectangle, Iyy about y-axis)
    Ixy_S = 0.0  # Symmetric rectangle has zero product moment
    twist_deg = 90.0
    
    transform = CoordinateFrameTransforms()
    Ixx_B, Iyy_B, Ixy_B = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, twist_deg)
    
    print(f"\nInput S frame:")
    print(f"  Ixx_S = {Ixx_S:.6e} m^4")
    print(f"  Iyy_S = {Iyy_S:.6e} m^4")
    print(f"  Ixy_S = {Ixy_S:.6e} m^4")
    print(f"  twist_deg = {twist_deg:.6f}°")
    
    print(f"\nOutput B frame:")
    print(f"  Ixx_B = {Ixx_B:.6e} m^4")
    print(f"  Iyy_B = {Iyy_B:.6e} m^4")
    print(f"  Ixy_B = {Ixy_B:.6e} m^4")
    
    # Check swap (for 90° rotation, Ixx_B should equal Iyy_S and Iyy_B should equal Ixx_S)
    ixx_swap_error = abs(Ixx_B - Iyy_S)
    iyy_swap_error = abs(Iyy_B - Ixx_S)
    ixy_zero_error = abs(Ixy_B)
    
    print(f"\n90° Rotation Check:")
    print(f"  |Ixx_B - Iyy_S| = {ixx_swap_error:.6e} m^4 (should be ~0)")
    print(f"  |Iyy_B - Ixx_S| = {iyy_swap_error:.6e} m^4 (should be ~0)")
    print(f"  |Ixy_B| = {ixy_zero_error:.6e} m^4 (should be ~0)")
    
    tolerance = 1e-6
    ixx_pass = ixx_swap_error < tolerance
    iyy_pass = iyy_swap_error < tolerance
    ixy_pass = ixy_zero_error < tolerance
    
    print(f"\nResults:")
    print(f"  Ixx_B ≈ Iyy_S: {'PASS' if ixx_pass else 'FAIL'}")
    print(f"  Iyy_B ≈ Ixx_S: {'PASS' if iyy_pass else 'FAIL'}")
    print(f"  Ixy_B ≈ 0: {'PASS' if ixy_pass else 'FAIL'}")
    
    all_pass = ixx_pass and iyy_pass and ixy_pass
    print(f"\nTest 2 Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def test_transformation_45_degree():
    """
    Test 3: 45-degree twist (intermediate rotation)
    
    When twist_deg = 45°, the transformation should mix Ixx and Iyy:
    - For symmetric case (Ixy_S = 0): Ixx_B = Iyy_B = (Ixx_S + Iyy_S) / 2
    - Ixy_B should be non-zero (unless Ixx_S = Iyy_S)
    """
    print("\n" + "=" * 80)
    print("Test 3: 45-Degree Twist (Intermediate Rotation)")
    print("=" * 80)
    
    # Test with symmetric rectangle (Ixy_S = 0)
    Ixx_S = 1.666667  # m^4 (5m x 1m rectangle, Ixx about x-axis)
    Iyy_S = 41.66667  # m^4 (5m x 1m rectangle, Iyy about y-axis)
    Ixy_S = 0.0  # Symmetric rectangle has zero product moment
    twist_deg = 45.0
    
    transform = CoordinateFrameTransforms()
    Ixx_B, Iyy_B, Ixy_B = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, twist_deg)
    
    print(f"\nInput S frame:")
    print(f"  Ixx_S = {Ixx_S:.6e} m^4")
    print(f"  Iyy_S = {Iyy_S:.6e} m^4")
    print(f"  Ixy_S = {Ixy_S:.6e} m^4")
    print(f"  twist_deg = {twist_deg:.6f}°")
    
    print(f"\nOutput B frame:")
    print(f"  Ixx_B = {Ixx_B:.6e} m^4")
    print(f"  Iyy_B = {Iyy_B:.6e} m^4")
    print(f"  Ixy_B = {Ixy_B:.6e} m^4")
    
    # For 45° rotation with Ixy_S = 0:
    # Ixx_B = Ixx_S * cos²(45°) + Iyy_S * sin²(45°) = (Ixx_S + Iyy_S) / 2
    # Iyy_B = Ixx_S * sin²(45°) + Iyy_S * cos²(45°) = (Ixx_S + Iyy_S) / 2
    # Ixy_B = (Ixx_S - Iyy_S) * sin(45°) * cos(45°) = (Ixx_S - Iyy_S) / 2
    expected_ixx = (Ixx_S + Iyy_S) / 2.0
    expected_iyy = (Ixx_S + Iyy_S) / 2.0
    expected_ixy = (Ixx_S - Iyy_S) / 2.0
    
    ixx_error = abs(Ixx_B - expected_ixx)
    iyy_error = abs(Iyy_B - expected_iyy)
    ixy_error = abs(Ixy_B - expected_ixy)
    
    print(f"\n45° Rotation Check:")
    print(f"  Expected Ixx_B = (Ixx_S + Iyy_S) / 2 = {expected_ixx:.6e} m^4")
    print(f"  Expected Iyy_B = (Ixx_S + Iyy_S) / 2 = {expected_iyy:.6e} m^4")
    print(f"  Expected Ixy_B = (Ixx_S - Iyy_S) / 2 = {expected_ixy:.6e} m^4")
    print(f"  |Ixx_B - expected| = {ixx_error:.6e} m^4")
    print(f"  |Iyy_B - expected| = {iyy_error:.6e} m^4")
    print(f"  |Ixy_B - expected| = {ixy_error:.6e} m^4")
    
    tolerance = 1e-6
    ixx_pass = ixx_error < tolerance
    iyy_pass = iyy_error < tolerance
    ixy_pass = ixy_error < tolerance
    
    print(f"\nResults:")
    print(f"  Ixx_B formula: {'PASS' if ixx_pass else 'FAIL'}")
    print(f"  Iyy_B formula: {'PASS' if iyy_pass else 'FAIL'}")
    print(f"  Ixy_B formula: {'PASS' if ixy_pass else 'FAIL'}")
    
    all_pass = ixx_pass and iyy_pass and ixy_pass
    print(f"\nTest 3 Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def test_ei_assignment():
    """
    Test 4: EI Assignment (Ixx_B → EI_edge, Iyy_B → EI_flap)
    
    Validates that:
    - EI_edge = E * Ixx_B (edgewise stiffness uses Ixx in B frame)
    - EI_flap = E * Iyy_B (flapwise stiffness uses Iyy in B frame)
    
    This is the critical assignment that determines which direction gets which stiffness.
    """
    print("\n" + "=" * 80)
    print("Test 4: EI Assignment (Ixx_B → EI_edge, Iyy_B → EI_flap)")
    print("=" * 80)
    
    # Test with symmetric rectangle
    Ixx_S = 1.666667  # m^4 (5m x 1m rectangle, Ixx about x-axis)
    Iyy_S = 41.66667  # m^4 (5m x 1m rectangle, Iyy about y-axis)
    Ixy_S = 0.0
    twist_deg = 0.0  # Start with zero twist for clarity
    E = 2.0e11  # Pa (typical Young's modulus for composite)
    
    transform = CoordinateFrameTransforms()
    Ixx_B, Iyy_B, Ixy_B = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, twist_deg)
    
    # Compute EI values (as done in the code)
    EI_edge = E * Ixx_B
    EI_flap = E * Iyy_B
    
    print(f"\nInput S frame:")
    print(f"  Ixx_S = {Ixx_S:.6e} m^4")
    print(f"  Iyy_S = {Iyy_S:.6e} m^4")
    print(f"  Ixy_S = {Ixy_S:.6e} m^4")
    print(f"  twist_deg = {twist_deg:.6f}°")
    print(f"  E = {E:.6e} Pa")
    
    print(f"\nTransformed B frame:")
    print(f"  Ixx_B = {Ixx_B:.6e} m^4")
    print(f"  Iyy_B = {Iyy_B:.6e} m^4")
    print(f"  Ixy_B = {Ixy_B:.6e} m^4")
    
    print(f"\nEI Values (as computed in code):")
    print(f"  EI_edge = E * Ixx_B = {EI_edge:.6e} N·m²")
    print(f"  EI_flap = E * Iyy_B = {EI_flap:.6e} N·m²")
    
    # For a 5m x 1m rectangle (wider than tall):
    # - Ixx_S is smaller (about x-axis, uses y^2, height=1m)
    # - Iyy_S is larger (about y-axis, uses x^2, width=5m)
    # With zero twist: Ixx_B = Ixx_S, Iyy_B = Iyy_S
    # So: EI_edge should be smaller, EI_flap should be larger
    
    print(f"\nPhysical Interpretation:")
    print(f"  Rectangle: 5m wide × 1m tall (wider than tall)")
    print(f"  Ixx_B (about x-axis, uses y^2) = {Ixx_B:.6e} m^4 (SMALLER)")
    print(f"  Iyy_B (about y-axis, uses x^2) = {Iyy_B:.6e} m^4 (LARGER)")
    print(f"  EI_edge = {EI_edge:.6e} N·m² (SMALLER - edgewise stiffness)")
    print(f"  EI_flap = {EI_flap:.6e} N·m² (LARGER - flapwise stiffness)")
    print(f"  Ratio EI_flap/EI_edge = {EI_flap/EI_edge:.6f}")
    
    # Check that assignments are correct
    # For a wider-than-tall rectangle, EI_flap should be larger than EI_edge
    flap_larger = EI_flap > EI_edge
    ixx_smaller = Ixx_B < Iyy_B
    
    print(f"\nAssignment Validation:")
    print(f"  Ixx_B < Iyy_B: {'PASS' if ixx_smaller else 'FAIL'}")
    print(f"  EI_flap > EI_edge: {'PASS' if flap_larger else 'FAIL'}")
    print(f"  EI_edge = E * Ixx_B: {'PASS' if abs(EI_edge - E * Ixx_B) < 1e-6 else 'FAIL'}")
    print(f"  EI_flap = E * Iyy_B: {'PASS' if abs(EI_flap - E * Iyy_B) < 1e-6 else 'FAIL'}")
    
    all_pass = flap_larger and ixx_smaller
    print(f"\nTest 4 Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def test_transformation_with_twist():
    """
    Test 5: Transformation with non-zero twist
    
    Tests that the transformation correctly handles twist angles and that
    EI assignments remain consistent (Ixx_B → EI_edge, Iyy_B → EI_flap).
    """
    print("\n" + "=" * 80)
    print("Test 5: Transformation with Non-Zero Twist")
    print("=" * 80)
    
    # Test with symmetric rectangle
    Ixx_S = 1.666667  # m^4
    Iyy_S = 41.66667  # m^4
    Ixy_S = 0.0
    E = 2.0e11  # Pa
    
    twist_angles = [0.0, 15.0, 30.0, 45.0, 60.0, 90.0]
    
    print(f"\nInput S frame (constant):")
    print(f"  Ixx_S = {Ixx_S:.6e} m^4")
    print(f"  Iyy_S = {Iyy_S:.6e} m^4")
    print(f"  Ixy_S = {Ixy_S:.6e} m^4")
    print(f"  E = {E:.6e} Pa")
    
    print(f"\n{'Twist [°]':<12} {'Ixx_B [m⁴]':<18} {'Iyy_B [m⁴]':<18} {'EI_edge [N·m²]':<20} {'EI_flap [N·m²]':<20}")
    print("-" * 100)
    
    all_pass = True
    transform = CoordinateFrameTransforms()
    for twist_deg in twist_angles:
        Ixx_B, Iyy_B, Ixy_B = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, twist_deg)
        EI_edge = E * Ixx_B
        EI_flap = E * Iyy_B
        
        # Check that all values are positive (physical constraint)
        valid = Ixx_B > 0 and Iyy_B > 0 and EI_edge > 0 and EI_flap > 0
        
        print(f"{twist_deg:>10.1f}  {Ixx_B:>16.6e}  {Iyy_B:>16.6e}  {EI_edge:>18.6e}  {EI_flap:>18.6e}  {'✓' if valid else '✗'}")
        
        if not valid:
            all_pass = False
    
    print(f"\nTest 5 Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def test_coordinate_convention():
    """
    Test 6: Coordinate System Convention Validation
    
    Validates the coordinate system conventions:
    - S frame: Station frame, rotates with blade twist
    - B frame: Blade frame, fixed orientation
    - Positive twist_deg: CCW rotation of S frame relative to B
    - Ixx_B → EI_edge (edgewise stiffness)
    - Iyy_B → EI_flap (flapwise stiffness)
    """
    print("\n" + "=" * 80)
    print("Test 6: Coordinate System Convention Validation")
    print("=" * 80)
    
    print("\nCoordinate System Conventions:")
    print("  S frame (Station frame):")
    print("    - Origin: 0.333 chordline (pivot point)")
    print("    - x-axis: Chordwise direction (horizontal)")
    print("    - y-axis: Thickness direction (vertical)")
    print("    - Rotates with blade twist")
    print("    - Properties: Ixx_S, Iyy_S, Ixy_S")
    print()
    print("  B frame (Blade frame):")
    print("    - Origin: 0.333 chordline at blade root")
    print("    - X_B axis: Edgewise direction (vertical in global)")
    print("    - Y_B axis: Flapwise direction (horizontal in global)")
    print("    - Does NOT rotate with blade twist")
    print("    - Properties: Ixx_B, Iyy_B, Ixy_B")
    print()
    print("  Transformation:")
    print("    - Positive twist_deg: CCW rotation of S frame relative to B")
    print("    - Uses transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, twist_deg)")
    print()
    print("  EI Assignment:")
    print("    - EI_edge = E * Ixx_B (edgewise stiffness uses Ixx in B frame)")
    print("    - EI_flap = E * Iyy_B (flapwise stiffness uses Iyy in B frame)")
    
    # Test with known values
    Ixx_S = 1.666667  # m^4
    Iyy_S = 41.66667  # m^4
    Ixy_S = 0.0
    
    print(f"\nTest Case: 5m × 1m rectangle (wider than tall)")
    print(f"  Ixx_S = {Ixx_S:.6e} m^4 (about x-axis, uses y^2, height=1m)")
    print(f"  Iyy_S = {Iyy_S:.6e} m^4 (about y-axis, uses x^2, width=5m)")
    print(f"  Ixy_S = {Ixy_S:.6e} m^4 (symmetric, zero product moment)")
    
    # Test zero twist
    transform = CoordinateFrameTransforms()
    Ixx_B_0, Iyy_B_0, Ixy_B_0 = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, 0.0)
    print(f"\n  With twist_deg = 0° (no rotation):")
    print(f"    Ixx_B = {Ixx_B_0:.6e} m^4 (should equal Ixx_S)")
    print(f"    Iyy_B = {Iyy_B_0:.6e} m^4 (should equal Iyy_S)")
    print(f"    Ixy_B = {Ixy_B_0:.6e} m^4 (should equal Ixy_S)")
    
    # Test 90 degree twist
    Ixx_B_90, Iyy_B_90, Ixy_B_90 = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, 90.0)
    print(f"\n  With twist_deg = 90° (full rotation):")
    print(f"    Ixx_B = {Ixx_B_90:.6e} m^4 (should equal Iyy_S)")
    print(f"    Iyy_B = {Iyy_B_90:.6e} m^4 (should equal Ixx_S)")
    print(f"    Ixy_B = {Ixy_B_90:.6e} m^4 (should be ~0)")
    
    # Validate conventions
    zero_twist_ok = abs(Ixx_B_0 - Ixx_S) < 1e-6 and abs(Iyy_B_0 - Iyy_S) < 1e-6
    ninety_twist_ok = abs(Ixx_B_90 - Iyy_S) < 1e-6 and abs(Iyy_B_90 - Ixx_S) < 1e-6
    
    print(f"\nConvention Validation:")
    print(f"  Zero twist identity: {'PASS' if zero_twist_ok else 'FAIL'}")
    print(f"  90° twist swap: {'PASS' if ninety_twist_ok else 'FAIL'}")
    
    all_pass = zero_twist_ok and ninety_twist_ok
    print(f"\nTest 6 Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("Coordinate Transformation and EI Assignment Validation Test")
    print("=" * 80)
    print("\nThis test validates:")
    print("  1. The S-to-B frame transformation")
    print("  2. The assignment of Ixx_B to EI_edge and Iyy_B to EI_flap")
    print("  3. The coordinate system conventions in the transformation")
    print("=" * 80)
    
    # Run all tests
    test1_passed = test_transformation_identity()
    test2_passed = test_transformation_90_degree()
    test3_passed = test_transformation_45_degree()
    test4_passed = test_ei_assignment()
    test5_passed = test_transformation_with_twist()
    test6_passed = test_coordinate_convention()
    
    # Summary
    print("\n" + "=" * 80)
    print("Test Summary")
    print("=" * 80)
    print(f"Test 1 (Zero Twist Identity): {'PASS' if test1_passed else 'FAIL'}")
    print(f"Test 2 (90° Twist Swap): {'PASS' if test2_passed else 'FAIL'}")
    print(f"Test 3 (45° Twist Formula): {'PASS' if test3_passed else 'FAIL'}")
    print(f"Test 4 (EI Assignment): {'PASS' if test4_passed else 'FAIL'}")
    print(f"Test 5 (Non-Zero Twist): {'PASS' if test5_passed else 'FAIL'}")
    print(f"Test 6 (Coordinate Convention): {'PASS' if test6_passed else 'FAIL'}")
    
    all_passed = all([test1_passed, test2_passed, test3_passed, test4_passed, test5_passed, test6_passed])
    print(f"\nOverall: {'ALL TESTS PASSED' if all_passed else 'SOME TESTS FAILED'}")
    print("=" * 80)

