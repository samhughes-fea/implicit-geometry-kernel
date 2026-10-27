"""
Test to verify coordinate system interpretation of inertia formulas.

This test validates that the standard inertia formulas (Ixx = ∫ y² dA, Iyy = ∫ x² dA)
are correctly interpreted in our coordinate systems:
- S frame: X_S = vertical (chord normal), Y_S = horizontal (chord parallel)
- B frame: X_B = vertical (flapwise), Y_B = horizontal (edgewise)

NOTE: This test creates polygons directly with coordinates already in S frame convention
(X=vertical, Y=horizontal). It does NOT test NACA coordinate mapping. For NACA coordinate
mapping tests, see test_naca_coordinate_mapping.py.

Key validation:
- For a wide rectangle (width > height), Iyy_S should be > Ixx_S (flapwise stiffer)
- After transformation to B frame, Iyy_B should be > Ixx_B
- EI_flap = E × Iyy_B should be > EI_edge = E × Ixx_B
"""

import sys
import numpy as np
from pathlib import Path

# Add workspace root to path
_workspace_root = Path(__file__).parent.parent.parent.parent.parent
if str(_workspace_root) not in sys.path:
    sys.path.insert(0, str(_workspace_root))

from shapely.geometry import Polygon as ShapelyPolygon
from blade_geometry_model.shoelace_integration import section_properties_S
from blade_geometry_model.coordinate_frames import CoordinateFrameTransforms


def create_rectangle_polygon(width, height, x0=0.0, y0=0.0):
    """
    Create a rectangle polygon with specified dimensions.
    
    Parameters
    ----------
    width : float
        Width in x-direction (horizontal/chord) [m]
    height : float
        Height in y-direction (vertical/chord normal) [m]
    x0 : float, optional
        Bottom-left corner x-coordinate. Default 0.0
    y0 : float, optional
        Bottom-left corner y-coordinate. Default 0.0
    
    Returns
    -------
    ShapelyPolygon
        Rectangle polygon (counter-clockwise)
    
    Note
    ----
    In standard notation:
    - x = horizontal (chord direction) → Y_S in S frame
    - y = vertical (chord normal) → X_S in S frame
    """
    # Create rectangle vertices (counter-clockwise)
    coords = [
        (x0, y0),                    # Bottom-left
        (x0 + width, y0),            # Bottom-right
        (x0 + width, y0 + height),   # Top-right
        (x0, y0 + height),            # Top-left
        (x0, y0),                     # Close polygon
    ]
    return ShapelyPolygon(coords)


def test_coordinate_system_interpretation():
    """
    Test coordinate system interpretation with wide rectangle.
    
    Creates a rectangle with width (horizontal/chord) > height (vertical/chord normal).
    Verifies:
    1. Ixx_S (about Y_S, horizontal/edgewise) < Iyy_S (about X_S, vertical/flapwise)
    2. After transformation to B frame: Ixx_B (about Y_B, edgewise) < Iyy_B (about X_B, flapwise)
    3. EI_flap = E × Iyy_B > EI_edge = E × Ixx_B
    """
    print("=" * 80)
    print("Coordinate System Interpretation Test")
    print("=" * 80)
    
    # Create wide rectangle: width (horizontal/chord) > height (vertical/chord normal)
    width = 5.0   # m (horizontal/chord direction → Y_S)
    height = 1.0  # m (vertical/chord normal → X_S)
    
    print(f"\nRectangle: width={width} m (horizontal/chord), height={height} m (vertical/chord normal)")
    print(f"In standard notation: x=horizontal, y=vertical")
    print(f"In S frame: X_S=vertical (chord normal), Y_S=horizontal (chord parallel)")
    
    # Create polygon
    rect_poly = create_rectangle_polygon(width, height, x0=0.0, y0=0.0)
    
    # Compute properties in S frame
    A_S, (Cx_S, Cy_S), Ixx_S, Iyy_S, Ixy_S = section_properties_S(rect_poly)
    
    print(f"\nS Frame Properties (computed from polygon):")
    print(f"  A_S = {A_S:.6e} m²")
    print(f"  Cx_S = {Cx_S:.6e} m (centroid X_S, vertical)")
    print(f"  Cy_S = {Cy_S:.6e} m (centroid Y_S, horizontal)")
    print(f"  Ixx_S = {Ixx_S:.6e} m⁴ (about Y_S axis, horizontal/edgewise)")
    print(f"  Iyy_S = {Iyy_S:.6e} m⁴ (about X_S axis, vertical/flapwise)")
    print(f"  Ixy_S = {Ixy_S:.6e} m⁴")
    
    # Verify interpretation
    print(f"\n" + "-" * 80)
    print("Verification 1: S Frame Interpretation")
    print("-" * 80)
    print(f"Ixx_S (about Y_S, horizontal/edgewise) = {Ixx_S:.6e} m⁴")
    print(f"Iyy_S (about X_S, vertical/flapwise) = {Iyy_S:.6e} m⁴")
    print(f"Ratio Iyy_S / Ixx_S = {Iyy_S/Ixx_S:.6f}")
    
    if Iyy_S > Ixx_S:
        print("✓ PASS: Iyy_S > Ixx_S (flapwise stiffer than edgewise for wide rectangle)")
    else:
        print("✗ FAIL: Iyy_S <= Ixx_S (INCORRECT INTERPRETATION!)")
        return False
    
    # Transform to B frame (test with zero twist first)
    transform = CoordinateFrameTransforms()
    twist_deg = 0.0
    
    Ixx_B, Iyy_B, Ixy_B = transform.transform_second_moments_S_to_B(
        Ixx_S, Iyy_S, Ixy_S, twist_deg
    )
    
    print(f"\n" + "-" * 80)
    print("Verification 2: B Frame Transformation (zero twist)")
    print("-" * 80)
    print(f"Ixx_B (about Y_B, horizontal/edgewise) = {Ixx_B:.6e} m⁴")
    print(f"Iyy_B (about X_B, vertical/flapwise) = {Iyy_B:.6e} m⁴")
    print(f"Ratio Iyy_B / Ixx_B = {Iyy_B/Ixx_B:.6f}")
    
    if Iyy_B > Ixx_B:
        print("✓ PASS: Iyy_B > Ixx_B (flapwise stiffer than edgewise)")
    else:
        print("✗ FAIL: Iyy_B <= Ixx_B (TRANSFORMATION ERROR!)")
        return False
    
    # Test with non-zero twist
    twist_deg = 15.0
    Ixx_B_twist, Iyy_B_twist, Ixy_B_twist = transform.transform_second_moments_S_to_B(
        Ixx_S, Iyy_S, Ixy_S, twist_deg
    )
    
    print(f"\n" + "-" * 80)
    print("Verification 3: B Frame Transformation (twist = 15°)")
    print("-" * 80)
    print(f"Ixx_B (about Y_B, horizontal/edgewise) = {Ixx_B_twist:.6e} m⁴")
    print(f"Iyy_B (about X_B, vertical/flapwise) = {Iyy_B_twist:.6e} m⁴")
    print(f"Ratio Iyy_B / Ixx_B = {Iyy_B_twist/Ixx_B_twist:.6f}")
    
    if Iyy_B_twist > Ixx_B_twist:
        print("✓ PASS: Iyy_B > Ixx_B (flapwise stiffer than edgewise with twist)")
    else:
        print("✗ FAIL: Iyy_B <= Ixx_B (TRANSFORMATION ERROR WITH TWIST!)")
        return False
    
    # Verify EI assignments
    E = 1.0e9  # Pa (example modulus)
    EI_flap = E * Iyy_B
    EI_edge = E * Ixx_B
    
    print(f"\n" + "-" * 80)
    print("Verification 4: EI Stiffness Assignments")
    print("-" * 80)
    print(f"EI_flap = E × Iyy_B = {EI_flap:.6e} N·m² (bending about X_B, vertical/flapwise)")
    print(f"EI_edge = E × Ixx_B = {EI_edge:.6e} N·m² (bending about Y_B, horizontal/edgewise)")
    print(f"Ratio EI_flap / EI_edge = {EI_flap/EI_edge:.6f}")
    
    if EI_flap > EI_edge:
        print("✓ PASS: EI_flap > EI_edge (correct for wide rectangle)")
    else:
        print("✗ FAIL: EI_flap <= EI_edge (INCORRECT ASSIGNMENT!)")
        return False
    
    print(f"\n" + "=" * 80)
    print("Test Result: PASS")
    print("=" * 80)
    print("\nSummary:")
    print("  - Standard formulas correctly interpreted in S frame")
    print("  - Ixx_S (about Y_S, edgewise) < Iyy_S (about X_S, flapwise) ✓")
    print("  - Transformation to B frame preserves relationship ✓")
    print("  - EI_flap = E × Iyy_B > EI_edge = E × Ixx_B ✓")
    print("=" * 80)
    
    return True


if __name__ == "__main__":
    success = test_coordinate_system_interpretation()
    sys.exit(0 if success else 1)

