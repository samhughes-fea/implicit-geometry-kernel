"""
Verify that the fix (using -theta) works correctly.
"""

import numpy as np
import sys
from pathlib import Path

# Add workspace root to path
_workspace_root = Path(__file__).parent.parent.parent.parent.parent
if str(_workspace_root) not in sys.path:
    sys.path.insert(0, str(_workspace_root))

from blade_geometry_model.coordinate_frames import CoordinateFrameTransforms


print("=" * 80)
print("Verification: Fixed Transformation (using -theta)")
print("=" * 80)

# Test case with asymmetric shape
Ixx_S = 0.100
Iyy_S = 0.066
Ixy_S = 0.020
theta_deg = 30.0

print(f"\nInput (S frame): Ixx_S={Ixx_S:.6f}, Iyy_S={Iyy_S:.6f}, Ixy_S={Ixy_S:.6f}")
print(f"Twist angle: {theta_deg} deg")
print(f"Expected: Transformation should use -{theta_deg} deg internally")

# Test with fixed transformation
transform = CoordinateFrameTransforms()
Ixx_B, Iyy_B, Ixy_B = transform.transform_second_moments_S_to_B(Ixx_S, Iyy_S, Ixy_S, theta_deg)

print(f"\nResult (B frame):")
print(f"  Ixx_B = {Ixx_B:.6f} {'[NEGATIVE!]' if Ixx_B < 0 else ''}")
print(f"  Iyy_B = {Iyy_B:.6f} {'[NEGATIVE!]' if Iyy_B < 0 else ''}")
print(f"  Ixy_B = {Ixy_B:.6f}")
print(f"  Sum = {Ixx_B + Iyy_B:.6f}")
print(f"  Both non-negative? {Ixx_B >= 0 and Iyy_B >= 0}")

# Expected from matrix approach: Ixx_B = 0.108821, Iyy_B = 0.057179
expected_Ixx_B = 0.108821
expected_Iyy_B = 0.057179

print(f"\nExpected (from matrix approach):")
print(f"  Ixx_B = {expected_Ixx_B:.6f}")
print(f"  Iyy_B = {expected_Iyy_B:.6f}")

print(f"\nMatch? {abs(Ixx_B - expected_Ixx_B) < 1e-4 and abs(Iyy_B - expected_Iyy_B) < 1e-4}")

# Test at 0 deg (should be identity)
print("\n" + "=" * 80)
print("Test: 0 deg (should be identity)")
print("=" * 80)

Ixx_S_0 = 0.1667
Iyy_S_0 = 0.6667
Ixy_S_0 = 0.0

Ixx_B_0, Iyy_B_0, Ixy_B_0 = transform.transform_second_moments_S_to_B(Ixx_S_0, Iyy_S_0, Ixy_S_0, 0.0)

print(f"Input (S frame): Ixx_S={Ixx_S_0:.6f}, Iyy_S={Iyy_S_0:.6f}, Ixy_S={Ixy_S_0:.6f}")
print(f"Output (B frame): Ixx_B={Ixx_B_0:.6f}, Iyy_B={Iyy_B_0:.6f}, Ixy_B={Ixy_B_0:.6f}")
print(f"Identity? {abs(Ixx_B_0 - Ixx_S_0) < 1e-6 and abs(Iyy_B_0 - Iyy_S_0) < 1e-6}")

# Test at 45 deg with symmetric shape
print("\n" + "=" * 80)
print("Test: 45 deg with symmetric shape (rectangle)")
print("=" * 80)

# Rectangle rotated 45°: Ixx_S = Iyy_S = (Ixx_B + Iyy_B) / 2
Ixx_S_45 = 0.4167
Iyy_S_45 = 0.4167
Ixy_S_45 = 0.0

Ixx_B_45, Iyy_B_45, Ixy_B_45 = transform.transform_second_moments_S_to_B(Ixx_S_45, Iyy_S_45, Ixy_S_45, 45.0)

print(f"Input (S frame): Ixx_S={Ixx_S_45:.6f}, Iyy_S={Iyy_S_45:.6f}, Ixy_S={Ixy_S_45:.6f}")
print(f"Output (B frame): Ixx_B={Ixx_B_45:.6f}, Iyy_B={Iyy_B_45:.6f}, Ixy_B={Ixy_B_45:.6f}")
print(f"Expected: Ixx_B ≈ Iyy_B ≈ 0.4167 (for 45° back-transformation)")
print(f"Match? {abs(Ixx_B_45 - 0.4167) < 0.01 and abs(Iyy_B_45 - 0.4167) < 0.01}")

print("\n" + "=" * 80)
print("Fix Verification Complete")
print("=" * 80)


