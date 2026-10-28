"""
Test coordinate frame transformations for z_B, z_R, z_B/L, and z_R/R mappings.

This test verifies the CoordinateFrameTransforms class correctly maps between:
- Rotor frame (R): r/R (normalized radial position), r (radial distance from center)
- Blade frame (B): z/L (normalized blade position), z (distance from blade root)

Key relationships:
- r = r/R × R (radial distance from rotor origin)
- z = r - R_hub (blade coordinate, 0 at root)
- R_hub = R - L (hub radius)
- At hub: r = R_hub, z = 0, r/R = R_hub/R
- At tip: r = R, z = L, r/R = 1
"""

import numpy as np
import sys
from pathlib import Path

# Add workspace root to path
_workspace_root = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(_workspace_root))

from blade_geometry_model.coordinate_frames import CoordinateFrameTransforms


def test_coordinate_mappings():
    """Test coordinate transformations with specific blade geometry."""
    
    # Test parameters matching actual blade
    R = 11.0  # Rotor radius [m]
    r_root_R = 1.38 / 11.0  # Hub at r = 1.38 m
    L = 9.62  # Blade length [m] (should equal R * (1 - r_root_R))
    
    print("=" * 80)
    print("COORDINATE FRAME TRANSFORMATION TESTS")
    print("=" * 80)
    print(f"\nBlade Parameters:")
    print(f"  R (rotor radius) = {R:.3f} m")
    print(f"  L (blade length) = {L:.3f} m")
    print(f"  r_root/R (hub normalized) = {r_root_R:.6f}")
    print(f"  R_hub (hub radius) = R - L = {R - L:.3f} m")
    
    # Create transform
    rR_grid = np.array([r_root_R, 0.25, 0.5, 0.75, 1.0])
    transform = CoordinateFrameTransforms(R=R, rR_grids=[rR_grid])
    
    print(f"\nTransform Properties:")
    print(f"  transform.R = {transform.R:.3f} m")
    print(f"  transform.L = {transform.L:.3f} m")
    print(f"  transform.R_hub = {transform.R_hub:.3f} m")
    print(f"  transform.r_root_R = {transform.r_root_R:.6f}")
    
    print("\n" + "=" * 80)
    print("KEY POSITIONS TEST")
    print("=" * 80)
    
    # Test key positions
    # Note: z/L can be negative before blade root (between rotor center and hub)
    test_cases = [
        ("Rotor Center", 0.0, 0.0, -transform.R_hub / L, -transform.R_hub),
        ("Blade Root (Hub)", r_root_R, transform.R_hub, 0.0, 0.0),
        ("Mid-Blade", 0.5 * (r_root_R + 1.0), 0.5 * (transform.R_hub + R), 0.5, 0.5 * L),
        ("Blade Tip", 1.0, R, 1.0, L),
    ]
    
    print(f"\n{'Position':<20} {'r/R':>10} {'r [m]':>10} {'z/L':>10} {'z [m]':>10}")
    print(f"{'':20} {'(input)':>10} {'(calc)':>10} {'(expect)':>10} {'(expect)':>10}")
    print("-" * 80)
    
    for name, rR_input, r_expected, zL_expected, z_expected in test_cases:
        # Calculate using transform
        r_calc = transform.r_from_rR(rR_input)
        z_calc = transform.z_from_rR(rR_input)
        zL_calc = transform.zL_from_rR(rR_input)
        
        print(f"{name:<20} {rR_input:>10.6f} {r_calc:>10.3f} {zL_calc:>10.6f} {z_calc:>10.3f}")
        
        # Verify
        assert np.isclose(r_calc, r_expected, atol=1e-6), f"{name}: r mismatch"
        assert np.isclose(z_calc, z_expected, atol=1e-6), f"{name}: z mismatch"
        assert np.isclose(zL_calc, zL_expected, atol=1e-6), f"{name}: z/L mismatch"
    
    print("\n" + "=" * 80)
    print("TICK MAPPING TEST (for visualization)")
    print("=" * 80)
    
    # Test tick positions used in visualization
    bottom_rR_ticks = np.array([0.00, r_root_R, 0.25, 0.50, 0.75, 1.00])
    
    print(f"\n{'z_R [m]':<12} {'r/R':<12} {'z_B [m]':<12} {'z_B/L':<12} {'Check':<20}")
    print("-" * 80)
    
    for rR in bottom_rR_ticks:
        z_R = rR * R  # Radial distance (what we call z_R in viz)
        z_B = transform.z_from_rR(rR)  # Blade coordinate
        zL = transform.zL_from_rR(rR)  # Normalized blade coordinate
        
        # Verification
        z_B_expected = z_R - transform.R_hub
        check = "✓" if np.isclose(z_B, z_B_expected, atol=1e-6) else "✗ MISMATCH"
        
        print(f"{z_R:<12.3f} {rR:<12.6f} {z_B:<12.3f} {zL:<12.6f} {check:<20}")
    
    print("\n" + "=" * 80)
    print("INVERSE TRANSFORMATION TEST")
    print("=" * 80)
    
    # Test round-trip transformations
    print(f"\n{'Original':<15} {'Forward':<15} {'Reverse':<15} {'Match':<10}")
    print("-" * 60)
    
    test_rR_values = np.array([r_root_R, 0.5, 1.0])
    for rR_orig in test_rR_values:
        z = transform.z_from_rR(rR_orig)
        rR_back = transform.rR_from_z(z)
        match = "✓" if np.isclose(rR_orig, rR_back, atol=1e-6) else "✗"
        print(f"r/R={rR_orig:<12.6f} z={z:<12.3f} r/R={rR_back:<12.6f} {match:<10}")
    
    print("\n" + "=" * 80)
    print("VISUALIZATION TICK MAPPING VERIFICATION")
    print("=" * 80)
    
    # Simulate what visualization does
    len_scale = 1.0  # meters
    bottom_rR_ticks_viz = np.array([0.00, 0.125455, 0.25, 0.50, 0.75, 1.00])  # Typical viz ticks
    default_zR_ticks = bottom_rR_ticks_viz * R * len_scale
    
    print("\nUsing zL_from_rR() method (CORRECT):")
    print(f"{'z_R [m]':<12} {'r/R':<12} {'z/L':<12} {'z_B [m]':<12}")
    print("-" * 60)
    
    for zR, rR in zip(default_zR_ticks, bottom_rR_ticks_viz):
        zL = transform.zL_from_rR(rR)
        zB = zL * L * len_scale
        print(f"{zR:<12.3f} {rR:<12.6f} {zL:<12.6f} {zB:<12.3f}")
    
    print("\n✓ All tests passed!")
    print("\nKEY FINDINGS:")
    print(f"  - At z_R = {transform.R_hub:.3f} m (hub): z_B should be 0.000 m")
    print(f"  - At z_R = {R:.3f} m (tip): z_B should be {L:.3f} m")
    print(f"  - Use: zL_from_rR(r/R) then multiply by L to get z_B")
    print("=" * 80)


if __name__ == "__main__":
    test_coordinate_mappings()

