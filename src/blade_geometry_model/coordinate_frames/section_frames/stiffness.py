"""B → S′ stiffness orchestration (depends on shoelace section properties)."""

from typing import Optional, Tuple, Union

import numpy as np

from blade_geometry_model.shoelace_integration import section_properties_S

from .tensors import (
    structural_product_moment,
    transform_second_moments_B_to_S,
    transform_second_moments_S_to_B,
)


def transform_stiffness_B_to_Sprime(
    E: float,
    EI_flap_B: Union[float, np.ndarray],
    scaled_dimensions,
    *,
    EI_edge_B: Optional[Union[float, np.ndarray]] = None,
    station_indices: Optional[np.ndarray] = None,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Transform generated EI curves from B frame to S' frame (centroid-referenced).
    
    This method maps the generated EI curves (structural stiffness distributions)
    from B frame to their equivalent values in S' frame. The generated EI curves
    are computed in B frame and represent the actual structural stiffness of the
    blade, not geometric properties.
    
    This is a convenience method that combines two transformations:
    1. B frame → S frame (pivot): Rotation transformation using twist angle with full Ixy support
    2. S frame (pivot) → S' frame (centroid): Parallel axis theorem using geometry
    
    This method uses full Ixy (product moment) in all transformations for accurate
    handling of asymmetric airfoils. It computes Ixy_B_structural by scaling geometric
    Ixy_B proportionally to the Iyy_B scaling factor (from generated EI_flap_B).
    
    This method is used when you have generated EI curves in B frame (from stiffness
    selection) and need centroid-referenced EI values in S' frame for stress computation.
    
    Parameters
    ----------
    E : float
        Elastic modulus [Pa]
    EI_flap_B : float or np.ndarray
        Generated flapwise bending stiffness in B frame [N·m²] (E * Iyy_B), shape (n_stations,)
        This is the generated EI curve from stiffness selection, not geometric properties.
    scaled_dimensions
        ScaledDimensions instance with:
        - outer_loop_geometry_S: S frame geometry (for computing centroid)
        - twist: twist angles [deg] for B→S transformation
        - outer_loop_point_counts: point counts per station
        - zL_grid: station grid (for indexing)
    EI_edge_B : float or np.ndarray, optional
        Generated edgewise bending stiffness in B frame [N·m²] (E * Ixx_B), shape (n_stations,)
        If not provided, derived from outline second moments transformed S→B.
        Degenerate loops raise (CF-002); no 1e-6 invent.
    station_indices : np.ndarray, optional
        Station indices to evaluate. If None, uses all stations.
        Shape (n_eval,) with values in [0, n_stations-1]
    
    Returns
    -------
    Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]
        (EI_xx_Sprime, EI_yy_Sprime) where:
        - EI_xx_Sprime: Equivalent bending stiffness about X_S' axis (centroid-referenced) [N·m²]
          in S' frame, transformed from EI_edge_B
        - EI_yy_Sprime: Equivalent bending stiffness about Y_S' axis (centroid-referenced) [N·m²]
          in S' frame, transformed from EI_flap_B
        Shape: (n_eval,) matching station_indices or all stations
        
    Notes
    -----
    This method transforms the generated EI curves (structural stiffness) from B frame
    to S' frame. It does NOT compute geometric properties from airfoil geometry.
    The input EI values should be the generated stiffness distributions from the
    stiffness selection process, which are B frame-consistent.
    """
    from shapely.geometry import Polygon as ShapelyPolygon
    
    EI_flap_B = np.asarray(EI_flap_B, dtype=np.float64)
    
    # Get geometry and twist data
    outer_loop_geometry_S = np.asarray(scaled_dimensions.outer_loop_geometry_S, dtype=np.float64)
    outer_loop_point_counts = np.asarray(scaled_dimensions.outer_loop_point_counts, dtype=np.int64)
    zL_grid = np.asarray(scaled_dimensions.zL_grid, dtype=np.float64)
    twist_deg = np.asarray(scaled_dimensions.twist, dtype=np.float64)
    
    # Handle station indices
    if station_indices is not None:
        station_indices = np.asarray(station_indices, dtype=np.int64)
        n_eval = len(station_indices)
        EI_flap_B = EI_flap_B[station_indices]
        twist_deg = twist_deg[station_indices]
        if EI_edge_B is not None:
            EI_edge_B = np.asarray(EI_edge_B, dtype=np.float64)
            # Only index if EI_edge_B is full size (not already sliced)
            if EI_edge_B.shape[0] == len(zL_grid):
                EI_edge_B = EI_edge_B[station_indices]
            # If already the right size, assume it's already been indexed
            elif EI_edge_B.shape[0] != n_eval:
                raise ValueError(
                    f"EI_edge_B.shape={EI_edge_B.shape} must match either "
                    f"full size ({len(zL_grid)}) or station_indices size ({n_eval})"
                )
    else:
        station_indices = np.arange(len(zL_grid), dtype=np.int64)
        n_eval = len(zL_grid)
        if EI_edge_B is not None:
            EI_edge_B = np.asarray(EI_edge_B, dtype=np.float64)
    
    # Missing EI_edge_B: derive from outline Ixx in B (not S, and never invent 1e-6).
    if EI_edge_B is None:
        EI_edge_B = np.zeros(n_eval, dtype=np.float64)
        for i, idx in enumerate(station_indices):
            n_points = int(outer_loop_point_counts[idx])
            if n_points < 3:
                raise ValueError(
                    f"station {int(idx)} outer loop has {n_points} points; "
                    "cannot invent EI_edge_B (CF-002)"
                )
            coords_S = outer_loop_geometry_S[idx, :n_points, :]
            poly_S = ShapelyPolygon(coords_S)
            _, _, Ixx_S, Iyy_S, Ixy_S = section_properties_S(poly_S)
            twist_deg_val = twist_deg[i] if hasattr(twist_deg, "__getitem__") else twist_deg
            Ixx_B, _Iyy_B, _Ixy_B = transform_second_moments_S_to_B(
                Ixx_S=Ixx_S,
                Iyy_S=Iyy_S,
                Ixy_S=Ixy_S,
                twist_deg=twist_deg_val,
            )
            EI_edge_B[i] = E * float(Ixx_B)
    
    # Step 1: Compute structural I values in B frame with full Ixy support
    # Get geometric properties from geometry (S frame, pivot-referenced)
    # Transform S→B with full Ixy, then scale to structural values
    Ixx_S_geom = np.zeros(n_eval, dtype=np.float64)
    Iyy_S_geom = np.zeros(n_eval, dtype=np.float64)
    Ixy_S_geom = np.zeros(n_eval, dtype=np.float64)
    
    for i, idx in enumerate(station_indices):
        n_points = int(outer_loop_point_counts[idx])
        if n_points < 3:
            raise ValueError(
                f"station {int(idx)} outer loop has {n_points} points; "
                "refusing 1e-12/E invent for geometric I (CF-002)"
            )
        coords_S = outer_loop_geometry_S[idx, :n_points, :]
        poly_S = ShapelyPolygon(coords_S)
        _, _, Ixx_S_val, Iyy_S_val, Ixy_S_val = section_properties_S(poly_S)
        Ixx_S_geom[i] = Ixx_S_val
        Iyy_S_geom[i] = Iyy_S_val
        Ixy_S_geom[i] = Ixy_S_val
    
    # Transform geometric properties S→B with full Ixy
    Ixx_B_geom = np.zeros(n_eval, dtype=np.float64)
    Iyy_B_geom = np.zeros(n_eval, dtype=np.float64)
    Ixy_B_geom = np.zeros(n_eval, dtype=np.float64)
    
    for i in range(n_eval):
        twist_deg_val = twist_deg[i] if hasattr(twist_deg, '__getitem__') else twist_deg
        Ixx_B_val, Iyy_B_val, Ixy_B_val = transform_second_moments_S_to_B(
            Ixx_S=Ixx_S_geom[i],
            Iyy_S=Iyy_S_geom[i],
            Ixy_S=Ixy_S_geom[i],  # Use full Ixy instead of 0.0 for accuracy
            twist_deg=twist_deg_val,
        )
        Ixx_B_geom[i] = Ixx_B_val
        Iyy_B_geom[i] = Iyy_B_val
        Ixy_B_geom[i] = Ixy_B_val
    
    # Get structural I values in B frame
    # Iyy_B_structural comes directly from generated EI_flap_B (from optimization)
    Iyy_B_structural = EI_flap_B / E  # Flapwise: Iyy_B from EI_flap_B
    # Ixx_B_structural comes from EI_edge_B (computed using ratio method)
    Ixx_B_structural = EI_edge_B / E  # Edgewise: Ixx_B from EI_edge_B
    
    # Ixy is closed from the outer geometry -- nobody supplies it. The closure
    # preserves the section's correlation coefficient, which keeps the tensor
    # positive-definite for any positive (Ixx, Iyy) and reproduces the historical
    # flapwise-only scaling exactly whenever EI_edge_B came from the ratio
    # fallback. `compute_stress_at_touch_points_Sprime` forwards a caller-supplied
    # EI_edge_B here, so this site carried the same defect as the S-frame one.
    Ixy_B_structural = structural_product_moment(
        Ixy_B_geom=Ixy_B_geom,
        Ixx_B_geom=Ixx_B_geom,
        Iyy_B_geom=Iyy_B_geom,
        Ixx_B_structural=Ixx_B_structural,
        Iyy_B_structural=Iyy_B_structural,
    )
    
    # Transform from B to S frame (about pivot) with full Ixy
    Ixx_S_pivot = np.zeros(n_eval, dtype=np.float64)
    Iyy_S_pivot = np.zeros(n_eval, dtype=np.float64)
    Ixy_S_pivot = np.zeros(n_eval, dtype=np.float64)
    
    for i in range(n_eval):
        twist_deg_val = twist_deg[i] if hasattr(twist_deg, '__getitem__') else twist_deg
        Ixx_S_val, Iyy_S_val, Ixy_S_val = transform_second_moments_B_to_S(
            Ixx_B=Ixx_B_structural[i],
            Iyy_B=Iyy_B_structural[i],
            Ixy_B=Ixy_B_structural[i],  # Use full Ixy instead of 0.0 for accuracy
            twist_deg=twist_deg_val,
        )
        Ixx_S_pivot[i] = Ixx_S_val
        Iyy_S_pivot[i] = Iyy_S_val
        Ixy_S_pivot[i] = Ixy_S_val
    
    # Step 2: Transform from S frame (pivot) to S' frame (centroid) using parallel axis theorem
    EI_xx_Sprime = np.zeros(n_eval, dtype=np.float64)
    EI_yy_Sprime = np.zeros(n_eval, dtype=np.float64)
    
    for i, idx in enumerate(station_indices):
        n_points = int(outer_loop_point_counts[idx])
        if n_points < 3:
            raise ValueError(
                f"station {int(idx)} outer loop has {n_points} points; "
                "cannot form S' centroid EI (CF-002)"
            )
        coords_S = outer_loop_geometry_S[idx, :n_points, :]
        poly_S = ShapelyPolygon(coords_S)
        A, (Cx, Cy), _, _, _ = section_properties_S(poly_S)
        Ixx_centroid = Ixx_S_pivot[i] - A * Cy**2
        Iyy_centroid = Iyy_S_pivot[i] - A * Cx**2
        if Ixx_centroid <= 0.0 or Iyy_centroid <= 0.0:
            raise ValueError(
                f"station {int(idx)} centroid second moment is non-positive "
                f"(Ixx={Ixx_centroid:.3e}, Iyy={Iyy_centroid:.3e}); "
                "refusing 1e-12 clamp (CF-002)"
            )
        EI_xx_Sprime[i] = E * Ixx_centroid
        EI_yy_Sprime[i] = E * Iyy_centroid
    
    # Return scalars if inputs were scalars
    if isinstance(EI_flap_B, (int, float)) or (hasattr(EI_flap_B, 'ndim') and EI_flap_B.ndim == 0):
        return float(EI_xx_Sprime), float(EI_yy_Sprime)
    return EI_xx_Sprime, EI_yy_Sprime
