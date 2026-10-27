"""R↔B span grid converters and locii."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union

import numpy as np


@dataclass
class SpanFrameTransforms:
    """
    Transforms between Radial (R), Blade (B), Station (S), and S′ frames.

    - **R**: \(r\), \(r/R\) — external turbine grid / visualization.
    - **B**: \(X_B\) flap (vertical), \(Y_B\) edge (horizontal), \(Z_B\) span.
      Section outline in B **includes geometric twist** about \(Z_B\).
      \(\mathrm{EI}_{flap}=E\,I_{yy,B}\), \(\mathrm{EI}_{edge}=E\,I_{xx,B}\).
    - **S**: \(X_S\) chord-normal, \(Y_S\) chord-parallel, \(Z_S\) span —
      **zero-twist** chord-aligned section view (origin at 0.333 chordline).
    - **S′**: same orientation as S; origin at geometric centroid.

    \(Y_B\)/\(Y_S\) sense: \(+Y\) toward TE (LE at negative \(Y\)); standard plots
    invert the horizontal axis so the LE still appears on the RHS.
    Always use frame names R/B/S/S′ — never bare local/global.

    These are coordinate-system conversions, not airfoil scaling. For profile
    geometry transforms see Scale / section_geometry.

    Parameters
    ----------
    R : float, optional
        Rotor radius [m], measured from rotor center to blade tip.
        Required for R↔B coordinate conversions. Default 0.0.
        Can be 0.0 if only using S↔B second moment transformations.
    rR_grids : Union[np.ndarray, List[np.ndarray]], optional
        r/R grid(s) from which the minimum r/R value will be found to compute L.
        Can be a single array or a list of arrays (e.g., chord source, twist source, sampling grid).
        Required for R↔B coordinate conversions. Default None.
        Can be None if only using S↔B second moment transformations.

    Notes
    -----
    For R↔B coordinate conversions:
    - R and rR_grids must be provided
    - L is computed as L = R * (1 - r_root_R) where r_root_R is the minimum value across all rR_grids.
    - Must satisfy 0 < L < R, which implies 0 < r_root_R < 1.
    - Downstream modules (InternalForce) reuse the r/R and z/L grids generated through this transform.
    
    For S↔B second moment transformations:
    - R and rR_grids are not required
    - Can instantiate with default values: CoordinateFrameTransforms()
    
    Always use explicit frame names (R, B, S) - avoid relative terms like "local/global".
    """

    R: float = 0.0  # Optional, only needed for R↔B conversions
    rR_grids: Optional[Union[np.ndarray, List[np.ndarray]]] = None  # Optional, only needed for R↔B conversions

    # Computed values
    _r_root_R: float = field(init=False, default=0.0)
    
    # Locii
    r_locii: Dict[str, float] = field(init=False, default_factory=dict)
    z_locii: Dict[str, float] = field(init=False, default_factory=dict)
    rR_locii: Dict[str, float] = field(init=False, default_factory=dict)
    zL_locii: Dict[str, float] = field(init=False, default_factory=dict)

    def __post_init__(self) -> None:
        # Only validate and compute R↔B values if R and rR_grids are provided
        if self.R > 0 and self.rR_grids is not None:
            # Validate R
            if self.R <= 0:
                raise ValueError("R must be positive")
            
            # Find minimum r/R from grids and compute all derived values
            self._find_min_rR()
        else:
            # S↔B only usage - set default values
            self._r_root_R = 0.0
            self.r_locii = {}
            self.z_locii = {}
            self.rR_locii = {}
            self.zL_locii = {}

    def _find_min_rR(self) -> None:
        """
        Find minimum r/R value from all provided grids and compute L and locii.
        
        This method:
        1. Finds the smallest r/R value across all grids
        2. Sets _r_root_R
        3. Computes L (via property access)
        4. Computes all locii
        """
        # Collect all grids into a list
        grids_to_check = []
        if isinstance(self.rR_grids, np.ndarray):
            grids_to_check.append(np.asarray(self.rR_grids, dtype=np.float64))
        else:
            for grid in self.rR_grids:
                grids_to_check.append(np.asarray(grid, dtype=np.float64))
        
        # Find minimum across all grids
        min_values = [float(np.min(grid)) for grid in grids_to_check if grid.size > 0]
        if not min_values:
            raise ValueError("At least one r/R grid must contain values")
        
        self._r_root_R = float(min(min_values))
        if not 0.0 <= self._r_root_R < 1.0:
            raise ValueError(f"Smallest r/R value must be in [0, 1), got {self._r_root_R}")
        
        # Compute L (triggers property computation)
        _ = self.L  # Access property to trigger computation
        
        # Compute all locii
        self._compute_fixed_locii()
    
    @property
    def r_root_R(self) -> float:
        """Normalized radial position of blade root (smallest r/R value)."""
        return self._r_root_R
    
    @property
    def L(self) -> float:
        """
        Blade length [m], computed from r_root_R as L = R * (1 - r_root_R).
        
        This is a computed property - L is always derived from the minimum r/R value.
        """
        L_computed = self.R * (1.0 - self._r_root_R)
        if L_computed <= 0 or L_computed > self.R:
            raise ValueError(f"Computed L ({L_computed}) must satisfy 0 < L <= R ({self.R})")
        return L_computed

    @property
    def R_hub(self) -> float:
        """Distance from rotor origin to blade root [m]."""
        return self.R - self.L
    
    def _compute_fixed_locii(self) -> None:
        """Compute fixed locii from current R and L values."""
        R = self.R
        L = self.L
        R_hub = self.R_hub

        self.rR_locii = {
            "rotor_origin": 0.0,
            "rotor_hub": R_hub / R,
            "rotor_edge": 1.0,
        }

        self.zL_locii = {
            "blade_root": 0.0,
            "blade_tip": 1.0,
        }

        self.r_locii = {key: val * R for key, val in self.rR_locii.items()}
        self.z_locii = {key: val * L for key, val in self.zL_locii.items()}

    # ──────────────────────────────────────────────
    # R↔B coordinate conversion utilities
    # ──────────────────────────────────────────────

    # ------ finding radial r coordinates from rR, z, zL ------
    def r_from_rR(self, rR: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert r/R to r [m]."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return rR * self.R

    def r_from_z(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert z [m] to r [m]."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return z + self.R_hub

    def r_from_zL(self, zL: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert z/L to r [m]."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return zL * self.L + self.R_hub

    # ------ finding spanwise z coordinates from r, rR, zL ------
    def z_from_r(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert r [m] to z [m]."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return r - self.R_hub

    def z_from_rR(self, rR: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert r/R to z [m]."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return self.r_from_rR(rR) - self.R_hub

    def z_from_zL(self, zL: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert z/L to z [m]."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return zL * self.L

    # ------ finding non-dimensional radial rR coordinates from r, z, zL ------
    def rR_from_r(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert r [m] to r/R."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return r / self.R
    
    def rR_from_z(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert z [m] to r/R."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return self.r_from_z(z) / self.R

    def rR_from_zL(self, zL: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert z/L to r/R."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return self.r_from_zL(zL) / self.R

    # ------ finding non-dimensional spanwise zL coordinates from r, rR, z ------
    def zL_from_r(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert r [m] to z/L."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return self.z_from_r(r) / self.L

    def zL_from_rR(self, rR: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert r/R to z/L."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return self.z_from_rR(rR) / self.L

    def zL_from_z(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Convert z [m] to z/L."""
        if self.R <= 0 or self.rR_grids is None:
            raise ValueError("R and rR_grids must be provided for R↔B transformations. Use CoordinateFrameTransforms(R=..., rR_grids=...) to create instance with R↔B support.")
        return z / self.L

