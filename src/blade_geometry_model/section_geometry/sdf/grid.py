"""
Structured background grid for SDF evaluation and post-processing.

SDFGrid wraps a 2-D Cartesian grid and provides:
  - Vectorised SDF evaluation for any callable
  - Gradient computation (∇φ via central finite differences)
  - Gradient magnitude |∇φ|
  - Zero-level-set extraction as a polyline
  - Masked interior / exterior arrays
  - 2-D quadrature (area, centroid, moments of inertia) for phi < 0
"""

import numpy as np

try:
    import contourpy
    _HAS_CONTOURPY = True
except Exception:
    contourpy = None
    _HAS_CONTOURPY = False


class SDFGrid:
    """Structured 2-D grid for SDF operations.

    Parameters
    ----------
    X, Y : ndarray, shape (ny, nx)
        Meshgrid coordinate arrays (from np.meshgrid with indexing='xy').
    """

    def __init__(self, X, Y):
        if X.shape != Y.shape:
            raise ValueError("X and Y must have identical shape.")
        self.X  = X
        self.Y  = Y
        self.nx = X.shape[1]
        self.ny = X.shape[0]
        self.dx = float(X[0, 1] - X[0, 0]) if self.nx > 1 else 1.0
        self.dy = float(Y[1, 0] - Y[0, 0]) if self.ny > 1 else 1.0
        self._rotated_coords_cache: dict[tuple[float, float, float], tuple[np.ndarray, np.ndarray]] = {}

    @classmethod
    def from_bbox(cls, x_min, x_max, y_min, y_max, nx=400, ny=200):
        """Create a grid covering a bounding box."""
        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)
        X, Y = np.meshgrid(x, y)
        return cls(X, Y)

    @classmethod
    def from_airfoil(cls, airfoil_sdf, padding=0.1, nx=400, ny=200):
        """Automatically size the grid to wrap an AirfoilSDF with padding."""
        verts = airfoil_sdf.vertices
        x_min = verts[:, 0].min() - padding
        x_max = verts[:, 0].max() + padding
        y_min = verts[:, 1].min() - padding
        y_max = verts[:, 1].max() + padding
        return cls.from_bbox(x_min, x_max, y_min, y_max, nx=nx, ny=ny)

    def eval(self, sdf_callable):
        """Evaluate an SDF callable on the full grid."""
        try:
            from blade_geometry_model.section_geometry.sdf.profiling import profiling_enabled, record_grid_eval

            if profiling_enabled():
                record_grid_eval(n_points=int(self.X.size))
        except Exception:
            pass
        return np.asarray(sdf_callable(self.X, self.Y), dtype=float)

    def rotated_coords(self, angle, cx=0.0, cy=0.0):
        """Return query coordinates inverse-rotated by ``angle`` around ``(cx, cy)``."""
        key = (float(angle), float(cx), float(cy))
        cached = self._rotated_coords_cache.get(key)
        if cached is not None:
            return cached
        cos_a = np.cos(-float(angle))
        sin_a = np.sin(-float(angle))
        xr = self.X - float(cx)
        yr = self.Y - float(cy)
        xp = cos_a * xr - sin_a * yr + float(cx)
        yp = sin_a * xr + cos_a * yr + float(cy)
        self._rotated_coords_cache[key] = (xp, yp)
        return xp, yp

    def rotated_eval(self, sdf_callable, angle, cx=0.0, cy=0.0):
        """Evaluate ``sdf_callable`` on inverse-rotated query coordinates."""
        xp, yp = self.rotated_coords(angle, cx=cx, cy=cy)
        return np.asarray(sdf_callable(xp, yp), dtype=float)

    def gradient(self, phi):
        """Central-difference gradient → (gx, gy), shape (ny, nx) each."""
        gx = np.gradient(phi, self.dx, axis=1)
        gy = np.gradient(phi, self.dy, axis=0)
        return gx, gy

    def grad_magnitude(self, phi):
        """Return |∇φ| on the grid."""
        gx, gy = self.gradient(phi)
        return np.sqrt(gx**2 + gy**2)

    def eikonal_error(self, phi):
        """Pointwise Eikonal residual  | |∇φ| − 1 |."""
        return np.abs(self.grad_magnitude(phi) - 1.0)

    def zero_contour(self, phi):
        """Extract the zero-level-set as a list of (N, 2) polyline arrays."""
        return self.level_set(phi, level=0.0)

    def level_set(self, phi, level=0.0):
        """Extract an arbitrary level-set contour."""
        if _HAS_CONTOURPY:
            x = np.asarray(self.X[0, :], dtype=float)
            y = np.asarray(self.Y[:, 0], dtype=float)
            z = np.asarray(phi, dtype=float)
            cg = contourpy.contour_generator(x=x, y=y, z=z, name="serial")
            raw_segments = cg.lines(float(level))
            return [np.asarray(seg, dtype=float).copy() for seg in raw_segments if len(seg) >= 2]
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        cs = ax.contour(self.X, self.Y, phi, levels=[level])
        segments = [path.vertices.copy() for path in cs.get_paths()]
        plt.close(fig)
        return segments

    def interior_mask(self, phi):
        """Boolean mask: True where phi < 0 (inside the shape)."""
        return phi < 0.0

    def boundary_mask(self, phi, tol=None):
        """Boolean mask near the zero-level-set."""
        tol = tol or max(self.dx, self.dy)
        return np.abs(phi) < tol

    def area(self, phi):
        """Enclosed area (phi < 0) by cell-counting quadrature."""
        return float(self.section_properties_fused(phi)["area"])

    def centroid(self, phi):
        """Centroid of the region phi < 0 → (cx, cy)."""
        props = self.section_properties_fused(phi)
        return float(props["cx"]), float(props["cy"])

    def second_moments(self, phi):
        """Second moments of area (Ixx, Iyy, Ixy) about the centroid."""
        props = self.section_properties_fused(phi)
        return float(props["Ixx"]), float(props["Iyy"]), float(props["Ixy"])

    def section_properties_fused(self, phi):
        """Return section properties for ``phi < 0`` in a single pass.

        Keys: area, cx, cy, Ixx, Iyy, Ixy, r_gyr_x, r_gyr_y
        """
        mask = np.asarray(phi < 0.0, dtype=float)
        dA = self.dx * self.dy
        area_cells = float(mask.sum())
        A = area_cells * dA
        if area_cells == 0.0:
            return {
                "area": 0.0, "cx": float("nan"), "cy": float("nan"),
                "Ixx": float("nan"), "Iyy": float("nan"), "Ixy": float("nan"),
                "r_gyr_x": float("nan"), "r_gyr_y": float("nan"),
            }
        sum_x = float((mask * self.X).sum())
        sum_y = float((mask * self.Y).sum())
        cx = sum_x / area_cells
        cy = sum_y / area_cells
        xr = self.X - cx
        yr = self.Y - cy
        Ixx = float((mask * yr**2).sum() * dA)
        Iyy = float((mask * xr**2).sum() * dA)
        Ixy = float((mask * xr * yr).sum() * dA)
        return {
            "area": A, "cx": float(cx), "cy": float(cy),
            "Ixx": Ixx, "Iyy": Iyy, "Ixy": Ixy,
            "r_gyr_x": float(np.sqrt(Ixx / A)) if A > 0 else float("nan"),
            "r_gyr_y": float(np.sqrt(Iyy / A)) if A > 0 else float("nan"),
        }

    def section_properties(self, phi):
        """Alias for section_properties_fused."""
        return self.section_properties_fused(phi)

    def __repr__(self):
        return (
            f"SDFGrid(nx={self.nx}, ny={self.ny}, "
            f"dx={self.dx:.4g}, dy={self.dy:.4g}, "
            f"x=[{self.X.min():.3g},{self.X.max():.3g}], "
            f"y=[{self.Y.min():.3g},{self.Y.max():.3g}])"
        )
