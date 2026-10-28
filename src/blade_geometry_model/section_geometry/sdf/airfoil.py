"""
Convert a discrete airfoil polyline into a robust SDF callable.

Approach
--------
1.  Accept (x, y) point arrays that trace the airfoil contour.
2.  Build an exact polyhedral SDF via ``sdf_polygon`` (winding-number sign).
3.  Expose convenience methods for chord normalisation, camber-line extraction,
    thickness distribution, and NACA 4/5/6-digit generation.
"""

import numpy as np

from .primitives import sdf_polygon


class AirfoilSDF:
    """Implicit SDF representation of an airfoil cross-section.

    Parameters
    ----------
    vertices : ndarray, shape (N, 2)
        Ordered (x, y) coordinates tracing the airfoil contour.
    chord : float, optional
        Reference chord length. If None, inferred as max(x) - min(x).
    """

    def __init__(self, vertices, chord=None):
        verts = np.asarray(vertices, dtype=float)
        if verts.ndim != 2 or verts.shape[1] != 2:
            raise ValueError(
                f"vertices must have shape (N, 2); got array with shape {verts.shape}."
            )
        if len(verts) < 3:
            raise ValueError(
                f"At least 3 vertices are required; got {len(verts)}."
            )
        if np.allclose(verts[0], verts[-1]):
            verts = verts[:-1]
        if len(verts) < 3:
            raise ValueError(
                "Airfoil polygon has fewer than 3 unique vertices after removing duplicated closing point."
            )
        self._verts = verts
        self.chord = float(chord) if chord is not None else float(verts[:, 0].max() - verts[:, 0].min())
        if not np.isfinite(self.chord) or self.chord <= 0.0:
            raise ValueError(f"chord must be a positive finite value; got {self.chord!r}.")
        self._le_idx = int(np.argmin(verts[:, 0]))

    def __call__(self, x, y):
        """Evaluate the SDF at arbitrary coordinates (phi < 0 inside)."""
        return sdf_polygon(x, y, self._verts)

    @classmethod
    def from_dat(cls, filepath, chord=None):
        """Load from a Selig-format .dat file."""
        pts = []
        with open(filepath, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        pts.append([float(parts[0]), float(parts[1])])
                    except ValueError:
                        continue
        if not pts:
            raise ValueError(f"No valid airfoil coordinate rows found in DAT file: {filepath}")
        return cls(np.array(pts), chord=chord)

    @classmethod
    def from_array(cls, xy, chord=None):
        """Construct from an (N, 2) NumPy array."""
        return cls(xy, chord=chord)

    @classmethod
    def from_naca(cls, naca_code, n_points=200, chord=1.0, closed_te=True):
        """Generate a NACA 4-digit airfoil.

        The returned contour polygon *is* the section geometry for this package
        (SG-008): downstream offsets and junctions are exact w.r.t. these
        vertices. Default ``n_points=200`` is ~62 µm from the analytic LE on
        NACA 2412 at 1 m chord — discretisation of the generator, not an error
        to close in the junction path.
        """
        code = str(naca_code).strip().zfill(4)
        if len(code) != 4:
            raise ValueError(f"Expected 4-digit NACA code, got '{naca_code}'.")
        m = int(code[0]) / 100.0
        p = int(code[1]) / 10.0
        tt = int(code[2:]) / 100.0
        from .naca import naca_four_digit_vertices
        verts = naca_four_digit_vertices(m, p, tt, int(n_points), float(chord), closed_te=closed_te)
        return cls(verts, chord=chord)

    @classmethod
    def from_naca_series(
        cls,
        series: int,
        m: float,
        p: float,
        xx: float,
        n_points: int = 200,
        chord: float = 1.0,
        closed_te: bool = True,
    ):
        """Build from precompute spanwise convention (4-, 5-, or 6-series)."""
        from .naca import airfoil_vertices_from_spanwise
        verts = airfoil_vertices_from_spanwise(
            int(series), float(m), float(p), float(xx), int(n_points), float(chord), closed_te=closed_te
        )
        return cls(verts, chord=float(chord))

    @property
    def vertices(self):
        """Airfoil contour vertices, shape (N, 2)."""
        return self._verts.copy()

    @property
    def leading_edge(self):
        """Approximate leading-edge point (min x)."""
        return self._verts[self._le_idx].copy()

    @property
    def trailing_edge(self):
        """Approximate trailing-edge point (mean of first and last vertex)."""
        return 0.5 * (self._verts[0] + self._verts[-1])

    def upper_surface(self):
        return self._verts[:self._le_idx + 1]

    def lower_surface(self):
        return self._verts[self._le_idx:]

    def camber_line(self, n_points=100):
        """Estimate the camber line as the midpoint between upper/lower surfaces."""
        upper = self.upper_surface()
        lower = self.lower_surface()
        upper = upper[np.argsort(upper[:, 0])]
        lower = lower[np.argsort(lower[:, 0])]
        x_min = max(upper[:, 0].min(), lower[:, 0].min())
        x_max = min(upper[:, 0].max(), lower[:, 0].max())
        xc = np.linspace(x_min, x_max, n_points)
        yu = np.interp(xc, upper[:, 0], upper[:, 1])
        yl = np.interp(xc, lower[:, 0], lower[:, 1])
        return xc, 0.5 * (yu + yl)

    def thickness_distribution(self, n_points=100):
        """Chord-normal thickness as a function of chordwise position."""
        upper = self.upper_surface()
        lower = self.lower_surface()
        upper = upper[np.argsort(upper[:, 0])]
        lower = lower[np.argsort(lower[:, 0])]
        x_min = max(upper[:, 0].min(), lower[:, 0].min())
        x_max = min(upper[:, 0].max(), lower[:, 0].max())
        xc = np.linspace(x_min, x_max, n_points)
        yu = np.interp(xc, upper[:, 0], upper[:, 1])
        yl = np.interp(xc, lower[:, 0], lower[:, 1])
        return xc, yu - yl

    def normalise(self):
        """Return a new AirfoilSDF with chord = 1, LE at origin."""
        le = self.leading_edge
        verts = (self._verts - le) / self.chord
        return AirfoilSDF(verts, chord=1.0)

    def scale(self, new_chord):
        """Return a new AirfoilSDF scaled to a different chord length."""
        return AirfoilSDF(self._verts * (new_chord / self.chord), chord=new_chord)

    def translate(self, dx, dy):
        """Return a translated copy."""
        return AirfoilSDF(self._verts + np.array([dx, dy]), chord=self.chord)

    def rotate(self, angle):
        """Return a new AirfoilSDF with vertices rotated CCW by angle radians about the origin."""
        c = np.cos(float(angle))
        s = np.sin(float(angle))
        v = self._verts
        xr = c * v[:, 0] - s * v[:, 1]
        yr = s * v[:, 0] + c * v[:, 1]
        return AirfoilSDF(np.column_stack([xr, yr]), chord=self.chord)

    def __repr__(self):
        return (
            f"AirfoilSDF(n_verts={len(self._verts)}, "
            f"chord={self.chord:.4g}, "
            f"LE={self.leading_edge}, "
            f"TE={self.trailing_edge})"
        )
