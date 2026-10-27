"""blade_geometry_model - blade and wing geometry, independent of materials and meshes.

Public API overview
-------------------
Airfoils        Airfoil (from coordinates, ``Airfoil.naca4``, ``Airfoil.ellipse``)
Layout          Station (span position, chord, twist, pitch axis, offset),
                SkinRegion (chordwise skin bands), ShearWeb (webs at chord fractions)

This package depends on numpy (and scipy, optionally, for spline-resampled
airfoils) and on nothing else in the suite.
"""

from .blade_geometry import Airfoil, ShearWeb, SkinRegion, Station

__version__ = "0.1.0"

__all__ = ["Airfoil", "ShearWeb", "SkinRegion", "Station"]
