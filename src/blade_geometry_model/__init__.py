"""blade_geometry_model - blade and wing geometry, independent of materials and meshes.

Public API overview
-------------------
Airfoils        Airfoil (from coordinates, ``Airfoil.naca4``, ``Airfoil.ellipse``)
Layout          Station (span position, chord, twist, pitch axis, offset),
                SkinRegion (chordwise skin bands), ShearWeb (webs at chord fractions)

Implicit geometry (forked from blade-structure, see FORK.md)
    section_geometry   signed-distance (SDF / CSG) cross-sections: 40 structural
                       layouts, wall midlines with thickness (ShellMidlineStrip) and
                       closed cell-wall loops (CellWallLoop)
    section_contract   the hand-off types those come in
    coordinate_frames  section and shell frame conventions
    area_2d, shoelace_integration   area integration

The subpackages are imported on their own (``from blade_geometry_model.section_geometry
import ...``). This package depends on numpy, shapely and numba (scipy, optionally,
for spline-resampled airfoils), and on nothing else in the suite.
"""

from .blade_geometry import Airfoil, ShearWeb, SkinRegion, Station

__version__ = "0.1.0"

__all__ = ["Airfoil", "ShearWeb", "SkinRegion", "Station"]
