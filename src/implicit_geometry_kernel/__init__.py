"""implicit_geometry_kernel - geometry of blade and wing cross-sections.

The kernel describes the geometry of blade and wing sections and lofts it along the
span. It holds nothing about materials, plies, stiffness, loads or solvers; the
laminate and section models build those on top of it.

Public API
----------
Airfoils   Airfoil (from coordinates, ``Airfoil.naca4``, ``Airfoil.ellipse``)
Layout     Station (span position, chord, twist, pitch axis, offset),
           SkinRegion (chordwise skin bands), ShearWeb (webs at chord fractions),
           SYSTEM_TYPE_KEYS and resolve_system_type (the 40 structural layouts)
Loft       BladeGeometry (a layout lofted between stations),
           BladeGeometry.section_at(z) / section_at -> SectionCut (wall midlines with
           thickness and labels, cell loops, the merged mid-surface, in the blade frame)

Subpackages
-----------
blade         airfoils, stations, lofting and cut sections (the API above)
section       one 2-D section: the 40 layouts, placement, midlines, junctions,
              cell-wall loops, area and second moments
sdf           2-D signed-distance modelling: primitives, CSG, offsets, Numba kernels
frames        geometric coordinate frames and point transforms
integration   quadtree and shoelace area and moment integrals
contract      the hand-off types (ShellMidlineStrip, CellWallLoop)
plotting      Matplotlib figures (needs the ``[plots]`` extra)

Dependencies: numpy, shapely and numba; scipy (``[splines]``) for spline-resampled
airfoils, matplotlib (``[plots]``) for ``plotting``. Nothing else in the suite.
"""

from implicit_geometry_kernel.blade import (
    Airfoil,
    BladeGeometry,
    CellLoop,
    MidSurface,
    SectionCut,
    SectionWall,
    ShearWeb,
    SkinRegion,
    Station,
    section_at,
)
from implicit_geometry_kernel.section.layout.system_layout import SYSTEM_TYPE_KEYS, resolve_system_type

__version__ = "0.2.0"

__all__ = ["Airfoil", "BladeGeometry", "CellLoop", "MidSurface", "SYSTEM_TYPE_KEYS", "SectionCut", "SectionWall",
           "ShearWeb", "SkinRegion", "Station", "resolve_system_type", "section_at"]
