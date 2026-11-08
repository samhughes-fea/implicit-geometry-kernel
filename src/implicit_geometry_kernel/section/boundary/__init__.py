"""The section's one gateway to kernel polygons and GEOS (see :mod:`.polygons`).

Every section module that needs polygons, offsets or Boolean clean-up imports them
from ``implicit_geometry_kernel.section.boundary``. Overlay is GEOS, not Shapely.
"""

from implicit_geometry_kernel.section.boundary.polygons import *  # noqa: F401,F403
from implicit_geometry_kernel.section.boundary.polygons import __all__  # noqa: F401
