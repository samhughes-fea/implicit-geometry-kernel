"""Deprecated alias of ``implicit_geometry_kernel.blade`` (removed in 0.3.0).

Import ``Airfoil``, ``Station``, ``SkinRegion`` and ``ShearWeb`` from
``implicit_geometry_kernel`` instead.
"""

import warnings

from implicit_geometry_kernel.blade.airfoil import Airfoil
from implicit_geometry_kernel.blade.stations import ShearWeb, SkinRegion, Station

warnings.warn(
    "implicit_geometry_kernel.blade_geometry is deprecated and will be removed in 0.3.0; "
    "import Airfoil, Station, SkinRegion and ShearWeb from implicit_geometry_kernel",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["Airfoil", "ShearWeb", "SkinRegion", "Station"]
