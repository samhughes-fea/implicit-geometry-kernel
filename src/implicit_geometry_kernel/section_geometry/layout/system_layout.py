"""Deprecated alias of ``implicit_geometry_kernel.section.layout.system_layout`` (removed in 0.3.0)."""

import warnings

from implicit_geometry_kernel.section.layout.system_layout import *  # noqa: F401,F403
from implicit_geometry_kernel.section.layout.system_layout import SYSTEM_TYPE_KEYS, __all__  # noqa: F401

warnings.warn(
    "implicit_geometry_kernel.section_geometry.layout.system_layout is deprecated and will be removed "
    "in 0.3.0; import SYSTEM_TYPE_KEYS from implicit_geometry_kernel",
    DeprecationWarning,
    stacklevel=2,
)
