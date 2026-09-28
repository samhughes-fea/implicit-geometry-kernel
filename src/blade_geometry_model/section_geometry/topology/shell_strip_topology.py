"""Compat alias module — use ``junction_preparation.shell_midline_junction_preparation``.

Expires 2027-03-31.
"""

from __future__ import annotations

import warnings

warnings.warn(
    "Import from "
    "blade_structure_utilities.geometry_utilities.section_geometry."
    "junction_preparation.shell_midline_junction_preparation "
    "(shell_strip_topology is a reverse-compatibility alias until 2027-03-31).",
    DeprecationWarning,
    stacklevel=2,
)

from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import *  # noqa: E402,F401,F403
from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.shell_midline_junction_preparation import (  # noqa: E402,F401
    __all__,
)
