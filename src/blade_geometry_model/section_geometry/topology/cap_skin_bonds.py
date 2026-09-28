"""Compat alias — use ``junction_preparation.cap_skin_bonds`` (expires 2027-03-31)."""

from __future__ import annotations

import warnings

warnings.warn(
    "Import from "
    "blade_structure_utilities.geometry_utilities.section_geometry."
    "junction_preparation.cap_skin_bonds "
    "(topology.cap_skin_bonds is a reverse-compatibility alias until 2027-03-31).",
    DeprecationWarning,
    stacklevel=2,
)

from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation.cap_skin_bonds import *  # noqa: E402,F401,F403
