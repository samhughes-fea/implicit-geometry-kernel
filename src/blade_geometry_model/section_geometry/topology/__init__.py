"""Compat alias — import junction preparation from ``junction_preparation``.

Logic-free warning aliases for one migration cycle. Canonical home is
``section_geometry.junction_preparation``. Expires 2027-03-31.
"""

from __future__ import annotations

import warnings

warnings.warn(
    "Import shell junction preparation from "
    "blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation "
    "(topology is a reverse-compatibility alias until 2027-03-31).",
    DeprecationWarning,
    stacklevel=2,
)

from blade_structure_utilities.geometry_utilities.section_geometry.junction_preparation import (  # noqa: E402,F401
    ShellMidlineJunctionPreparationOptions,
    ShellStripTopologySpec,
    build_shell_strip_topology,
    prepare_shell_midline_junctions,
)

__all__ = [
    "ShellMidlineJunctionPreparationOptions",
    "prepare_shell_midline_junctions",
    "ShellStripTopologySpec",
    "build_shell_strip_topology",
]
