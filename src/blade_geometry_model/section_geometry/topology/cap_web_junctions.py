"""Compat alias — use ``junction_preparation.cap_web_junctions`` (expires 2027-03-31)."""

from __future__ import annotations

import warnings

warnings.warn(
    "Import from "
    "blade_geometry_model.section_geometry."
    "junction_preparation.cap_web_junctions "
    "(topology.cap_web_junctions is a reverse-compatibility alias until 2027-03-31).",
    DeprecationWarning,
    stacklevel=2,
)

from blade_geometry_model.section_geometry.junction_preparation.cap_web_junctions import *  # noqa: E402,F401,F403
