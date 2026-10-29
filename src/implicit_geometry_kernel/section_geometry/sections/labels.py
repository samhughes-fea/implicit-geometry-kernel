"""Canonical component label enum."""

import enum


class ComponentLabel(str, enum.Enum):
    """Canonical string labels for section subcomponents."""

    OUTER_SKIN = "outer_skin"  # legacy full-annulus alias
    SKIN_SUCTION_SURFACE = "skin_suction_surface"
    SKIN_PRESSURE_SURFACE = "skin_pressure_surface"
    SPAR_CAP = "spar_cap"
    SHEAR_WEB = "shear_web"
    SANDWICH_CORE = "sandwich_core"
    TRAILING_EDGE_INSERT = "trailing_edge_insert"
    LEADING_EDGE_INSERT = "leading_edge_insert"
    SPAR_CAP_SUCTION_SURFACE = "spar_cap_suction_surface"
    SPAR_CAP_PRESSURE_SURFACE = "spar_cap_pressure_surface"
