"""
interface.export
================
Data export utilities for section geometry.

Thin re-export surface over :mod:`properties_report` and :mod:`export_section_json`.
Standalone medial-axis CSV export was removed from the public API (G3).
"""

from .properties_report import SectionPropertiesReport
from .export_section_json import (
    _multicell_host,
    _shared_grid_phi,
    export_section_json,
)

__all__ = [
    "SectionPropertiesReport",
    "export_section_json",
    "_shared_grid_phi",
    "_multicell_host",
]
