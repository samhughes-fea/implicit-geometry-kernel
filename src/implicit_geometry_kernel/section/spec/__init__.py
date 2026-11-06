"""Section input types: ``SectionSpec`` (geometry) and ``SectionOptions`` (compute)."""

from implicit_geometry_kernel.section.spec.options import SectionOptions
from implicit_geometry_kernel.section.spec.section_spec import Gauges, SectionSpec

__all__ = ["Gauges", "SectionOptions", "SectionSpec"]
