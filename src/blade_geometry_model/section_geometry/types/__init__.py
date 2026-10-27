"""Geometry worker product DTOs / aliases."""

from ..api.worker import SectionGeometryResult, SectionPropertiesMeta
from ..layout.types import SDFCallable
from ..props_emit.export import SectionPropertiesReport

__all__ = [
    "SectionGeometryResult",
    "SectionPropertiesMeta",
    "SDFCallable",
    "SectionPropertiesReport",
]
