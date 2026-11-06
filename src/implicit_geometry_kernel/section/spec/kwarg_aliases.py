"""Historical keyword aliases; canonical resolution lives on ``SectionSpec``.

Re-exports until 0.6.0. Prefer ``SectionSpec.from_kwargs`` for the public
geometry input.
"""

from implicit_geometry_kernel.section.spec.section_spec import (
    SkinGaugeArg,
    WebThicknessArg,
    resolve_geometry_mode,
    resolve_shear_web_thickness,
    resolve_skin_gauges,
    resolve_skin_thickness,
    resolve_spar_cap_thickness,
    resolve_twist_rad,
)

__all__ = [
    "SkinGaugeArg",
    "WebThicknessArg",
    "resolve_skin_gauges",
    "resolve_skin_thickness",
    "resolve_spar_cap_thickness",
    "resolve_shear_web_thickness",
    "resolve_twist_rad",
    "resolve_geometry_mode",
]
