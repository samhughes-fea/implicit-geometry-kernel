"""Section factory: build section views from layout specs for property reporting."""
from __future__ import annotations

from blade_geometry_model.section_geometry.api.kwarg_aliases import (
    resolve_shear_web_thickness,
    resolve_skin_thickness,
    resolve_spar_cap_thickness,
    resolve_twist_rad,
)
from blade_geometry_model.section_geometry.layout.section_layout import (
    SectionLayout,
)
from blade_geometry_model.section_geometry.layout.system_layout import (
    SYSTEM_TYPE_KEYS,
    SystemLayoutSpec,
    resolve_system_type,
)
from blade_geometry_model.section_geometry.midlines.section_midline_host import (
    SectionMidlineHost,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)
from blade_geometry_model.section_geometry.sections.from_layout import (
    multicell_from_layout,
)

LayoutInput = SystemLayoutSpec | SectionLayout

# Nominal skin thickness passed to section views (pre-clamp).
_DEFAULT_SECTION_SKIN_THICKNESS_M = 0.003


def build_section_view(
    airfoil_sdf: AirfoilSDF,
    layout: LayoutInput,
    *,
    twist_rad: float | None = None,
    twist_angle_rad: float | None = None,
    twist_angle: float | None = None,
    skin_thickness: float | None = None,
    spar_cap_thickness: float | None = None,
    shear_web_thickness: float | None = None,
    web_thickness: float | None = None,
    cap_height: float | None = None,
) -> SectionMidlineHost:
    """Return a section object usable by :class:`SectionPropertiesReport`.

    Accepts catalogue :class:`SystemLayoutSpec` or parametric
    :class:`SectionLayout`. Always assembles via :class:`MultiCellSection`
    (including zero-web ``0A`` / ``0B``). Same host path as
    :func:`build_section_with_thicknesses`.

    Prefer ``skin_thickness`` / ``spar_cap_thickness`` / ``shear_web_thickness``
    and ``twist_rad``. Secondary aliases ``twist_angle_rad`` / ``cap_height`` /
    ``web_thickness`` remain accepted.
    """
    skin = resolve_skin_thickness(
        skin_thickness=skin_thickness,
        default=_DEFAULT_SECTION_SKIN_THICKNESS_M,
    )
    web = resolve_shear_web_thickness(
        shear_web_thickness=shear_web_thickness,
        web_thickness=web_thickness,
        default=0.004,
    )
    spar = resolve_spar_cap_thickness(
        spar_cap_thickness=spar_cap_thickness,
        cap_height=cap_height,
        default=0.012,
    )
    twist = resolve_twist_rad(
        twist_rad=twist_rad,
        twist_angle_rad=twist_angle_rad,
        twist_angle=twist_angle,
    )

    return multicell_from_layout(
        airfoil_sdf,
        layout,
        twist_rad=twist,
        skin_thickness=skin,
        spar_cap_thickness=spar,
        shear_web_thickness=web,
    )


__all__ = [
    "SYSTEM_TYPE_KEYS",
    "build_section_view",
    "resolve_system_type",
    "_DEFAULT_SECTION_SKIN_THICKNESS_M",
]
