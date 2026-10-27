"""Build :class:`MultiCellSection` from parametric :class:`SectionLayout`.

Single production host path for every registered layout key (including
``0A`` / ``0B``). Catalogue :class:`SystemLayoutSpec` is a **preset factory** —
:func:`multicell_from_layout` normalises it to :class:`SectionLayout` before
assemble. Cap placement is :class:`CapPlacementRule`, not a family-letter
assemble switch.
"""

from __future__ import annotations

from collections.abc import Sequence

from blade_geometry_model.section_geometry.api.kwarg_aliases import (
    resolve_shear_web_thickness,
    resolve_skin_gauges,
    resolve_spar_cap_thickness,
    resolve_twist_rad,
)
from blade_geometry_model.section_geometry.layout.cap_placement import (
    CapPlacementRule,
)
from blade_geometry_model.section_geometry.layout.section_layout import (
    SectionLayout,
)
from blade_geometry_model.section_geometry.layout.system_layout import (
    SystemLayoutSpec,
)
from blade_geometry_model.section_geometry.midlines.section_midline_host import (
    SectionMidlineHost,
)
from blade_geometry_model.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)

# Matches the former OuterSkinOnlySection band half-width for zero-web family B.
_ZERO_WEB_CAP_HALF_BAND_FRAC = 0.10

LayoutInput = SystemLayoutSpec | SectionLayout
"""Public entry may accept a catalogue preset; assemble normalises to SectionLayout."""


def fixed_cap_anchor_from_layout(layout: LayoutInput) -> str:
    """Map layout spar-cap positioning onto MultiCell ``fixed_cap_anchor``."""
    sl = (
        layout
        if isinstance(layout, SectionLayout)
        else SectionLayout.from_system_layout_spec(layout)
    )
    if sl.cap_placement.kind != "fixed_band":
        return "pitching"
    anchor = sl.cap_placement.anchor or "max_thickness"
    return "max_thickness" if anchor == "max_thickness" else "pitching"


def _family_from_cap_rule(rule: CapPlacementRule) -> str:
    """Legacy family letter for host metadata / catalogue preset re-attach."""
    if rule.kind == "none":
        return "A"
    if rule.kind == "fixed_band":
        return "B"
    if rule.kind == "per_web_pads":
        return "C"
    if rule.kind == "continuous_box":
        return "D"
    raise ValueError(f"Unknown CapPlacementRule.kind {rule.kind!r}")


def as_section_layout(layout: LayoutInput) -> SectionLayout:
    """Normalise catalogue specs to the parametric layout SoT."""
    if isinstance(layout, SectionLayout):
        return layout
    return SectionLayout.from_system_layout_spec(layout)


def multicell_from_section_layout(
    airfoil: AirfoilSDF,
    layout: SectionLayout,
    *,
    twist_rad: float | None = None,
    twist_angle_rad: float | None = None,
    twist_angle: float | None = None,
    skin_thickness: float | None = None,
    spar_cap_thickness: float | None = None,
    shear_web_thickness: float | None = None,
    web_thickness: float | None = None,
    cap_height: float | None = None,
    web_x_positions_m: Sequence[float] | None = None,
    web_y_coords: Sequence[tuple[float, float]] | None = None,
    stage_sink: dict[str, float] | None = None,
    spar_cap_ply_pitch_m: float | None = None,
) -> SectionMidlineHost:
    """Assemble from parametric :class:`SectionLayout` (catalogue key optional)."""
    from blade_geometry_model.section_geometry.sections.multicell import (
        MultiCellSection,
    )

    skin = resolve_skin_gauges(skin_thickness=skin_thickness, required=True)
    spar = resolve_spar_cap_thickness(
        spar_cap_thickness=spar_cap_thickness,
        cap_height=cap_height,
        required=True,
    )
    web = resolve_shear_web_thickness(
        shear_web_thickness=shear_web_thickness,
        web_thickness=web_thickness,
        required=True,
    )
    twist = resolve_twist_rad(
        twist_rad=twist_rad,
        twist_angle_rad=twist_angle_rad,
        twist_angle=twist_angle,
    )

    if web_x_positions_m is not None:
        xs = [float(x) for x in web_x_positions_m]
    elif layout.web_stations_m is not None:
        xs = [float(x) for x in layout.web_stations_m]
    else:
        xs = [float(f) * float(airfoil.chord) for f in layout.web_chord_fracs]

    web_endpoints = None
    if layout.web_axes_m is not None:
        web_endpoints = [
            (float(top[0]), float(top[1]), float(bot[0]), float(bot[1]))
            for top, bot in layout.web_axes_m
        ]
        xs = [0.5 * (e[0] + e[2]) for e in web_endpoints]

    align = "flapwise" if layout.orientation == "flapwise" else "chord_normal"
    y_coords = None
    if web_y_coords is not None:
        y_coords = [(float(yt), float(yb)) for yt, yb in web_y_coords]
        if len(y_coords) != len(xs):
            raise ValueError(
                f"web_y_coords length ({len(y_coords)}) must match "
                f"web_x_positions length ({len(xs)})."
            )

    cap_rule = layout.cap_placement
    fam = str(layout.structural_family or _family_from_cap_rule(cap_rule)).upper()
    fhw: float | None = cap_rule.half_width
    if len(xs) == 0 and cap_rule.kind == "fixed_band" and fhw is None:
        fhw = _ZERO_WEB_CAP_HALF_BAND_FRAC * float(airfoil.chord)
        cap_rule = CapPlacementRule.fixed_band(
            anchor=cap_rule.anchor or "max_thickness",
            half_width=fhw,
        )

    return MultiCellSection(
        airfoil,
        web_x_positions=list(xs),
        web_thickness=web,
        web_alignment=align,
        cap_height=float(spar),
        skin_thickness=skin,
        twist_angle=float(twist),
        core_enabled=True,
        structural_family=fam,
        fixed_cap_anchor=fixed_cap_anchor_from_layout(layout),
        fixed_cap_chord_half_width=fhw,
        web_y_coords=y_coords,
        stage_sink=stage_sink,
        cap_placement=cap_rule,
        web_endpoints=web_endpoints,
        spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
    )


def multicell_from_layout(
    airfoil: AirfoilSDF,
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
    web_x_positions_m: Sequence[float] | None = None,
    web_y_coords: Sequence[tuple[float, float]] | None = None,
    stage_sink: dict[str, float] | None = None,
    spar_cap_ply_pitch_m: float | None = None,
) -> SectionMidlineHost:
    """Assemble a sized :class:`MultiCellSection`.

    Catalogue :class:`SystemLayoutSpec` is converted to :class:`SectionLayout`
    (the assemble SoT) before host construction.
    """
    return multicell_from_section_layout(
        airfoil,
        as_section_layout(layout),
        twist_rad=twist_rad,
        twist_angle_rad=twist_angle_rad,
        twist_angle=twist_angle,
        skin_thickness=skin_thickness,
        spar_cap_thickness=spar_cap_thickness,
        shear_web_thickness=shear_web_thickness,
        web_thickness=web_thickness,
        cap_height=cap_height,
        web_x_positions_m=web_x_positions_m,
        web_y_coords=web_y_coords,
        stage_sink=stage_sink,
        spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
    )


__all__ = [
    "as_section_layout",
    "fixed_cap_anchor_from_layout",
    "multicell_from_layout",
    "multicell_from_section_layout",
]
