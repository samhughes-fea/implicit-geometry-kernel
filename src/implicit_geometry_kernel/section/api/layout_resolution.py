"""Station layout resolution: drop webs that do not fit, then step the layout down."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from implicit_geometry_kernel.sdf.airfoil import AirfoilSDF
from implicit_geometry_kernel.section.api.results import DroppedWeb, LayoutInput
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
from implicit_geometry_kernel.section.layout.system_layout import SystemLayoutSpec


def _layout_web_fracs(layout: LayoutInput) -> tuple[float, ...]:
    return tuple(float(f) for f in (layout.web_chord_fracs or ()))


def _layout_n_webs(layout: LayoutInput) -> int:
    if isinstance(layout, SectionLayout):
        return int(layout.n_webs)
    return int(layout.n_webs or 0)


def _step_down_section_layout(layout: SectionLayout, n_webs: int) -> SectionLayout:
    """Reduce a keyless parametric layout to ``n_webs`` survivors."""
    from implicit_geometry_kernel.section.layout.cap_placement import (
        CapPlacementRule,
    )

    n = int(n_webs)
    if n < 0:
        raise ValueError(f"n_webs must be >= 0; got {n_webs!r}.")
    fracs = tuple(float(f) for f in layout.web_chord_fracs[:n])
    stations = (
        None
        if layout.web_stations_m is None
        else tuple(float(x) for x in layout.web_stations_m[:n])
    )
    axes = (
        None
        if layout.web_axes_m is None
        else tuple(layout.web_axes_m[:n])
    )
    cap = layout.cap_placement
    if n < 2 and cap.kind == "continuous_box":
        cap = CapPlacementRule.per_web_pads()
    fam = layout.structural_family
    if n < 2 and str(fam or "").upper() == "D":
        fam = "C"
    return SectionLayout(
        web_chord_fracs=fracs,
        orientation=layout.orientation,
        cap_placement=cap,
        catalogue_key=None if n != _layout_n_webs(layout) else layout.catalogue_key,
        web_stations_m=stations,
        web_axes_m=axes,
        structural_family=fam,
        flags=layout.flags,
    )


def _as_section_layout(layout: LayoutInput) -> SectionLayout:
    """Normalise catalogue specs to the parametric layout SoT."""
    from implicit_geometry_kernel.section.assembly.from_layout import (
        as_section_layout,
    )

    return as_section_layout(layout)


def _reattach_catalogue_preset(layout: SectionLayout) -> SectionLayout:
    """Best-effort ``N{family}-{orient}`` preset re-attach after a parametric drop."""
    from implicit_geometry_kernel.section.layout.system_layout import (
        resolve_system_type,
    )

    n = int(layout.n_webs)
    fam = str(layout.structural_family or "").upper()
    if not fam:
        return layout
    orient = "F" if str(layout.orientation) == "flapwise" else "CN"
    key = f"{n}{fam}-{orient}"
    try:
        preset = resolve_system_type(key)
    except KeyError:
        return layout
    reattached = SectionLayout.from_system_layout_spec(preset)
    return SectionLayout(
        web_chord_fracs=layout.web_chord_fracs or reattached.web_chord_fracs,
        orientation=reattached.orientation,
        cap_placement=reattached.cap_placement,
        catalogue_key=reattached.catalogue_key,
        web_stations_m=layout.web_stations_m,
        web_axes_m=layout.web_axes_m,
        structural_family=reattached.structural_family,
        flags=reattached.flags,
    )


def resolve_web_x_positions(
    layout: LayoutInput,
    airfoil: AirfoilSDF,
    web_x_positions_m: Sequence[float] | None = None,
) -> list[float]:
    """Web chordwise positions in metres for one section."""
    if web_x_positions_m is not None:
        return [float(x) for x in web_x_positions_m]
    if isinstance(layout, SectionLayout) and layout.web_stations_m is not None:
        return [float(x) for x in layout.web_stations_m]
    return [float(f) * float(airfoil.chord) for f in _layout_web_fracs(layout)]


@dataclass(frozen=True)
class StationLayoutResolution:
    """What a station actually builds, once webs that do not fit are dropped."""

    layout: SectionLayout
    web_x_positions_m: tuple[float, ...]
    dropped: tuple[DroppedWeb, ...]
    nominal_layout: SectionLayout
    web_y_coords: tuple[tuple[float, float], ...] = ()

    @property
    def stepped_down(self) -> bool:
        eff = self.layout
        nom = self.nominal_layout
        return (
            _layout_n_webs(eff) != _layout_n_webs(nom)
            or eff.structural_family != nom.structural_family
            or eff.cap_placement != nom.cap_placement
            or eff.catalogue_key != nom.catalogue_key
        )


def resolve_station_layout(
    layout: LayoutInput,
    airfoil: AirfoilSDF,
    *,
    skin_thickness: float | None = None,
    web_x_positions_m: Sequence[float] | None = None,
    min_cavity_depth_m: float = 0.0,
    y_search: Any | None = None,
    phi_inner: Any | None = None,
) -> StationLayoutResolution:
    """Resolve the layout a station can actually carry, before building anything."""
    from implicit_geometry_kernel.section.api.kwarg_aliases import (
        resolve_skin_thickness,
    )
    from implicit_geometry_kernel.section.layout.system_layout import (
        InfeasibleSectionGeometry,
    )
    from implicit_geometry_kernel.section.assembly.web_layout import (
        resolve_web_anchors,
    )

    t_skin_v = resolve_skin_thickness(
        skin_thickness=skin_thickness, required=True
    )

    wanted_preset = isinstance(layout, SystemLayoutSpec) or (
        isinstance(layout, SectionLayout) and bool(layout.catalogue_key)
    )
    nominal_sl = _as_section_layout(layout)

    requested = resolve_web_x_positions(layout, airfoil, web_x_positions_m)
    if not requested:
        return StationLayoutResolution(
            nominal_sl, (), (), nominal_sl, web_y_coords=()
        )

    resolved, dropped = resolve_web_anchors(
        airfoil,
        float(t_skin_v),
        requested,
        y_search=y_search,
        phi_inner=phi_inner,
        min_cavity_depth_m=float(min_cavity_depth_m),
        return_dropped=True,
    )
    survivors = tuple(float(a[0]) for a in resolved)
    y_coords = tuple((float(a[1]), float(a[2])) for a in resolved)

    if len(survivors) == 0 and _layout_n_webs(nominal_sl) > 0:
        label = nominal_sl.catalogue_key or "SectionLayout"
        raise InfeasibleSectionGeometry(
            f"{label} cannot anchor a single web at this station. The "
            "drop cascade floors at one web, because a section with no web "
            "has no shear path; 0A/0B are selectable layouts, not degradation "
            "targets. Report this station as over-constrained."
        )

    effective_sl = _step_down_section_layout(nominal_sl, len(survivors))
    if wanted_preset:
        effective_sl = _reattach_catalogue_preset(effective_sl)

    return StationLayoutResolution(
        effective_sl,
        survivors,
        tuple(dropped),
        nominal_sl,
        web_y_coords=y_coords,
    )


__all__ = [
    "StationLayoutResolution",
    "resolve_station_layout",
    "resolve_web_x_positions",
]
