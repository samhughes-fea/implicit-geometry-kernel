"""Use-case-agnostic section geometry worker facade.

Mirrors the FE worker contract: no sizing/fatigue modes; callers choose
output gates (compute-cost knobs) and interpret results in orchestration.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Collection, Iterable, Literal, Sequence

from blade_structure_utilities.geometry_utilities.section_contract.cell_wall_loop import (
    CellWallLoop,
)
from blade_structure_utilities.geometry_utilities.section_contract.shell_midline_strip import (
    ShellMidlineStrip,
)
from blade_structure_utilities.geometry_utilities.section_geometry.layout.section_layout import (
    SectionLayout,
)
from blade_structure_utilities.geometry_utilities.section_geometry.layout.system_layout import (
    SystemLayoutSpec,
)
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.section_midline_host import (
    SectionMidlineHost,
)
from blade_structure_utilities.geometry_utilities.section_geometry.props_emit.props_compute import (
    PropsBackendArg,
)
from blade_structure_utilities.geometry_utilities.section_geometry.sdf.airfoil import (
    AirfoilSDF,
)
from blade_structure_utilities.mesh_utilities.area_2d.quadtree import QuadTreeMesh
from blade_structure_utilities.mesh_utilities.area_2d.spec import QuadtreeSpec

LayoutInput = SystemLayoutSpec | SectionLayout


def _compat_kwargs_from_zone_gauges(zone_gauges_m: Mapping[str, float]) -> dict[str, Any]:
    """One-way projection: group zone-id thicknesses into geometry wall gauges.

    Geometry stays gauge-agnostic and does not import ``section_gauge``.
    Callers that already hold a realization should pass ``zone_gauges_m``
    rather than collapsing to component maxima first.
    """
    items = {str(k): float(v) for k, v in dict(zone_gauges_m).items()}
    skins = [v for k, v in items.items() if k.startswith("skin:")]
    caps = [v for k, v in items.items() if k.startswith("cap:") or k.startswith("spar:")]
    webs = sorted(
        ((k, v) for k, v in items.items() if k.startswith("web:")),
        key=lambda kv: kv[0],
    )
    ss = float(items.get("skin:suction", max(skins) if skins else 0.0))
    ps = float(items.get("skin:pressure", ss))
    web_vals = [float(v) for _, v in webs]
    skin_kw: float | tuple[float, float] = ss if ss == ps else (ss, ps)
    if not web_vals:
        web_kw: float | list[float] = 0.0
    elif all(v == web_vals[0] for v in web_vals):
        web_kw = web_vals[0]
    else:
        web_kw = web_vals
    return {
        "skin_thickness": skin_kw,
        "spar_cap_thickness": float(max(caps) if caps else 0.0),
        "shear_web_thickness": web_kw,
    }


def _layout_web_fracs(layout: LayoutInput) -> tuple[float, ...]:
    return tuple(float(f) for f in (layout.web_chord_fracs or ()))


def _layout_n_webs(layout: LayoutInput) -> int:
    if isinstance(layout, SectionLayout):
        return int(layout.n_webs)
    return int(layout.n_webs or 0)


def _as_system_layout_if_possible(layout: LayoutInput) -> SystemLayoutSpec | None:
    if isinstance(layout, SystemLayoutSpec):
        return layout
    if layout.catalogue_key:
        return layout.to_system_layout_spec()
    return None


def _step_down_section_layout(layout: SectionLayout, n_webs: int) -> SectionLayout:
    """Reduce a keyless parametric layout to ``n_webs`` survivors."""
    from blade_structure_utilities.geometry_utilities.section_geometry.layout.cap_placement import (
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
    # Mirror catalogue D→C: a continuous box needs ≥2 webs.
    if n < 2 and cap.kind == "continuous_box":
        cap = CapPlacementRule.per_web_pads()
    fam = layout.structural_family
    if n < 2 and str(fam or "").upper() == "D":
        fam = "C"
    return SectionLayout(
        web_chord_fracs=fracs,
        orientation=layout.orientation,
        cap_placement=cap,
        # Survivors are parametric until a catalogue preset is re-attached.
        catalogue_key=None if n != _layout_n_webs(layout) else layout.catalogue_key,
        web_stations_m=stations,
        web_axes_m=axes,
        structural_family=fam,
        flags=layout.flags,
    )


OutputGate = Literal["properties", "midlines", "cell_wall_loops"]
ALL_OUTPUT_GATES: frozenset[str] = frozenset(
    {"properties", "midlines", "cell_wall_loops"}
)


def normalize_output_gates(outputs: Collection[str]) -> frozenset[str]:
    """Validate ``outputs`` gate names.

    Unknown names raise ``ValueError``.
    """
    raw = {str(o) for o in outputs}
    unknown = raw - ALL_OUTPUT_GATES
    if unknown:
        raise ValueError(f"Unknown section geometry outputs: {sorted(unknown)}")
    return frozenset(raw)


@dataclass(frozen=True)
class DroppedWeb:
    """Web dropped while resolving spanwise anchors at a station."""

    x_m: float
    reason: str


@dataclass(frozen=True)
class SectionPropertiesMeta:
    """Frozen bag of property-backend metadata from :func:`compute_section_properties`."""

    backend: str
    grid: dict[str, int] | None = None
    quadtree_spec: dict[str, Any] | None = None
    constructive_labels: tuple[str, ...] | None = None
    wall_solid_labels: tuple[str, ...] | None = None
    residual_labels: tuple[str, ...] | None = None
    residual_maturity: dict[str, str] | None = None
    constructive_s: float | None = None
    polygon_count: int | None = None
    residual_quadtree_s: float | None = None
    quadtree_stats: dict[str, Any] | None = None
    constructive_polygons: dict[str, str] | None = None
    profiling: dict[str, Any] | None = None

    @classmethod
    def from_mapping(cls, meta: Mapping[str, Any]) -> SectionPropertiesMeta:
        def _str_tuple(key: str) -> tuple[str, ...] | None:
            raw = meta.get(key)
            if raw is None:
                return None
            return tuple(str(x) for x in raw)

        grid_raw = meta.get("grid")
        grid: dict[str, int] | None = None
        if isinstance(grid_raw, Mapping):
            grid = {str(k): int(v) for k, v in grid_raw.items()}

        maturity_raw = meta.get("residual_maturity")
        maturity: dict[str, str] | None = None
        if isinstance(maturity_raw, Mapping):
            maturity = {str(k): str(v) for k, v in maturity_raw.items()}

        polys_raw = meta.get("constructive_polygons")
        polys: dict[str, str] | None = None
        if isinstance(polys_raw, Mapping):
            polys = {str(k): str(v) for k, v in polys_raw.items()}

        qt_stats = meta.get("quadtree_stats")
        profiling = meta.get("profiling")
        return cls(
            backend=str(meta.get("backend", "")),
            grid=grid,
            quadtree_spec=(
                dict(meta["quadtree_spec"])
                if isinstance(meta.get("quadtree_spec"), Mapping)
                else None
            ),
            constructive_labels=_str_tuple("constructive_labels"),
            wall_solid_labels=_str_tuple("wall_solid_labels"),
            residual_labels=_str_tuple("residual_labels"),
            residual_maturity=maturity,
            constructive_s=(
                float(meta["constructive_s"])
                if meta.get("constructive_s") is not None
                else None
            ),
            polygon_count=(
                int(meta["polygon_count"])
                if meta.get("polygon_count") is not None
                else None
            ),
            residual_quadtree_s=(
                float(meta["residual_quadtree_s"])
                if meta.get("residual_quadtree_s") is not None
                else None
            ),
            quadtree_stats=dict(qt_stats) if isinstance(qt_stats, Mapping) else None,
            constructive_polygons=polys,
            profiling=dict(profiling) if isinstance(profiling, Mapping) else None,
        )


@dataclass(frozen=True)
class SectionGeometryResult:
    """Stable worker output; gated carriers are ``None`` when not requested.

    When midlines are gated, also fills shell-handoff scalars (``chord_m``,
    ``twist_rad``, ``leading_edge_twisted_pitch_axis_blade_frame``, ``trailing_edge_twisted_pitch_axis_blade_frame``) so FE can build
    ``ShellMeshInputs`` without re-evaluating SDF.

    ``section`` implements :class:`SectionMidlineHost` for midline / mesh export.

    ``properties_map`` (when gated) carries area / first / second moments in the
    **solids' authoring frame as built** — with host twist that is LE-origin
    rotated-NACA, **not** frame **B** (`twisted_pitch_axis_blade_frame`). Do not
    silently re-twist the map to match emit midlines.

    ``quadtree_mesh`` is **G-PROPS debug only** (populated when ``retain_mesh``
    is true). It is not a shell-mesh product; density / Q4 / +x3 stay in
    ``section_mesh``.

    ``stage_timings`` is a small ``{stage: seconds}`` dict populated when
    ``SECTION_GEOMETRY_PROFILE=1`` (assemble / assemble_* / properties /
    properties_* / emit_* / midlines / topology / topology_* / cell_wall_loops).
    """

    section: SectionMidlineHost
    prep_key: tuple[Any, ...]
    properties_map: dict[str, dict[str, float]] | None = None
    properties_meta: SectionPropertiesMeta | None = None
    midline_strips: tuple[ShellMidlineStrip, ...] | None = None
    cell_wall_loops: tuple[CellWallLoop, ...] | None = None
    quadtree_mesh: QuadTreeMesh | None = None
    chord_m: float | None = None
    twist_rad: float | None = None
    leading_edge_twisted_pitch_axis_blade_frame: tuple[float, float] | None = None
    trailing_edge_twisted_pitch_axis_blade_frame: tuple[float, float] | None = None
    stage_timings: dict[str, float] | None = None


def _thickness_cache_token(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return tuple(round(float(x), 12) for x in value)
    return round(float(value), 12)


def section_geometry_prep_key(
    layout: LayoutInput,
    airfoil: AirfoilSDF,
    *,
    skin_thickness: float | None = None,
    spar_cap_thickness: float | None = None,
    shear_web_thickness: float | None = None,
    twist_rad: float = 0.0,
    props_backend: str | None = None,
    outputs: Collection[str] = ALL_OUTPUT_GATES,
    web_x_positions_m: Sequence[float] | None = None,
    zone_gauges_m: Mapping[str, float] | None = None,
    placement_digest_by_kind: Mapping[str, str] | None = None,
) -> tuple[Any, ...]:
    """Stable cache key for a geometry worker call (layout + airfoil + DVs + gates).

    ``web_x_positions_m`` participates: anchored positions change the built
    section without changing the layout key or any thickness, so omitting it
    would let two different sections collide on one cache entry.
    ``placement_digest_by_kind`` participates so a surface-aware policy change
    invalidates geometry caches (WP3 D3.3 / WP5 D5.3).
    """
    from blade_structure_utilities.geometry_utilities.section_geometry.api.kwarg_aliases import (
        resolve_shear_web_thickness,
        resolve_skin_gauges,
        resolve_spar_cap_thickness,
    )

    if zone_gauges_m:
        derived = _compat_kwargs_from_zone_gauges(zone_gauges_m)
        skin_thickness = derived["skin_thickness"]
        spar_cap_thickness = derived["spar_cap_thickness"]
        shear_web_thickness = derived["shear_web_thickness"]

    t_skin_gauges = resolve_skin_gauges(
        skin_thickness=skin_thickness, required=True
    )
    t_spar_v = resolve_spar_cap_thickness(
        spar_cap_thickness=spar_cap_thickness, required=True
    )
    t_web_v = resolve_shear_web_thickness(
        shear_web_thickness=shear_web_thickness, required=True
    )
    chord = float(airfoil.chord or 0.0)
    naca = getattr(airfoil, "naca_code", None) or getattr(airfoil, "code", None)
    airfoil_id: Any
    if naca is not None:
        airfoil_id = ("naca", str(naca), chord)
    else:
        verts = airfoil.vertices
        if verts is not None:
            import numpy as np

            arr = np.asarray(verts, dtype=float)
            airfoil_id = ("verts", chord, int(arr.shape[0]), float(arr[0, 0]), float(arr[0, 1]))
        else:
            airfoil_id = ("airfoil", chord, id(airfoil))

    # Always key on parametric SectionLayout (catalogue specs are presets).
    layout_sot = layout if isinstance(layout, SectionLayout) else _as_section_layout(layout)
    cap = layout_sot.cap_placement
    layout_key = (
        "section_layout",
        layout_sot.catalogue_key,
        _layout_web_fracs(layout_sot),
        None
        if layout_sot.web_stations_m is None
        else tuple(float(x) for x in layout_sot.web_stations_m),
        str(layout_sot.orientation),
        str(cap.kind),
        str(cap.anchor),
        None if cap.half_width is None else float(cap.half_width),
        str(layout_sot.structural_family),
    )
    out_key = tuple(sorted(normalize_output_gates(outputs)))
    web_key = (
        None if web_x_positions_m is None
        else tuple(round(float(x), 12) for x in web_x_positions_m)
    )
    zone_key = (
        None
        if not zone_gauges_m
        else tuple(sorted((str(k), round(float(v), 12)) for k, v in zone_gauges_m.items()))
    )
    placement_key = (
        None
        if not placement_digest_by_kind
        else tuple(sorted((str(k), str(v)) for k, v in placement_digest_by_kind.items()))
    )
    return (
        layout_key,
        airfoil_id,
        _thickness_cache_token(t_skin_gauges),
        _thickness_cache_token(t_spar_v),
        _thickness_cache_token(t_web_v),
        float(twist_rad),
        props_backend,
        out_key,
        web_key,
        zone_key,
        placement_key,
    )


def _as_section_layout(layout: LayoutInput) -> SectionLayout:
    """Normalise catalogue specs to the parametric layout SoT."""
    from blade_structure_utilities.geometry_utilities.section_geometry.sections.from_layout import (
        as_section_layout,
    )

    return as_section_layout(layout)


def _reattach_catalogue_preset(layout: SectionLayout) -> SectionLayout:
    """Best-effort ``N{family}-{orient}`` preset re-attach after a parametric drop.

    Returns a :class:`SectionLayout` whose ``catalogue_key`` matches the
    registry entry when one exists; otherwise returns ``layout`` unchanged.
    """
    from blade_structure_utilities.geometry_utilities.section_geometry.layout.system_layout import (
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
    # Preserve survivor stations / axes from the parametric drop; only
    # re-attach the catalogue key and registry cap rule / orientation.
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
    """Web chordwise positions in metres for one section.

    ``web_x_positions_m`` is the **spanwise-anchored** form: absolute offsets
    from the leading edge, resolved once for the blade (conventionally at the
    root) and reused unchanged at every station, so the web centroid does not
    migrate as the chord tapers. Webs that no longer fit the local cavity are
    dropped downstream by ``resolve_web_anchors``, not slid forward.

    Omitting it falls back to the legacy per-station form ``frac × local_chord``,
    under which every web shrinks with the section and none is ever dropped.
    Absolute ``web_stations_m`` on a :class:`SectionLayout` win over fractions
    when anchors are not supplied.
    """
    if web_x_positions_m is not None:
        return [float(x) for x in web_x_positions_m]
    if isinstance(layout, SectionLayout) and layout.web_stations_m is not None:
        return [float(x) for x in layout.web_stations_m]
    return [float(f) * float(airfoil.chord) for f in _layout_web_fracs(layout)]


@dataclass(frozen=True)
class StationLayoutResolution:
    """What a station actually builds, once webs that do not fit are dropped.

    ``layout`` and ``nominal_layout`` are always :class:`SectionLayout`
    (catalogue :class:`SystemLayoutSpec` is normalised at resolve). For keyed
    presets, survivors re-attach the matching ``N{family}-{orient}`` catalogue
    key when one exists (``2D`` may hand over to ``1C``). For keyless layouts,
    survivors stay parametric: ``N`` webs + :class:`CapPlacementRule`
    (continuous box demotes to per-web pads when ``N < 2``).
    ``web_x_positions_m`` are the survivors, still at their anchored offsets.
    ``web_y_coords`` are the matching ``(y_top, y_bot)`` pairs from the same
    :func:`resolve_web_anchors` call — pass them into
    :func:`build_section_geometry` so ``MultiCellSection`` does not re-resolve.
    ``dropped`` records each lost web with a reason.
    """

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
    """Resolve the layout a station can actually carry, before building anything.

    One drop interpreter: survivors → parametric ``N`` webs +
    :class:`CapPlacementRule` via :func:`_step_down_section_layout`. Catalogue
    keys are presets — when the nominal layout carried a key, the effective
    result is re-attached to the matching ``N{family}-{orient}`` registry entry
    when one exists (so callers still see ``1C-F`` after a ``5D-F`` drop).

    Raises :class:`InfeasibleSectionGeometry` when no web can be anchored.
    """
    from blade_structure_utilities.geometry_utilities.section_geometry.api.kwarg_aliases import (
        resolve_skin_thickness,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.layout.system_layout import (
        InfeasibleSectionGeometry,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.sections.web_layout import (
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
        # Zero-web layouts (0A / 0B) are selected, never stepped into.
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


def build_section_with_thicknesses(
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
    """Build sized implicit section geometry for a layout and free DVs.

    Accepts catalogue :class:`SystemLayoutSpec` or parametric
    :class:`SectionLayout`. Every catalogue key (including ``0A`` / ``0B``)
    builds a :class:`MultiCellSection`. Prefer long thickness / twist names;
    ``cap_height`` / ``web_thickness`` remain accepted aliases.
    """
    from blade_structure_utilities.geometry_utilities.section_geometry.api.kwarg_aliases import (
        resolve_shear_web_thickness,
        resolve_skin_gauges,
        resolve_spar_cap_thickness,
        resolve_twist_rad,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.sections.from_layout import (
        multicell_from_layout,
    )

    return multicell_from_layout(
        airfoil,
        layout,
        twist_rad=resolve_twist_rad(
            twist_rad=twist_rad,
            twist_angle_rad=twist_angle_rad,
            twist_angle=twist_angle,
        ),
        skin_thickness=resolve_skin_gauges(
            skin_thickness=skin_thickness, required=True
        ),
        spar_cap_thickness=resolve_spar_cap_thickness(
            spar_cap_thickness=spar_cap_thickness,
            cap_height=cap_height,
            required=True,
        ),
        shear_web_thickness=resolve_shear_web_thickness(
            shear_web_thickness=shear_web_thickness,
            web_thickness=web_thickness,
            required=True,
        ),
        web_x_positions_m=web_x_positions_m,
        web_y_coords=web_y_coords,
        stage_sink=stage_sink,
        spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
    )


def build_section_geometry(
    layout: LayoutInput,
    airfoil: AirfoilSDF,
    *,
    skin_thickness: float | None = None,
    spar_cap_thickness: float | None = None,
    shear_web_thickness: float | None = None,
    twist_rad: float | None = None,
    twist_angle_rad: float | None = None,
    twist_angle: float | None = None,
    outputs: Iterable[str] = ("properties", "midlines", "cell_wall_loops"),
    props_backend: PropsBackendArg | None = None,
    n_web_samples: int = 20,
    n_cap_samples: int = 80,
    retain_mesh: bool = False,
    quadtree_spec: QuadtreeSpec | None = None,
    airfoil_for_spec: AirfoilSDF | None = None,
    web_x_positions_m: Sequence[float] | None = None,
    web_y_coords: Sequence[tuple[float, float]] | None = None,
    spar_cap_ply_pitch_m: float | None = None,
    zone_gauges_m: Mapping[str, float] | None = None,
) -> SectionGeometryResult:
    """Geometry worker entry: assemble section and emit gated carriers.

    ``layout`` may be a catalogue :class:`SystemLayoutSpec` or a parametric
    :class:`SectionLayout` (catalogue key optional).

    Parameters
    ----------
    outputs
        Compute-cost gates (not use-case modes): ``properties``, ``midlines``,
        ``cell_wall_loops``. ``cell_wall_loops`` implies midline extraction.
    web_y_coords
        Optional ``(y_top, y_bot)`` pairs aligned with ``web_x_positions_m``.
        When set (e.g. from :func:`resolve_station_layout`), skips a second
        :func:`resolve_web_anchors` inside ``MultiCellSection``.
    """
    from blade_structure_utilities.geometry_utilities.section_geometry.api.kwarg_aliases import (
        resolve_shear_web_thickness,
        resolve_skin_gauges,
        resolve_spar_cap_thickness,
        resolve_twist_rad,
    )

    if zone_gauges_m:
        if any(
            v is not None
            for v in (skin_thickness, spar_cap_thickness, shear_web_thickness)
        ):
            raise ValueError(
                "geometry received both zone_gauges_m and component thickness kwargs"
            )
        derived = _compat_kwargs_from_zone_gauges(zone_gauges_m)
        skin_thickness = derived["skin_thickness"]
        spar_cap_thickness = derived["spar_cap_thickness"]
        shear_web_thickness = derived["shear_web_thickness"]

    skin = resolve_skin_gauges(
        skin_thickness=skin_thickness, required=True
    )
    spar = resolve_spar_cap_thickness(
        spar_cap_thickness=spar_cap_thickness, required=True
    )
    web = resolve_shear_web_thickness(
        shear_web_thickness=shear_web_thickness, required=True
    )
    twist_rad = resolve_twist_rad(
        twist_rad=twist_rad,
        twist_angle_rad=twist_angle_rad,
        twist_angle=twist_angle,
    )

    out_set = set(normalize_output_gates(outputs))

    need_props = "properties" in out_set
    need_midlines = "midlines" in out_set or "cell_wall_loops" in out_set
    need_void = "cell_wall_loops" in out_set

    from blade_structure_utilities.geometry_utilities.section_geometry.props_emit.props_backend import (
        props_backend_from_env,
    )
    from blade_structure_utilities.geometry_utilities.section_geometry.sdf.profiling import (
        profiling_enabled,
        timed_stage,
    )

    backend = props_backend if props_backend is not None else (
        props_backend_from_env() if need_props else None
    )
    backend_key = str(getattr(backend, "value", backend)) if backend is not None else None
    prep_key = section_geometry_prep_key(
        layout,
        airfoil,
        skin_thickness=skin,
        spar_cap_thickness=spar,
        shear_web_thickness=web,
        twist_rad=twist_rad,
        props_backend=backend_key,
        outputs=out_set,
        web_x_positions_m=web_x_positions_m,
        zone_gauges_m=zone_gauges_m,
    )

    # Stage sink is always collected; attached to the result only when profiling.
    stage_s: dict[str, float] = {}

    with timed_stage("assemble", stage_s):
        section = build_section_with_thicknesses(
            airfoil,
            layout,
            twist_rad=float(twist_rad),
            skin_thickness=skin,
            spar_cap_thickness=float(spar),
            shear_web_thickness=web,
            web_x_positions_m=web_x_positions_m,
            web_y_coords=web_y_coords,
            stage_sink=stage_s,
            spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
        )

    properties_map: dict[str, dict[str, float]] | None = None
    properties_meta: SectionPropertiesMeta | None = None
    quadtree_mesh: QuadTreeMesh | None = None
    if need_props:
        from blade_structure_utilities.geometry_utilities.section_geometry.props_emit.props_compute import (
            compute_section_properties,
        )

        af_spec = airfoil_for_spec
        if af_spec is None:
            af_spec = airfoil
            if abs(float(twist_rad)) > 1e-12:
                af_spec = airfoil.rotate(float(twist_rad))
        with timed_stage("properties", stage_s):
            props_out = compute_section_properties(
                section,
                backend=backend,
                chord=float(airfoil.chord),
                twist_rad=float(twist_rad),
                airfoil_for_spec=af_spec,
                quadtree_spec=quadtree_spec,
                with_mesh=bool(retain_mesh),
            )
        if retain_mesh:
            properties_map, meta_raw, quadtree_mesh = props_out
        else:
            properties_map, meta_raw = props_out
            quadtree_mesh = None
        properties_meta = SectionPropertiesMeta.from_mapping(meta_raw)
        # Fold hybrid props / emit substage timers into stage_timings.
        if properties_meta.constructive_s is not None:
            stage_s["properties_constructive"] = float(
                stage_s.get("properties_constructive", 0.0)
            ) + float(properties_meta.constructive_s)
        if properties_meta.residual_quadtree_s is not None:
            stage_s["properties_residual"] = float(
                stage_s.get("properties_residual", 0.0)
            ) + float(properties_meta.residual_quadtree_s)
        for emit_key in (
            "emit_wall_solids",
            "emit_unions",
            "emit_cavity",
            "emit_cores",
            "properties_shoelace",
            "properties_emit",
            "properties_filter",
        ):
            raw = meta_raw.get(emit_key) if isinstance(meta_raw, Mapping) else None
            if raw is not None:
                stage_s[emit_key] = float(stage_s.get(emit_key, 0.0)) + float(raw)

    midline_strips: tuple[ShellMidlineStrip, ...] | None = None
    cell_wall_loops: tuple[CellWallLoop, ...] | None = None
    chord_m: float | None = None
    leading_edge_twisted_pitch_axis_blade_frame: tuple[float, float] | None = None
    trailing_edge_twisted_pitch_axis_blade_frame: tuple[float, float] | None = None
    if need_midlines:
        import numpy as np

        from blade_structure_utilities.geometry_utilities.section_geometry.midlines.shell_midline_export import (
            build_shell_midline_strips,
        )
        from blade_structure_utilities.geometry_utilities.section_geometry.topology.shell_strip_topology import (
            build_shell_strip_topology,
        )

        with timed_stage("midlines", stage_s):
            raw = build_shell_midline_strips(
                section,
                twist_rad=float(twist_rad),
                n_web_samples=int(n_web_samples),
                n_cap_samples=int(n_cap_samples),
            )
        with timed_stage("topology", stage_s):
            # Offset-locus interiors require the same frame as the strips.
            # Midlines are emitted in B; the airfoil locus is authored in S —
            # only share them when twist is ~0.
            locus_af = airfoil if abs(float(twist_rad)) < 1e-10 else None
            cap_depth = None
            if locus_af is not None:
                skin_max = max(
                    float(x)
                    for x in (skin if isinstance(skin, (list, tuple)) else (skin,))
                )
                cap_depth = float(skin_max) + 0.5 * float(spar)
            midline_strips = tuple(
                build_shell_strip_topology(
                    raw,
                    stage_sink=stage_s,
                    airfoil=locus_af,
                    cap_mid_depth_m=cap_depth,
                )
            )
        if need_void:
            from blade_structure_utilities.geometry_utilities.section_geometry.cell_wall_loops.cell_wall_loops import (
                build_cell_wall_loops,
            )

            with timed_stage("cell_wall_loops", stage_s):
                cell_wall_loops = tuple(build_cell_wall_loops(midline_strips))

        chord_m = float(airfoil.chord or 1.0)
        twist_out = float(twist_rad)
        # Full chord→B map (swap + pivot shift + twist), matching the midline
        # strips above. A bare rotation would leave LE/TE in the raw-NACA frame
        # while the strips were canonical B.
        from blade_structure_utilities.geometry_utilities.coordinate_frames import (
            naca_to_twisted_pitch_axis_blade_frame,
        )

        le_S = np.asarray(airfoil.leading_edge, dtype=float).reshape(1, 2)
        te_S = np.asarray(airfoil.trailing_edge, dtype=float).reshape(1, 2)
        leading_edge_twisted_pitch_axis_blade_frame = tuple(
            float(v)
            for v in naca_to_twisted_pitch_axis_blade_frame(
                le_S, float(twist_rad), chord_m
            )[0]
        )
        trailing_edge_twisted_pitch_axis_blade_frame = tuple(
            float(v)
            for v in naca_to_twisted_pitch_axis_blade_frame(
                te_S, float(twist_rad), chord_m
            )[0]
        )
    else:
        twist_out = None

    return SectionGeometryResult(
        section=section,
        prep_key=prep_key,
        properties_map=properties_map,
        properties_meta=properties_meta,
        midline_strips=midline_strips,
        cell_wall_loops=cell_wall_loops,
        quadtree_mesh=quadtree_mesh,
        chord_m=chord_m,
        twist_rad=twist_out,
        leading_edge_twisted_pitch_axis_blade_frame=leading_edge_twisted_pitch_axis_blade_frame,  # type: ignore[arg-type]
        trailing_edge_twisted_pitch_axis_blade_frame=trailing_edge_twisted_pitch_axis_blade_frame,  # type: ignore[arg-type]
        stage_timings=dict(stage_s) if profiling_enabled() else None,
    )


__all__ = [
    "ALL_OUTPUT_GATES",
    "DroppedWeb",
    "LayoutInput",
    "OutputGate",
    "SectionGeometryResult",
    "SectionPropertiesMeta",
    "StationLayoutResolution",
    "build_section_geometry",
    "build_section_with_thicknesses",
    "normalize_output_gates",
    "resolve_station_layout",
    "resolve_web_x_positions",
    "section_geometry_prep_key",
]
