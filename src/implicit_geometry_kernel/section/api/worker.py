"""Use-case-agnostic section geometry worker facade.

Mirrors the FE worker contract: no sizing/fatigue modes; callers choose
output gates (compute-cost knobs) and interpret results in orchestration.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping
from typing import Any, Iterable, Sequence

from implicit_geometry_kernel.contract.cell_wall_loop import (
    CellWallLoop,
)
from implicit_geometry_kernel.contract.shell_midline_strip import (
    ShellMidlineStrip,
)
from implicit_geometry_kernel.section.api.layout_resolution import (
    StationLayoutResolution,
    _as_section_layout,
    _layout_web_fracs,
    resolve_station_layout,
    resolve_web_x_positions,
)
from implicit_geometry_kernel.section.api.results import (
    ALL_OUTPUT_GATES,
    DroppedWeb,
    LayoutInput,
    OutputGate,
    SectionGeometryResult,
    SectionPropertiesMeta,
    normalize_output_gates,
)
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
from implicit_geometry_kernel.section.midlines.host import (
    SectionMidlineHost,
)
from implicit_geometry_kernel.section.props.compute import (
    PropsBackendArg,
)
from implicit_geometry_kernel.sdf.airfoil import (
    AirfoilSDF,
)
from implicit_geometry_kernel.integration.quadtree import QuadTreeMesh
from implicit_geometry_kernel.integration.spec import QuadtreeSpec


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
    from implicit_geometry_kernel.section.api.kwarg_aliases import (
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
    from implicit_geometry_kernel.section.api.kwarg_aliases import (
        resolve_shear_web_thickness,
        resolve_skin_gauges,
        resolve_spar_cap_thickness,
        resolve_twist_rad,
    )
    from implicit_geometry_kernel.section.assembly.from_layout import (
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
    web_offsets: Sequence[float] | None = None,
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
    from implicit_geometry_kernel.section.api.kwarg_aliases import (
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
    if web_offsets is not None and web_x_positions_m is not None:
        raise ValueError("pass web_offsets or web_x_positions_m, not both")
    if web_offsets is not None:
        web_x_positions_m = web_offsets

    out_set = set(normalize_output_gates(outputs))

    need_props = "properties" in out_set
    need_midlines = "midlines" in out_set or "cell_wall_loops" in out_set
    need_void = "cell_wall_loops" in out_set

    from implicit_geometry_kernel.section.props.backend import (
        props_backend_from_env,
    )
    from implicit_geometry_kernel.profiling import (
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
        from implicit_geometry_kernel.section.props.compute import (
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

        from implicit_geometry_kernel.section.midlines.export import (
            build_shell_midline_strips,
        )
        from implicit_geometry_kernel.section.junctions.prepare import (
            prepare_shell_midline_junctions,
        )

        with timed_stage("midlines", stage_s):
            raw = build_shell_midline_strips(
                section,
                twist_rad=float(twist_rad),
                n_web_samples=int(n_web_samples),
                n_cap_samples=int(n_cap_samples),
            )
        with timed_stage("junction_preparation", stage_s):
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
                prepare_shell_midline_junctions(
                    raw,
                    stage_sink=stage_s,
                    airfoil=locus_af,
                    cap_mid_depth_m=cap_depth,
                )
            )
        if need_void:
            from implicit_geometry_kernel.section.cells.loops import (
                build_cell_wall_loops,
            )

            with timed_stage("cell_wall_loops", stage_s):
                cell_wall_loops = tuple(build_cell_wall_loops(midline_strips))

        chord_m = float(airfoil.chord or 1.0)
        twist_out = float(twist_rad)
        # Full chord→B map (swap + pivot shift + twist), matching the midline
        # strips above. A bare rotation would leave LE/TE in the raw-NACA frame
        # while the strips were canonical B.
        from implicit_geometry_kernel.frames import (
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
