"""Geometry worker entry: assemble a section and emit gated carriers."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from implicit_geometry_kernel.contract.cell_wall_loop import CellWallLoop
from implicit_geometry_kernel.contract.shell_midline_strip import ShellMidlineStrip
from implicit_geometry_kernel.integration.quadtree.tree import QuadTreeMesh
from implicit_geometry_kernel.section.spec.results import (
    SectionGeometryResult,
    SectionPropertiesMeta,
    normalize_output_gates,
)
from implicit_geometry_kernel.section.spec.prep import (
    build_section_with_thicknesses,
    section_geometry_prep_key,
)


def build_section_geometry(
    spec: Any,
    options: Any | None = None,
) -> SectionGeometryResult:
    """Geometry worker entry: assemble section and emit gated carriers.

    ``spec`` is a :class:`SectionSpec`; ``options`` is a :class:`SectionOptions`.
    """
    from implicit_geometry_kernel.section.spec.options import SectionOptions
    from implicit_geometry_kernel.section.spec.section_spec import SectionSpec

    if not isinstance(spec, SectionSpec):
        raise TypeError("build_section_geometry expects a SectionSpec")
    return _emit_section(spec, options or SectionOptions())


def _emit_section(spec: Any, options: Any) -> SectionGeometryResult:
    """Assemble and emit gated carriers from a validated spec."""
    layout = spec.layout
    airfoil = spec.airfoil
    skin = (spec.gauges.skin_ss, spec.gauges.skin_ps)
    spar = spec.gauges.cap
    web = spec.gauges.web_arg
    twist_rad = spec.twist_rad
    web_x_positions_m = spec.webs
    web_indices = spec.web_indices
    web_y_coords = spec.web_y_anchors
    spar_cap_ply_pitch_m = spec.spar_cap_ply_pitch_m
    zone_gauges_m = dict(spec.zone_gauges_m) if spec.zone_gauges_m else None
    outputs = options.outputs
    props_backend = options.props_backend
    n_web_samples = options.n_web_samples
    n_cap_samples = options.n_cap_samples
    retain_mesh = options.retain_mesh
    quadtree_spec = options.quadtree_spec
    airfoil_for_spec = options.airfoil_for_spec

    out_set = set(normalize_output_gates(outputs))

    need_props = "properties" in out_set
    need_midlines = "midlines" in out_set or "cell_wall_loops" in out_set
    need_void = "cell_wall_loops" in out_set

    from implicit_geometry_kernel.section.properties.backend import (
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
        skin_thickness=skin,  # type: ignore
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
            skin_thickness=skin,  # type: ignore
            spar_cap_thickness=spar,
            shear_web_thickness=web,
            web_x_positions_m=web_x_positions_m,
            web_indices=web_indices,
            web_y_coords=web_y_coords,
            stage_sink=stage_s,
            spar_cap_ply_pitch_m=spar_cap_ply_pitch_m,
        )

    properties_map: dict[str, dict[str, float]] | None = None
    properties_meta: SectionPropertiesMeta | None = None
    quadtree_mesh: QuadTreeMesh | None = None
    if need_props:
        from implicit_geometry_kernel.section.properties.compute import (
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
            properties_map, meta_raw, quadtree_mesh = props_out  # type: ignore
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
                pivot_fraction=spec.pivot_fraction,
            )
        with timed_stage("junction_preparation", stage_s):
            midline_strips = tuple(
                prepare_shell_midline_junctions(raw, stage_sink=stage_s)
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
            PIVOT_FRACTION_DEFAULT,
            naca_to_twisted_pitch_axis_blade_frame,
        )

        f_pivot = (
            PIVOT_FRACTION_DEFAULT
            if spec.pivot_fraction is None
            else float(spec.pivot_fraction)
        )
        le_S = np.asarray(airfoil.leading_edge, dtype=float).reshape(1, 2)
        te_S = np.asarray(airfoil.trailing_edge, dtype=float).reshape(1, 2)
        leading_edge_twisted_pitch_axis_blade_frame = tuple(  # type: ignore
            float(v)
            for v in naca_to_twisted_pitch_axis_blade_frame(
                le_S, float(twist_rad), chord_m, f_pivot
            )[0]
        )
        trailing_edge_twisted_pitch_axis_blade_frame = tuple(  # type: ignore
            float(v)
            for v in naca_to_twisted_pitch_axis_blade_frame(
                te_S, float(twist_rad), chord_m, f_pivot
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


__all__ = ["_emit_section", "build_section_geometry"]
