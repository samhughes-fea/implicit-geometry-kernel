"""Shell midline junction preparation for strips (B-frame).

Thin orchestrator over focused helpers.
"""

from __future__ import annotations

from dataclasses import dataclass

from implicit_geometry_kernel.contract.shell_midline_strip import ShellMidlineStrip
from implicit_geometry_kernel.profiling import timed_stage

from implicit_geometry_kernel.section.midlines.cap_skin_ray_resample import (
    concat_skin_ring_vertices,
    pick_inward_ray_direction,
    resample_cap_interiors_from_skin_rays,
    skin_ring_vertex_tangents,
)
from implicit_geometry_kernel.section.junctions.cap_skin_bonds import bond_caps_to_skin
from implicit_geometry_kernel.section.junctions.cap_web_junctions import insert_cap_web_interior_junctions
from implicit_geometry_kernel.section.midlines.ring_geometry import normal_2d_ccw, point_to_open_polyline_nearest
from implicit_geometry_kernel.section.midlines.skin_split import split_skin_at_junctions



@dataclass(frozen=True)
class ShellMidlineJunctionPreparationOptions:
    cap_resample_from_skin_rays: bool = True
    #: ``None`` derives the tolerance from the arc being split, which scales
    #: with the section; a number pins it absolutely.
    merge_s_tol_m: float | None = None


def prepare_shell_midline_junctions(
    strips: tuple[ShellMidlineStrip, ...],
    spec: ShellMidlineJunctionPreparationOptions | None = None,
    *,
    stage_sink: dict[str, float] | None = None,
    airfoil=None,
    cap_mid_depth_m: float | None = None,
) -> tuple[ShellMidlineStrip, ...]:
    """Prepare strip junctions for shell meshing.

    When ``stage_sink`` is provided (or ``SECTION_GEOMETRY_PROFILE=1``), records
    ``junction_preparation_split`` / ``junction_preparation_cap_ray`` /
    ``junction_preparation_cap_web`` / ``junction_preparation_cap_bond`` substages.

    Optional ``airfoil`` + ``cap_mid_depth_m`` (same frame as strips) let interior
    junctions prefer axis × offset locus before falling back to the final cap
    polyline.
    """
    cfg = spec if spec is not None else ShellMidlineJunctionPreparationOptions()
    has_webs = any(s.kind == "web" for s in strips)

    # N = 0 (skin ± caps): do not invent web feet / LE–bay–TE panelisation.
    # Cap ray resample still applies for family-B caps; cap–web inserts no-op.
    if has_webs:
        with timed_stage("junction_preparation_split", stage_sink):
            out = split_skin_at_junctions(strips, merge_s_tol_m=cfg.merge_s_tol_m)
    else:
        with timed_stage("junction_preparation_split", stage_sink):
            out = strips

    if cfg.cap_resample_from_skin_rays:
        with timed_stage("junction_preparation_cap_ray", stage_sink):
            out = resample_cap_interiors_from_skin_rays(out)

    # Cap–web interior junctions are a no-op without webs; skip the stage.
    if has_webs:
        with timed_stage("junction_preparation_cap_web", stage_sink):
            out = insert_cap_web_interior_junctions(
                out,
                airfoil=airfoil,
                cap_mid_depth_m=cap_mid_depth_m,
            )

    # Last: every later pass must leave cap and skin vertices alone, because the
    # bonds name them (thin_wall_section TW-001).
    if any(s.kind == "cap" for s in out):
        with timed_stage("junction_preparation_cap_bond", stage_sink):
            out = bond_caps_to_skin(tuple(out))
    return out


__all__ = [
    "ShellMidlineJunctionPreparationOptions",
    "prepare_shell_midline_junctions",
    "split_skin_at_junctions",
    "resample_cap_interiors_from_skin_rays",
    "insert_cap_web_interior_junctions",
    "bond_caps_to_skin",
    "concat_skin_ring_vertices",
    "skin_ring_vertex_tangents",
    "normal_2d_ccw",
    "pick_inward_ray_direction",
    "point_to_open_polyline_nearest",
]
