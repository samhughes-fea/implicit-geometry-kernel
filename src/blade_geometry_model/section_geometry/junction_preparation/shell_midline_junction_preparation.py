"""Shell midline junction preparation for strips (B-frame).

Thin orchestrator over focused helpers. Canonical public names live here;
historic ``ShellStripTopologySpec`` / ``build_shell_strip_topology`` names are
warning aliases until 2027-03-31 (also via the ``topology`` package shim).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

from blade_geometry_model.section_contract.shell_midline_strip import ShellMidlineStrip
from blade_geometry_model.section_geometry.profiling import timed_stage

from blade_geometry_model.section_geometry.midlines.cap_skin_ray_resample import (
    concat_skin_ring_vertices,
    pick_inward_ray_direction,
    resample_cap_interiors_from_skin_rays,
    skin_ring_vertex_tangents,
)
from .cap_skin_bonds import bond_caps_to_skin
from .cap_web_junctions import insert_cap_web_interior_junctions
from ..midlines.ring_geometry import normal_2d_ccw, point_to_open_polyline_nearest
from ..midlines.skin_split import split_skin_at_junctions

# Compatibility aliases (historic private names used by tests / callers).
_split_skin_at_junctions = split_skin_at_junctions
_resample_cap_interiors_from_skin_rays = resample_cap_interiors_from_skin_rays
_insert_cap_web_interior_junctions = insert_cap_web_interior_junctions
_concat_skin_ring_vertices = concat_skin_ring_vertices
_skin_ring_vertex_tangents = skin_ring_vertex_tangents
_normal_2d_ccw = normal_2d_ccw
_pick_inward_ray_direction = pick_inward_ray_direction
_point_to_open_polyline_nearest = point_to_open_polyline_nearest

_ALIAS_EXPIRY = "2027-03-31"


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


def build_shell_strip_topology(
    strips: tuple[ShellMidlineStrip, ...],
    spec: ShellMidlineJunctionPreparationOptions | None = None,
    *,
    stage_sink: dict[str, float] | None = None,
    airfoil=None,
    cap_mid_depth_m: float | None = None,
) -> tuple[ShellMidlineStrip, ...]:
    """Deprecated alias of :func:`prepare_shell_midline_junctions` (expires 2027-03-31)."""
    warnings.warn(
        "build_shell_strip_topology is deprecated; use "
        f"prepare_shell_midline_junctions (expires {_ALIAS_EXPIRY}).",
        DeprecationWarning,
        stacklevel=2,
    )
    return prepare_shell_midline_junctions(
        strips,
        spec,
        stage_sink=stage_sink,
        airfoil=airfoil,
        cap_mid_depth_m=cap_mid_depth_m,
    )


# Historic type name — prefer ShellMidlineJunctionPreparationOptions.
ShellStripTopologySpec = ShellMidlineJunctionPreparationOptions


__all__ = [
    "ShellMidlineJunctionPreparationOptions",
    "prepare_shell_midline_junctions",
    "ShellStripTopologySpec",
    "build_shell_strip_topology",
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
