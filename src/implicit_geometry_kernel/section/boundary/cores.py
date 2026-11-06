"""Exact core and insert polygons, and the combined constructive polygon set."""

from __future__ import annotations

from typing import Any

from implicit_geometry_kernel.section.boundary import MultiPolygon, Polygon
from implicit_geometry_kernel.section.boundary import unary_union

from implicit_geometry_kernel.profiling import timed_stage

from implicit_geometry_kernel.section.boundary.emitters import (
    emit_inner_cavity_polygon,
    emit_outer_skin_polygon,
)
from implicit_geometry_kernel.section.boundary import airfoil_polygon, ensure_polygon, rotate_geom
from implicit_geometry_kernel.section.boundary.labels import (
    is_core_label,
    is_insert_label,
    is_skin_label,
    is_wall_solid_label,
)
from implicit_geometry_kernel.section.boundary.walls import (
    _cap_from_recipes,
    _component_map,
    _core_poly_via_recipe,
    _host_cap_recipes,
    _host_recipes,
    _insert_poly_via_recipe,
    _is_web_key,
    _label_key,
    _resolve_host,
    _shared_deeper_by_surface,
    _skin_via_recipe,
    _spar_cap_label_for_surface,
    _web_poly_via_recipe,
    _web_union,
    emit_wall_solid_polygons,
)


def emit_core_polygons(
    section: Any,
    wall_solid_polys: dict[str, Polygon | MultiPolygon] | None = None,
) -> dict[str, Polygon | MultiPolygon]:
    """Emit constructive ``core_*`` bay polygons (chord frame then section twist)."""
    host = _resolve_host(section)
    comps = _component_map(section)
    twist = float(getattr(host, "twist_rad", 0.0) or 0.0)

    if wall_solid_polys is None:
        # Wall solids needed in chord frame for exclusions â€” emit then un-rotate if twisted.
        wall_solid_polys = emit_wall_solid_polygons(section)
        if abs(twist) > 1e-10:
            wall_solid_polys = {
                k: rotate_geom(v, -twist, cx=0.0, cy=0.0) for k, v in wall_solid_polys.items()
            }
    elif abs(twist) > 1e-10:
        # Caller passed section-frame wall solids; undo twist for chord-frame Booleans.
        wall_solid_polys = {
            k: rotate_geom(v, -twist, cx=0.0, cy=0.0) for k, v in wall_solid_polys.items()
        }

    exclusion_parts = [
        v
        for k, v in wall_solid_polys.items()
        if not is_skin_label(k) and not v.is_empty
    ]
    exclusion = unary_union(exclusion_parts) if exclusion_parts else None

    # Shared cavity + one cavity−wall-solids difference for all bays.
    cavity = None
    open_region: Polygon | MultiPolygon | None = None
    out: dict[str, Polygon | MultiPolygon] = {}
    for label, comp in comps.items():
        key = _label_key(label)
        if not is_core_label(key):
            continue
        if cavity is None:
            cavity = emit_inner_cavity_polygon(comp._af, float(comp.skin_thickness))
            open_region = cavity
            if exclusion is not None and not exclusion.is_empty:
                open_region = cavity.difference(exclusion)
        poly = _core_poly_via_recipe(
            section, key, comp, cavity=open_region, wall_exclusion_union=None
        )
        if poly.is_empty or float(poly.area) <= 0.0:
            out[key] = Polygon()
        else:
            out[key] = poly

    if abs(twist) > 1e-10:
        out = {
            k: (rotate_geom(v, twist, cx=0.0, cy=0.0) if not v.is_empty else v)
            for k, v in out.items()
        }
    return out


def _keep_poly(poly: Polygon | MultiPolygon | None) -> bool:
    if poly is None or poly.is_empty:
        return False
    try:
        return float(poly.area) > 1e-16
    except Exception:
        return False


def emit_constructive_polygons(
    section: Any,
    *,
    stage_sink: dict[str, float] | None = None,
) -> dict[str, Polygon | MultiPolygon]:
    """Emit all constructive property carriers: wall solids + sandwich cores.

    Empty / near-zero Boolean leftovers are omitted so shoelace never sees
    degenerate rings (Shapely area can be tiny while shoelace A â†’ 0).

    When ``stage_sink`` is provided (or ``SECTION_GEOMETRY_PROFILE=1``), records
    ``emit_wall_solids`` / ``emit_unions`` / ``emit_cavity`` / ``emit_cores``.
    """
    host = _resolve_host(section)
    twist = float(getattr(host, "twist_rad", 0.0) or 0.0)

    # Emit wall solids in chord frame first so cores can subtract without twist games.
    comps = _component_map(section)
    discrete = getattr(host, "spar_cap_components_unrotated", None) or {}
    wall_chord: dict[str, Polygon | MultiPolygon] = {}
    web_polys: list[Polygon | MultiPolygon] = []
    core_items = [
        (_label_key(label), comp)
        for label, comp in comps.items()
        if is_core_label(_label_key(label))
    ]

    # Shared airfoil / inner-cavity polygon for skin, webs, caps, and cores.
    af_ref = None
    skin_t_ref: float | None = None
    for _lbl, comp in comps.items():
        af_ref = getattr(comp, "_af", None)
        if af_ref is None:
            continue
        if hasattr(comp, "skin_thickness"):
            skin_t_ref = float(comp.skin_thickness)
        elif hasattr(comp, "thickness") and is_wall_solid_label(_label_key(_lbl)):
            # OuterSkin stores thickness as the skin wall thickness.
            if is_skin_label(_label_key(_lbl)):
                skin_t_ref = float(comp.thickness)
        if skin_t_ref is not None:
            break
    if af_ref is None and core_items:
        af_ref = getattr(core_items[0][1], "_af", None)
        skin_t_ref = float(core_items[0][1].skin_thickness)

    af_poly: Polygon | MultiPolygon | None = None
    cavity: Polygon | MultiPolygon = Polygon()
    with timed_stage("emit_cavity", stage_sink):
        if af_ref is not None and skin_t_ref is not None:
            try:
                af_poly = airfoil_polygon(af_ref)
                cavity = emit_inner_cavity_polygon(
                    af_ref, float(skin_t_ref), airfoil_poly=af_poly
                )
            except ValueError:
                cavity = Polygon()

    raw_caps: dict[str, Polygon | MultiPolygon] = {}
    with timed_stage("emit_wall_solids", stage_sink):
        for label, comp in comps.items():
            key = _label_key(label)
            if not is_wall_solid_label(key):
                continue
            if is_skin_label(key):
                poly = emit_outer_skin_polygon(
                    _skin_via_recipe(section, key, comp, af_ref),
                    airfoil_poly=af_poly,
                    inner_cavity=cavity,
                )
                if _keep_poly(poly):
                    wall_chord[key] = poly
            elif _is_web_key(key):
                wp = _web_poly_via_recipe(section, key, comp, cavity=cavity)
                if _keep_poly(wp):
                    wall_chord[key] = wp
                    web_polys.append(wp)
            elif is_insert_label(key):
                ip = _insert_poly_via_recipe(
                    section, key, comp, airfoil_poly=af_poly, inner_cavity=cavity
                )
                if _keep_poly(ip):
                    wall_chord[key] = ip

        deeper_by_surface = _shared_deeper_by_surface(
            section, comps, discrete, airfoil_poly=af_poly
        )
        for surface in ("upper", "lower"):
            key = _spar_cap_label_for_surface(surface)
            keys = {_label_key(k) for k in comps}
            recipe_keys = set(_host_recipes(section)) | set(_host_cap_recipes(section))
            if (
                key not in keys
                and f"spar_cap_{surface}" not in keys
                and key not in recipe_keys
                and f"spar_cap_{surface}" not in recipe_keys
            ):
                continue
            cap = _cap_from_recipes(
                section,
                comps,
                discrete,
                surface=surface,
                airfoil_poly=af_poly,
                inner_cavity=cavity,
                deeper_cavity=deeper_by_surface.get(surface),
            )
            if cap is not None and _keep_poly(cap):
                raw_caps[key] = ensure_polygon(cap, allow_empty=True)

    exclusion = None
    open_region: Polygon | MultiPolygon | None = cavity
    with timed_stage("emit_unions", stage_sink):
        # One web union for both upper/lower cap differences (was 2Ã— unary_union).
        webs_u = _web_union(web_polys)
        for key, cap in raw_caps.items():
            if webs_u is not None:
                cap = cap.difference(webs_u)
            if _keep_poly(cap):
                wall_chord[key] = ensure_polygon(cap, allow_empty=True)
        exclusion_parts = [
            v for k, v in wall_chord.items() if not is_skin_label(k) and _keep_poly(v)
        ]
        exclusion = unary_union(exclusion_parts) if exclusion_parts else None
        # One cavity − wall solids difference shared by every bay.
        open_region = cavity
        if exclusion is not None and not exclusion.is_empty and not cavity.is_empty:
            open_region = cavity.difference(exclusion)

    cores: dict[str, Polygon | MultiPolygon] = {}
    if core_items:
        with timed_stage("emit_cores", stage_sink):
            for key, comp in core_items:
                poly = _core_poly_via_recipe(
                    section,
                    key,
                    comp,
                    cavity=open_region,
                    wall_exclusion_union=None,
                )
                if _keep_poly(poly):
                    cores[key] = poly

    out = dict(wall_chord)
    out.update(cores)
    if abs(twist) > 1e-10:
        twisted: dict[str, Polygon | MultiPolygon] = {}
        for k, v in out.items():
            tv = rotate_geom(v, twist, cx=0.0, cy=0.0)
            if _keep_poly(tv):
                twisted[k] = ensure_polygon(tv, allow_empty=True)
        out = twisted
    return out

