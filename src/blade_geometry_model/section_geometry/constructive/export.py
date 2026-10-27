"""Section-level constructive wall-solid and core polygon export."""

from __future__ import annotations

from typing import Any

from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

from blade_geometry_model.section_geometry.midlines.section_midline_host import (
    resolve_section_host,
)
from blade_geometry_model.section_geometry.profiling import timed_stage

from .emitters import (
    emit_inner_cavity_polygon,
    emit_outer_skin_polygon,
)
from .offset_locus import airfoil_polygon, ensure_polygon, offset_polygon_inward, rotate_geom
from .labels import (
    is_core_label,
    is_insert_label,
    is_skin_label,
    is_wall_solid_label,
)


def _label_key(label: Any) -> str:
    return str(getattr(label, "value", label))


def _resolve_host(section: Any) -> Any:
    return resolve_section_host(section)


def _component_map(section: Any) -> dict[str, Any]:
    host = _resolve_host(section)
    unrot = getattr(host, "components_unrotated", None)
    if isinstance(unrot, dict) and unrot:
        return dict(unrot)
    comps = getattr(host, "components", None)
    if isinstance(comps, dict):
        return dict(comps)
    out = {}
    for lbl in getattr(host, "labels", ()) or list(host):
        key = _label_key(lbl)
        out[key] = host[key]
    return out


def _union_discrete_caps(
    section: Any,
    discrete: dict[str, Any],
    *,
    surface: str,
    airfoil_poly: Polygon | MultiPolygon | None = None,
    inner_cavity: Polygon | MultiPolygon | None = None,
    deeper_cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon | None:
    long_prefix = f"{_spar_cap_label_for_surface(surface)}_"
    legacy_prefix = f"spar_cap_{surface}_"
    parts = []
    # One deeper offset for all discrete bands on this surface (same skin+cap_h).
    shared_deeper = deeper_cavity
    for key, comp in discrete.items():
        lk = str(_label_key(key))
        if not (lk.startswith(long_prefix) or lk.startswith(legacy_prefix)):
            continue
        if shared_deeper is None and airfoil_poly is not None:
            depth = float(comp.skin_thickness) + float(comp.cap_height)
            shared_deeper = offset_polygon_inward(airfoil_poly, depth)
        parts.append(
            _cap_poly_via_recipe(
                section,
                lk,
                comp,
                airfoil_poly=airfoil_poly,
                inner_cavity=inner_cavity,
                deeper_cavity=shared_deeper,
            )
        )
    if not parts:
        return None
    return ensure_polygon(unary_union(parts), allow_empty=True)


def _cap_from_recipes(
    section: Any,
    comps: dict[str, Any],
    discrete: dict[str, Any],
    *,
    surface: str,
    airfoil_poly: Polygon | MultiPolygon | None = None,
    inner_cavity: Polygon | MultiPolygon | None = None,
    deeper_cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon | None:
    """Emit a surface spar-cap polygon from assemble-time CapWall recipes."""
    key = _spar_cap_label_for_surface(surface)
    legacy = f"spar_cap_{surface}"
    merged = _union_discrete_caps(
        section,
        discrete,
        surface=surface,
        airfoil_poly=airfoil_poly,
        inner_cavity=inner_cavity,
        deeper_cavity=deeper_cavity,
    )
    if merged is not None:
        return merged
    # Prefer recipe keys on the host; component maps only supply airfoil refs.
    for lbl, wall in _host_cap_recipes(section).items():
        lk = _label_key(lbl)
        if lk in (key, legacy) or str(lk).startswith(f"{key}_") or str(lk).startswith(f"{legacy}_"):
            af = getattr(_resolve_host(section), "airfoil", None)
            if af is None:
                for src in (discrete, comps):
                    for _, comp in src.items():
                        af = getattr(comp, "_af", None)
                        if af is not None:
                            break
                    if af is not None:
                        break
            return wall.to_boundary_polygon(
                af,
                airfoil_poly=airfoil_poly,
                inner_cavity=inner_cavity,
                deeper_cavity=deeper_cavity,
            )
    return None


    return None


def _shared_deeper_by_surface(
    section: Any,
    comps: dict[str, Any],
    discrete: dict[str, Any],
    *,
    airfoil_poly: Polygon | MultiPolygon | None,
) -> dict[str, Polygon | MultiPolygon | None]:
    """Build deeper-offset polygons once per distinct (skin, cap_h)."""
    if airfoil_poly is None:
        return {"upper": None, "lower": None}

    def _height_for(surface: str) -> tuple[float, float] | None:
        key = _spar_cap_label_for_surface(surface)
        legacy = f"spar_cap_{surface}"
        prefix = f"{key}_"
        legacy_prefix = f"{legacy}_"
        for lbl, wall in _host_cap_recipes(section).items():
            lk = _label_key(lbl)
            if lk in (key, legacy) or str(lk).startswith(prefix) or str(lk).startswith(legacy_prefix):
                return float(wall.skin_thickness), float(wall.cap_height)
        for lbl, wall in _host_recipes(section).items():
            lk = _label_key(lbl)
            if lk in (key, legacy) or str(lk).startswith(prefix) or str(lk).startswith(legacy_prefix):
                return float(wall.skin_thickness), float(wall.cap_height)
        for lbl, comp in discrete.items():
            lk = _label_key(lbl)
            if lk in (key, legacy) or str(lk).startswith(prefix) or str(lk).startswith(legacy_prefix):
                return float(comp.skin_thickness), float(comp.cap_height)
        for lbl, comp in comps.items():
            if _label_key(lbl) in (key, legacy):
                return float(comp.skin_thickness), float(comp.cap_height)
        return None

    cache: dict[tuple[float, float], Polygon | MultiPolygon] = {}
    out: dict[str, Polygon | MultiPolygon | None] = {"upper": None, "lower": None}
    for surface in ("upper", "lower"):
        pair = _height_for(surface)
        if pair is None:
            continue
        if pair not in cache:
            skin_t, cap_h = pair
            cache[pair] = offset_polygon_inward(airfoil_poly, skin_t + cap_h)
        out[surface] = cache[pair]
    return out


def _web_union(web_polys: list[Polygon | MultiPolygon]) -> Polygon | MultiPolygon | None:
    if not web_polys:
        return None
    if len(web_polys) == 1:
        return web_polys[0]
    return unary_union(web_polys)


def _spar_cap_label_for_surface(surface: str) -> str:
    """Map camber ``upper``/``lower`` to canonical spar-cap label."""
    if surface in ("upper", "suction_surface"):
        return "spar_cap_suction_surface"
    if surface in ("lower", "pressure_surface"):
        return "spar_cap_pressure_surface"
    return f"spar_cap_{surface}"


def _host_recipes(section: Any) -> dict[str, Any]:
    """Prefer assemble-time recipes on the host; empty if absent."""
    host = _resolve_host(section)
    recipes = getattr(host, "wall_recipes", None)
    if isinstance(recipes, dict) and recipes:
        return dict(recipes)
    return {}


def _host_cap_recipes(section: Any) -> dict[str, Any]:
    host = _resolve_host(section)
    recipes = getattr(host, "spar_cap_wall_recipes", None)
    if isinstance(recipes, dict) and recipes:
        return dict(recipes)
    return {}


def _recipe_for(section: Any, key: str, *, kind: str):
    """Look up a stored assemble-time recipe (required on production hosts)."""
    recipes = _host_recipes(section)
    if key in recipes:
        return recipes[key]
    if kind == "cap":
        cap_recipes = _host_cap_recipes(section)
        if key in cap_recipes:
            return cap_recipes[key]
        for rk, wall in cap_recipes.items():
            if rk.startswith(f"{key}_") or key.startswith(f"{rk}_"):
                return wall
    raise KeyError(
        f"No {kind} wall recipe for {key!r} on section host; "
        "MultiCellSection must author recipes at assemble."
    )


def _skin_via_recipe(section: Any, key: str, comp: Any, af_ref: Any):
    """Route outer-skin emission through :class:`SkinWall`."""
    if not hasattr(comp, "thickness"):
        return comp
    wall = _recipe_for(section, key, kind="skin")
    return wall.to_sdf_component(
        getattr(comp, "_af", None) or af_ref,
        camber_xy=getattr(comp, "_camber_xy", None),
    )


def _web_poly_via_recipe(section: Any, key: str, comp: Any, *, cavity=None):
    wall = _recipe_for(section, key, kind="web")
    af = getattr(comp, "_af", None)
    return wall.to_boundary_polygon(af, cavity=cavity)


def _cap_poly_via_recipe(section: Any, key: str, comp: Any, **emit_kwargs):
    wall = _recipe_for(section, key, kind="cap")
    af = getattr(comp, "_af", None)
    return wall.to_boundary_polygon(af, **emit_kwargs)


def _core_poly_via_recipe(
    section: Any, key: str, comp: Any, *, cavity=None, wall_exclusion_union=None
):
    wall = _recipe_for(section, key, kind="core")
    return wall.to_boundary_polygon(
        getattr(comp, "_af", None),
        cavity=cavity,
        wall_exclusion_union=wall_exclusion_union,
    )


def _insert_poly_via_recipe(
    section: Any, key: str, comp: Any, *, airfoil_poly=None, inner_cavity=None
):
    wall = _recipe_for(section, key, kind="insert")
    af = getattr(comp, "_af", None)
    return wall.to_boundary_polygon(
        af, airfoil_poly=airfoil_poly, inner_cavity=inner_cavity
    )


def _is_skin_key(key: str) -> bool:
    return (
        key == "outer_skin"
        or key.startswith("skin_suction_surface")
        or key.startswith("skin_pressure_surface")
        or key in {"skin_suction_surface", "skin_pressure_surface"}
    )


def _is_web_key(key: str) -> bool:
    return key.startswith("shear_web_")


def emit_wall_solid_polygons(section: Any) -> dict[str, Polygon | MultiPolygon]:
    """Emit constructive polygons for wall-solid labels (skin / webs / caps).

    Polygons are returned in the **section frame** (B-frame when twisted).
    """
    host = _resolve_host(section)
    comps = _component_map(section)
    twist = float(getattr(host, "twist_rad", 0.0) or 0.0)
    discrete = getattr(host, "spar_cap_components_unrotated", None) or {}

    out: dict[str, Polygon | MultiPolygon] = {}
    web_polys: list[Polygon | MultiPolygon] = []

    af_ref = None
    skin_t_ref: float | None = None
    for lbl, comp in comps.items():
        af_ref = getattr(comp, "_af", None)
        if af_ref is None:
            continue
        key = _label_key(lbl)
        if _is_skin_key(key):
            skin_t_ref = float(comp.thickness)
            break
        if hasattr(comp, "skin_thickness"):
            skin_t_ref = float(comp.skin_thickness)
            break
    af_poly = None
    cavity = None
    if af_ref is not None and skin_t_ref is not None:
        try:
            af_poly = airfoil_polygon(af_ref)
            cavity = emit_inner_cavity_polygon(
                af_ref, float(skin_t_ref), airfoil_poly=af_poly
            )
        except ValueError:
            cavity = Polygon()

    for label, comp in comps.items():
        key = _label_key(label)
        if not is_wall_solid_label(key):
            continue
        if _is_skin_key(key):
            out[key] = emit_outer_skin_polygon(
                _skin_via_recipe(section, key, comp, af_ref),
                airfoil_poly=af_poly,
                inner_cavity=cavity,
            )
        elif _is_web_key(key):
            wp = _web_poly_via_recipe(section, key, comp, cavity=cavity)
            out[key] = wp
            web_polys.append(wp)
        elif is_insert_label(key):
            out[key] = _insert_poly_via_recipe(
                section, key, comp, airfoil_poly=af_poly, inner_cavity=cavity
            )

    webs_u = _web_union(web_polys)
    deeper_by_surface = _shared_deeper_by_surface(
        section, comps, discrete, airfoil_poly=af_poly
    )
    for surface in ("upper", "lower"):
        key = _spar_cap_label_for_surface(surface)
        keys = {_label_key(k) for k in comps}
        if key not in keys and f"spar_cap_{surface}" not in keys:
            # Also accept recipes-only hosts.
            recipe_keys = set(_host_recipes(section)) | set(_host_cap_recipes(section))
            if key not in recipe_keys and f"spar_cap_{surface}" not in recipe_keys:
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
        if cap is None:
            continue
        if webs_u is not None:
            cap = cap.difference(webs_u)
        if cap.is_empty or float(cap.area) <= 0.0:
            continue
        out[key] = ensure_polygon(cap, allow_empty=False)

    if abs(twist) > 1e-10:
        out = {k: rotate_geom(v, twist, cx=0.0, cy=0.0) for k, v in out.items()}
    return out


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

