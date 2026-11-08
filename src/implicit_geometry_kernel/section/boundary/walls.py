"""Exact wall-solid polygons (skin, webs, caps) for a section."""

from __future__ import annotations

from typing import Any

from implicit_geometry_kernel.section.boundary import MultiPolygon, Polygon
from implicit_geometry_kernel.section.boundary import OPS, unary_union

from implicit_geometry_kernel.section.midlines.host import (
    resolve_section_host,
)

from implicit_geometry_kernel.section.boundary.emitters import (
    emit_inner_cavity_polygon,
    emit_outer_skin_polygon,
)
from implicit_geometry_kernel.section.boundary import airfoil_polygon, ensure_polygon, offset_polygon_inward, rotate_geom
from implicit_geometry_kernel.section.assembly.members import (
    CapBuilder,
    wall_type_for_label,
)
from implicit_geometry_kernel.section.boundary.labels import (
    is_insert_label,
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
    parts = []
    # One deeper offset for all discrete bands on this surface (same skin+cap_h).
    # Indexed pads only (``*_0``, ``*_1``, …). The canonical surface label is the
    # merged C-family alias or the D-family continuous box, not a pad.
    shared_deeper = deeper_cavity
    canon = CapBuilder.canonical_surface_label(surface)
    legacy = f"spar_cap_{surface}"
    for key, comp in discrete.items():
        lk = str(_label_key(key))
        if lk in (canon, legacy):
            continue
        if not CapBuilder.matches_surface(lk, surface):
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
    key = CapBuilder.canonical_surface_label(surface)
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
        if CapBuilder.matches_surface(lk, surface) or lk == key:
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
        for lbl, wall in _host_cap_recipes(section).items():
            if CapBuilder.matches_surface(_label_key(lbl), surface):
                return float(wall.skin_thickness), float(wall.cap_height)
        for lbl, wall in _host_recipes(section).items():
            if CapBuilder.matches_surface(_label_key(lbl), surface):
                return float(wall.skin_thickness), float(wall.cap_height)
        for lbl, comp in discrete.items():
            if CapBuilder.matches_surface(_label_key(lbl), surface):
                return float(comp.skin_thickness), float(comp.cap_height)
        for lbl, comp in comps.items():
            if CapBuilder.matches_surface(_label_key(lbl), surface):
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
    return CapBuilder.canonical_surface_label(surface)


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
    wt = wall_type_for_label(key)
    if wt is not None and wt.polygon_by_surface:
        cap_recipes = _host_cap_recipes(section)
        if key in cap_recipes:
            return cap_recipes[key]
        for rk, wall in cap_recipes.items():
            if CapBuilder.matches_label(rk) and (
                rk == key or key in rk or rk in key
            ):
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
    wt = wall_type_for_label(key)
    return wt is not None and wt.kind == "skin"


def _is_web_key(key: str) -> bool:
    wt = wall_type_for_label(key)
    return wt is not None and wt.cell_face_kind == "web"


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
        else:
            wt = wall_type_for_label(key)
            if wt is None or wt.polygon_by_surface:
                continue
            recipes = dict(getattr(section, "wall_recipes", None) or {})
            wall = recipes.get(key)
            if wall is not None and hasattr(wall, "to_boundary_polygon"):
                out[key] = wall.to_boundary_polygon(af_ref, cavity=cavity)

    webs_u = _web_union(web_polys)
    deeper_by_surface = _shared_deeper_by_surface(
        section, comps, discrete, airfoil_poly=af_poly
    )
    for surface in ("upper", "lower"):
        key = CapBuilder.canonical_surface_label(surface)
        keys = {_label_key(k) for k in comps}
        recipe_keys = set(_host_recipes(section)) | set(_host_cap_recipes(section))
        present = any(CapBuilder.matches_surface(k, surface) for k in keys | recipe_keys)
        if not present:
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
            cap = OPS.difference(cap, webs_u)
        if cap.is_empty or float(cap.area) <= 0.0:
            continue
        out[key] = ensure_polygon(cap, allow_empty=False)

    if abs(twist) > 1e-10:
        out = {k: rotate_geom(v, twist, cx=0.0, cy=0.0) for k, v in out.items()}
    return out
