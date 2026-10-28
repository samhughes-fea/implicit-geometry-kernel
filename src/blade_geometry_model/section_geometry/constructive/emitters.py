"""Recipe-faithful constructive polygons for wall-solid subcomponents."""

from __future__ import annotations

from typing import Any

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from .helpers import (
    ensure_polygon,
    rotate_geom,
    stadium_polygon,
)
from .offset_locus import airfoil_polygon, offset_polygon_inward


def _clip_algebra():
    from ..sections.clip_algebra import (
        CamberSide,
        ChordBand,
        FlapwiseStrip,
        apply_shapely,
    )

    return CamberSide, ChordBand, FlapwiseStrip, apply_shapely


def _af_from_component(comp: Any):
    af = getattr(comp, "_af", None)
    if af is not None:
        return af
    raise ValueError(f"{type(comp).__name__} lacks airfoil reference (_af).")


def emit_outer_skin_polygon(
    comp: Any,
    *,
    airfoil_poly: Polygon | MultiPolygon | None = None,
    inner_cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """Skin annulus grown inward: airfoil polygon minus inner(−t) offset.

    When ``comp.skin_side`` is ``SS``/``PS``, the annulus is clipped to the
    suction (camber-upper) or pressure (camber-lower) half.
    """
    af = _af_from_component(comp)
    t = float(getattr(comp, "thickness"))
    if t <= 0:
        raise ValueError("OuterSkin thickness must be positive.")
    poly = airfoil_poly if airfoil_poly is not None else airfoil_polygon(af)
    inner = (
        inner_cavity
        if inner_cavity is not None
        else offset_polygon_inward(poly, t)
    )
    if inner.is_empty:
        annulus = ensure_polygon(poly)
    else:
        annulus = ensure_polygon(poly.difference(inner))

    side = getattr(comp, "skin_side", None)
    if side is None:
        return annulus
    s = str(side).strip().lower()
    if s in ("ss", "suction_surface", "upper"):
        surface = "upper"
    elif s in ("ps", "pressure_surface", "lower"):
        surface = "lower"
    else:
        return annulus
    camber_xy = getattr(comp, "_camber_xy", None)
    CamberSide, _ChordBand, _FlapwiseStrip, apply_shapely = _clip_algebra()
    return apply_shapely(
        CamberSide(surface=surface),  # type: ignore[arg-type]
        annulus,
        camber_xy=camber_xy,
    )


def emit_shear_web_polygon(
    comp: Any,
    *,
    resolution: int = 16,
    cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """Stadium capsule clipped to the interior of the inner skin."""
    # Prefer the extended capsule axis (round ends land outside the cavity so
    # the clipped web is a clean constant-width strip flush against the skin).
    x_top = float(getattr(comp, "x_top_ext", comp.x_top))
    y_top = float(getattr(comp, "y_top_ext", comp.y_top))
    x_bot = float(getattr(comp, "x_bot_ext", comp.x_bot))
    y_bot = float(getattr(comp, "y_bot_ext", comp.y_bot))
    r = 0.5 * float(comp.thickness)
    stadium = stadium_polygon(x_top, y_top, x_bot, y_bot, r, resolution=resolution)

    # Inner-skin clip: Offset(airfoil, -skin).
    if cavity is not None:
        if not cavity.is_empty:
            stadium = stadium.intersection(cavity)
    else:
        af = getattr(comp, "_af", None)
        skin_thickness = getattr(comp, "skin_thickness", None)
        if af is not None and skin_thickness is not None:
            poly = airfoil_polygon(af)
            cav = offset_polygon_inward(poly, float(skin_thickness))
            if not cav.is_empty:
                stadium = stadium.intersection(cav)
    # Thick skin / shallow sections can empty the cavity ∩ stadium Boolean.
    return ensure_polygon(stadium, allow_empty=True)


def _camber_clip(geom, camber_xy, *, surface: str):
    """Clip to the cap's own side of the sampled cavity midline polyline."""
    cam = np.asarray(camber_xy, dtype=float)
    minx, miny, maxx, maxy = geom.bounds
    pad = max(maxx - minx, maxy - miny, 1.0) + 1.0
    xs = np.concatenate(([minx - pad], cam[:, 0], [maxx + pad]))
    ys = np.concatenate(([cam[0, 1]], cam[:, 1], [cam[-1, 1]]))
    y_far = (maxy + pad) if surface == "upper" else (miny - pad)
    ring = [(x, y) for x, y in zip(xs, ys)] + [(xs[-1], y_far), (xs[0], y_far)]
    return geom.intersection(Polygon(ring))


def emit_spar_cap_polygon(
    comp: Any,
    *,
    airfoil_poly: Polygon | MultiPolygon | None = None,
    inner_cavity: Polygon | MultiPolygon | None = None,
    deeper_cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """Slab against inner skin, camber/half-plane + chord/strip clipped, local twist applied."""
    af = _af_from_component(comp)
    skin_t = float(comp.skin_thickness)
    cap_h = float(comp.cap_height)
    surface_raw = str(comp.surface)
    s = surface_raw.strip().lower()
    if s in ("ss", "suction_surface", "upper"):
        surface = "upper"
    elif s in ("ps", "pressure_surface", "lower"):
        surface = "lower"
    else:
        surface = surface_raw
    poly = airfoil_poly if airfoil_poly is not None else airfoil_polygon(af)
    inner = (
        inner_cavity
        if inner_cavity is not None
        else offset_polygon_inward(poly, skin_t)
    )
    deeper = (
        deeper_cavity
        if deeper_cavity is not None
        else offset_polygon_inward(poly, skin_t + cap_h)
    )
    if inner.is_empty:
        raise ValueError("SparCap: inner-skin offset is empty.")
    slab = inner if deeper.is_empty else inner.difference(deeper)
    CamberSide, ChordBand, FlapwiseStrip, apply_shapely = _clip_algebra()
    camber_xy = getattr(comp, "_camber_xy", None)
    slab = apply_shapely(
        CamberSide(surface=surface),  # type: ignore[arg-type]
        slab,
        camber_xy=camber_xy,
    )
    if slab.is_empty:
        # The side clip consumed the whole band: this cap has no material on
        # this surface. Both clip branches below build their clip region from
        # ``slab.bounds``, and Shapely raises an opaque LinearRing error on the
        # empty bounds rather than returning an empty geometry.
        return Polygon()
    # Prefer CapWall flapwise_* (same contract as SandwichCore) so recipe-driven
    # emission stays flush to web midplanes under twist (SG-011). Explicit
    # ``_strip_clip`` still requires the flapwise params for the Shapely path.
    strip_clip = getattr(comp, "_strip_clip", None)
    fw_nx = getattr(comp, "flapwise_nx", None)
    fw_ny = getattr(comp, "flapwise_ny", None)
    fw_lo = getattr(comp, "flapwise_c_lo", None)
    fw_hi = getattr(comp, "flapwise_c_hi", None)
    mode = str(getattr(comp, "flapwise_clip_mode", "strip") or "strip")
    use_flapwise = (
        fw_nx is not None
        and fw_ny is not None
        and fw_lo is not None
        and fw_hi is not None
        and mode == "strip"
    )
    if strip_clip is not None or use_flapwise:
        if fw_nx is None or fw_ny is None or fw_lo is None or fw_hi is None:
            raise ValueError(
                "SparCap flapwise clip needs flapwise_nx/ny/c_lo/c_hi/clip_mode; "
                "construct via CapWall / spar_family or pass them explicitly."
            )
        slab = apply_shapely(
            FlapwiseStrip(
                nx=float(fw_nx),
                ny=float(fw_ny),
                c_lo=fw_lo,
                c_hi=fw_hi,
                clip_mode=mode,  # type: ignore[arg-type]
            ),
            slab,
        )
        if slab.is_empty:
            return Polygon()
    else:
        slab = apply_shapely(
            ChordBand(x_start=float(comp.x_start), x_end=float(comp.x_end)),
            slab,
        )

    twist = float(getattr(comp, "twist_angle", 0.0) or 0.0)
    if abs(twist) > 1e-10:
        cx = 0.5 * (float(comp.x_start) + float(comp.x_end))
        slab = rotate_geom(slab, twist, cx=cx, cy=0.0)
    # Cap/cavity clips can collapse on heavily thickened tip stations.
    return ensure_polygon(slab, allow_empty=True)


def emit_inner_cavity_polygon(
    airfoil_sdf,
    skin_thickness: float,
    *,
    airfoil_poly: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """Interior of the inner skin: Offset(airfoil, -skin)."""
    poly = airfoil_poly if airfoil_poly is not None else airfoil_polygon(airfoil_sdf)
    cavity = offset_polygon_inward(poly, float(skin_thickness))
    # Fully closed sections (skin ≥ local thickness) yield an empty cavity.
    return ensure_polygon(cavity, allow_empty=True)


def emit_sandwich_core_polygon(
    comp: Any,
    wall_exclusion_union: Polygon | MultiPolygon | None,
    *,
    cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """Core bay: (inner cavity − wall solids) ∩ bay clip.

    Matches SDF ``SandwichCore``: Offset(airfoil,-skin/2) minus wall-solid
    exclusions, then chordwise or flapwise bay clip.
    """
    af = _af_from_component(comp)
    skin_t = float(comp.skin_thickness)
    if cavity is None:
        cavity = emit_inner_cavity_polygon(af, skin_t)

    region = cavity
    if wall_exclusion_union is not None and not wall_exclusion_union.is_empty:
        region = region.difference(wall_exclusion_union)

    bay_mode = getattr(comp, "bay_mode", "chordwise")
    _CamberSide, ChordBand, FlapwiseStrip, apply_shapely = _clip_algebra()
    target = region if not region.is_empty else cavity
    if bay_mode == "flapwise":
        clipped = apply_shapely(
            FlapwiseStrip(
                nx=float(comp.flapwise_nx),
                ny=float(comp.flapwise_ny),
                c_lo=comp.flapwise_c_lo,
                c_hi=comp.flapwise_c_hi,
                clip_mode=str(comp.flapwise_clip_mode or "strip"),  # type: ignore[arg-type]
            ),
            target,
        )
    else:
        x_start = comp.x_start
        x_end = comp.x_end
        # Open ends: omit chord clip when both bounds are None.
        if x_start is None and x_end is None:
            clipped = target
        else:
            clipped = apply_shapely(
                ChordBand(
                    x_start=float(x_start) if x_start is not None else -1e9,
                    x_end=float(x_end) if x_end is not None else 1e9,
                ),
                target,
            )
        if getattr(comp, "y_min", None) is not None or getattr(comp, "y_max", None) is not None:
            if clipped.is_empty:
                return Polygon()
            minx, miny, maxx, maxy = clipped.bounds
            pad = max(maxx - minx, maxy - miny, 1.0) + 1.0
            y0 = float(comp.y_min) if comp.y_min is not None else miny - pad
            y1 = float(comp.y_max) if comp.y_max is not None else maxy + pad
            ybox = Polygon(
                [(minx - pad, y0), (maxx + pad, y0), (maxx + pad, y1), (minx - pad, y1)]
            )
            clipped = clipped.intersection(ybox)

    if clipped.is_empty or float(getattr(clipped, "area", 0.0) or 0.0) <= 0.0:
        return Polygon()
    return ensure_polygon(clipped, allow_empty=True)


def emit_trailing_edge_insert_polygon(
    comp: Any,
    *,
    airfoil_poly: Polygon | MultiPolygon | None = None,
    inner_cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """TE insert: Offset(airfoil, -skin) ∩ {x >= x_start}."""
    af = _af_from_component(comp)
    skin_t = float(comp.skin_thickness)
    x_start = float(comp.x_start)
    cavity = (
        inner_cavity
        if inner_cavity is not None
        else emit_inner_cavity_polygon(af, skin_t, airfoil_poly=airfoil_poly)
    )
    if cavity.is_empty:
        return Polygon()
    # Half-plane x >= x_start via ChordBand from x_start to far TE.
    _CamberSide, ChordBand, _FlapwiseStrip, apply_shapely = _clip_algebra()
    minx, miny, maxx, maxy = cavity.bounds
    pad = max(maxx - minx, maxy - miny, 1.0) + 1.0
    region = apply_shapely(
        ChordBand(x_start=x_start, x_end=maxx + pad),
        cavity,
    )
    return ensure_polygon(region, allow_empty=True)


def emit_leading_edge_insert_polygon(
    comp: Any,
    *,
    airfoil_poly: Polygon | MultiPolygon | None = None,
    inner_cavity: Polygon | MultiPolygon | None = None,
) -> Polygon | MultiPolygon:
    """LE insert: Offset(airfoil, -skin) ∩ {x <= x_end} ∩ Circle(le, r)."""
    from shapely.geometry import Point

    af = _af_from_component(comp)
    skin_t = float(comp.skin_thickness)
    x_end = float(comp.x_end)
    cavity = (
        inner_cavity
        if inner_cavity is not None
        else emit_inner_cavity_polygon(af, skin_t, airfoil_poly=airfoil_poly)
    )
    if cavity.is_empty:
        return Polygon()
    _CamberSide, ChordBand, _FlapwiseStrip, apply_shapely = _clip_algebra()
    minx, miny, maxx, maxy = cavity.bounds
    pad = max(maxx - minx, maxy - miny, 1.0) + 1.0
    region = apply_shapely(
        ChordBand(x_start=minx - pad, x_end=x_end),
        cavity,
    )
    # Circle clip matches IR Circle(cx, cy, r).
    le_x = float(getattr(comp, "le_x", 0.0) or 0.0)
    le_y = float(getattr(comp, "le_y", 0.0) or 0.0)
    # LEInsert stores radius only via IR; recover from x_end default if absent.
    radius = getattr(comp, "radius", None)
    if radius is None:
        radius = x_end
    disk = Point(le_x, le_y).buffer(float(radius), quad_segs=32)
    region = region.intersection(disk)
    return ensure_polygon(region, allow_empty=True)


def emit_component_polygon(comp: Any) -> Polygon | MultiPolygon:
    """Dispatch constructive polygon emission by component type / attributes."""
    name = type(comp).__name__
    if name == "OuterSkin" or (
        hasattr(comp, "thickness") and hasattr(comp, "midline_polyline") and hasattr(comp, "_af")
        and not hasattr(comp, "x_top") and not hasattr(comp, "cap_height")
    ):
        if hasattr(comp, "cap_height"):
            return emit_spar_cap_polygon(comp)
        if hasattr(comp, "x_top"):
            return emit_shear_web_polygon(comp)
        return emit_outer_skin_polygon(comp)
    if name == "ShearWeb" or hasattr(comp, "x_top"):
        return emit_shear_web_polygon(comp)
    if name == "SparCap" or hasattr(comp, "cap_height"):
        return emit_spar_cap_polygon(comp)
    if name == "TEInsert" or (
        hasattr(comp, "x_start") and hasattr(comp, "skin_thickness") and name.endswith("Insert")
        and "TE" in name
    ):
        return emit_trailing_edge_insert_polygon(comp)
    if name == "LEInsert" or (
        hasattr(comp, "x_end") and hasattr(comp, "skin_thickness") and name.endswith("Insert")
        and "LE" in name
    ):
        return emit_leading_edge_insert_polygon(comp)
    if name == "SandwichCore" or hasattr(comp, "bay_mode"):
        raise TypeError(
            "SandwichCore requires emit_sandwich_core_polygon with wall-solid exclusions."
        )
    raise TypeError(f"No constructive polygon emitter for {name}.")
