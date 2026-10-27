"""Full section JSON export (properties, boundaries, optional medial axes)."""

from __future__ import annotations

import json

import numpy as np

from blade_geometry_model.section_geometry.sdf.cache import SectionEvalCache
from blade_geometry_model.section_geometry.sdf.ir import eval_expr_batch

from .properties_report import SectionPropertiesReport


def _multicell_host(section_geometry):
    """Resolve the concrete section graph host (e.g. MultiCellSection)."""
    from blade_geometry_model.section_geometry.midlines.section_midline_host import (
        resolve_section_host,
    )

    return resolve_section_host(section_geometry)


def _zero_contour_to_nested(
    grid,
    sdf_callable,
    *,
    label: str,
    eval_cache: SectionEvalCache | None = None,
):
    """Extract φ=0 as list of polylines (each Nx2), JSON-ready."""
    if eval_cache is None:
        phi = grid.eval(sdf_callable)
    else:
        phi = eval_cache.get_or_eval_with_owner(label, sdf_callable, grid, owner=None)
    segs = grid.zero_contour(phi)
    return [np.asarray(s, dtype=float).tolist() for s in segs]


def _shared_grid_phi(section_geometry, grid) -> dict:
    """Evaluate every export field's φ on the grid with one shared IR memo.

    All component and skin-boundary SDFs embed the same airfoil polygon subtree,
    which dominates the per-field ``grid.eval`` cost. Evaluating them through a
    single shared ``eval_expr_batch`` memo computes that subtree once instead of
    once per field. ``SectionComponent.__call__`` is itself
    ``eval_expr_batch(self.expr, ...)``, so the φ is bit-identical to
    ``grid.eval`` -- only the evaluation memo is shared.

    Returns ``{id(field_callable): phi_grid}`` for fields exposing an IR expr;
    opaque callables are omitted so the caller falls back to ``grid.eval``.
    """
    x = np.asarray(grid.X, dtype=float).ravel()
    y = np.asarray(grid.Y, dtype=float).ravel()
    shape = np.asarray(grid.X).shape
    memo: dict = {}
    out: dict[int, np.ndarray] = {}

    # (id_key, expr) pairs. id_key is the object the caller looks up by id().
    pairs: list[tuple[int, object]] = []
    for label in section_geometry.labels:
        fn = section_geometry[label]
        expr = getattr(fn, "expr", getattr(fn, "_expr", None))
        if expr is not None:
            pairs.append((id(fn), expr))

    # Skin boundaries: the rotate_field-wrapped callables carry no .expr in the
    # twisted case, but the host exposes the equivalent IR nodes (Offset over the
    # shared Rotate(airfoil) subtree). Evaluating those through the shared memo is
    # bit-identical to grid.eval of the callable and reuses the rotated airfoil
    # polygon already computed for the twisted components.
    host = _multicell_host(section_geometry)
    for sdf_attr, expr_attr in (
        ("_skin_outer_boundary_sdf", "_skin_outer_boundary_expr"),
        ("_skin_inner_boundary_sdf", "_skin_inner_boundary_expr"),
    ):
        sdf = getattr(host, sdf_attr, None)
        expr = getattr(host, expr_attr, None)
        if sdf is not None and expr is not None:
            pairs.append((id(sdf), expr))

    for id_key, expr in pairs:
        phi = np.asarray(eval_expr_batch(expr, x, y, memo), dtype=float)
        out[id_key] = phi.reshape(shape)
    return out


def _segments_to_nested(segs: np.ndarray) -> list:
    """Convert ``(m, 2, 2)`` marching-squares segments to JSON polyline lists."""
    if segs is None or len(segs) == 0:
        return []
    arr = np.asarray(segs, dtype=float)
    return [arr[i].tolist() for i in range(arr.shape[0])]


def _mesh_label_boundary(mesh, label: str) -> list:
    """φ=0 boundary polylines for ``label`` from shared quadtree corner φ."""
    from blade_geometry_model.area_2d.quadrature import label_boundary_segments

    def _as_key(lbl) -> str:
        return str(getattr(lbl, "value", lbl))

    phi_map = getattr(mesh, "phi_by_label", None) or {}
    keys = {_as_key(k): k for k in phi_map}
    key = keys.get(_as_key(label))
    if key is None:
        return []
    return _segments_to_nested(label_boundary_segments(mesh, key))


def export_section_json(
    section_geometry,
    grid,
    midline_dict,
    filepath,
    *,
    include_component_zero_contours=False,
    include_geometry_detail=True,
    eval_cache: SectionEvalCache | None = None,
    properties_report: SectionPropertiesReport | None = None,
):
    """Export a complete section description to JSON.

    Includes:
      - Reference airfoil contour vertices (chord-frame ``AirfoilSDF`` polyline)
      - Per-component section properties
      - Medial axis polylines (flat ``medial_axes`` dict, legacy — pass an empty
        ``midline_dict`` unless using deprecated grid-medial extraction)
      - When ``include_geometry_detail`` is True, a nested ``geometry`` block:
        * ``skin``: outer mold line, inner mold line, and ``outer_skin`` medial polylines
        * ``components``: per label, boundary loops (φ=0 of that solid) and medial polylines
      - Optional legacy flat ``component_zero_contours`` (same loops as
        ``geometry["components"][label]["boundary"]`` when both are enabled).

    Parameters
    ----------
    section_geometry : MultiCellSection or SectionMidlineHost
    grid : SDFGrid
        Used for skin-boundary fallback contours and as bbox for property rebuild
        when ``properties_report`` is omitted.
    midline_dict : dict  {label: list of ndarray}
    filepath : str
    include_component_zero_contours : bool
        If True, also write flat ``component_zero_contours`` (matplotlib contours).
    include_geometry_detail : bool
        If True, populate ``geometry`` with skin mold lines and per-component
        boundaries plus medials. If False, medial polylines are still written
        to the flat ``medial_axes`` key when ``midline_dict`` is non-empty.
    properties_report : SectionPropertiesReport, optional
        Reuse an already-computed property report (and its shared quadtree mesh)
        to skip the property rebuild and derive component boundaries from the mesh
        instead of a dense contour grid.
    """
    # Airfoil (reference section shape; not the twisted boundary SDFs)
    af_verts = section_geometry.airfoil.vertices.tolist()

    # Section properties — reuse caller report when provided
    report = properties_report or SectionPropertiesReport(
        section_geometry, grid, eval_cache=eval_cache
    )
    mesh = report.mesh

    def _arr(a):
        return np.asarray(a).tolist()

    props_out = {}
    for label, props in report.to_dict().items():
        props_out[label] = {k: float(v) if k != "label" else v
                            for k, v in props.items()}

    # Medial axes (flat, backward compatible)
    midlines_out = {}
    for label, polylines in midline_dict.items():
        midlines_out[label] = [_arr(p) for p in polylines]

    payload = {
        "schema_version": 2,
        "airfoil_vertices": af_verts,
        "reference_airfoil_vertices": af_verts,
        "section_properties": props_out,
        "medial_axes": midlines_out,
    }

    # Prefer shared-quadtree boundaries for components (no dense 384×180 φ eval).
    # Skin mold-line SDFs are not mesh labels — fall back to a shared-memo grid eval
    # only for those two fields when geometry detail is requested.
    use_mesh_components = mesh is not None and getattr(mesh, "phi_by_label", None)

    _shared_phi = None

    def _ensure_shared_phi():
        nonlocal _shared_phi
        if _shared_phi is None:
            _shared_phi = _shared_grid_phi(section_geometry, grid)
        return _shared_phi

    def _field_contour(sdf_callable, label):
        shared = _ensure_shared_phi()
        phi = shared.get(id(sdf_callable))
        if phi is None:
            return _zero_contour_to_nested(
                grid, sdf_callable, label=label, eval_cache=eval_cache
            )
        segs = grid.zero_contour(phi)
        return [np.asarray(s, dtype=float).tolist() for s in segs]

    _component_boundaries: dict = {}

    def _constructive_boundary(lab: str):
        """Exterior (+ holes) rings from constructive wall-solid polygons when available."""
        from blade_geometry_model.section_geometry.constructive import (
            emit_constructive_polygons,
            is_constructive_label,
        )

        if not is_constructive_label(lab):
            return None
        try:
            polys = emit_constructive_polygons(section_geometry)
        except Exception:
            return None
        geom = polys.get(lab)
        if geom is None or geom.is_empty:
            return None
        rings = []
        geoms = getattr(geom, "geoms", None) or (geom,)
        for g in geoms:
            if g.is_empty:
                continue
            rings.append(np.asarray(g.exterior.coords, dtype=float).tolist())
            for hole in g.interiors:
                rings.append(np.asarray(hole.coords, dtype=float).tolist())
        return rings or None

    def _component_boundary(label):
        lab = str(getattr(label, "value", label))
        if lab not in _component_boundaries:
            constructive = _constructive_boundary(lab)
            if constructive is not None and report.backend_meta.get("backend") == "hybrid":
                _component_boundaries[lab] = constructive
            elif use_mesh_components:
                nested = _mesh_label_boundary(mesh, label)
                if nested:
                    _component_boundaries[lab] = nested
                else:
                    _component_boundaries[lab] = _field_contour(
                        section_geometry[label], lab
                    )
            else:
                _component_boundaries[lab] = _field_contour(
                    section_geometry[label], lab
                )
        return _component_boundaries[lab]

    if include_geometry_detail:
        geometry: dict = {"components": {}}
        host = _multicell_host(section_geometry)
        o_sdf = getattr(host, "_skin_outer_boundary_sdf", None)
        i_sdf = getattr(host, "_skin_inner_boundary_sdf", None)

        skin_block = {}
        if o_sdf is not None and i_sdf is not None:
            skin_block["outer_boundary"] = _field_contour(o_sdf, "_skin_outer_boundary")
            skin_block["inner_boundary"] = _field_contour(i_sdf, "_skin_inner_boundary")
        else:
            skin_block["outer_boundary"] = None
            skin_block["inner_boundary"] = None
        skin_axes = []
        for key in ("skin_suction_surface", "skin_pressure_surface", "outer_skin"):
            skin_axes.extend(midlines_out.get(key, []))
        skin_block["medial_axes"] = skin_axes
        geometry["skin"] = skin_block

        for label in section_geometry.labels:
            lab = str(getattr(label, "value", label))
            geometry["components"][lab] = {
                "boundary": _component_boundary(label),
                "medial_axes": midlines_out.get(lab, midlines_out.get(label, [])),
            }
        payload["geometry"] = geometry

    if include_component_zero_contours:
        payload["component_zero_contours"] = {
            str(getattr(label, "value", label)): _component_boundary(label)
            for label in section_geometry.labels
        }

    with open(filepath, "w") as fh:
        json.dump(payload, fh, indent=2)
    return filepath


__all__ = ["export_section_json", "_shared_grid_phi"]
