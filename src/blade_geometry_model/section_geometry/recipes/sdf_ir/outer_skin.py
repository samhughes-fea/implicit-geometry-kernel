"""Outer skin subcomponent (full annulus or SS/PS half-annulus)."""

from __future__ import annotations

import numpy as np

from blade_geometry_model.section_geometry.sdf.ir import (
    CallableField,
    Intersect,
    Offset,
    Slab,
)

from blade_geometry_model.section_geometry.sections.clip_algebra import (
    airfoil_expr,
    camber_side_field,
    lower_half_plane_expr,
    upper_half_plane_expr,
)
from blade_geometry_model.section_geometry.sections.labels import (
    ComponentLabel,
)
from .base import SectionComponent, expr_callable

# Blade-side vocabulary: suction = camber-upper, pressure = camber-lower.
_SKIN_SIDE_TO_SURFACE = {
    "SS": "upper",
    "PS": "lower",
    "SUCTION_SURFACE": "upper",
    "PRESSURE_SURFACE": "lower",
    "suction_surface": "upper",
    "pressure_surface": "lower",
}
_SURFACE_TO_SKIN_SIDE = {
    "upper": "suction_surface",
    "lower": "pressure_surface",
}


def _normalise_skin_side(skin_side: str | None) -> str | None:
    if skin_side is None:
        return None
    s = str(skin_side).strip()
    key = s.upper() if s.upper() in ("SS", "PS") else s
    if key not in _SKIN_SIDE_TO_SURFACE and s not in _SKIN_SIDE_TO_SURFACE:
        raise ValueError(
            "skin_side must be 'suction_surface', 'pressure_surface', "
            "'SS', 'PS', or None."
        )
    surface = _SKIN_SIDE_TO_SURFACE.get(s, _SKIN_SIDE_TO_SURFACE.get(key))
    return "suction_surface" if surface == "upper" else "pressure_surface"


def _closed_mid_ring(airfoil_sdf, thickness: float, *, mid_depth: float | None = None) -> np.ndarray:
    """Closed mid-surface ring from ``buffer(-mid_depth)`` (may omit aero tip)."""
    from shapely.geometry import MultiPolygon, Polygon

    from blade_geometry_model.section_geometry.constructive.offset_locus import (
        airfoil_polygon,
        offset_polygon_inward,
    )

    verts = getattr(airfoil_sdf, "vertices", None)
    if verts is None:
        raise AttributeError(
            "OuterSkin midline requires the underlying airfoil SDF to expose "
            "a 'vertices' attribute (AirfoilSDF-compatible)."
        )
    v = np.asarray(verts, dtype=float)
    if v.ndim != 2 or v.shape[0] < 3 or v.shape[1] != 2:
        raise ValueError(f"Invalid airfoil vertices shape {v.shape}.")

    t = float(thickness)
    if t <= 0.0:
        raise ValueError("OuterSkin thickness must be positive.")
    depth = 0.5 * t if mid_depth is None else float(mid_depth)
    if depth <= 0.0:
        raise ValueError(f"OuterSkin mid-surface depth must be positive, got {depth!r}.")

    poly = airfoil_polygon(airfoil_sdf)
    mid_geom = offset_polygon_inward(poly, depth)
    if mid_geom.is_empty:
        raise ValueError(
            f"OuterSkin mid-surface buffer(-{depth:g}) is empty for thickness={t:g}."
        )
    if isinstance(mid_geom, MultiPolygon):
        mid_geom = max(mid_geom.geoms, key=lambda g: g.area)
    if not isinstance(mid_geom, Polygon) or mid_geom.exterior is None:
        raise ValueError(
            f"OuterSkin mid-surface buffer(-{depth:g}) did not yield a Polygon "
            f"(got {type(mid_geom).__name__})."
        )
    mid = np.asarray(mid_geom.exterior.coords, dtype=float)
    if mid.shape[0] < 3:
        raise ValueError("OuterSkin mid-surface buffer produced too few vertices.")
    if float(np.hypot(*(mid[0] - mid[-1]))) <= 1e-15:
        mid = mid[:-1]
    if mid.shape[0] < 3:
        raise ValueError("OuterSkin mid-surface buffer produced too few unique vertices.")
    return mid


def split_closed_ring_ss_ps(ring: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split a closed mid-ring into open SS (TE→LE) and PS (LE→TE) polylines."""
    pts = np.asarray(ring, dtype=float)
    if pts.shape[0] < 3:
        raise ValueError("Need at least 3 ring vertices to split SS/PS.")
    if float(np.hypot(*(pts[0] - pts[-1]))) <= 1e-15:
        pts = pts[:-1]
    n = int(pts.shape[0])
    i_le = int(np.argmin(pts[:, 0]))
    i_te = int(np.argmax(pts[:, 0]))
    if i_le == i_te:
        raise ValueError("Degenerate mid-ring: LE and TE indices coincide.")

    def _arc(i0: int, i1: int) -> np.ndarray:
        if i0 == i1:
            return pts[i0 : i0 + 1].copy()
        out = [pts[i0]]
        i = i0
        while i != i1:
            i = (i + 1) % n
            out.append(pts[i])
        return np.asarray(out, dtype=float)

    a = _arc(i_te, i_le)
    b = _arc(i_le, i_te)
    # SS = suction = camber-upper path TE→LE (higher mean y).
    if float(np.mean(a[:, 1])) >= float(np.mean(b[:, 1])):
        ss, ps = a, b
    else:
        ss, ps = b[::-1].copy(), a[::-1].copy()
        # Ensure SS is TE→LE (increasing then decreasing x toward LE) and PS LE→TE.
        if ss[0, 0] < ss[-1, 0]:
            ss = ss[::-1].copy()
        if ps[0, 0] > ps[-1, 0]:
            ps = ps[::-1].copy()
    # Normalise orientation: SS TE→LE (start x >= end x), PS LE→TE (start x <= end x).
    if ss[0, 0] < ss[-1, 0]:
        ss = ss[::-1].copy()
    if ps[0, 0] > ps[-1, 0]:
        ps = ps[::-1].copy()
    return ss, ps


class OuterSkin(SectionComponent):
    """Skin wall of uniform thickness grown from a fixed reference face.

    Default placement (``outer_mould_line`` / inward): the supplied airfoil
    contour is the aerodynamic *outer* surface; the skin occupies the band
    ``-t <= phi <= 0`` so the outer shape never grows with skin thickness, and
    the midline sits at ``buffer(-t/2)``.

    Named non-default reference faces (WP3):
    * ``inner_mould_line`` — IML fixed at ``-reference_iml_depth_m`` from the
      authoring OML; midline at ``iml - t/2`` (toward the fluid);
    * ``midsurface`` — midline held at the authoring mid depth while both faces
      move (solid still authored from the OML for SDF compatibility).

    When ``skin_side`` is ``"SS"`` / ``"PS"``, the solid is the half-annulus on
    the suction (camber-upper) or pressure (camber-lower) side.
    """

    label = ComponentLabel.OUTER_SKIN

    def __init__(
        self,
        airfoil_sdf,
        thickness,
        *,
        airfoil_expr_node=None,
        skin_side: str | None = None,
        camber_xy=None,
        reference_face: str = "outer_mould_line",
        reference_iml_depth_m: float | None = None,
    ):
        from blade_geometry_model.section_geometry.placement.anchoring import (
            airfoil_inward_depth_for_midsurface,
        )

        self._af = airfoil_sdf
        self.thickness = float(thickness)
        self.reference_face = str(reference_face).strip().lower()
        self.reference_iml_depth_m = (
            None if reference_iml_depth_m is None else float(reference_iml_depth_m)
        )
        self._mid_depth = airfoil_inward_depth_for_midsurface(
            self.reference_face,
            self.thickness,
            reference_iml_depth_m=float(self.reference_iml_depth_m or 0.0),
        )
        side = _normalise_skin_side(skin_side)
        self.skin_side = side
        if side == "suction_surface":
            self.label = ComponentLabel.SKIN_SUCTION_SURFACE
        elif side == "pressure_surface":
            self.label = ComponentLabel.SKIN_PRESSURE_SURFACE
        else:
            self.label = ComponentLabel.OUTER_SKIN

        poly = airfoil_expr_node if airfoil_expr_node is not None else airfoil_expr(airfoil_sdf)
        if poly is None:
            raise ValueError("OuterSkin requires an AirfoilSDF with vertices for IR evaluation.")

        slab = Slab(poly, depth=self.thickness)
        if side is None:
            expr = slab
            self._camber_xy = None
        else:
            surface = _SKIN_SIDE_TO_SURFACE[side]
            cam = None if camber_xy is None else np.asarray(camber_xy, dtype=float)
            if cam is None:
                try:
                    xc, yc = airfoil_sdf.camber_line(n_points=120)
                    cam = np.column_stack([np.asarray(xc, dtype=float), np.asarray(yc, dtype=float)])
                except Exception:
                    cam = None
            self._camber_xy = cam
            if cam is not None and cam.shape[0] >= 2:
                side_expr = CallableField(
                    key=f"skin_side_{side}_{id(self)}",
                    fn=camber_side_field(cam, surface),
                )
            else:
                side_expr = (
                    upper_half_plane_expr() if surface == "upper" else lower_half_plane_expr()
                )
            expr = Intersect(slab, side_expr)

        super().__init__(expr)
        self._poly = poly
        self.inner_sdf = expr_callable(Offset(poly, amount=-self.thickness))
        self.outer_sdf = expr_callable(Offset(poly, amount=0.0))

    def midline_polyline(self):
        """Return the skin mid-surface polyline.

        * Full skin (``skin_side is None``): closed ring from the placement mid depth.
        * Suction/pressure: open half-ring (suction: TE→LE, pressure: LE→TE).

        Where local airfoil thickness drops below ``t``, no constant-thickness
        mid-surface exists, so the locus ends short of the aerodynamic tip by
        construction — not by post-hoc trim.
        """
        mid = _closed_mid_ring(self._af, self.thickness, mid_depth=self._mid_depth)
        if self.skin_side is None:
            return np.vstack([mid, mid[0:1]])
        ss, ps = split_closed_ring_ss_ps(mid)
        return ss if self.skin_side == "suction_surface" else ps

    def __repr__(self):
        side = self.skin_side or "full"
        return (
            f"OuterSkin(thickness={self.thickness:.4g}, skin_side={side!r}, "
            f"reference_face={self.reference_face!r})"
        )


__all__ = [
    "OuterSkin",
    "split_closed_ring_ss_ps",
    "_SKIN_SIDE_TO_SURFACE",
    "_SURFACE_TO_SKIN_SIDE",
]
