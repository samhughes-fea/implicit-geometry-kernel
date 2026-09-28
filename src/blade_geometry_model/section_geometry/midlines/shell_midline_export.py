"""
shell_midline_export
====================
Export ``section_geometry`` subcomponent midlines as FE-agnostic
:class:`ShellMidlineStrip` records in the blade (B) frame.

**Dependency rule:** this module may import from ``blade_structure_utilities.geometry_utilities.section_contract``
and from ``blade_structure_utilities.geometry_utilities.section_geometry``, but must **never** import from
``blade_structure_utilities.fem_utilities.section_shell_fea``.

**Strip order:** ``build_shell_midline_strips`` yields skins and webs in
``components_unrotated`` dict-insertion order, then spar caps in
``spar_cap_components_unrotated`` dict-insertion order.  Any consumer that
cares about a canonical ordering (e.g. ``skin → caps → webs`` as required by
``bind_shell_material_zones``) must sort or reorder independently.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from blade_structure_utilities.geometry_utilities.section_contract.shell_midline_strip import ShellMidlineStrip
from blade_structure_utilities.geometry_utilities.section_geometry.midlines.section_midline_host import (
    SectionMidlineHost,
    SkinShellComponentLike,
    WebShellComponentLike,
)


def _skin_side_from_label(label: object) -> str | None:
    """``suction_surface`` / ``pressure_surface`` from a skin label, or ``None``.

    Accepts long labels and legacy ``SS`` / ``PS`` / ``outer_skin_*`` tails.
    """
    key = str(getattr(label, "value", label)).lower()
    if "suction" in key or key.endswith("_ss") or key == "ss":
        return "suction_surface"
    if "pressure" in key or key.endswith("_ps") or key == "ps":
        return "pressure_surface"
    return None


def build_shell_midline_strips(
    section: SectionMidlineHost,
    *,
    twist_rad: float,
    n_web_samples: int = 20,
    n_cap_samples: int = 80,
    chord_m: float | None = None,
    pivot_fraction: float | None = None,
) -> tuple[ShellMidlineStrip, ...]:
    """Extract per-subcomponent midline strips from a :class:`SectionMidlineHost`.

    Iterates ``components_unrotated`` (skin + webs) then
    ``spar_cap_components_unrotated`` (caps), calls ``midline_polyline()``
    on each subcomponent in its chord frame, and maps each polyline into the
    canonical **B frame** via
    :func:`coordinate_frames.naca_to_twisted_pitch_axis_blade_frame`
    (short alias ``chord_to_blade``).

    The map is a component **swap** plus a **pivot shift** plus the twist
    rotation, and all three matter. The SDF interior stays in the raw-NACA
    chord frame by design (``AirfoilSDF`` is built on
    ``(x_chord, y_thickness)`` with the LE at ``x = 0``); this function is the
    authoring→handoff boundary. Output arrays are ``(X_B, Y_B)`` with
    ``Y_LE < 0 < Y_TE`` and the pivot at the origin.

    A bare S→B rotation (``coordinate_frames.rotate_S_to_B``) is **not**
    sufficient alone — it leaves the array chordwise-first with its origin at
    the LE, which is the rotated-NACA adapter frame, not B. Everything
    downstream (``+x3`` orientation, the flap/edge load-channel map, the
    section plot contract) reads canonical B.

    Parameters
    ----------
    section
        Section implementing :class:`SectionMidlineHost`.
    twist_rad
        Section twist in radians applied to every chord-frame midline.
    n_web_samples, n_cap_samples
        Number of sample points for web / cap polylines.  Skin midline
        resolution is set by the underlying airfoil tessellation.
    chord_m
        Section chord [m] for the pivot shift. Defaults to
        ``section.airfoil.chord``.
    pivot_fraction
        Pitch-axis chordline as a fraction of chord from the LE. Defaults to
        ``coordinate_frames.PIVOT_FRACTION_DEFAULT``.

    Returns
    -------
    tuple[ShellMidlineStrip, ...]
        Strips in dict-insertion order: skins + webs first, then caps.
    """
    from blade_structure_utilities.geometry_utilities.coordinate_frames import (
        PIVOT_FRACTION_DEFAULT,
        naca_to_twisted_pitch_axis_blade_frame,
    )

    try:
        components_S = section.components_unrotated
    except AttributeError:
        components_S = {}

    if chord_m is None:
        airfoil = getattr(section, "airfoil", None)
        chord_m = float(getattr(airfoil, "chord", 0.0) or 0.0)
    chord_m = float(chord_m)
    if chord_m <= 0.0:
        raise ValueError(
            "chord_m is required for the chord→B pivot shift and could not be read "
            "from section.airfoil.chord; pass chord_m explicitly."
        )
    f_pivot = (
        PIVOT_FRACTION_DEFAULT if pivot_fraction is None else float(pivot_fraction)
    )

    def _to_twisted_pitch_axis_blade_frame(mid_zero_twist_pitch_axis_station_frame):
        return naca_to_twisted_pitch_axis_blade_frame(
            mid_zero_twist_pitch_axis_station_frame, twist_rad, chord_m, f_pivot
        )

    # Web before skin: webs also satisfy SkinShellComponentLike (thickness + midline).
    def _label_key(label: object) -> str:
        """Stable string id; prefer Enum.value over ``ComponentLabel.NAME`` repr."""
        value = getattr(label, "value", None)
        if isinstance(value, str) and value:
            return value
        return str(label)

    strips: list[ShellMidlineStrip] = []

    recipes = dict(getattr(section, "wall_recipes", None) or {})
    cap_recipes = dict(getattr(section, "spar_cap_wall_recipes", None) or {})
    airfoil = getattr(section, "airfoil", None)
    if airfoil is None:
        raise ValueError("section must expose airfoil for recipe midline projection.")

    def _is_web_key(key: str) -> bool:
        return key.startswith("shear_web_") or key.startswith("web_")

    def _is_skin_key(key: str) -> bool:
        return (
            key == "outer_skin"
            or key.startswith("skin_suction_surface")
            or key.startswith("skin_pressure_surface")
            or key in {"skin_suction_surface", "skin_pressure_surface"}
        )

    def _is_cap_key(key: str) -> bool:
        return key.startswith("spar_cap")

    # Prefer recipe map; fall back to component iteration only for thickness
    # when a host lacks wall_recipes (should not happen for MultiCellSection).
    skin_web_items: list[tuple[str, Any, str]] = []
    for key, wall in recipes.items():
        if _is_web_key(key):
            skin_web_items.append((key, wall, "web"))
        elif _is_skin_key(key):
            skin_web_items.append((key, wall, "skin"))

    if not skin_web_items:
        for label, comp in components_S.items():
            key = _label_key(label)
            if isinstance(comp, WebShellComponentLike):
                skin_web_items.append((key, None, "web"))
            elif isinstance(comp, SkinShellComponentLike):
                skin_web_items.append((key, None, "skin"))

    # Web before skin (legacy strip order; webs also look skin-like by attrs).
    skin_web_items.sort(key=lambda item: 0 if item[2] == "web" else 1)

    for key, wall, kind in skin_web_items:
        if wall is None:
            raise KeyError(
                f"No {kind} wall recipe for {key!r}; "
                "MultiCellSection must author recipes at assemble."
            )
        if kind == "web":
            mid_S = wall.to_midline(airfoil, n=n_web_samples)
            strips.append(
                ShellMidlineStrip(
                    label=key,
                    kind="web",
                    midline_twisted_pitch_axis_blade_frame=_to_twisted_pitch_axis_blade_frame(mid_S),
                    thickness_m=float(wall.thickness),
                    closed=False,
                    alignment=str(wall.alignment),
                )
            )
        else:
            mid_S = wall.to_midline(airfoil)
            side = _skin_side_from_label(key)
            strips.append(
                ShellMidlineStrip(
                    label=key,
                    kind="skin",
                    midline_twisted_pitch_axis_blade_frame=_to_twisted_pitch_axis_blade_frame(mid_S),
                    thickness_m=float(wall.thickness),
                    closed=side is None,
                    skin_side=side,
                )
            )

    cap_items: list[tuple[str, Any]] = []
    try:
        cap_components = dict(section.spar_cap_components_unrotated or {})
    except AttributeError:
        cap_components = {}
    for label in cap_components:
        key = _label_key(label)
        wall = cap_recipes.get(key) or recipes.get(key)
        if wall is None:
            raise KeyError(
                f"No cap wall recipe for {key!r}; "
                "MultiCellSection must author CapWall recipes at assemble."
            )
        cap_items.append((key, wall))

    for key, wall in cap_items:
        mid_S = wall.to_midline(airfoil, n=n_cap_samples)
        strips.append(
            ShellMidlineStrip(
                label=key,
                kind="cap",
                midline_twisted_pitch_axis_blade_frame=_to_twisted_pitch_axis_blade_frame(mid_S),
                thickness_m=float(wall.cap_height),
                closed=False,
                surface=str(wall.surface),
            )
        )

    return tuple(_insert_web_cap_junction_nodes(strips))


def _insert_web_cap_junction_nodes(
    strips: list[ShellMidlineStrip],
) -> list[ShellMidlineStrip]:
    """Put a node on **both** strips wherever a web crosses a cap.

    Webs and caps are sampled independently, so a junction between them lands
    part-way along a segment of one or both and no shared node exists. Before
    this pass the mesh held together only because ``section_mesh`` clusters
    endpoints within ``ENDPOINT_TOL_FRAC`` — 0.6 % of chord, 6 mm at 1 m —
    with measured gaps of 1.88 mm at cap ends and 1.95 mm at interior
    crossings. That is a tolerance scaled by *chord*, not by mesh density, so
    refining the mesh erodes the margin instead of improving it.

    Both junction kinds are the same computation: intersect the web's axis with
    the cap polyline. Nothing is relocated and no membership tolerance is
    needed — the crossing is solved for rather than searched for.

    Only the **web** is given nodes here. The cap's interior is rebuilt later by
    ``resample_cap_interiors_from_skin_rays``, which would discard anything
    inserted at this stage; the cap side is closed afterwards by
    ``insert_cap_web_interior_junctions``, once its polyline is final.
    Inserting into the web is exact because a web midline is straight.
    """
    web_indices = [i for i, s in enumerate(strips) if s.kind == "web"]
    cap_indices = [i for i, s in enumerate(strips) if s.kind == "cap"]
    if not web_indices or not cap_indices:
        return strips

    points_along_web: dict[int, list[float]] = {}

    for web_index in web_indices:
        web_points = np.asarray(
            strips[web_index].midline_twisted_pitch_axis_blade_frame, dtype=float
        )
        if web_points.shape[0] < 2:
            continue
        start = web_points[0]
        axis = web_points[-1] - start
        axis_length = float(np.linalg.norm(axis))
        if axis_length <= 1e-12:
            continue
        direction = axis / axis_length

        for cap_index in cap_indices:
            cap_points = np.asarray(
                strips[cap_index].midline_twisted_pitch_axis_blade_frame, dtype=float
            )
            if cap_points.shape[0] < 2:
                continue
            segment_start = cap_points[:-1]
            segment = cap_points[1:] - segment_start

            denominator = (
                segment[:, 0] * direction[1] - segment[:, 1] * direction[0]
            )
            usable = np.abs(denominator) > 1e-14
            if not np.any(usable):
                continue
            safe = np.where(usable, denominator, 1.0)
            offset = start - segment_start
            along_segment = (
                offset[:, 0] * direction[1] - offset[:, 1] * direction[0]
            ) / safe
            crossing = segment_start + along_segment[:, None] * segment
            along_axis = (crossing - start) @ direction

            hit = (
                usable
                & (along_segment >= -1e-12)
                & (along_segment <= 1.0 + 1e-12)
                & (along_axis >= -1e-12)
                & (along_axis <= axis_length + 1e-12)
            )
            for segment_index in np.flatnonzero(hit):
                t_axis = float(along_axis[segment_index])
                u_segment = float(np.clip(along_segment[segment_index], 0.0, 1.0))
                crossing_point = (
                    segment_start[segment_index]
                    + u_segment * segment[segment_index]
                )
                # Only junctions at a cap **end** are handled here. Interior
                # crossings must not be: inserting a web vertex exactly on the
                # cap makes the two polylines touch rather than cross, and
                # `insert_cap_web_interior_junctions` — which finds interior
                # junctions by polyline-polyline intersection after the cap has
                # been resampled — then stops finding them. Measured: doing it
                # here moved the web gap to 0 but pushed the cap gap from 0 to
                # 1056 um on 4D-F. See docs/bug_log.md.
                at_cap_end = (
                    min(
                        float(np.linalg.norm(crossing_point - cap_points[0])),
                        float(np.linalg.norm(crossing_point - cap_points[-1])),
                    )
                    <= 1e-9
                )
                if at_cap_end and 1e-9 < t_axis < axis_length - 1e-9:
                    points_along_web.setdefault(web_index, []).append(t_axis)

    if not points_along_web:
        return strips

    out = list(strips)

    for web_index, parameters in points_along_web.items():
        web_points = np.asarray(
            strips[web_index].midline_twisted_pitch_axis_blade_frame, dtype=float
        )
        start = web_points[0]
        direction = web_points[-1] - start
        direction = direction / float(np.linalg.norm(direction))
        existing = (web_points - start) @ direction
        merged = np.unique(
            np.round(np.concatenate([existing, np.asarray(parameters)]), 12)
        )
        out[web_index] = replace(
            out[web_index],
            midline_twisted_pitch_axis_blade_frame=start
            + merged[:, None] * direction,
        )

    return out


__all__ = [
    "build_shell_midline_strips",
    "SectionMidlineHost",
    "ShellComponentLike",
    "SkinShellComponentLike",
    "WebShellComponentLike",
    "CapShellComponentLike",
]
