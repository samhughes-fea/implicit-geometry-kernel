"""
shell_midline_export
====================
Export ``section_geometry`` subcomponent midlines as FE-agnostic
:class:`ShellMidlineStrip` records in the blade (B) frame.

**Dependency rule:** this module may import from ``implicit_geometry_kernel.contract``
and from ``implicit_geometry_kernel.section``, but must **never** import a mesh
or FE package.

**Strip order:** ``build_shell_midline_strips`` yields skins and webs in
``components_unrotated`` dict-insertion order. Caps are partitions on the
skin-mid loop (``shell_strips_from_graph`` / ``offset_regions``), not a
second polyline.
"""

from __future__ import annotations

from typing import Any

from implicit_geometry_kernel.contract.shell_midline_strip import ShellMidlineStrip
from implicit_geometry_kernel.section.midlines.host import (
    SectionMidlineHost,
    SkinShellComponentLike,
    WebShellComponentLike,
)

def _label_key(label: object) -> str:
    """Stable string id; prefer Enum.value over ``ComponentLabel.NAME`` repr."""
    value = getattr(label, "value", None)
    if isinstance(value, str) and value:
        return value
    return str(label)


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

    Iterates ``components_unrotated`` (skin + webs), calls ``midline_polyline()``
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
    n_web_samples
        Number of sample points for web polylines. Skin midline resolution
        is set by the underlying airfoil tessellation.
    n_cap_samples
        Ignored. Cap footprints are offset regions on the skin mid, not strips.
    chord_m
        Section chord [m] for the pivot shift. Defaults to
        ``section.airfoil.chord``.
    pivot_fraction
        Pitch-axis chordline as a fraction of chord from the LE. Defaults to
        ``coordinate_frames.PIVOT_FRACTION_DEFAULT``.

    Returns
    -------
    tuple[ShellMidlineStrip, ...]
        Skin and web strips (web before skin). No extra cap polyline.
    """
    from implicit_geometry_kernel.frames import (
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
    strips: list[ShellMidlineStrip] = []

    _ = n_cap_samples
    recipes = dict(getattr(section, "wall_recipes", None) or {})
    airfoil = getattr(section, "airfoil", None)
    if airfoil is None:
        raise ValueError("section must expose airfoil for recipe midline projection.")

    from implicit_geometry_kernel.section.assembly.members import wall_type_for_label

    # Prefer recipe map; fall back to component iteration only for thickness
    # when a host lacks wall_recipes (should not happen for MultiCellSection).
    skin_web_items: list[tuple[str, Any, str]] = []
    for key, wall in recipes.items():
        wt = wall_type_for_label(key)
        if wt is None or not wt.has_midline or not wt.on_shell_product:
            continue
        role = wt.cell_face_kind or wt.kind
        skin_web_items.append((key, wall, role))

    if not skin_web_items:
        for label, _comp in components_S.items():
            key = _label_key(label)
            wt = wall_type_for_label(key)
            if wt is None or not wt.has_midline or not wt.on_shell_product:
                continue
            skin_web_items.append((key, None, wt.cell_face_kind or wt.kind))

    # Web before skin (legacy strip order; webs also look skin-like by attrs).
    skin_web_items.sort(key=lambda item: 0 if item[2] == "web" else 1)

    for key, wall, kind in skin_web_items:
        if wall is None:
            raise KeyError(
                f"No {kind} wall recipe for {key!r}; "
                "MultiCellSection must author recipes at assemble."
            )
        if kind == "skin":
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
        else:
            mid_S = wall.to_midline(airfoil, n=n_web_samples)
            strips.append(
                ShellMidlineStrip(
                    label=key,
                    kind=kind,
                    midline_twisted_pitch_axis_blade_frame=_to_twisted_pitch_axis_blade_frame(mid_S),
                    thickness_m=float(wall.thickness),
                    closed=False,
                    alignment=str(getattr(wall, "alignment", "chord_normal")),
                )
            )

    return tuple(strips)


def build_solid_cap_strips(
    section: SectionMidlineHost,
    *,
    twist_rad: float,
    n_cap_samples: int = 80,
    chord_m: float | None = None,
    pivot_fraction: float | None = None,
) -> tuple[ShellMidlineStrip, ...]:
    """Cap recipe midlines for stacked **solids** on ``cut.walls``.

    Not the conformal mesh product. ``result.midline_strips`` stays skin+web;
    eccentricity for constitutives is ``offset_regions()``, not this polyline.
    """
    from implicit_geometry_kernel.frames import (
        PIVOT_FRACTION_DEFAULT,
        naca_to_twisted_pitch_axis_blade_frame,
    )

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
    airfoil = getattr(section, "airfoil", None)
    if airfoil is None:
        raise ValueError("section must expose airfoil for recipe midline projection.")
    recipes = dict(getattr(section, "wall_recipes", None) or {})
    cap_recipes = dict(getattr(section, "spar_cap_wall_recipes", None) or {})
    try:
        cap_components = dict(section.spar_cap_components_unrotated or {})
    except AttributeError:
        cap_components = {}
    strips: list[ShellMidlineStrip] = []
    for label in cap_components:
        key = _label_key(label)
        wall = cap_recipes.get(key) or recipes.get(key)
        if wall is None:
            raise KeyError(
                f"No cap wall recipe for {key!r}; "
                "MultiCellSection must author CapWall recipes at assemble."
            )
        mid_S = wall.to_midline(airfoil, n=n_cap_samples)
        strips.append(
            ShellMidlineStrip(
                label=key,
                kind="cap",
                midline_twisted_pitch_axis_blade_frame=naca_to_twisted_pitch_axis_blade_frame(
                    mid_S, twist_rad, chord_m, f_pivot
                ),
                thickness_m=float(wall.cap_height),
                closed=False,
                surface=str(wall.surface),
            )
        )
    return tuple(strips)


__all__ = [
    "build_shell_midline_strips",
    "build_solid_cap_strips",
    "SectionMidlineHost",
    "SkinShellComponentLike",
    "WebShellComponentLike",
]
