"""Spanwise blade geometry: the implicit section geometry lofted between stations.

``BladeGeometry`` takes the ``Station`` list a ``BladeBuilder`` takes (chord, twist,
pitch axis, offset and airfoil at each span position) and a structural layout from
``section`` (one of its 40 system types, such as ``"2D-F-LC-0R"``: two webs and a
box spar cap between them). ``interpolate(z)`` is the only interpolator: airfoils morph with
``Airfoil.blend`` / ``cylinder_to_airfoil`` (circle roots) or
``interpolate_by_thickness`` (foil to foil). Chord, twist, pitch axis and offset
come from the stations or optional span tables. ``section_at`` and ``oml_polyline``
both call it. Scalar tables may be linear or PCHIP; airfoils are not 3-D-morphed.

``section_at(z)`` returns the wall midlines with their thickness, region labels and the
closed cell-wall loops, placed in the blade frame. ``oml_polyline(z)`` is the closed
outer mould line in the same frame.

Blade frame (the ``BladeBuilder`` axes): z along the span, x along the chord from the
pitch axis towards the trailing edge, y through the thickness (suction side positive).
Positive twist is nose up: the leading edge moves to +y. Required
``pivot_fraction`` is the S-frame origin — the chord fraction the spanwise axis
pierces. ``offset`` then shifts that point (prebend, sweep). The section geometry is built at the station's twist, so
a layout's flapwise webs (``"-F"``) stay normal to the rotor plane (along y) while its
chord-normal webs (``"-CN"``) turn with the chord; the two are the same untwisted.

Wall thicknesses are geometry, not material: they set where each wall's midline lies
(skins grow inward from the outer mould line, caps sit inside the skin, webs end on the
skin midline). Give each as a number or as a function of z. Laminates are bound to the
walls by label downstream (``ThinWalledSection.from_geometry``).

Web positions: ``LC`` uses the layout's chord fractions at every station. ``LE``
holds metres from the leading edge (root chord as the metre stick). ``PA`` holds
metres from the pitch axis on a root lattice ``k·c_root/n`` (``k = 1…n``); the
last point sits on the root TE and may miss the inner skin until the chord
grows. Two-web PA keeps the legacy ``(0.15, 0.50)`` pair. Give ``web_positions``
(chord fractions) or ``web_offsets`` (metres from the leading edge) to override.
Webs that do not intercept a station are omitted (same index may seed later)
rather than slid, and ``SectionCut.dropped_webs`` lists them. Loft identity
follows the same metre laws (not catalogue LC fractions), so a PA web walks
with the pitch axis until that intercept fails.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import replace
from typing import Callable, Optional, Sequence, Union

import numpy as np

from implicit_geometry_kernel.blade.cut import CellPolygon, SectionCut, SectionWall, _twist_offset
from implicit_geometry_kernel.blade.stations import PlyDrop, Station, Stringer, TransverseMember
from implicit_geometry_kernel.section.layout.thickness_limits import (
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M as _MIN_SKIN,
    MIN_REALISTIC_SPAR_WALL_THICKNESS_M as _MIN_CAP,
    MIN_REALISTIC_WEB_WALL_THICKNESS_M as _MIN_WEB,
)


def _catalogue_ribs(stations: Sequence[Station], n_ribs: int) -> tuple[TransverseMember, ...]:
    n = int(n_ribs)
    if n <= 0:
        return ()
    z0, z1 = float(stations[0].z), float(stations[-1].z)
    if abs(z1 - z0) < 1e-12:
        return (TransverseMember("tip_closure", z1, 0.012),)
    zs = np.linspace(z0, z1, n + 1)[1:]
    out = []
    for i, z in enumerate(zs):
        name = "tip_closure" if i == n - 1 else f"rib_{i + 1}"
        out.append(TransverseMember(name, float(z), 0.012))
    return tuple(out)

Thickness = Union[float, Callable[[float], float]]
# A span table: f(z), or (z, values) with values shape (n,) or (n, 2) for offset.
SpanTable = Union[Callable, Sequence]

# section_geometry turns twist the other way (its positive twist is nose down). The
# section is built twisted, so that flapwise ("-F") webs stay normal to the rotor
# plane while chord-normal ("-CN") webs turn with the chord.
_FORK_TWIST_SIGN = -1.0


def _as_function(value: Thickness, name: str) -> Callable[[float], float]:
    if callable(value):
        return value
    v = float(value)
    if v < 0.0:
        raise ValueError(f"{name} must be >= 0")
    return lambda z: v


def _as_pivot_fraction(value: Union[float, SpanTable]):
    """Required piercing-point fraction: a constant in (0, 1) or a span table."""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        f = float(value)
        if not 0.0 < f < 1.0:
            raise ValueError("pivot_fraction must be a chord fraction in (0, 1)")
        return lambda z, _f=f: _f
    table = _as_span_table(value, "pivot_fraction")
    if table is None:
        raise TypeError("pivot_fraction is required (a chord fraction or span table)")
    return table


def _as_span_table(value: Optional[SpanTable], name: str, *, cols: int = 1):
    """Linear lookup ``f(z)``. ``None`` means use the bounding stations."""
    if value is None:
        return None
    if callable(value):
        return value
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError(f"{name} must be a callable or a (z, values) pair")
    z = np.asarray(value[0], float).reshape(-1)
    v = np.asarray(value[1], float)
    if z.size < 2 or np.any(np.diff(z) <= 0.0):
        raise ValueError(f"{name} needs strictly increasing z of length >= 2")
    if v.shape[0] != z.size:
        raise ValueError(f"{name} values must match z")
    if cols == 1:
        v = np.asarray(v, float).reshape(z.size)

        def lookup(zz: float) -> float:
            return float(np.interp(zz, z, v))

        return lookup
    if v.ndim != 2 or v.shape[1] != cols:
        raise ValueError(f"{name} values must have shape (n, {cols})")

    def lookup_vec(zz: float):
        return tuple(float(np.interp(zz, z, v[:, i])) for i in range(cols))

    return lookup_vec


def _smoothstep(w: float) -> float:
    w = float(np.clip(w, 0.0, 1.0))
    return w * w * (3.0 - 2.0 * w)


def _near_circle(foil) -> bool:
    name = str(getattr(foil, "name", "")).lower()
    return float(foil.max_thickness) >= 0.6 or "ellipse" in name


def _morph_surfaces(fa, fb, w: float, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """TE-preserving morph: smoothstep blend on circle roots, else thickness blend."""
    from implicit_geometry_kernel.blade.airfoil import Airfoil

    w = float(np.clip(w, 0.0, 1.0))
    n = max(int(len(x)), 5)
    if w <= 1e-15:
        return np.asarray(fa.surface("upper", x), float), np.asarray(fa.surface("lower", x), float)
    if w >= 1.0 - 1e-15:
        return np.asarray(fb.surface("upper", x), float), np.asarray(fb.surface("lower", x), float)
    if _near_circle(fa) or _near_circle(fb):
        s = _smoothstep(w) if _near_circle(fa) != _near_circle(fb) else w
        foil = Airfoil.blend(fa, fb, s, n=n)
    else:
        t_c = (1.0 - w) * float(fa.max_thickness) + w * float(fb.max_thickness)
        foil = Airfoil.interpolate_by_thickness(fa, fb, t_c, n=n)
    upper = np.asarray(foil.surface("upper", x), float)
    lower = np.asarray(foil.surface("lower", x), float)
    gap_a = abs(float(np.asarray(fa.surface("upper", 1.0)) - np.asarray(fa.surface("lower", 1.0))))
    gap_b = abs(float(np.asarray(fb.surface("upper", 1.0)) - np.asarray(fb.surface("lower", 1.0))))
    gap = (1.0 - w) * gap_a + w * gap_b
    mid = 0.5 * (float(upper[-1]) + float(lower[-1]))
    upper[-1] = mid + 0.5 * gap
    lower[-1] = mid - 0.5 * gap
    return upper, lower


class BladeGeometry:
    """A blade or wing lofted from stations, cut into thin-walled sections on demand.

    Parameters
    ----------
    stations:
        Two or more ``Station`` objects in increasing ``z``.
    layout:
        A catalogue key (``"2D-F-LC-0R"``, …), a ``SystemLayoutSpec``, or a parametric
        ``SectionLayout`` (web fractions and cap rule, as for reference blades).
    skin_thickness, spar_cap_thickness, shear_web_thickness:
        Wall thicknesses in metres, as numbers or functions of z.
    web_positions:
        Web chord fractions, LE to TE (default: the layout's). A callable ``f(z)``
        or ``(z, values)`` table with shape ``(n, n_webs)`` varies stations along
        the span.
    cap_windows:
        Optional chord windows ``((lo, hi), …)`` or a span table of the same
        (fractions from the LE). Applied at each ``z`` as that interval, not a
        pitch-centred pad, even when webs are present.
    web_offsets:
        Web positions in metres from the leading edge, fixed along the span. Exclusive
        with ``web_positions``.
    n_airfoil_points:
        Points per surface used to hand the interpolated airfoil to the section
        geometry.
        chord, twist, offset:
        Optional span tables that override the station scalars at ``z`` (callable or
        ``(z, values)``; ``offset`` values are shape ``(n, 2)``). Airfoil graphs morph
        between the bounding ``Station`` foils (circle roots use a smoothstep blend).
    pivot_fraction:
        Required pitch-axis chord fraction from the LE (S-frame origin and XY
        stacking point). A float is constant along the span; a callable or
        ``(z, values)`` table varies it. ``pitch_axis`` is an alias for the same
        table. There is no 0.25 / 0.333 fallback.
    span_kind:
        ``"linear"`` (default) or ``"pchip"`` for station chord, twist, pitch axis
        and offset. Airfoil graphs stay linear blends.
    transverse:
        Ribs, bulkheads and frames at fixed ``z``. A member appears as a closed
        rib face on ``section_at(z).shell_graph().faces`` at that station only.
    stringers:
        Partial-span or free-end branches. Present on ``shell_graph().branches``
        inside ``extent``; ends are topology events.
        ply_drops:
        Geometry-only extra thickness windows in ``(z, s)`` with linear tapers.
    merge_cap_and_skin:
        If True, ``section_at`` walls follow ``shell_graph`` runs (caps as
        partitions, combined thickness). Default False keeps stacked cap solids
        (golden numbers).
    """

    def __init__(
        self,
        stations: Sequence[Station],
        layout,
        skin_thickness: Thickness,
        spar_cap_thickness: Thickness = 0.0,
        shear_web_thickness: Thickness = 0.0,
        web_positions: Optional[SpanTable] = None,
        web_offsets: Optional[Sequence[float]] = None,
        n_airfoil_points: int = 161,
        chord: Optional[SpanTable] = None,
        twist: Optional[SpanTable] = None,
        pitch_axis: Optional[SpanTable] = None,
        pivot_fraction: Optional[Union[float, SpanTable]] = None,
        offset: Optional[SpanTable] = None,
        cap_windows: Optional[SpanTable] = None,
        transverse: Sequence[TransverseMember] = (),
        stringers: Sequence[Stringer] = (),
        span_kind: str = "linear",
        ply_drops: Sequence[PlyDrop] = (),
        merge_cap_and_skin: bool = False,
        region_graph=None,
    ):
        from implicit_geometry_kernel.section.layout.system_layout import resolve_system_type

        self.stations = tuple(sorted(stations, key=lambda s: s.z))
        if len(self.stations) < 2:
            raise ValueError("a blade needs at least two stations")
        if np.any(np.diff([s.z for s in self.stations]) <= 0.0):
            raise ValueError("stations need distinct z")
        if isinstance(layout, str):
            spec = resolve_system_type(layout)
        else:
            spec = layout
        if web_positions is not None and web_offsets is not None:
            raise ValueError("give web_positions or web_offsets, not both")
        self._web_positions_fn = None
        if web_positions is not None and callable(web_positions):
            self._web_positions_fn = web_positions
        elif web_positions is not None and isinstance(web_positions, (tuple, list)) and len(web_positions) == 2 and not _is_frac_sequence(web_positions):
            n_cols = int(np.asarray(web_positions[1], float).reshape(len(web_positions[0]), -1).shape[1])
            self._web_positions_fn = _as_span_table(web_positions, "web_positions", cols=n_cols)
        elif web_positions is not None:
            fr = tuple(float(f) for f in web_positions)
            if len(fr) != spec.n_webs:
                raise ValueError(f"layout has {spec.n_webs} webs; got {len(fr)} positions")
            if any(not 0.0 < f < 1.0 for f in fr) or np.any(np.diff(fr) <= 0.0):
                raise ValueError("web_positions must be increasing chord fractions inside (0, 1)")
            spec = _with_web_fractions(spec, fr)
        self.spec = spec
        self.layout = getattr(spec, "key", None) or getattr(spec, "catalogue_key", None)
        self.web_offsets = None if web_offsets is None else tuple(float(v) for v in web_offsets)
        self._web_x_from_pitch = None
        n_webs = int(getattr(spec, "n_webs", 0) or 0)
        anchoring = getattr(spec, "web_anchoring", "local_chord")
        if pitch_axis is not None and pivot_fraction is not None:
            raise ValueError("give pitch_axis or pivot_fraction, not both")
        src = pitch_axis if pivot_fraction is None else pivot_fraction
        if src is None:
            raise TypeError("pivot_fraction is required (a chord fraction or span table)")
        self._pitch = _as_pivot_fraction(src)
        if self.web_offsets is None and web_positions is None and n_webs > 0:
            c0 = float(self.stations[0].chord)
            p0 = float(self._pitch(self.stations[0].z))
            fracs = spec.web_chord_fracs
            if anchoring == "leading_edge":
                self.web_offsets = tuple(float(f) * c0 for f in fracs)
            elif anchoring == "pitch_axis":
                if n_webs <= 2:
                    self._web_x_from_pitch = tuple((float(f) - p0) * c0 for f in fracs)
                else:
                    self._web_x_from_pitch = tuple(
                        (k / n_webs - p0) * c0 for k in range(1, n_webs + 1)
                    )
        if self.web_offsets is not None and len(self.web_offsets) != spec.n_webs:
            raise ValueError(f"layout has {spec.n_webs} webs; got {len(self.web_offsets)} offsets")
        if self._web_x_from_pitch is not None and len(self._web_x_from_pitch) != spec.n_webs:
            raise ValueError(f"layout has {spec.n_webs} webs; got {len(self._web_x_from_pitch)} pitch offsets")
        self._cap_windows_fn = None
        if cap_windows is not None:
            if callable(cap_windows):
                self._cap_windows_fn = cap_windows
            elif isinstance(cap_windows, (tuple, list)) and cap_windows and isinstance(cap_windows[0], (tuple, list)):
                frozen = tuple((float(a), float(b)) for a, b in cap_windows)
                self._cap_windows_fn = lambda z, _w=frozen: _w
            else:
                self._cap_windows_fn = _as_span_table(cap_windows, "cap_windows", cols=2)
        self._skin = _as_function(skin_thickness, "skin_thickness")
        self._cap = _as_function(spar_cap_thickness, "spar_cap_thickness")
        self._web = _as_function(shear_web_thickness, "shear_web_thickness")
        self._chord = _as_span_table(chord, "chord")
        self._twist = _as_span_table(twist, "twist")
        self._offset = _as_span_table(offset, "offset", cols=2)
        self.n_airfoil_points = int(n_airfoil_points)
        self.transverse = tuple(transverse)
        if int(getattr(spec, "n_ribs", 0) or 0) > 0 and not self.transverse:
            self.transverse = _catalogue_ribs(self.stations, spec.n_ribs)
        self.stringers = tuple(stringers)
        if span_kind not in ("linear", "pchip"):
            raise ValueError("span_kind must be 'linear' or 'pchip'")
        self.span_kind = span_kind
        self.ply_drops = tuple(ply_drops)
        self.merge_cap_and_skin = bool(merge_cap_and_skin)
        if region_graph is None:
            from implicit_geometry_kernel.blade.topology import WallRegionGraph
            self.region_graph = self._with_anchored_web_laws(
                WallRegionGraph.from_section_layout(self.spec)
            )
        else:
            self.region_graph = region_graph
        # Live fields (wall_meta / stack T) key stations by round(z, 9). Cache the
        # cut itself so those fields do not rebuild geometry on a second hit.
        self._section_cache: OrderedDict[float, SectionCut] = OrderedDict()
        self._section_cache_size = 256

    @classmethod
    def from_concept(cls, concept, stations: Sequence[Station], skin_thickness: Thickness, **kwargs):
        """Build from a :class:`~implicit_geometry_kernel.blade.concept.BladeConcept`.

        Extra kwargs overlay the compiled scheme (thicknesses, web positions, …).
        ``pivot_fraction`` is required (or the ``pitch_axis`` alias).
        """
        compiled = dict(concept.to_geometry_kwargs())
        compiled.update(kwargs)
        layout = compiled.pop("layout")
        return cls(stations, layout, skin_thickness, **compiled)

    @property
    def concept(self):
        """Layout plus the spanwise extras already stored on this geometry."""
        from implicit_geometry_kernel.blade.concept import BladeConcept, SpanwiseScheme

        layout = self.layout if self.layout is not None else self.spec
        return BladeConcept(
            layout,
            SpanwiseScheme(
                transverse=self.transverse,
                stringers=self.stringers,
                ply_drops=self.ply_drops,
            ),
        )

    def forced_stations(self) -> tuple[float, ...]:
        """Sorted unique member stations (ribs and stringer ends; mesh hint)."""
        by_key: dict[float, float] = {}
        for member in self.transverse:
            by_key[round(float(member.z), 9)] = float(member.z)
        for stringer in self.stringers:
            if stringer.extent is None:
                continue
            for z in stringer.extent:
                by_key[round(float(z), 9)] = float(z)
        for drop in self.ply_drops:
            for z in (drop.z0, drop.z1):
                by_key[round(float(z), 9)] = float(z)
        return tuple(by_key[key] for key in sorted(by_key))

    def _web_chord_frac_law(self, index: int):
        """Chord-fraction law for loft identity: LC catalogue, or PA/LE metres."""
        if self._web_x_from_pitch is not None:
            x_pa = float(self._web_x_from_pitch[int(index)])

            def frac(z, _x=x_pa):
                zz = float(z)
                chord = max(self._span_scalar(self._chord, [s.chord for s in self.stations], zz), 1e-12)
                return float(self._pitch(zz)) + _x / chord

            return frac
        if self.web_offsets is not None:
            x_le = float(self.web_offsets[int(index)])

            def frac(z, _x=x_le):
                zz = float(z)
                chord = max(self._span_scalar(self._chord, [s.chord for s in self.stations], zz), 1e-12)
                return _x / chord

            return frac
        fracs = tuple(float(f) for f in getattr(self.spec, "web_chord_fracs", ()) or ())
        f = float(fracs[int(index)])
        return lambda z, _f=f: _f

    def _with_anchored_web_laws(self, graph):
        """PA/LE loft anchors follow the metre freeze, not catalogue LC fractions."""
        if self._web_x_from_pitch is None and self.web_offsets is None:
            return graph
        webs = tuple(replace(w, frac=self._web_chord_frac_law(w.index)) for w in graph.webs)
        anchors = []
        for a in graph.anchors:
            if a.name.startswith("web_"):
                try:
                    i = int(a.name.split("_", 1)[1])
                except ValueError:
                    anchors.append(a)
                    continue
                anchors.append(replace(a, frac=self._web_chord_frac_law(i)))
            else:
                anchors.append(a)
        return replace(graph, webs=webs, anchors=tuple(anchors))

    def _transverse_at(self, z: float) -> tuple[TransverseMember, ...]:
        key = round(float(z), 9)
        return tuple(m for m in self.transverse if round(float(m.z), 9) == key)

    def _stringers_at(self, z: float) -> tuple[Stringer, ...]:
        key = round(float(z), 9)
        out = []
        for stringer in self.stringers:
            if stringer.extent is None:
                out.append(stringer)
                continue
            z0, z1 = round(float(stringer.extent[0]), 9), round(float(stringer.extent[1]), 9)
            if z0 <= key <= z1:
                out.append(stringer)
        return tuple(out)

    def _layout_at(self, z: float):
        from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
        from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
        from implicit_geometry_kernel.section.layout.system_layout import SystemLayoutSpec

        spec = self.spec
        if self._web_positions_fn is None and self._cap_windows_fn is None:
            return spec
        if isinstance(spec, SystemLayoutSpec):
            layout = SectionLayout.from_system_layout_spec(spec)
        else:
            layout = spec
        fracs = None
        if self._web_positions_fn is not None:
            raw = self._web_positions_fn(float(z))
            fracs = tuple(float(f) for f in np.atleast_1d(raw).reshape(-1))
        if fracs is not None:
            layout = _with_web_fractions(layout, fracs)
        if self._cap_windows_fn is not None:
            windows = self._cap_windows_fn(float(z))
            if windows is not None:
                raw = tuple(windows) if not isinstance(windows, tuple) else windows
                if raw and not isinstance(raw[0], (tuple, list)):
                    raw = (tuple(float(x) for x in raw),)
                layout = SectionLayout(
                    web_chord_fracs=layout.web_chord_fracs,
                    orientation=layout.orientation,
                    cap_placement=CapPlacementRule.from_chord_windows(raw),
                    catalogue_key=layout.catalogue_key,
                    web_stations_m=layout.web_stations_m,
                    web_axes_m=layout.web_axes_m,
                    structural_family=layout.structural_family,
                    flags=layout.flags,
                )
        return layout

    # ------------------------------------------------------------ stations
    def interpolate(self, z: float):
        """(chord, twist, pitch_axis, offset, (x, upper, lower)) interpolated at z.

        Circle-like parents morph with a smoothstep blend; otherwise surfaces follow
        ``Airfoil.interpolate_by_thickness``. Chord, twist and offset come from span
        tables when given, otherwise from the stations. Pitch axis is always
        ``pivot_fraction``.
        """
        zs = np.array([s.z for s in self.stations])
        if z < zs[0] - 1e-12 * max(1.0, abs(zs[0])) or z > zs[-1] + 1e-12 * max(1.0, abs(zs[-1])):
            raise ValueError(f"z = {z} lies outside the stations ({zs[0]} to {zs[-1]})")
        k = int(np.clip(np.searchsorted(zs, z, side="right") - 1, 0, len(zs) - 2))
        a, b = self.stations[k], self.stations[k + 1]
        w = float(np.clip((z - a.z) / (b.z - a.z), 0.0, 1.0))
        x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, self.n_airfoil_points)))
        upper, lower = _morph_surfaces(a.airfoil, b.airfoil, w, x)
        zz = float(z)
        chord = self._span_scalar(self._chord, [s.chord for s in self.stations], zz)
        twist = self._span_scalar(self._twist, [s.twist for s in self.stations], zz)
        pitch = float(self._pitch(zz))
        if self._offset is not None:
            off = self._offset(zz)
            offset = (float(off[0]), float(off[1]))
        else:
            ox = self._span_scalar(None, [s.offset[0] for s in self.stations], zz)
            oy = self._span_scalar(None, [s.offset[1] for s in self.stations], zz)
            offset = (ox, oy)
        return (chord, twist, pitch, offset, (x, upper, lower))

    def _span_scalar(self, table, values, z: float) -> float:
        if table is not None:
            return float(table(z))
        if self.span_kind == "linear":
            zs = np.array([s.z for s in self.stations], float)
            return float(np.interp(z, zs, np.asarray(values, float)))
        from implicit_geometry_kernel.structure.laws import TableLaw

        knots = [float(s.z) for s in self.stations]
        return float(TableLaw(knots, [float(v) for v in values], kind="pchip")(z))

    def oml_polyline(self, z: float) -> np.ndarray:
        """Closed outer mould line in the blade frame at ``z``, shape ``(N, 3)``.

        Same interpolator and placement as ``section_at`` (pitch axis, twist, offset).
        """
        z = float(z)
        chord, twist, pitch, offset, (x, upper, lower) = self.interpolate(z)
        naca = np.vstack([np.column_stack([x[::-1], upper[::-1]]), np.column_stack([x[1:], lower[1:]])]) * chord
        if np.allclose(naca[0], naca[-1]):
            naca = naca[:-1]
        local = np.column_stack([naca[:, 0] - pitch * chord, naca[:, 1]])
        xy = _twist_offset(local, twist, offset)
        pts = np.column_stack([xy, np.full(len(xy), z)])
        return np.vstack([pts, pts[0]])

    # ------------------------------------------------------------ cutting
    def section_at(self, z: float) -> SectionCut:
        """The wall midlines, thicknesses, labels and cell loops of the section at z.

        Cuts are cached per ``round(z, 9)`` (LRU 256). ``CutWallMetaSDF`` and
        ``StackThickness`` use the same key; this is the G7 cost fix, not a BVH.
        """
        from implicit_geometry_kernel.frames import twisted_pitch_axis_blade_frame_to_naca
        from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF
        from implicit_geometry_kernel.section import SectionGeometry

        z = round(float(z), 9)
        hit = self._section_cache.get(z)
        if hit is not None:
            self._section_cache.move_to_end(z)
            return hit
        chord, twist, pitch, offset, (x, upper, lower) = self.interpolate(z)
        layout = self._layout_at(z)
        t_skin, t_cap, t_web = self._skin(z), (self._cap(z) if layout.has_spar_caps else 0.0), self._web(z)
        # section_geometry floors wall thicknesses at wind-turbine manufacturing minima
        # (2.5 mm skin, 4 mm web, 5 mm cap). The geometry is similar under scaling, so a
        # thinner-walled section is built k times larger and scaled back: exact, and the
        # frozen geometry code is used as it is.
        k = max(1.0, *(lo / t for lo, t in ((_MIN_SKIN, t_skin), (_MIN_WEB, t_web), (_MIN_CAP, t_cap)) if t > 0.0))
        if k > 1.0:
            k *= 1.0 + 1e-9
        kwargs = dict(
            layout=layout,
            airfoil=AirfoilSDF.from_surface_graphs(x, upper, lower, chord * k),
            skin_thickness=t_skin * k,
            spar_cap_thickness=t_cap * k,
            shear_web_thickness=t_web * k,
            twist_rad=_FORK_TWIST_SIGN * np.radians(twist),
            pivot_fraction=pitch,
            outputs=("midlines", "cell_wall_loops"),
        )
        if self.web_offsets is not None:
            kwargs["web_offsets"] = tuple(v * k for v in self.web_offsets)
        elif self._web_x_from_pitch is not None:
            kwargs["web_offsets"] = tuple((x + pitch * chord) * k for x in self._web_x_from_pitch)
        runner = SectionGeometry(**kwargs).execute()
        result = runner.get_results()
        resolution = runner.get_layout_resolution()

        def place(pts_B):
            pts = np.asarray(pts_B, float) / k
            naca = twisted_pitch_axis_blade_frame_to_naca(
                pts, _FORK_TWIST_SIGN * np.radians(twist), chord, pitch,
            )
            local = np.column_stack([naca[:, 0] - pitch * chord, naca[:, 1]])
            xy = _twist_offset(local, twist, offset)
            return np.column_stack([xy, np.full(len(xy), z)])

        from implicit_geometry_kernel.section.junctions.cap_skin_bonds import (
            bond_caps_to_skin,
        )
        from implicit_geometry_kernel.section.midlines.export import build_solid_cap_strips
        from implicit_geometry_kernel.section.midlines.skin_split import (
            split_skin_at_cap_footprints,
        )

        mesh_strips = tuple(result.midline_strips or ())
        cap_strips = build_solid_cap_strips(
            result.section,
            twist_rad=float(kwargs["twist_rad"]),
            chord_m=float(chord * k),
            pivot_fraction=pitch,
        )
        # Solids path: cap recipe midlines + skin/web mesh strips, bonded for
        # graph partitions. Split host skins at footprint ends so cap colouring
        # cannot emit two loop runs with one wall name.
        if cap_strips:
            strips = list(
                split_skin_at_cap_footprints(
                    bond_caps_to_skin(mesh_strips + cap_strips)
                )
            )
        else:
            strips = list(mesh_strips)
        by_label = {s.label: s for s in strips}
        walls = []
        for s in strips:
            pts = place(s.midline_twisted_pitch_axis_blade_frame)
            side = s.skin_side or s.surface
            side = None if side is None else side.replace("_surface", "")
            bonds = {}
            for bond in s.bonds:
                host = by_label.get(bond.host_label)
                if host is None:
                    continue
                bonds[bond.host_label] = (
                    _match(s.midline_twisted_pitch_axis_blade_frame, bond.points),
                    _match(host.midline_twisted_pitch_axis_blade_frame, bond.host_points),
                )
            walls.append(SectionWall(s.label, s.kind, pts, float(s.thickness_m) / k, side, bonds))
        cells = tuple(
            CellPolygon(int(c.id), place(c.polyline_twisted_pitch_axis_blade_frame), float(c.area_m2) / k**2,
                     tuple(c.wall_labels))
            for c in (result.cell_wall_loops or ())
        )
        dropped = tuple(f"web at {d.x_m / k:.6g} m from the LE: {d.reason}" for d in resolution.dropped)
        dropped_idx = tuple(int(d.index) for d in resolution.dropped if getattr(d, "index", None) is not None)
        key = getattr(resolution.layout, "catalogue_key", None) or self.layout
        cut = SectionCut(
            z, chord, twist, pitch, (float(offset[0]), float(offset[1])), str(key), tuple(walls), cells,
            dropped, self._transverse_at(z), self._stringers_at(z), result.section, k,
            dropped_web_indices=dropped_idx,
        )
        if self.merge_cap_and_skin:
            from implicit_geometry_kernel.blade.shell_graph import shell_strips_from_graph
            from dataclasses import replace as _replace

            merged = []
            for strip in shell_strips_from_graph(cut):
                xy = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, float)
                if xy.shape[1] == 2:
                    pts = np.column_stack([xy, np.full(len(xy), z)])
                else:
                    pts = xy
                stack_t = (
                    float(sum(strip.layer_thicknesses_m))
                    if strip.layer_thicknesses_m
                    else float(strip.thickness_m)
                )
                merged.append(SectionWall(strip.label, strip.kind, pts, stack_t, None, {}))
            cut = _replace(cut, walls=tuple(merged))
        self._section_cache[z] = cut
        while len(self._section_cache) > self._section_cache_size:
            self._section_cache.popitem(last=False)
        return cut

    def oml_sdf(self, meta: int = 0):
        """Exact 2-D polygon SDF of ``oml_polyline(z)`` at each point's own ``z``."""
        from implicit_geometry_kernel.sdf.nodes.lofted import LoftedPolygonSDF

        return LoftedPolygonSDF(lambda z: self.oml_polyline(float(z))[:, :2], meta=meta)

    def cavity_sdf(self, meta: int = 0):
        """Inner mould line: ``oml + T``. ``T`` is a step (skin, plus cap inside footprints); no taper."""
        from implicit_geometry_kernel.blade.fields import StackThickness
        from implicit_geometry_kernel.sdf.nodes.tagged import FieldOffsetNode

        return FieldOffsetNode(self.oml_sdf(meta), StackThickness(self), scale=1.0)

    def wall_mid_sdf(self, meta: int = 0):
        """Signed polygon SDF of the ``shell_graph`` skin loop at each point's ``z``.

        The graph loop is the spec (run interiors and partition vertices on
        ``phi = 0``). Webs are branches, not this contour. Colour is ``wall_meta``.
        """
        from implicit_geometry_kernel.contract import InfeasibleSectionGeometry
        from implicit_geometry_kernel.sdf.nodes.lofted import LoftedPolygonSDF

        def provider(z: float):
            try:
                return self.section_at(float(z)).shell_graph().loop_polygon()
            except InfeasibleSectionGeometry:
                return self.oml_polyline(float(z))[:, :2]

        return LoftedPolygonSDF(provider, meta=meta)

    def wall_meta(self, table=None):
        """Tagged wall-strip DAG: ``evaluate((x,y,z)).meta`` is the region at that ``z``."""
        from implicit_geometry_kernel.blade.fields import CutWallMetaSDF
        from implicit_geometry_kernel.sdf.nodes.base import FeatureTable

        return CutWallMetaSDF(self, table if table is not None else FeatureTable())

    def thickness_law(self):
        """:class:`~implicit_geometry_kernel.structure.laws.ThicknessLaw` over ``(z, s)``."""
        from implicit_geometry_kernel.structure.laws import ThicknessLaw

        geo = self

        class _BladeThickness(ThicknessLaw):
            def __call__(self, z, s) -> np.ndarray:
                return geo.thickness_at(z, s)

        return _BladeThickness()

    def named_edges(self, z: float):
        """Partition edges and ply-drop ``s`` bounds at ``z`` (mesher feature curves)."""
        from implicit_geometry_kernel.blade.shell_graph import PartitionEdge

        z = float(z)
        edges = list(self.section_at(z).shell_graph().partition_edges())
        chord, _, pitch, offset, (x, upper, lower) = self.interpolate(z)
        for drop in self.ply_drops:
            if not (float(drop.z0) - 1e-12 <= z <= float(drop.z1) + 1e-12):
                continue
            for s_frac, tag in ((drop.s0, "s0"), (drop.s1, "s1")):
                s_frac = float(np.clip(s_frac, 0.0, 1.0))
                i = int(np.clip(round(s_frac * (len(x) - 1)), 0, len(x) - 1))
                local = np.array([[x[i] - pitch * chord, 0.5 * (upper[i] + lower[i])]], float)
                xy = _twist_offset(local, self.interpolate(z)[1], offset)[0]
                pt = np.array([xy[0], xy[1], z], float)
                name = f"{drop.name}:{tag}"
                edges.append(PartitionEdge(z, name, name, pt))
        return tuple(edges)

    def size_fields(self, points) -> dict[str, np.ndarray]:
        """Distance, thickness and a polyline curvature estimate at 3-D points."""
        pts = np.asarray(points, float).reshape(-1, 3)
        out = self.query(pts)
        curv = np.zeros(len(pts))
        for i, p in enumerate(pts):
            poly = self.oml_polyline(float(p[2]))[:, :2]
            j = int(np.argmin(np.linalg.norm(poly - p[:2], axis=1)))
            a = poly[(j - 1) % len(poly)]
            b = poly[j]
            c = poly[(j + 1) % len(poly)]
            ab, bc = b - a, c - b
            lab = float(np.linalg.norm(ab))
            lbc = float(np.linalg.norm(bc))
            if lab < 1e-15 or lbc < 1e-15:
                continue
            cross = ab[0] * bc[1] - ab[1] * bc[0]
            curv[i] = abs(cross) / max(lab * lbc * 0.5 * (lab + lbc), 1e-15)
        out["curvature"] = curv
        return out

    def thickness_at(self, z, s) -> np.ndarray:
        """Skin thickness at span ``z`` and chord fraction ``s``, plus ply-drop ramps."""
        z = np.asarray(z, float)
        s = np.asarray(s, float)
        z, s = np.broadcast_arrays(z, s)
        out = np.empty(z.shape, dtype=float)
        for idx in np.ndindex(z.shape):
            t = float(self._skin(float(z[idx])))
            for drop in self.ply_drops:
                t += drop.extra(float(z[idx]), float(s[idx]))
            out[idx] = t
        return out

    def shell_graph(self, z: float, include_solids: bool = False):
        """Loftable shell graph at ``z``. Solids (cores/inserts) are opt-in."""
        return self.section_at(z).shell_graph(include_solids=include_solids)

    def key_points(self, z: float) -> dict[str, np.ndarray]:
        """Named points at ``z``: LE, TE, cap edges, web feet."""
        return self.section_at(z).key_points()

    def query(self, points) -> dict[str, np.ndarray]:
        """Batched thickness, owning member, local frame and OML distance at 3-D points."""
        pts = np.asarray(points, float).reshape(-1, 3)
        oml = self.oml_sdf().evaluate(pts)
        field = self.wall_meta()
        sample = field.evaluate(pts)
        table = field.table
        names, kinds = [], []
        for m in sample.meta.astype(int):
            if m <= 0:
                names.append("")
                kinds.append("")
            else:
                spec = table[int(m)]
                names.append(spec.name)
                kinds.append(spec.kind)
        s = []
        for p in pts:
            chord, _, pitch, _, _ = self.interpolate(float(p[2]))
            s.append(float(np.clip((p[0] + pitch * chord) / max(chord, 1e-15), 0.0, 1.0)))
        T = np.array([float(np.ravel(self.thickness_at(p[2], ss))[0]) for p, ss in zip(pts, s)])
        frames = []
        for p in pts:
            _, twist, _, _, _ = self.interpolate(float(p[2]))
            c, si = np.cos(np.radians(twist)), np.sin(np.radians(twist))
            frames.append(np.array([[c, -si, 0.0], [si, c, 0.0], [0.0, 0.0, 1.0]]))
        return {
            "thickness": T,
            "member": np.array(names, dtype=object),
            "kind": np.array(kinds, dtype=object),
            "distance": oml.phi,
            "frame": np.stack(frames, axis=0),
        }

    def write_oml(self, path, n_span: int = 41, n_chord: int = 81) -> str:
        """Watertight OML as STL or OBJ from a quad grid of span rings."""
        from pathlib import Path as _Path

        path = _Path(path)
        zs = np.linspace(self.stations[0].z, self.stations[-1].z, int(n_span))
        rings = []
        for z in zs:
            poly = self.oml_polyline(z)
            if np.allclose(poly[0], poly[-1]):
                poly = poly[:-1]
            idx = np.linspace(0, len(poly) - 1, int(n_chord), endpoint=False)
            rings.append(np.vstack([poly[int(i) % len(poly)] for i in np.round(idx)]))
        verts = np.vstack(rings)
        ns, nc = len(rings), int(n_chord)
        faces = []
        for i in range(ns - 1):
            for j in range(nc):
                a = i * nc + j
                b = i * nc + (j + 1) % nc
                c = (i + 1) * nc + (j + 1) % nc
                d = (i + 1) * nc + j
                faces.append((a, b, c))
                faces.append((a, c, d))
        # close root and tip
        for base, sign in ((0, 1), ((ns - 1) * nc, -1)):
            hub = verts[base:base + nc].mean(axis=0)
            verts = np.vstack([verts, hub])
            h = len(verts) - 1
            for j in range(nc):
                if sign > 0:
                    faces.append((h, base + j, base + (j + 1) % nc))
                else:
                    faces.append((h, base + (j + 1) % nc, base + j))
        _write_tri_mesh(path, verts, faces, solid_name="oml")
        return str(path)

    def write_mid_surfaces(
        self,
        path,
        n_span: int = 41,
        n_chord: int = 81,
        members=None,
        include_solids: bool = False,
    ) -> str:
        """One lofted patch per authored region / web, as STL or OBJ (groups).

        Skins and webs always; cores/inserts only if ``include_solids`` (G9).
        ``members`` is an optional iterable of region names. STEP uses extra ``[cad]``.
        Stations join with ``map_vertex``; infeasible cuts are skipped. Rib faces
        are n-gons at member ``z`` from loft skin/web nodes (OBJ SoT). STL fans
        those n-gons from the first loft vertex.
        """
        from pathlib import Path as _Path

        from implicit_geometry_kernel.blade.topology import (
            _template_from_cut,
            loft_mapped_strips,
            rib_cell_polygons_from_samples,
            sample_region_polylines,
            section_cut_near,
        )
        from implicit_geometry_kernel.contract import InfeasibleSectionGeometry

        path = _Path(path)
        want = None if members is None else set(members)
        z0, z1 = float(self.stations[0].z), float(self.stations[-1].z)
        zs = list(np.linspace(z0, z1, int(n_span)))
        for z in self.forced_stations():
            if z0 - 1e-12 <= float(z) <= z1 + 1e-12:
                zs.append(float(z))
        zs = sorted(set(round(float(z), 9) for z in zs))
        samples: list[tuple[float, object, dict]] = []
        for z in zs:
            cut = section_cut_near(self, float(z))
            if cut is None:
                continue
            sg = cut.shell_graph(include_solids=include_solids)
            tmpl = _template_from_cut(self.region_graph, cut, float(z), int(n_chord))
            samples.append((float(z), tmpl, sample_region_polylines(cut, sg, tmpl, self.region_graph)))
        groups: dict[str, list[tuple[np.ndarray, bool]]] = {}
        for (_za, t0, s0), (_zb, t1, s1) in zip(samples[:-1], samples[1:]):
            for name, face, collapsing in loft_mapped_strips(s0, s1, t0, t1):
                if want is not None and name not in want:
                    continue
                groups.setdefault(name, []).append((face, collapsing))
        if path.suffix.lower() in {".step", ".stp"}:
            from implicit_geometry_kernel.blade.cad import write_mid_surface_step
            from implicit_geometry_kernel.blade.locked_loft import locked_mid_surface_patches

            graphs = []
            for z in zs:
                try:
                    graphs.append(self.shell_graph(float(z), include_solids=include_solids))
                except InfeasibleSectionGeometry:
                    continue
            patches = locked_mid_surface_patches(graphs, n_chord=int(n_chord), members=members)
            return write_mid_surface_step(path, [(name, rings) for name, rings, _surface in patches])
        chunks: list[tuple[str, np.ndarray, list]] = []
        for name, faces in groups.items():
            verts: list[np.ndarray] = []
            idx: list[tuple[int, ...]] = []
            for face, _collapsing in faces:
                base = len(verts)
                for p in face:
                    verts.append(np.asarray(p, float).reshape(3))
                n = len(face)
                if n == 3:
                    idx.append((base, base + 1, base + 2))
                elif n >= 4:
                    idx.append(tuple(base + i for i in range(n)))
            if verts:
                chunks.append((name, np.vstack(verts), idx))
        for z, tmpl, samp in samples:
            for member in self.transverse:
                if abs(float(member.z) - float(z)) > 1e-8:
                    continue
                for name, ring in rib_cell_polygons_from_samples(tmpl, samp, member):
                    if want is not None and name not in want and "rib" not in want:
                        continue
                    n = len(ring)
                    if n < 3:
                        continue
                    chunks.append((name, np.asarray(ring, float), [tuple(range(n))]))
        if path.suffix.lower() == ".obj":
            lines: list[str] = []
            base = 0
            for name, verts, faces in chunks:
                lines.append(f"g {name}")
                lines += [f"v {x:g} {y:g} {z:g}" for x, y, z in verts]
                for face in faces:
                    lines.append("f " + " ".join(str(base + i + 1) for i in face))
                base += len(verts)
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        else:
            with path.open("w", encoding="utf-8") as fh:
                for name, verts, faces in chunks:
                    tris = []
                    for face in faces:
                        if len(face) == 3:
                            tris.append(face)
                        else:
                            tris.extend((face[0], face[i], face[i + 1]) for i in range(1, len(face) - 1))
                    _write_stl_solid(fh, name, verts, tris)
        return str(path)

    def validate(self) -> list[str]:
        """Named station/member checks: span order, thin walls, members outside the OML."""
        issues = []
        zs = [s.z for s in self.stations]
        if any(b <= a for a, b in zip(zs, zs[1:])):
            issues.append("stations must increase in z")
        for i, st in enumerate(self.stations):
            if st.chord <= 0:
                issues.append(f"station {i} at z={st.z:g}: non-positive chord")
        for member in self.transverse:
            if member.z < zs[0] or member.z > zs[-1]:
                issues.append(f"transverse {member.name!r} at z={member.z:g} leaves the blade")
            if member.thickness < 1e-6:
                issues.append(f"transverse {member.name!r}: too thin")
        for stringer in self.stringers:
            if stringer.extent is not None:
                if stringer.extent[0] < zs[0] - 1e-9 or stringer.extent[1] > zs[-1] + 1e-9:
                    issues.append(f"stringer {stringer.name!r} extent leaves the blade")
            if stringer.thickness < 1e-6:
                issues.append(f"stringer {stringer.name!r}: too thin")
        for z in (zs[0], 0.5 * (zs[0] + zs[-1]), zs[-1]):
            try:
                cut = self.section_at(z)
            except Exception as exc:
                issues.append(f"station z={z:g}: {exc}")
                continue
            oml = self.oml_polyline(z)
            if not np.allclose(oml[0], oml[-1]):
                issues.append(f"station z={z:g}: OML is not closed")
            for wall in cut.walls:
                if wall.thickness < 1e-6:
                    issues.append(f"station z={z:g} member {wall.label}: too thin")
        return issues

    def snap_oml(self, points, solver=None):
        """Project points onto the OML in each point's own ``z``-plane."""
        from implicit_geometry_kernel.sdf.solve import GeometrySolver

        solver = solver or GeometrySolver(tolerance=1e-8, max_step=0.05)
        return solver.project(self.oml_sdf(), np.asarray(points, float), plane=True)

    def snap_wall_mid(self, points, solver=None):
        """Project onto the ``shell_graph`` loop field in each point's ``z``-plane."""
        from implicit_geometry_kernel.sdf.solve import GeometrySolver

        solver = solver or GeometrySolver(tolerance=1e-8, max_step=0.05)
        return solver.project(self.wall_mid_sdf(), np.asarray(points, float), plane=True)

    def __repr__(self) -> str:
        return f"BladeGeometry({self.layout!r}, z = {self.stations[0].z:g} to {self.stations[-1].z:g})"


def _match(polyline: np.ndarray, points: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Index of each of ``points`` among the vertices of ``polyline``."""
    P = np.asarray(polyline, float)
    out = []
    for p in np.asarray(points, float).reshape(-1, 2):
        d = np.linalg.norm(P - p, axis=1)
        k = int(np.argmin(d))
        if d[k] > tol * max(1.0, float(np.abs(P).max())):
            raise ValueError("a bond point is not a vertex of its strip")
        out.append(k)
    return np.array(out, dtype=int)


def _is_frac_sequence(value) -> bool:
    try:
        vals = [float(v) for v in value]
    except (TypeError, ValueError):
        return False
    return bool(vals) and all(0.0 < v < 1.0 for v in vals)


def _with_web_fractions(spec, fractions):
    import dataclasses

    return dataclasses.replace(spec, web_chord_fracs=tuple(fractions))


def section_at(stations: Sequence[Station], z: float, layout, skin_thickness: Thickness, **options) -> SectionCut:
    """One-off cut: ``BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)``.

    ``pivot_fraction`` is required on ``options`` (or the ``pitch_axis`` alias).
    """
    return BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)


def _loft_rings(rings, cap: bool) -> tuple[np.ndarray, list]:
    verts = np.vstack(rings)
    ns, nc = len(rings), len(rings[0])
    faces = []
    for i in range(ns - 1):
        for j in range(nc):
            a = i * nc + j
            b = i * nc + (j + 1) % nc
            c = (i + 1) * nc + (j + 1) % nc
            d = (i + 1) * nc + j
            faces.append((a, b, c))
            faces.append((a, c, d))
    if cap:
        for base, sign in ((0, 1), ((ns - 1) * nc, -1)):
            hub = verts[base:base + nc].mean(axis=0)
            verts = np.vstack([verts, hub])
            h = len(verts) - 1
            for j in range(nc):
                if sign > 0:
                    faces.append((h, base + j, base + (j + 1) % nc))
                else:
                    faces.append((h, base + (j + 1) % nc, base + j))
    return verts, faces


def _write_stl_solid(fh, name: str, verts: np.ndarray, faces: list) -> None:
    fh.write(f"solid {name}\n")
    for i, j, k in faces:
        n = np.cross(verts[j] - verts[i], verts[k] - verts[i])
        ln = np.linalg.norm(n)
        if ln:
            n = n / ln
        fh.write(f" facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n  outer loop\n")
        for p in (verts[i], verts[j], verts[k]):
            fh.write(f"   vertex {p[0]:.6e} {p[1]:.6e} {p[2]:.6e}\n")
        fh.write("  endloop\n endfacet\n")
    fh.write(f"endsolid {name}\n")


def _write_tri_mesh(path, verts: np.ndarray, faces: list, solid_name: str = "oml") -> None:
    if path.suffix.lower() == ".obj":
        lines = [f"v {x:g} {y:g} {z:g}" for x, y, z in verts]
        lines += [f"f {i + 1} {j + 1} {k + 1}" for i, j, k in faces]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8") as fh:
        _write_stl_solid(fh, solid_name, verts, faces)


__all__ = ["BladeGeometry", "SpanTable", "Thickness", "section_at"]
