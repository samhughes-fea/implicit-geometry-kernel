"""Spanwise blade geometry: the implicit section geometry lofted between stations.

``BladeGeometry`` takes the ``Station`` list a ``BladeBuilder`` takes (chord, twist,
pitch axis, offset and airfoil at each span position) and a structural layout from
``section`` (one of its 40 system types, such as ``"2D-F-LC-0R"``: two webs and a
box spar cap between them). ``interpolate(z)`` is the only interpolator: airfoils morph with
``Airfoil.blend`` / ``cylinder_to_airfoil`` (circle roots) or
``interpolate_by_thickness`` (foil to foil). Chord, twist, pitch axis and offset
come from the stations or optional span tables. ``section_at`` and ``oml_polyline``
both call it. Scalars interpolate linearly; airfoils are not 3-D-morphed.

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
walls by label downstream (``mid_surface`` / ``shell_graph`` ``layers``; the
thin-wall consumer is not in this package).

Web positions are hard-coded in each mode's native unit (see
``blade.concepts.web_anchoring``). ``LC`` is authored chord fractions at every
station. ``LE`` is authored metres from the leading edge (catalogue fractions
times root chord). ``PA`` is authored blade-x from the pitch axis: up to three
primary slots are uniform interiors of the root ``[x_LE, x_TE]`` band; extra
slots continue that lattice (first extra at root TE). Give
``web_positions`` (chord fractions) or ``web_offsets`` (metres from the leading
edge) to override. Webs that do not intercept a station are omitted (same
index may seed later) rather than slid, and ``SectionCut.dropped_webs`` lists
them. ``SectionCut.layout`` stays the authored concept; ``station_webs`` /
``station_family`` are the realised W/X slice. Loft identity follows the same
metre laws (not catalogue LC fractions),
so a PA web walks with the pitch axis until that intercept fails.
"""

from __future__ import annotations

import warnings

from collections import OrderedDict
from typing import TYPE_CHECKING, Optional, Sequence, Union

import numpy as np

from implicit_geometry_kernel.blade.concepts.concept import catalogue_ribs
from implicit_geometry_kernel.blade.cut import SectionCut
from implicit_geometry_kernel.blade.stations import PlyDrop, Station, Stringer, TipClosure, TransverseMember
from implicit_geometry_kernel.blade.loft.loft_tables import (
    SpanTable,
    Thickness,
    _as_function,
    _as_pivot_fraction,
    _as_span_table,
    _is_frac_sequence,
    _metres_law,
    _with_web_fractions,
)
from implicit_geometry_kernel.blade.loft.loft_cut import LoftCutMixin
from implicit_geometry_kernel.blade.loft.loft_fields import LoftFieldsMixin
from implicit_geometry_kernel.blade.loft.loft_layout import LoftLayoutMixin

if TYPE_CHECKING:
    from implicit_geometry_kernel.catalogue.system_layout import SystemLayoutSpec
    from implicit_geometry_kernel.section.layout.section_layout import SectionLayout


class BladeGeometry(LoftCutMixin, LoftFieldsMixin, LoftLayoutMixin):
    """A blade or wing lofted from stations, cut into thin-walled sections on demand.

    This class stores interpolants and layout, not a snapshot of every wall along
    the span. Public attributes after construction: ``stations``, ``spec``,
    ``layout``, ``web_offsets``, ``transverse``, ``stringers``, ``ply_drops``,
    ``region_graph``, ``n_airfoil_points``, ``merge_cap_and_skin``.
    ``section_at(z)`` returns a frozen ``SectionCut`` (``walls``, ``cells``);
    ``cut.shell_graph()`` derives the loftable graph. That is not
    ``SectionGeometryResult`` (the gated station product). Object graph:
    ``docs/blade/README.md`` (In-memory objects; Loft cut versus station result).

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
        with ``web_positions`` and ``web_x_from_pitch``.
    web_x_from_pitch:
        Web positions in metres from the pitch axis, toward the trailing edge.
        A float per web, nested ``(suction, pressure)`` pairs, or a callable of
        ``z``. Converted with this station's chord and pivot, so the fraction is
        of that station's own chord. Exclusive with ``web_positions`` and
        ``web_offsets``.
    cap_x_from_pitch:
        Cap windows in the same metres. A ``(x0, x1)`` pair, a sequence of
        pairs, or a callable of ``z``. Exclusive with ``cap_windows``.
    twist_sense:
        ``"nose_up"`` (default) places a positive station twist with the leading
        edge toward +y. ``"nose_down"`` uses the other sense for the same number.
        The default is the one the golden sections were built with.
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
    transverse:
        Ribs, bulkheads and frames at fixed ``z``. A member appears as a closed
        rib face on ``section_at(z).shell_graph().faces`` at that station only.
        A layout's ``nR`` fills this with ``nR`` interior ribs when it is empty.
    stringers:
        Partial-span or free-end branches. Present on ``shell_graph().branches``
        inside ``extent``; ends are topology events.
        ply_drops:
        Geometry-only extra thickness windows in ``(z, s)`` with linear tapers.
    merge_cap_and_skin:
        If True, ``section_at`` walls follow ``shell_graph`` runs (caps as
        partitions, combined thickness). Default False keeps stacked cap solids
        (golden numbers).
    tip_closure:
        The face that closes the skin at the tip station. Every blade has one
        (default ``TipClosure()``: skin partition, skin thickness at the tip). It
        is not a rib and not in ``transverse``.
    """

    def __init__(
        self,
        stations: Sequence[Station],
        layout: str | SystemLayoutSpec | SectionLayout,
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
        ply_drops: Sequence[PlyDrop] = (),
        merge_cap_and_skin: bool = False,
        region_graph=None,
        tip_closure: Optional[TipClosure] = None,
        airfoil_morph: str = "auto",
        web_x_from_pitch: Optional[SpanTable] = None,
        cap_x_from_pitch: Optional[SpanTable] = None,
        twist_sense: str = "nose_up",
    ):
        from implicit_geometry_kernel.catalogue import resolve_system_type

        self.stations = tuple(sorted(stations, key=lambda s: s.z))
        if len(self.stations) < 2:
            raise ValueError("a blade needs at least two stations")
        if np.any(np.diff([s.z for s in self.stations]) <= 0.0):
            raise ValueError("stations need distinct z")
        spec: Union["SystemLayoutSpec", "SectionLayout"]
        if isinstance(layout, str):
            spec = resolve_system_type(layout)
        else:
            spec = layout
        if sum(v is not None for v in (web_positions, web_offsets, web_x_from_pitch)) > 1:
            raise ValueError("give only one of web_positions, web_offsets, web_x_from_pitch")
        if cap_windows is not None and cap_x_from_pitch is not None:
            raise ValueError("give cap_windows or cap_x_from_pitch, not both")
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
        if self.web_offsets is None and web_positions is None and web_x_from_pitch is None and n_webs > 0:
            from implicit_geometry_kernel.blade.concepts.web_anchoring import (
                lc_web_fractions,
                le_web_offsets_m,
                pa_web_x_from_pitch,
            )

            c0 = float(self.stations[0].chord)
            fracs = lc_web_fractions(spec.web_chord_fracs)
            if anchoring == "leading_edge":
                self.web_offsets = le_web_offsets_m(fracs, c0)
            elif anchoring == "pitch_axis":
                self._web_x_from_pitch = pa_web_x_from_pitch(
                    self.stations, self._pitch, n_webs
                )
        if self.web_offsets is not None and len(self.web_offsets) != spec.n_webs:
            raise ValueError(f"layout has {spec.n_webs} webs; got {len(self.web_offsets)} offsets")
        if self._web_x_from_pitch is not None and len(self._web_x_from_pitch) != spec.n_webs:
            raise ValueError(f"layout has {spec.n_webs} webs; got {len(self._web_x_from_pitch)} pitch offsets")
        self._web_x_from_pitch_fn = _metres_law(web_x_from_pitch, "web_x_from_pitch")
        self._cap_x_from_pitch_fn = _metres_law(cap_x_from_pitch, "cap_x_from_pitch", pairs=True)
        if self._web_x_from_pitch_fn is not None:
            sample = self._pitch_metres_at(float(self.stations[0].z))
            if sample is not None and len(sample) != n_webs:
                raise ValueError(
                    f"layout has {n_webs} webs; got {len(sample)} pitch-axis positions"
                )
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
        self.transverse, self.tip_closure = _split_tip_closure(
            tuple(transverse), tip_closure, float(self.stations[-1].z)
        )
        n_ribs = int(getattr(spec, "n_ribs", 0) or 0)
        if n_ribs > 0 and not self.transverse:
            self.transverse = catalogue_ribs(self.stations, n_ribs)
        self.stringers = tuple(stringers)
        self.ply_drops = tuple(ply_drops)
        self.merge_cap_and_skin = bool(merge_cap_and_skin)
        from implicit_geometry_kernel.blade.loft.loft_tables import _MORPH_RULES

        if airfoil_morph not in _MORPH_RULES:
            raise ValueError(
                f"airfoil_morph must be one of {', '.join(_MORPH_RULES)}; got {airfoil_morph!r}"
            )
        self.airfoil_morph = str(airfoil_morph)
        if twist_sense not in ("nose_up", "nose_down"):
            raise ValueError("twist_sense must be 'nose_up' or 'nose_down'")
        self.twist_sense = str(twist_sense)
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
        # Set by from_section: already-blended samples, parallel to self.stations.
        # interpolate lerps these and does not call _morph_surfaces.
        self._section_samples = None
        # interpolate / section_at / query / write / validate live on mixins

    @classmethod
    def from_section(
        cls,
        section,
        z,
        tolerance: float,
        layout: str | SystemLayoutSpec | SectionLayout,
        skin_thickness: Thickness,
        *,
        knots: Sequence[float] = (),
        n_airfoil_points: int = 161,
        **kwargs,
    ):
        """Loft ``section(z)`` until the linear loft is within ``tolerance`` metres.

        ``section(z)`` is the consumer's blend at that span: a
        :class:`~implicit_geometry_kernel.blade.loft.section_callable.SectionSample`,
        a mapping, or ``(chord, twist, pivot_fraction, offset, airfoil)``. Stations
        are added, including ``knots``, until a straight loft between neighbours lies
        within ``tolerance`` of ``section(z)`` on the placed outer mould line. The
        60 % smoothstep morph is not applied on this path; it stays the rule for
        blades built from a station table.

        ``z`` is ``(z_root, z_tip)`` or a longer list of spans that must be stations.
        Chord, twist, pitch axis and offset come from the callable. ``pivot_fraction``
        is still required by the constructor: pass it, or it is taken from the first
        sample (a varying pitch axis is then carried by the samples, not the table).
        """
        from implicit_geometry_kernel.blade.loft.section_callable import (
            airfoil_from_surfaces,
            cosine_x,
            refine_section_samples,
        )

        grid = cosine_x(n_airfoil_points)
        samples = refine_section_samples(section, z, tolerance, grid, knots=knots)
        stations = tuple(
            Station(
                s.z,
                s.chord,
                airfoil_from_surfaces(s.x, s.upper, s.lower, name=f"section-z{s.z:g}"),
                twist=s.twist,
                pitch_axis=s.pivot_fraction,
                offset=s.offset,
            )
            for s in samples
        )
        if "pivot_fraction" not in kwargs and "pitch_axis" not in kwargs:
            kwargs["pivot_fraction"] = float(samples[0].pivot_fraction)
        kwargs["n_airfoil_points"] = n_airfoil_points
        geo = cls(stations, layout, skin_thickness, **kwargs)
        by_z = {s.z: s for s in samples}
        geo._section_samples = tuple(by_z[st.z] for st in geo.stations)
        geo._section_cache.clear()
        return geo

    @classmethod
    def from_concept(cls, concept, stations: Sequence[Station], skin_thickness: Thickness, **kwargs):
        """Build from a :class:`~implicit_geometry_kernel.blade.concepts.BladeConcept`.

        Extra kwargs overlay the compiled scheme (thicknesses, web positions, …).
        ``pivot_fraction`` is required (or the ``pitch_axis`` alias).
        """
        compiled = dict(concept.to_geometry_kwargs(stations=stations))
        compiled.update(kwargs)
        layout = compiled.pop("layout")
        return cls(stations, layout, skin_thickness, **compiled)

    @property
    def concept(self):
        """Layout plus the spanwise extras already stored on this geometry."""
        from implicit_geometry_kernel.blade.concepts import BladeConcept, SpanwiseScheme

        layout = self.layout if self.layout is not None else self.spec
        return BladeConcept(
            layout,
            SpanwiseScheme(
                transverse=self.transverse,
                stringers=self.stringers,
                ply_drops=self.ply_drops,
                tip_closure=self.tip_closure,
            ),
        )

    @property
    def tip_closure_z(self) -> float:
        """Span position of the tip closure: the last station."""
        return float(self.stations[-1].z)

    def tip_closure_thickness(self) -> float:
        """The tip closure's thickness: its own, or the skin thickness law at the tip."""
        if self.tip_closure.thickness is not None:
            return float(self.tip_closure.thickness)
        return float(self._skin(self.tip_closure_z))

    def forced_stations(self) -> tuple[float, ...]:
        """Sorted unique member stations (ribs, stringer ends, ply drops and the ends of
        partial-span webs; mesh hint)."""
        by_key: dict[float, float] = {}
        z_lo, z_hi = float(self.stations[0].z), float(self.stations[-1].z)
        for web in getattr(self.region_graph, "webs", ()) or ():
            for z in (web.z0, web.z1):
                if z is not None and z_lo < float(z) < z_hi:
                    by_key[round(float(z), 9)] = float(z)
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

    def __repr__(self) -> str:
        return f"BladeGeometry({self.layout!r}, z = {self.stations[0].z:g} to {self.stations[-1].z:g})"



def _split_tip_closure(
    transverse: tuple[TransverseMember, ...], tip_closure: Optional[TipClosure], z_tip: float
) -> tuple[tuple[TransverseMember, ...], TipClosure]:
    """Separate a legacy ``TransverseMember("tip_closure")`` from the ribs.

    Before 0.9.2 the tip closure was a transverse member. A member of that name at
    the tip becomes the :class:`TipClosure` (keeping its thickness), with a
    ``DeprecationWarning``; one anywhere else is an error.
    """
    legacy = [m for m in transverse if str(m.name) == "tip_closure"]
    if not legacy:
        return transverse, tip_closure if tip_closure is not None else TipClosure()
    if tip_closure is not None:
        raise ValueError("give tip_closure= or a transverse member named 'tip_closure', not both")
    if len(legacy) > 1 or abs(float(legacy[0].z) - z_tip) > 1e-9:
        raise ValueError(
            f"a transverse member named 'tip_closure' must sit at the tip (z = {z_tip:g}); "
            "the tip closure is now BladeGeometry(tip_closure=TipClosure(...))"
        )
    warnings.warn(
        "TransverseMember('tip_closure') is deprecated: every BladeGeometry has a tip closure; "
        "pass tip_closure=TipClosure(thickness=...) to set its thickness",
        DeprecationWarning,
        stacklevel=3,
    )
    ribs = tuple(m for m in transverse if str(m.name) != "tip_closure")
    return ribs, TipClosure(thickness=float(legacy[0].thickness))


for _mixin in (LoftCutMixin, LoftFieldsMixin):
    for _name, _val in vars(_mixin).items():
        if _name.startswith("_"):
            continue
        if callable(_val) or isinstance(_val, (property, classmethod, staticmethod)):
            setattr(BladeGeometry, _name, _val)


def section_at(stations: Sequence[Station], z: float, layout, skin_thickness: Thickness, **options) -> SectionCut:
    """One-off cut: ``BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)``.

    ``pivot_fraction`` is required on ``options`` (or the ``pitch_axis`` alias).
    """
    return BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)




__all__ = ["BladeGeometry", "SpanTable", "Thickness", "section_at"]
