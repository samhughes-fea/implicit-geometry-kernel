"""Spanwise blade geometry: the implicit section geometry lofted between stations.

``BladeGeometry`` takes the ``Station`` list a ``BladeBuilder`` takes (chord, twist,
pitch axis, offset and airfoil at each span position) and a structural layout from
``section`` (one of its 40 system types, such as ``"2D-F"``: two webs and a
box spar cap between them). ``interpolate(z)`` is the only interpolator: it blends the
airfoil graphs linearly between the bounding stations and reads chord, twist, pitch
axis and offset either from those stations or from optional span tables. ``section_at``
and ``oml_polyline`` both call it. Tables and station scalars are linear in this
release; a later PCHIP may apply to scalars only (``[splines]``), never to 3-D-morph
airfoils.

``section_at(z)`` returns the wall midlines with their thickness, region labels and the
closed cell-wall loops, placed in the blade frame. ``oml_polyline(z)`` is the closed
outer mould line in the same frame.

Blade frame (the ``BladeBuilder`` axes): z along the span, x along the chord from the
pitch axis towards the trailing edge, y through the thickness (suction side positive).
Positive twist is nose up: the leading edge moves to +y. ``offset`` shifts the pitch
axis point (prebend, sweep). The section geometry is built at the station's twist, so
a layout's flapwise webs (``"-F"``) stay normal to the rotor plane (along y) while its
chord-normal webs (``"-CN"``) turn with the chord; the two are the same untwisted.

Wall thicknesses are geometry, not material: they set where each wall's midline lies
(skins grow inward from the outer mould line, caps sit inside the skin, webs end on the
skin midline). Give each as a number or as a function of z. Laminates are bound to the
walls by label downstream (``ThinWalledSection.from_geometry``).

Web positions: by default the layout's chord fractions at every station. Give
``web_positions`` (chord fractions, the same at every station) to place them yourself,
or ``web_offsets`` (metres from the leading edge, fixed along the span, as on a tapered
blade whose webs are straight); with offsets, webs that no longer fit a thin station
are dropped rather than slid, and ``SectionCut.dropped_webs`` lists them.
"""

from __future__ import annotations

import warnings
from collections import OrderedDict
from typing import Callable, Optional, Sequence, Union

import numpy as np

from implicit_geometry_kernel.blade.cut import CellPolygon, SectionCut, SectionWall, _twist_offset
from implicit_geometry_kernel.blade.stations import Station, Stringer, TransverseMember
from implicit_geometry_kernel.section.layout.thickness_limits import (
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M as _MIN_SKIN,
    MIN_REALISTIC_SPAR_WALL_THICKNESS_M as _MIN_CAP,
    MIN_REALISTIC_WEB_WALL_THICKNESS_M as _MIN_WEB,
)

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


class BladeGeometry:
    """A blade or wing lofted from stations, cut into thin-walled sections on demand.

    Parameters
    ----------
    stations:
        Two or more ``Station`` objects in increasing ``z``.
    layout:
        A catalogue key (``"2D-F"``, …), a ``SystemLayoutSpec``, or a parametric
        ``SectionLayout`` (web fractions and cap rule, as for reference blades).
    skin_thickness, spar_cap_thickness, shear_web_thickness:
        Wall thicknesses in metres, as numbers or functions of z.
    web_positions:
        Web chord fractions, LE to TE (default: the layout's). A callable ``f(z)``
        or ``(z, values)`` table with shape ``(n, n_webs)`` varies stations along
        the span.
    cap_windows:
        Optional chord windows ``((lo, hi), …)`` or a span table of the same,
        used as a ``CapPlacementRule`` when ``layout`` does not already carry one.
    web_offsets:
        Web positions in metres from the leading edge, fixed along the span. Exclusive
        with ``web_positions``.
    n_airfoil_points:
        Points per surface used to hand the interpolated airfoil to the section
        geometry.
    chord, twist, pitch_axis, offset:
        Optional span tables that override the station scalars at ``z`` (callable or
        ``(z, values)``; ``offset`` values are shape ``(n, 2)``). Airfoil graphs still
        blend linearly between the bounding ``Station`` foils.
    span_kind:
        ``"linear"`` (default) or ``"pchip"`` for station chord, twist, pitch axis
        and offset. Airfoil graphs stay linear blends.
    transverse:
        Ribs, bulkheads and frames at fixed ``z``. A member appears as a closed
        rib face on ``section_at(z).shell_graph().faces`` at that station only.
    stringers:
        Partial-span or free-end branches. Present on ``shell_graph().branches``
        inside ``extent``; ends are topology events.
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
        offset: Optional[SpanTable] = None,
        cap_windows: Optional[SpanTable] = None,
        transverse: Sequence[TransverseMember] = (),
        stringers: Sequence[Stringer] = (),
        span_kind: str = "linear",
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
        if self.web_offsets is not None and len(self.web_offsets) != spec.n_webs:
            raise ValueError(f"layout has {spec.n_webs} webs; got {len(self.web_offsets)} offsets")
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
        self._pitch = _as_span_table(pitch_axis, "pitch_axis")
        self._offset = _as_span_table(offset, "offset", cols=2)
        self.n_airfoil_points = int(n_airfoil_points)
        self.transverse = tuple(transverse)
        self.stringers = tuple(stringers)
        if span_kind not in ("linear", "pchip"):
            raise ValueError("span_kind must be 'linear' or 'pchip'")
        self.span_kind = span_kind
        # Live fields (wall_meta / stack T) key stations by round(z, 9). Cache the
        # cut itself so those fields do not rebuild geometry on a second hit.
        self._section_cache: OrderedDict[float, SectionCut] = OrderedDict()
        self._section_cache_size = 256

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
        return tuple(by_key[key] for key in sorted(by_key))

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
        if self._cap_windows_fn is not None and layout.n_webs == 0:
            windows = self._cap_windows_fn(float(z))
            layout = SectionLayout(
                web_chord_fracs=layout.web_chord_fracs,
                orientation=layout.orientation,
                cap_placement=CapPlacementRule.from_chord_windows(windows),
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

        Airfoil graphs blend linearly between the bounding stations. Chord, twist,
        pitch axis and offset come from span tables when given, otherwise from those
        stations.
        """
        zs = np.array([s.z for s in self.stations])
        if z < zs[0] - 1e-12 * max(1.0, abs(zs[0])) or z > zs[-1] + 1e-12 * max(1.0, abs(zs[-1])):
            raise ValueError(f"z = {z} lies outside the stations ({zs[0]} to {zs[-1]})")
        k = int(np.clip(np.searchsorted(zs, z, side="right") - 1, 0, len(zs) - 2))
        a, b = self.stations[k], self.stations[k + 1]
        w = float(np.clip((z - a.z) / (b.z - a.z), 0.0, 1.0))

        def lerp(p, q):
            return (1.0 - w) * np.asarray(p, float) + w * np.asarray(q, float)

        x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, self.n_airfoil_points)))
        upper = lerp(a.airfoil.surface("upper", x), b.airfoil.surface("upper", x))
        lower = lerp(a.airfoil.surface("lower", x), b.airfoil.surface("lower", x))
        zz = float(z)
        chord = self._span_scalar(self._chord, [s.chord for s in self.stations], zz)
        twist = self._span_scalar(self._twist, [s.twist for s in self.stations], zz)
        pitch = self._span_scalar(self._pitch, [s.pitch_axis for s in self.stations], zz)
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

    def station_at(self, z: float):
        """Deprecated alias of :meth:`interpolate`. Removed in 0.6.0."""
        warnings.warn(
            "BladeGeometry.station_at is deprecated; use interpolate. Removed in 0.6.0.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.interpolate(z)

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
        from implicit_geometry_kernel.frames import PIVOT_FRACTION_DEFAULT, twisted_pitch_axis_blade_frame_to_naca
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
            outputs=("midlines", "cell_wall_loops"),
        )
        if self.web_offsets is not None:
            kwargs["web_offsets"] = tuple(v * k for v in self.web_offsets)
        runner = SectionGeometry(**kwargs).execute()
        result = runner.get_results()
        resolution = runner.get_layout_resolution()

        def place(pts_B):
            pts = np.asarray(pts_B, float) / k
            naca = twisted_pitch_axis_blade_frame_to_naca(pts, _FORK_TWIST_SIGN * np.radians(twist), chord,
                                                          PIVOT_FRACTION_DEFAULT)
            local = np.column_stack([naca[:, 0] - pitch * chord, naca[:, 1]])
            xy = _twist_offset(local, twist, offset)
            return np.column_stack([xy, np.full(len(xy), z)])

        strips = list(result.midline_strips or ())
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
        key = getattr(resolution.layout, "catalogue_key", None) or self.layout
        cut = SectionCut(z, chord, twist, pitch, (float(offset[0]), float(offset[1])), str(key), tuple(walls), cells,
                         dropped, self._transverse_at(z), self._stringers_at(z))
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
        """Stacked wall mid-surface: ``(oml + cavity) / 2``.

        Snap oracle: run interiors match ``shell_graph`` midpoints;
        partition vertices are not a Newton spec (``|phi|`` up to ``t_cap/2``).
        """
        from implicit_geometry_kernel.sdf.nodes.tagged import MidSurfaceNode

        return MidSurfaceNode(self.oml_sdf(meta), self.cavity_sdf(meta))

    def wall_meta(self, table=None):
        """Tagged wall-strip DAG: ``evaluate((x,y,z)).meta`` is the region at that ``z``."""
        from implicit_geometry_kernel.blade.fields import CutWallMetaSDF
        from implicit_geometry_kernel.sdf.nodes.base import FeatureTable

        return CutWallMetaSDF(self, table if table is not None else FeatureTable())

    def snap_oml(self, points, solver=None):
        """Project points onto the OML in each point's own ``z``-plane."""
        from implicit_geometry_kernel.sdf.solve import GeometrySolver

        solver = solver or GeometrySolver(tolerance=1e-8, max_step=0.05)
        return solver.project(self.oml_sdf(), np.asarray(points, float), plane=True)

    def snap_wall_mid(self, points, solver=None):
        """Project onto ``wall_mid_sdf`` in each point's ``z``-plane (not a spec at ``T``-step joints)."""
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
    """One-off cut: ``BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)``."""
    return BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)


__all__ = ["BladeGeometry", "SpanTable", "Thickness", "section_at"]
