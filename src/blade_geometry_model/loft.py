"""Spanwise blade geometry: the implicit section geometry lofted between stations.

``BladeGeometry`` takes the ``Station`` list a ``BladeBuilder`` takes (chord, twist,
pitch axis, offset and airfoil at each span position, interpolated linearly between
them) and a structural layout from ``section_geometry`` (one of its 40 system types,
such as ``"2D-F"``: two webs and a box spar cap between them), and cuts the section at
any span position: ``section_at(z)`` returns the wall midlines with their thickness,
region labels and the closed cell-wall loops, placed in the blade frame.

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

from dataclasses import dataclass, field
from typing import Callable, Dict, Optional, Sequence, Tuple, Union

import numpy as np

from .blade_geometry import Station

Thickness = Union[float, Callable[[float], float]]

# the section geometry's wall-thickness floors (layout/thickness_limits.py)
_MIN_SKIN, _MIN_WEB, _MIN_CAP = 0.0025, 0.0040, 0.0050
# section_geometry turns twist the other way (its positive twist is nose down). The
# section is built twisted, so that flapwise ("-F") webs stay normal to the rotor
# plane while chord-normal ("-CN") webs turn with the chord.
_FORK_TWIST_SIGN = -1.0


@dataclass(frozen=True)
class SectionWall:
    """One wall midline of a cut section, in the blade frame.

    ``kind`` is ``"skin"``, ``"web"`` or ``"cap"``. ``label`` is the section
    geometry's region label (``"skin_suction_surface:bay_0"``, ``"shear_web_1"``,
    ``"spar_cap_pressure_surface"``, ...). A cap lies inside the skin it is bonded to:
    ``bonds`` maps each host skin label to the (cap vertex, skin vertex) pairs that
    face each other across the bond, as index arrays into the two midlines.
    """

    label: str
    kind: str
    points: np.ndarray  # (n, 3): x, y, z
    thickness: float
    side: Optional[str] = None  # "suction" / "pressure" for skins and caps
    bonds: Dict[str, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)

    @property
    def xy(self) -> np.ndarray:
        return self.points[:, :2]

    @property
    def length(self) -> float:
        return float(np.sum(np.linalg.norm(np.diff(self.xy, axis=0), axis=1)))

    @property
    def region(self) -> str:
        """The label without its part suffix: ``"skin_suction_surface"`` for ``"skin_suction_surface:bay_0"``."""
        return self.label.split(":")[0]

    @property
    def lookup_keys(self) -> Tuple[str, ...]:
        """Names a property (a laminate, an element count) is looked up by, most specific first."""
        return tuple(dict.fromkeys((self.label, self.region, self.kind)))


@dataclass(frozen=True)
class CellLoop:
    """A closed cell of the section: its wall-midline polygon and the walls around it."""

    id: int
    points: np.ndarray  # (n, 3), first point repeated at the end
    area: float
    wall_labels: Tuple[str, ...]


@dataclass(frozen=True)
class SectionCut:
    """The section of a ``BladeGeometry`` at span position ``z``."""

    z: float
    chord: float
    twist: float
    pitch_axis: float
    offset: Tuple[float, float]
    layout: str
    walls: Tuple[SectionWall, ...]
    cells: Tuple[CellLoop, ...]
    dropped_webs: Tuple[str, ...] = ()

    def wall(self, label: str) -> SectionWall:
        for w in self.walls:
            if w.label == label:
                return w
        raise KeyError(f"no wall '{label}'. Walls: {', '.join(w.label for w in self.walls)}")

    def walls_of(self, kind: str) -> Tuple[SectionWall, ...]:
        return tuple(w for w in self.walls if w.kind == kind)

    def to_blade(self, x, y) -> np.ndarray:
        """Blade-frame (x, y) of untwisted section coordinates measured from the pitch axis."""
        return _twist_offset(np.column_stack([np.ravel(x), np.ravel(y)]), self.twist, self.offset)

    def mid_surface(self) -> "MidSurface":
        """The section's walls as one graph of straight segments on their mid-surfaces.

        Skins and webs become one segment per midline segment, joined where they meet.
        A spar cap is merged with the skin it is bonded to: along the bond the segment
        is the cap and skin stacked (``layers`` = (cap, skin), inside first), on the
        stack's mid-surface, and is named after the cap. Segments are oriented so that
        z x (end - start), the laminate normal of a shell or thin wall laid along them,
        points out of the section on the skins and towards the trailing edge on the
        webs. Both ``ThinWalledSection.from_geometry`` and the shell mesher build on it.
        """
        return _mid_surface(self)


@dataclass(frozen=True)
class MidSegment:
    """One straight piece of a wall mid-surface: point indices, name and wall layers (inside first)."""

    start: int
    end: int
    name: str
    kind: str  # "skin" or "web" ("skin" also for a cap merged into its skin)
    layers: Tuple[str, ...]  # SectionWall labels, inside first: (cap, skin) along a bond


@dataclass(frozen=True)
class MidSurface:
    """The walls of a ``SectionCut`` as a graph of oriented mid-surface segments."""

    points: np.ndarray  # (n, 2) blade-frame x, y
    segments: Tuple[MidSegment, ...]
    z: float = 0.0

    @property
    def length(self) -> np.ndarray:
        P = self.points
        return np.array([float(np.linalg.norm(P[s.end] - P[s.start])) for s in self.segments])


def _chain(strips, tol):
    """Order open polylines into one closed loop: [(label, xy, reversed)], or None."""
    if not strips:
        return None
    left = list(strips)
    label, xy = left.pop(0)
    out = [(label, xy, False)]
    end = xy[-1]
    while left:
        for k, (lb, p) in enumerate(left):
            if np.linalg.norm(p[0] - end) <= tol:
                out.append((lb, p, False))
                end = p[-1]
                break
            if np.linalg.norm(p[-1] - end) <= tol:
                out.append((lb, p, True))
                end = p[0]
                break
        else:
            return None
        left.pop(k)
    if np.linalg.norm(end - xy[0]) > tol:
        return None
    return out


def _mid_surface(cut: SectionCut) -> MidSurface:
    tol = 1e-9 * max(1.0, max(float(np.abs(w.xy).max()) for w in cut.walls))
    pts = []
    lookup: Dict[Tuple[int, int], list] = {}

    def index(p):
        key = (int(round(p[0] / (100 * tol))), int(round(p[1] / (100 * tol))))
        for dk in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)):
            for k in lookup.get((key[0] + dk[0], key[1] + dk[1]), ()):
                if np.linalg.norm(pts[k] - p) <= tol:
                    return k
        pts.append(np.array(p, float))
        lookup.setdefault(key, []).append(len(pts) - 1)
        return len(pts) - 1

    skins, webs, caps = cut.walls_of("skin"), cut.walls_of("web"), cut.walls_of("cap")
    # skins: orient every skin so that z x s points out of the section
    flip = {w.label: False for w in skins}
    loop = _chain([(w.label, w.xy) for w in skins], tol)
    if loop is not None:
        ring = np.vstack([xy[::-1] if rev else xy for _, xy, rev in loop])
        area = 0.5 * np.sum(ring[:, 0] * np.roll(ring[:, 1], -1) - np.roll(ring[:, 0], -1) * ring[:, 1])
        for label, _, rev in loop:
            flip[label] = rev != (area > 0.0)  # counter-clockwise: z x s points in
    # caps: which skin vertices carry which cap, and where the stack's mid-surface lies
    capped: Dict[Tuple[str, int], str] = {}
    moved: Dict[Tuple[str, int], np.ndarray] = {}
    skin_of = {w.label: w for w in skins}
    for cap in caps:
        for host_label, (ic, ih) in cap.bonds.items():
            host = skin_of.get(host_label)
            if host is None:
                continue
            if len(ih) < 2 or np.sum(np.linalg.norm(np.diff(host.xy[ih], axis=0), axis=1)) < cap.thickness:
                # the cap's edge reaching just past a web onto the next skin: whether that
                # sliver holds a vertex depends on sampling, so it is left to the skin
                continue
            f = cap.thickness / (cap.thickness + host.thickness)
            for a, b in zip(ic, ih):
                capped[(host_label, int(b))] = cap.label
                moved[(host_label, int(b))] = host.xy[b] + f * (cap.xy[a] - host.xy[b])
            # likewise a sliver of bare skin between the cap's edge and the end of its skin
            arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(host.xy, axis=0), axis=1))])
            lo, hi = int(np.min(ih)), int(np.max(ih))
            for edge, rest, gap in ((lo, range(0, lo), arc[lo]), (hi, range(hi + 1, len(host.xy)), arc[-1] - arc[hi])):
                if len(rest) and gap < cap.thickness:
                    step = moved[(host_label, edge)] - host.xy[edge]
                    for v in rest:
                        capped[(host_label, v)] = cap.label
                        moved[(host_label, v)] = host.xy[v] + step
    segments = []
    shift: Dict[int, np.ndarray] = {}
    for w in skins:
        ids = [index(p) for p in w.xy]
        for v, k in enumerate(ids):
            if (w.label, v) in moved:
                shift[k] = moved[(w.label, v)]
        run = []
        for v in range(len(ids) - 1):
            i, j = ids[v], ids[v + 1]
            if i == j:
                continue
            ca, cb = capped.get((w.label, v)), capped.get((w.label, v + 1))
            if ca is not None and ca == cb:
                seg = (i, j, ca, (ca, w.label))
            else:
                seg = (i, j, w.label, (w.label,))
            run.append(seg)
        if flip[w.label]:
            run = [(j, i, name, layers) for i, j, name, layers in reversed(run)]
        segments.extend(MidSegment(i, j, name, "skin", layers) for i, j, name, layers in run)
    t = np.radians(cut.twist)
    chord_dir = np.array([np.cos(t), -np.sin(t)])
    for w in webs:
        ids = [index(p) for p in w.xy]
        run = [(i, j) for i, j in zip(ids[:-1], ids[1:]) if i != j]
        if run:
            d = pts[run[-1][1]] - pts[run[0][0]]
            if np.array([-d[1], d[0]]) @ chord_dir < 0.0:
                run = [(j, i) for i, j in reversed(run)]
        segments.extend(MidSegment(i, j, w.label, "web", (w.label,)) for i, j in run)
    P = np.array(pts)
    for k, p in shift.items():
        P[k] = p
    return MidSurface(P, tuple(segments), cut.z)


def _twist_offset(xy: np.ndarray, twist_deg: float, offset) -> np.ndarray:
    t = np.radians(twist_deg)
    c, s = np.cos(t), np.sin(t)
    x, y = xy[:, 0], xy[:, 1]
    return np.column_stack([c * x + s * y + offset[0], -s * x + c * y + offset[1]])


def _as_function(value: Thickness, name: str) -> Callable[[float], float]:
    if callable(value):
        return value
    v = float(value)
    if v < 0.0:
        raise ValueError(f"{name} must be >= 0")
    return lambda z: v


class BladeGeometry:
    """A blade or wing lofted from stations, cut into thin-walled sections on demand.

    Parameters
    ----------
    stations:
        Two or more ``Station`` objects in increasing ``z``.
    layout:
        A ``section_geometry`` system type key (``"2D-F"``, ``"3D-F"``, ``"0B"``,
        ...; ``section_geometry.SYSTEM_TYPE_KEYS`` lists all 40) or a
        ``SystemLayoutSpec``.
    skin_thickness, spar_cap_thickness, shear_web_thickness:
        Wall thicknesses in metres, as numbers or functions of z.
    web_positions:
        Web chord fractions, LE to TE (default: the layout's).
    web_offsets:
        Web positions in metres from the leading edge, fixed along the span. Exclusive
        with ``web_positions``.
    n_airfoil_points:
        Points per surface used to hand the interpolated airfoil to the section
        geometry.
    """

    def __init__(
        self,
        stations: Sequence[Station],
        layout,
        skin_thickness: Thickness,
        spar_cap_thickness: Thickness = 0.0,
        shear_web_thickness: Thickness = 0.0,
        web_positions: Optional[Sequence[float]] = None,
        web_offsets: Optional[Sequence[float]] = None,
        n_airfoil_points: int = 161,
    ):
        from .section_geometry import resolve_system_type

        self.stations = tuple(sorted(stations, key=lambda s: s.z))
        if len(self.stations) < 2:
            raise ValueError("a blade needs at least two stations")
        if np.any(np.diff([s.z for s in self.stations]) <= 0.0):
            raise ValueError("stations need distinct z")
        spec = resolve_system_type(layout) if isinstance(layout, str) else layout
        if web_positions is not None and web_offsets is not None:
            raise ValueError("give web_positions or web_offsets, not both")
        if web_positions is not None:
            fr = tuple(float(f) for f in web_positions)
            if len(fr) != spec.n_webs:
                raise ValueError(f"layout {spec.key} has {spec.n_webs} webs; got {len(fr)} positions")
            if any(not 0.0 < f < 1.0 for f in fr) or np.any(np.diff(fr) <= 0.0):
                raise ValueError("web_positions must be increasing chord fractions inside (0, 1)")
            spec = _with_web_fractions(spec, fr)
        self.spec = spec
        self.layout = spec.key
        self.web_offsets = None if web_offsets is None else tuple(float(v) for v in web_offsets)
        if self.web_offsets is not None and len(self.web_offsets) != spec.n_webs:
            raise ValueError(f"layout {spec.key} has {spec.n_webs} webs; got {len(self.web_offsets)} offsets")
        self._skin = _as_function(skin_thickness, "skin_thickness")
        self._cap = _as_function(spar_cap_thickness, "spar_cap_thickness")
        self._web = _as_function(shear_web_thickness, "shear_web_thickness")
        self.n_airfoil_points = int(n_airfoil_points)

    # ------------------------------------------------------------ stations
    def station_at(self, z: float):
        """(chord, twist, pitch_axis, offset, (upper, lower) surfaces) interpolated at z."""
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
        return (float(lerp(a.chord, b.chord)), float(lerp(a.twist, b.twist)), float(lerp(a.pitch_axis, b.pitch_axis)),
                tuple(lerp(a.offset, b.offset)), (x, upper, lower))

    # ------------------------------------------------------------ cutting
    def section_at(self, z: float) -> SectionCut:
        """The wall midlines, thicknesses, labels and cell loops of the section at z."""
        from .coordinate_frames import PIVOT_FRACTION_DEFAULT, twisted_pitch_axis_blade_frame_to_naca
        from .section_geometry import AirfoilSDF, SectionGeometryRunner

        z = float(z)
        chord, twist, pitch, offset, (x, upper, lower) = self.station_at(z)
        t_skin, t_cap, t_web = self._skin(z), (self._cap(z) if self.spec.has_spar_caps else 0.0), self._web(z)
        # section_geometry floors wall thicknesses at wind-turbine manufacturing minima
        # (2.5 mm skin, 4 mm web, 5 mm cap). The geometry is similar under scaling, so a
        # thinner-walled section is built k times larger and scaled back: exact, and the
        # frozen geometry code is used as it is.
        k = max(1.0, *(lo / t for lo, t in ((_MIN_SKIN, t_skin), (_MIN_WEB, t_web), (_MIN_CAP, t_cap)) if t > 0.0))
        if k > 1.0:
            k *= 1.0 + 1e-9
        verts = np.vstack([np.column_stack([x[::-1], upper[::-1]]), np.column_stack([x[1:], lower[1:]])]) * chord * k
        if np.allclose(verts[0], verts[-1]):
            verts = verts[:-1]  # closed trailing edge: one TE point
        kwargs = dict(
            layout=self.spec,
            airfoil=AirfoilSDF(verts, chord=chord * k),
            skin_thickness=t_skin * k,
            spar_cap_thickness=t_cap * k,
            shear_web_thickness=t_web * k,
            twist_rad=_FORK_TWIST_SIGN * np.radians(twist),
            outputs=("midlines", "cell_wall_loops"),
        )
        if self.web_offsets is not None:
            kwargs["web_x_positions_m"] = tuple(v * k for v in self.web_offsets)
        runner = SectionGeometryRunner(**kwargs).execute()
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
            CellLoop(int(c.id), place(c.polyline_twisted_pitch_axis_blade_frame), float(c.area_m2) / k**2,
                     tuple(c.wall_labels))
            for c in (result.cell_wall_loops or ())
        )
        dropped = tuple(f"web at {d.x_m / k:.6g} m from the LE: {d.reason}" for d in resolution.dropped)
        key = getattr(resolution.layout, "catalogue_key", None) or self.layout
        return SectionCut(z, chord, twist, pitch, (float(offset[0]), float(offset[1])), str(key), tuple(walls), cells,
                          dropped)

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


def _with_web_fractions(spec, fractions):
    import dataclasses

    return dataclasses.replace(spec, web_chord_fracs=tuple(fractions))


def section_at(stations: Sequence[Station], z: float, layout, skin_thickness: Thickness, **options) -> SectionCut:
    """One-off cut: ``BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)``."""
    return BladeGeometry(stations, layout, skin_thickness, **options).section_at(z)


__all__ = ["BladeGeometry", "CellLoop", "MidSegment", "MidSurface", "SectionCut", "SectionWall", "section_at"]
