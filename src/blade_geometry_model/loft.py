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
axis point (prebend, sweep). The section geometry itself is built untwisted, with the
leading edge at the origin, and then placed.

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

from .blade_geometry import Airfoil, Station

Thickness = Union[float, Callable[[float], float]]

# the section geometry's wall-thickness floors (layout/thickness_limits.py)
_MIN_SKIN, _MIN_WEB, _MIN_CAP = 0.0025, 0.0040, 0.0050


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
            twist_rad=0.0,
            outputs=("midlines", "cell_wall_loops"),
        )
        if self.web_offsets is not None:
            kwargs["web_x_positions_m"] = tuple(v * k for v in self.web_offsets)
        runner = SectionGeometryRunner(**kwargs).execute()
        result = runner.get_results()
        resolution = runner.get_layout_resolution()

        def place(pts_B):
            naca = twisted_pitch_axis_blade_frame_to_naca(np.asarray(pts_B, float) / k, 0.0, chord, PIVOT_FRACTION_DEFAULT)
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


__all__ = ["BladeGeometry", "CellLoop", "SectionCut", "SectionWall", "section_at"]
