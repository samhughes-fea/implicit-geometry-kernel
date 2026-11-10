"""A cut section in the blade frame: wall midlines, cells and the merged mid-surface.

``BladeGeometry.section_at(z)`` returns a ``SectionCut``. Two views of the same cut:

* **Solids / bonds / cells** — ``walls`` always lists every cap as its own
  ``SectionWall`` (``kind="cap"``). Do not delete those rows. Cap–skin
  ``bonds`` are the solid footprint on the host skin.
* **Shell graph** — conformal **skin mid**. A cap is a *partition* on that
  loop (``kind="cap"`` as colour, xy stays the skin). Mesh thickness is
  ``t_skin``; ``eccentricity_m = t_cap/2`` for constitutive OFFSET. Stack
  ``T = t_cap + t_skin`` remains the thickness oracle (``StackThickness`` /
  cavity), not the mesh gauge. ``shell_graph()`` exports one loftable skin
  loop plus web branches. Sandwich cores and inserts stay off the loop by
  default (G9); ``shell_graph(include_solids=True)`` adds polygon-boundary
  runs with ``surface`` ``core`` / ``insert``.

``cells`` are the closed cell interiors (``CellPolygon``).

Graph types and builders live in :mod:`implicit_geometry_kernel.blade.shell_graph`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple

import numpy as np

from implicit_geometry_kernel.blade.stations import Stringer, TransverseMember

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.shell_graph import MidSurface, ShellGraph


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
class CellPolygon:
    """A closed cell of the section: its wall-midline polygon and the walls around it."""

    id: int
    points: np.ndarray  # (n, 3), first point repeated at the end
    area: float
    wall_labels: Tuple[str, ...]


@dataclass(frozen=True)
class SectionCut:
    """The section of a ``BladeGeometry`` at span position ``z``.

    ``layout`` is the authored **concept** token (W/X/Y, span Z, transverse R).
    Webs that no longer fit are listed on ``dropped_web_indices``; the realised
    station slice is ``station_webs`` (W) and ``station_family`` (X). Those two
    are not a new catalogue token.
    """

    z: float
    chord: float
    twist: float
    pitch_axis: float
    offset: Tuple[float, float]
    layout: str
    walls: Tuple[SectionWall, ...]
    cells: Tuple[CellPolygon, ...]
    dropped_webs: Tuple[str, ...] = ()
    transverse: Tuple[TransverseMember, ...] = ()
    stringers: Tuple[Stringer, ...] = ()
    host: Any = field(default=None, repr=False, compare=False)
    scale: float = 1.0
    dropped_web_indices: Tuple[int, ...] = ()
    station_webs: int = 0
    station_family: str = ""

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
        webs. Shell meshers and an external thin-wall consumer build on it;
        ``implicit_geometry_kernel`` has no ``thin_wall`` module.
        """
        from implicit_geometry_kernel.blade.shell_graph import _mid_surface

        return _mid_surface(self)

    def place_host_xy(self, pts_B) -> np.ndarray:
        """Blade-frame polyline of section-frame points on this cut (G9 solids)."""
        from implicit_geometry_kernel.frames import twisted_pitch_axis_blade_frame_to_naca

        k = float(self.scale) if float(self.scale) else 1.0
        pts = np.asarray(pts_B, float) / k
        if pts.ndim != 2 or pts.shape[1] < 2:
            return np.zeros((0, 3))
        naca = twisted_pitch_axis_blade_frame_to_naca(
            pts[:, :2], -np.radians(self.twist), self.chord, self.pitch_axis
        )
        local = np.column_stack([naca[:, 0] - self.pitch_axis * self.chord, naca[:, 1]])
        xy = _twist_offset(local, self.twist, self.offset)
        return np.column_stack([xy, np.full(len(xy), self.z)])

    def shell_graph(self, include_solids: bool = False) -> "ShellGraph":
        """Loftable shell: one wall loop (caps as partitions) plus web branches.

        ``include_solids=True`` adds core and insert runs from polygon boundaries (G9).
        """
        from implicit_geometry_kernel.blade.shell_graph import shell_graph as _shell_graph

        return _shell_graph(self, include_solids=include_solids)

    def key_points(self) -> dict[str, np.ndarray]:
        """LE, TE, cap-edge vertices and web feet at this station."""
        oml = np.vstack([w.points for w in self.walls_of("skin")])
        le = oml[int(np.argmin(oml[:, 0]))]
        te = oml[int(np.argmax(oml[:, 0]))]
        out: dict[str, np.ndarray] = {"LE": le, "TE": te}
        graph = self.shell_graph()
        from implicit_geometry_kernel.section.assembly.members import wall_type_for_label

        def _is_cap(name: str) -> bool:
            wt = wall_type_for_label(name)
            return wt is not None and wt.kind == "cap"

        for i, edge in enumerate(graph.partition_edges()):
            if _is_cap(edge.name_a) or _is_cap(edge.name_b):
                out[f"cap_edge_{i}"] = np.asarray(edge.point, float)
        for web in graph.webs:
            pts = np.asarray(web.points, float)
            out[f"{web.name}_foot_0"] = pts[0]
            out[f"{web.name}_foot_1"] = pts[-1]
        return out


def _twist_offset(xy: np.ndarray, twist_deg: float, offset) -> np.ndarray:
    t = np.radians(twist_deg)
    c, s = np.cos(t), np.sin(t)
    x, y = xy[:, 0], xy[:, 1]
    return np.column_stack([c * x + s * y + offset[0], -s * x + c * y + offset[1]])


__all__ = [
    "CellPolygon",
    "SectionCut",
    "SectionWall",
]
