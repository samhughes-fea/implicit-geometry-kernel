"""A cut section in the blade frame: wall midlines, cells and the merged mid-surface.

``BladeGeometry.section_at(z)`` returns a ``SectionCut``. Two views of the same cut:

* **Solids / bonds / cells** — ``walls`` always lists every cap as its own
  ``SectionWall`` (``kind="cap"``). Do not delete those rows.
* **Shell graph** — ``mid_surface()`` never has ``kind="cap"``. A cap is a
  *partition* on the skin wall: ``layers = (cap, skin)``, combined
  ``T = t_cap + t_skin``, and the cap edge stays a feature line. ``shell_graph()``
  exports that as one loftable skin loop plus web branches. Sandwich cores and
  trailing-edge inserts stay constructive solids (G9); they are not loop runs.

``cells`` are the closed cell interiors (``CellPolygon``; ``CellLoop`` is a
deprecated alias until 0.6.0).

Graph types and builders live in :mod:`implicit_geometry_kernel.blade.shell_graph`
and are re-exported here until 0.6.0.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np

from implicit_geometry_kernel.blade.stations import TransverseMember


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
    """The section of a ``BladeGeometry`` at span position ``z``."""

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
        from implicit_geometry_kernel.blade.shell_graph import _mid_surface

        return _mid_surface(self)

    def shell_graph(self) -> "ShellGraph":
        """Loftable shell: one wall loop (caps as partitions) plus web branches."""
        from implicit_geometry_kernel.blade.shell_graph import shell_graph as _shell_graph

        return _shell_graph(self)


def _twist_offset(xy: np.ndarray, twist_deg: float, offset) -> np.ndarray:
    t = np.radians(twist_deg)
    c, s = np.cos(t), np.sin(t)
    x, y = xy[:, 0], xy[:, 1]
    return np.column_stack([c * x + s * y + offset[0], -s * x + c * y + offset[1]])


from implicit_geometry_kernel.blade.shell_graph import (  # noqa: E402
    MidSegment,
    MidSurface,
    PartitionEdge,
    ShellGraph,
    ShellRun,
    TopologyDiff,
    TopologyEvent,
    partition_edge_curve,
    shared_loft_names,
    shell_graph,
    shell_strips_from_graph,
    topology_diff,
)

def __getattr__(name: str):
    if name == "CellLoop":
        import warnings

        warnings.warn(
            "CellLoop is deprecated; use CellPolygon. Removed in 0.6.0.",
            DeprecationWarning,
            stacklevel=2,
        )
        return CellPolygon
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [  # CellLoop via __getattr__ until 0.6.0
    "CellLoop",  # noqa: F822
    "CellPolygon",
    "MidSegment",
    "MidSurface",
    "PartitionEdge",
    "SectionCut",
    "SectionWall",
    "ShellGraph",
    "ShellRun",
    "TopologyDiff",
    "TopologyEvent",
    "partition_edge_curve",
    "shared_loft_names",
    "shell_graph",
    "shell_strips_from_graph",
    "topology_diff",
]
