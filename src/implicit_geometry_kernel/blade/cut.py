"""A cut section in the blade frame: wall midlines, cells and the merged mid-surface.

``BladeGeometry.section_at(z)`` returns a ``SectionCut``. Its ``walls`` are the wall
midlines (``SectionWall``) with thickness, region label, side and cap-to-skin bonds;
its ``cells`` the closed cell loops (``CellLoop``); ``mid_surface()`` merges the caps
into their skins to give one wall graph (``MidSurface``) for shell meshing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np


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


__all__ = ["CellLoop", "MidSegment", "MidSurface", "SectionCut", "SectionWall"]
