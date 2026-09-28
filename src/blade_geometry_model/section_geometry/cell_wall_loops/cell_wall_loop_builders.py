"""Build void-cell wall loops as the faces of the labelled wall graph (SG-014).

The skin is already split at every web foot (``split_skin_at_junctions``), so
the section is a planar graph whose **nodes** are the skin-arc junctions and
whose **edges** are the labelled strips: skin arcs around the ring, webs as
chords across it. The cells are that graph's bounded faces, found by one walk:
follow the ring, and turn onto a web wherever one starts.

Nothing here reads a coordinate axis. The loop builder this replaces sorted webs
on ``x`` and chose each web's top on ``y`` -- a chord-first convention applied to
canonical B midlines, where ``X_B`` is flapwise. It paired the wrong web ends
and emitted overlapping cells on most multi-web sections (SG-014). The only
geometric test left is that a web end **coincides** with a skin junction, which
holds to ~1e-13 m by construction (SG-002: walls emit at their true extent), so
a mismatch is a defect to report, never a tolerance to widen.

Each loop carries ``wall_labels``: the exact strips that bound it. That is the
wall-to-cell membership, and consumers read it rather than re-deriving it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from blade_structure_utilities.geometry_utilities.section_contract.cell_wall_loop import CellWallLoop
from blade_structure_utilities.geometry_utilities.section_contract.shell_midline_strip import ShellMidlineStrip

from ..midlines.skin_split import split_skin_at_junctions
from .cell_wall_geom import (
    _normalise_skin_arc,
    _normalise_skin_side,
    _skin_ring_order_key,
    drop_closing_vertex,
    shoelace_area,
)

#: A web end and the skin junction it was split at are the same point (SG-002);
#: measured <= 5e-13 m on 1 m and 0.2 m chords. 1 nm matches SG-007's coincide
#: tolerance and sits ~6 orders below the closest pair of junctions seen (4 mm).
WEB_FOOT_COINCIDENCE_TOL_M = 1e-9

#: Skin arcs from one split share bit-identical endpoints; a legacy
#: unlabelled skin may be chained with the historical stitch gap.
_SKIN_CHAIN_GAP_TOL_M = 1e-4

#: Σ cell area must equal the skin-ring area -- exact for a planar subdivision.
_AREA_IDENTITY_RTOL = 1e-9


@dataclass(frozen=True)
class _Arc:
    label: str
    points: NDArray[np.float64]
    """Oriented along the ring; ``points[0]`` is node ``index``."""
    skin_arc: str | None


@dataclass(frozen=True)
class _Web:
    label: str
    points: NDArray[np.float64]
    """Oriented from node ``start`` to node ``end``."""
    start: int
    end: int


def _is_hierarchical(skins: Sequence[ShellMidlineStrip]) -> bool:
    return bool(skins) and all(
        _normalise_skin_side(s.skin_side) is not None and _normalise_skin_arc(s.skin_arc) is not None
        for s in skins
    )


def _chain_ring(skins: Sequence[ShellMidlineStrip]) -> list[_Arc]:
    """Orient and order the skin strips head-to-tail around the ring."""
    items = [
        (s, drop_closing_vertex(np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=np.float64)))
        for s in skins
    ]
    items = [(s, p) for s, p in items if p.shape[0] >= 2]
    if not items:
        raise ValueError("no skin strip has a usable midline; there is no ring to bound a cell")
    if len(items) == 1:
        strip, points = items[0]
        closed = np.vstack([points, points[:1]])
        return [_Arc(label=str(strip.label), points=closed, skin_arc=_normalise_skin_arc(strip.skin_arc))]

    if _is_hierarchical([s for s, _ in items]):
        ordered = sorted(items, key=lambda item: _skin_ring_order_key(item[0]))
    else:
        ordered = [items.pop(0)]
        while items:
            tail = ordered[-1][1][-1]
            index = min(
                range(len(items)),
                key=lambda i: min(np.hypot(*(items[i][1][0] - tail)), np.hypot(*(items[i][1][-1] - tail))),
            )
            ordered.append(items.pop(index))

    arcs: list[_Arc] = []
    tail: NDArray[np.float64] | None = None
    for position, (strip, points) in enumerate(ordered):
        if tail is None:
            # Orient the first arc so that it meets the second one head-to-tail.
            following = ordered[1][1]
            gap_forward = min(np.hypot(*(following[0] - points[-1])), np.hypot(*(following[-1] - points[-1])))
            gap_reverse = min(np.hypot(*(following[0] - points[0])), np.hypot(*(following[-1] - points[0])))
            if gap_reverse < gap_forward:
                points = points[::-1].copy()
        else:
            if np.hypot(*(points[-1] - tail)) < np.hypot(*(points[0] - tail)):
                points = points[::-1].copy()
            gap = float(np.hypot(*(points[0] - tail)))
            if gap > _SKIN_CHAIN_GAP_TOL_M:
                raise ValueError(
                    f"skin strip {strip.label!r} starts {gap:.3e} m from the end of the previous "
                    "arc; the skin strips do not chain into one ring"
                )
        arcs.append(_Arc(label=str(strip.label), points=points, skin_arc=_normalise_skin_arc(strip.skin_arc)))
        tail = points[-1]
    closure = float(np.hypot(*(arcs[0].points[0] - arcs[-1].points[-1])))
    if closure > _SKIN_CHAIN_GAP_TOL_M:
        raise ValueError(f"the skin ring does not close (gap {closure:.3e} m)")
    return arcs


def _node_of(point: NDArray[np.float64], nodes: NDArray[np.float64], *, web: str) -> int:
    distances = np.hypot(*(nodes - point).T)
    index = int(np.argmin(distances))
    if float(distances[index]) > WEB_FOOT_COINCIDENCE_TOL_M:
        raise ValueError(
            f"web {web!r} ends {float(distances[index]):.3e} m from the nearest skin junction; "
            "a web must end on the skin mid-surface where the skin was split (SG-002). "
            "Emit feet via ShearWeb._extend_to_skin_mid_surface / solve_web_foot_on_skin_mid; "
            "prepare strips with build_shell_strip_topology (index-only — do not relocate feet); "
            "and do not pass frozen web_y that is stale relative to the current skin under twist."
        )
    return index


def _attach_webs(webs: Sequence[ShellMidlineStrip], nodes: NDArray[np.float64]) -> dict[int, _Web]:
    """``{node: web leaving that node}`` -- each web registered from both ends."""
    by_node: dict[int, _Web] = {}
    for strip in webs:
        points = drop_closing_vertex(np.asarray(strip.midline_twisted_pitch_axis_blade_frame, dtype=np.float64))
        if points.shape[0] < 2:
            raise ValueError(f"web {strip.label!r} midline has fewer than two points")
        start = _node_of(points[0], nodes, web=str(strip.label))
        end = _node_of(points[-1], nodes, web=str(strip.label))
        if start == end:
            raise ValueError(f"web {strip.label!r} starts and ends at the same skin junction")
        for node, web in (
            (start, _Web(label=str(strip.label), points=points, start=start, end=end)),
            (end, _Web(label=str(strip.label), points=points[::-1].copy(), start=end, end=start)),
        ):
            if node in by_node:
                raise ValueError(
                    f"webs {by_node[node].label!r} and {strip.label!r} meet the skin at the same junction"
                )
            by_node[node] = web
    return by_node


def _segments_meet(a0, a1, b0, b1) -> bool:
    """Closed segments ``a0a1`` and ``b0b1`` share a point (touching counts).

    Two webs never share an end -- each skin junction carries one web -- so any
    contact between them, even at a vertex, is a crossing.
    """

    def orient(p, q, r) -> float:
        return float((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))

    def within(p, q, r) -> bool:
        return min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and min(p[1], q[1]) <= r[1] <= max(p[1], q[1])

    d1, d2 = orient(b0, b1, a0), orient(b0, b1, a1)
    d3, d4 = orient(a0, a1, b0), orient(a0, a1, b1)
    if ((d1 > 0.0 and d2 < 0.0) or (d1 < 0.0 and d2 > 0.0)) and ((d3 > 0.0 and d4 < 0.0) or (d3 < 0.0 and d4 > 0.0)):
        return True
    return (
        (d1 == 0.0 and within(b0, b1, a0))
        or (d2 == 0.0 and within(b0, b1, a1))
        or (d3 == 0.0 and within(a0, a1, b0))
        or (d4 == 0.0 and within(a0, a1, b1))
    )


def _check_webs_do_not_cross(webs: Sequence[_Web]) -> None:
    for i, first in enumerate(webs):
        for second in webs[i + 1 :]:
            for a0, a1 in zip(first.points[:-1], first.points[1:]):
                for b0, b1 in zip(second.points[:-1], second.points[1:]):
                    if _segments_meet(a0, a1, b0, b1):
                        raise ValueError(f"webs {first.label!r} and {second.label!r} cross")


def _faces(arcs: Sequence[_Arc], web_at: dict[int, _Web]) -> list[list[tuple[str, str, NDArray[np.float64]]]]:
    """Bounded faces: walk the ring, turning onto a web wherever one starts.

    Every arc is used once (in ring direction) and every web twice (once each
    way), so each face is an interior cell and the exterior never appears.
    """
    count = len(arcs)
    used = [False] * count
    faces: list[list[tuple[str, str, NDArray[np.float64]]]] = []
    for first in range(count):
        if used[first]:
            continue
        face: list[tuple[str, str, NDArray[np.float64]]] = []
        arc = first
        for _ in range(2 * count + 2 * len(web_at) + 2):
            if used[arc]:
                raise ValueError("the wall graph walk revisited a skin arc; the cells are not a planar subdivision")
            used[arc] = True
            face.append(("skin", arcs[arc].label, arcs[arc].points))
            node = (arc + 1) % count
            web = web_at.get(node)
            if web is not None:
                face.append(("web", web.label, web.points))
                node = web.end
            if node == first:
                break
            arc = node
        else:
            raise ValueError("the wall graph walk did not close")
        faces.append(face)
    return faces


def _polyline(face: Sequence[tuple[str, str, NDArray[np.float64]]]) -> NDArray[np.float64]:
    points = [face[0][2]]
    for _kind, _label, piece in face[1:]:
        points.append(piece[1:])
    closed = np.vstack(points)
    closed[-1] = closed[0]
    return closed


def _order_le_to_te(
    faces: list[list[tuple[str, str, NDArray[np.float64]]]],
    arcs: Sequence[_Arc],
) -> list[list[tuple[str, str, NDArray[np.float64]]]]:
    """LE cell first, then across one shared web at a time to the TE."""
    if len(faces) == 1:
        return faces
    leading = {arc.label for arc in arcs if arc.skin_arc == "leading_edge"}
    starts = [i for i, face in enumerate(faces) if any(label in leading for kind, label, _ in face if kind == "skin")]
    if len(starts) != 1:
        raise ValueError(
            f"expected exactly one cell holding the leading-edge skin arcs, found {len(starts)}"
        )
    webs_of = [{label for kind, label, _ in face if kind == "web"} for face in faces]
    order = [starts[0]]
    while len(order) < len(faces):
        current = order[-1]
        following = [
            i for i in range(len(faces)) if i not in order and webs_of[i] & webs_of[current]
        ]
        if len(following) != 1:
            raise ValueError(
                f"cell after {current} is ambiguous ({len(following)} candidates); the webs do not "
                "divide the section into a single LE-to-TE chain of cells"
            )
        order.append(following[0])
    return [faces[i] for i in order]


def build_cell_wall_loops(
    strips: Sequence[ShellMidlineStrip],
) -> tuple[CellWallLoop, ...]:
    """Build the void-cell wall loops of a station from its midline strips.

    Parameters
    ----------
    strips
        Skin / web / cap midlines; caps bound no cell and are ignored.
        Topology-prepared strips (``build_shell_strip_topology``) are expected;
        an unsplit skin is split at the web feet here. Webs must end on the
        skin mid-surface.

    Returns
    -------
    tuple[CellWallLoop, ...]
        LE → TE ordered loops, ``n_webs + 1`` of them (one with no webs), each
        carrying the exact ``wall_labels`` that bound it.

    Raises
    ------
    ValueError
        When the walls do not form a planar subdivision of the skin ring: a web
        end off the skin, two webs on one junction, crossing webs, a skin that
        does not chain, or cell areas that do not tile the ring.
    """
    skins = [s for s in strips if s.kind == "skin"]
    webs = [s for s in strips if s.kind == "web"]
    if webs and not _is_hierarchical(skins):
        prepared = split_skin_at_junctions(tuple(s for s in strips if s.kind in ("skin", "web")))
        skins = [s for s in prepared if s.kind == "skin"]
        webs = [s for s in prepared if s.kind == "web"]

    arcs = _chain_ring(skins)
    nodes = np.asarray([arc.points[0] for arc in arcs], dtype=np.float64)
    web_at = _attach_webs(webs, nodes)
    _check_webs_do_not_cross(list({web.label: web for web in web_at.values()}.values()))

    faces = _order_le_to_te(_faces(arcs, web_at), arcs)
    if len(faces) != len(webs) + 1:
        raise ValueError(f"{len(webs)} web(s) produced {len(faces)} cells; expected {len(webs) + 1}")

    ring = np.vstack([arcs[0].points] + [arc.points[1:] for arc in arcs[1:]])
    ring_area = shoelace_area(drop_closing_vertex(ring))
    web_rank = {}
    for index in range(len(faces) - 1):
        shared = {l for k, l, _ in faces[index] if k == "web"} & {l for k, l, _ in faces[index + 1] if k == "web"}
        for label in shared:
            web_rank[label] = index

    loops: list[CellWallLoop] = []
    signed_total = 0.0
    for cell_id, face in enumerate(faces):
        polyline = _polyline(face)
        signed = shoelace_area(drop_closing_vertex(polyline))
        if signed == 0.0 or np.sign(signed) != np.sign(ring_area):
            raise ValueError(f"cell {cell_id} has zero or reversed area; the walls fold over themselves")
        signed_total += signed
        pieces: list[str] = []
        for kind, label, _ in face:
            tag = "skin" if kind == "skin" else label
            if not (tag == "skin" and pieces and pieces[-1] == "skin"):
                pieces.append(tag)
        loops.append(
            CellWallLoop(
                id=cell_id,
                polyline_twisted_pitch_axis_blade_frame=polyline,
                area_m2=float(abs(signed)),
                web_labels=tuple(sorted((l for k, l, _ in face if k == "web"), key=lambda l: web_rank[l])),
                piece_labels=tuple(pieces),
                wall_labels=tuple(label for _kind, label, _ in face),
            )
        )
    if not np.isclose(signed_total, ring_area, rtol=_AREA_IDENTITY_RTOL, atol=0.0):
        raise ValueError(
            f"cell areas sum to {abs(signed_total):.9e} m^2 but the skin ring encloses "
            f"{abs(ring_area):.9e} m^2; the cells do not tile the section"
        )
    return tuple(loops)


__all__ = ["WEB_FOOT_COINCIDENCE_TOL_M", "build_cell_wall_loops"]
