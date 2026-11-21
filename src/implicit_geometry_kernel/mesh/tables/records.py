"""Mesh table records: the vocabularies, the edge and corner dtypes, and the decision and seam dataclasses.

Part of :mod:`implicit_geometry_kernel.mesh.tables.tables`, split out of it in 0.9.4.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np



ORIGINS = ("station", "span_mid", "centroid", "station_as_midpoint", "interior")

CONSTRUCTORS = ("ruled", "midspan_split", "step_patch", "transition", "fan", "bowtie",
                "omitted", "inherited", "interior")

KINDS = ("station", "spanwise", "interior")

MEETS = ("skin", "skin_web", "rib", "tip", "open", "other")

EDGE_DTYPE = np.dtype([
    ("a", np.int64), ("b", np.int64), ("use", np.int64),
    ("e0", np.int64), ("e1", np.int64), ("e2", np.int64),
    ("length", np.float64), ("dx", np.float64), ("dy", np.float64), ("dz", np.float64),
    ("dihedral_deg", np.float64), ("chain_yaw_deg", np.float64), ("chain_pitch_deg", np.float64),
    ("kind", "U16"), ("meet", "U16"), ("sagitta", np.float64),
])

CORNER_DTYPE = np.dtype([
    ("element_id", np.int64), ("corner", np.int64), ("node_id", np.int64),
    ("x", np.float64), ("y", np.float64), ("z", np.float64),
    ("angle_deg", np.float64), ("yaw_deg", np.float64), ("pitch_deg", np.float64), ("roll_deg", np.float64),
    ("dx_in", np.float64), ("dy_in", np.float64), ("dz_in", np.float64),
    ("dx_out", np.float64), ("dy_out", np.float64), ("dz_out", np.float64),
    ("prev_id", np.int64), ("next_id", np.int64),
    ("constructor", "U32"), ("bay", np.int64), ("region", "U64"), ("surface", "U16"),
])


def _edge_users(row) -> list[int]:
    return [int(row[name]) for name in ("e0", "e1", "e2") if int(row[name]) >= 0]


@dataclass(frozen=True)
class StationDecision:
    """Why a span station is in the mesh, or why a grade refused to add one."""

    z: float
    origin: str
    refused: bool
    detail: str


@dataclass(frozen=True)
class StripDecision:
    """What one lofted strip did in one bay, including which strip made the bay an event."""

    bay: int
    z0: float
    z1: float
    region: str
    n0: int
    n1: int
    constructor: str
    trigger: str
    adv_a: int
    adv_b: int
    ties: int
    note: str
    chord_ratio: float = 0.0
    """``|c1 - c0| / max(c0, c1)`` of the bay. Zero when the chord does not change."""
    site0: int = -1
    """Coarse-edge index of the first 2:1 site. ``-1`` when the strip has no site."""
    site1: int = -1
    """Coarse-edge index of the second 2:1 site. ``-1`` when the strip has no site."""


@dataclass(frozen=True)
class MeshEdge:
    """One undirected pair of node ids."""

    a: int
    b: int
    use: int
    elements: tuple[int, ...]
    length: float
    dx: float
    dy: float
    dz: float
    """Edge vector from node ``a`` to node ``b``, in metres. ``a`` is the lower node id."""
    dihedral_deg: float
    """Signed fold from one face normal to the other about the edge vector, degrees.

    The pair is the two faces on the same surface when the edge meets more than
    two, otherwise the first two faces. The sign runs from the lower element id
    to the higher, about the vector from node ``a`` to node ``b``. Zero on a
    boundary edge and on a coplanar pair. The included angle where two skins
    meet is a steady value of this fold. A crinkle is ``chain_pitch_deg``.
    """
    chain_yaw_deg: float
    """Turn in the face plane from the previous spanwise edge of the same regions, degrees.

    Zero when the edge is not spanwise or is the first edge of its chain.
    """
    chain_pitch_deg: float
    """Turn out of the face plane from the previous spanwise edge, degrees.

    A trailing-edge crinkle is a jump in this angle from one spanwise edge to
    the next. Zero on the first edge of a chain.
    """
    kind: str
    meet: str
    sagitta: float


@dataclass(frozen=True)
class OpenLoop:
    """A chain of edges used by one element. The root ring is the expected one."""

    z: float
    n_edges: int
    is_root: bool
    nodes: tuple[int, ...]


@dataclass(frozen=True)
class Corner:
    """One interior angle of one element, named by the node it sits on.

    A node shared by four quads has four angles, so the angle is not a node
    attribute. ``x``, ``y``, ``z`` are that node's coordinates, copied onto the
    row so the corner can be plotted without a join.
    """

    element_id: int
    corner: int
    node_id: int
    x: float
    y: float
    z: float
    angle_deg: float
    """Interior angle in the plane of the two edges, degrees. This is the in-face corner."""
    yaw_deg: float
    """Turn of the outgoing edge in the element plane, degrees. Positive is a left turn about the normal."""
    pitch_deg: float
    """Turn of the outgoing edge out of the element plane, degrees. A crinkle lifts this off zero."""
    roll_deg: float
    """Twist of the surface normal about the outgoing edge, degrees, onto the face that shares that edge."""
    dx_in: float
    dy_in: float
    dz_in: float
    """Incoming edge vector, previous node to this node, metres."""
    dx_out: float
    dy_out: float
    dz_out: float
    """Outgoing edge vector, this node to the next, metres."""
    prev_id: int
    next_id: int
    constructor: str
    bay: int
    region: str
    surface: str


@dataclass(frozen=True)
class Seam:
    """Two node ids that did not merge, or a node lying in the interior of an edge."""

    kind: str
    a: int
    b: int
    node: int
    distance: float
