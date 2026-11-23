"""Authored sub-partitions: bands of a partition that step properties (decision 0016).

A :class:`SubPartition` is a band of one host partition of the shell graph (a
skin partition, a cap or a web), given as a window in the host's normalised arc
length ``t`` and, optionally, a span window in ``z``. Inside the band the
thickness and eccentricity are the sub-partition's own; outside, the host's.
Spar caps, bond flanges, web feet, wedges, promoted ply drops, leading-edge
protection and core extents are all sub-partitions.

Properties only: a sub-partition never moves a node or builds a second
reference surface. Its eccentricity is a vector ``(normal, tangential)`` in
metres, carried as data for a consumer to apply.

This module sits in ``structure`` and imports nothing above it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike

from implicit_geometry_kernel.structure.laws import SpanLaw, ThicknessField

#: Kinds of sub-partition named in decision 0016. Any other string is allowed.
SUB_PARTITION_KINDS = (
    "spar_cap",
    "bond_flange",
    "web_foot",
    "wedge",
    "ply_drop",
    "le_protection",
    "core",
)


@dataclass(frozen=True)
class SubPartition:
    """A band of one host partition with its own thickness and eccentricity.

    Attributes
    ----------
    name:
        Unique name of the band.
    host:
        Name of the shell-graph run it lies on (a loop partition or a web).
    t_range:
        ``(t0, t1)`` with ``0 <= t0 < t1 <= 1``, in the host's normalised arc
        length (``blade.skeleton.run_t``).
    z_range:
        ``(z0, z1)`` span window in metres, or ``None`` for wherever the host
        exists.
    thickness:
        Metres: a number, a :class:`SpanLaw` of ``z`` or a
        :class:`ThicknessField` of ``(t, z)`` (``t`` in the host's coordinate).
    eccentricity:
        ``(normal, tangential)`` offset of the band's stack mid-plane from the
        host mid-surface, metres. Normal is positive inward (along -normal of
        the skin, as ``ShellMesh.eccentricity``).
    kind:
        Free label, usually one of :data:`SUB_PARTITION_KINDS`.
    layers:
        Layer names, passed through for consumers.
    """

    name: str
    host: str
    t_range: Tuple[float, float]
    thickness: float | SpanLaw | ThicknessField
    z_range: Tuple[float, float] | None = None
    eccentricity: Tuple[float, float] = (0.0, 0.0)
    kind: str = "sub_partition"
    layers: Tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        t0, t1 = (float(v) for v in self.t_range)
        if not 0.0 <= t0 < t1 <= 1.0:
            raise ValueError(f"sub-partition {self.name!r}: t_range must satisfy 0 <= t0 < t1 <= 1, got {self.t_range}")
        object.__setattr__(self, "t_range", (t0, t1))
        if self.z_range is not None:
            z0, z1 = (float(v) for v in self.z_range)
            if not z0 < z1:
                raise ValueError(f"sub-partition {self.name!r}: z_range must satisfy z0 < z1, got {self.z_range}")
            object.__setattr__(self, "z_range", (z0, z1))
        e = tuple(float(v) for v in self.eccentricity)
        if len(e) != 2:
            raise ValueError(f"sub-partition {self.name!r}: eccentricity is (normal, tangential)")
        object.__setattr__(self, "eccentricity", e)
        if isinstance(self.thickness, (int, float)) and float(self.thickness) < 0.0:
            raise ValueError(f"sub-partition {self.name!r}: negative thickness")

    def contains(self, t: ArrayLike, z: ArrayLike) -> np.ndarray:
        """True where ``(t, z)`` lies in the band (edges included)."""
        tt, zz = np.broadcast_arrays(np.asarray(t, float), np.asarray(z, float))
        t0, t1 = self.t_range
        inside = (tt >= t0) & (tt <= t1)
        if self.z_range is not None:
            inside &= (zz >= self.z_range[0]) & (zz <= self.z_range[1])
        return inside

    def overlaps(self, other: "SubPartition") -> bool:
        """True when both bands lie on one host and share more than an edge in t and in z."""
        if self.host != other.host:
            return False
        if min(self.t_range[1], other.t_range[1]) - max(self.t_range[0], other.t_range[0]) <= 0.0:
            return False
        za = self.z_range or (-np.inf, np.inf)
        zb = other.z_range or (-np.inf, np.inf)
        return min(za[1], zb[1]) - max(za[0], zb[0]) > 0.0

    def thickness_at(self, t: ArrayLike, z: ArrayLike) -> np.ndarray:
        """Band thickness in metres at ``(t, z)`` (broadcast)."""
        tt, zz = np.broadcast_arrays(np.asarray(t, float), np.asarray(z, float))
        th = self.thickness
        if isinstance(th, ThicknessField):
            out = np.asarray(th(tt, zz), float)
        elif isinstance(th, SpanLaw):
            out = np.reshape(np.asarray(th(zz.ravel()), float), zz.shape)
        else:
            out = np.full(zz.shape, float(th))
        return np.broadcast_to(out, zz.shape).astype(float)


def check_sub_partitions(subs: Sequence[SubPartition]) -> list[str]:
    """Problems with a set of sub-partitions: duplicate names and overlaps on one host.

    Two bands on the same host overlap when their ``t`` windows overlap by more
    than a shared edge and their span windows overlap. Returns messages; empty
    when the set is consistent.
    """
    out: list[str] = []
    seen: set = set()
    for s in subs:
        if s.name in seen:
            out.append(f"duplicate sub-partition name {s.name!r}")
        seen.add(s.name)
    for i, a in enumerate(subs):
        for b in subs[i + 1:]:
            if a.overlaps(b):
                out.append(f"sub-partitions {a.name!r} and {b.name!r} overlap on {a.host!r}")
    return out


__all__ = ["SUB_PARTITION_KINDS", "SubPartition", "check_sub_partitions"]
