"""Node counts per region that follow region length and stay constant between events.

A zone is a run of consecutive stations with the same topology signature. In a
zone every region keeps one node count, so its strips are ruled quads; the
count follows the region's longest length in the zone at the target element
size ``h``. Counts change only across an event band, where the loft takes an
all-quad transition. Counts are odd (an even number of edges per region), as
the transition needs.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Mapping, Sequence

import numpy as np

from implicit_geometry_kernel.blade.topology.template import ShellTopologyTemplate


def _length(pts: np.ndarray) -> float:
    pts = np.asarray(pts, float)
    if len(pts) < 2:
        return 0.0
    return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))


def _odd_nodes(length: float, h: float, minimum: int = 3) -> int:
    edges = max(int(np.ceil(length / h - 1e-9)), minimum - 1)
    edges += edges % 2
    return edges + 1


def zones_of(templates: Sequence[ShellTopologyTemplate]) -> list[list[int]]:
    """Station indices grouped into runs of equal signature, in span order."""
    zones: list[list[int]] = []
    for i, tmpl in enumerate(templates):
        if zones and templates[zones[-1][-1]].signature == tmpl.signature:
            zones[-1].append(i)
        else:
            zones.append([i])
    return zones


def zone_sized_templates(templates: Sequence[ShellTopologyTemplate],
                         samples: Sequence[Mapping[str, np.ndarray]], h: float,
                         *, minimum: int = 3) -> list[ShellTopologyTemplate]:
    """``templates`` with each branch's node count set from its zone's longest length.

    ``samples`` are the station polylines sampled with ``templates`` (any
    count: only their lengths are used). ``h`` is the target element edge.
    """
    if h <= 0.0:
        raise ValueError("element size h must be positive")
    if len(templates) != len(samples):
        raise ValueError("one sample set per template")
    out = list(templates)
    for zone in zones_of(templates):
        longest: dict[str, float] = {}
        for i in zone:
            for name, pts in samples[i].items():
                longest[name] = max(longest.get(name, 0.0), _length(pts))
        for i in zone:
            tmpl = templates[i]
            loop = tuple(
                replace(b, n=_odd_nodes(longest.get(b.name, 0.0), h, minimum), sized=True)
                for b in tmpl.loop
            )
            webs = tuple(
                replace(w, n=_odd_nodes(longest.get(w.name, 0.0), h, minimum), sized=True)
                for w in tmpl.webs
            )
            out[i] = replace(tmpl, loop=loop, webs=webs)
    return out


__all__ = ["zone_sized_templates", "zones_of"]
