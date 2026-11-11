"""Node counts per region from its length at each station, at a target element size ``h``.

Every region's element size stays between ``h / 2`` and ``h``. Walking the
stations from root to tip, a region keeps its count while that holds and
takes the count from its own length when it does not (it has grown past ``h``
per element, or shrunk below ``h / 2``). Counts that follow length station by
station change in small steps, each taken by the loft's all-quad transition,
so a long tapering zone does not end in one band that drops from a fine root
count to a coarse one (bug log BG-034). Counts are odd (an even number of
edges per region), as that transition needs.
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


def _even_edges(length: float, h: float, minimum: int = 3) -> int:
    edges = max(int(np.ceil(length / h - 1e-9)), minimum - 1)
    return edges + edges % 2


def _odd_nodes(length: float, h: float, minimum: int = 3) -> int:
    return _even_edges(length, h, minimum) + 1


def _held_edges(length: float, previous: int | None, h: float, minimum: int) -> int:
    """``previous`` while ``length / previous`` is in ``[h / 2, h]``, else the count from ``length``."""
    edges = _even_edges(length, h, minimum)
    if previous is not None and previous >= edges and length >= 0.5 * h * previous:
        return previous
    return edges


def region_counts(templates: Sequence[ShellTopologyTemplate],
                  samples: Sequence[Mapping[str, np.ndarray]], h: float,
                  *, minimum: int = 3,
                  region_h: Mapping[str, float] | None = None) -> list[ShellTopologyTemplate]:
    """``templates`` (in span order) with each region's node count set from its length.

    ``samples`` are the station polylines sampled with ``templates`` (any
    count: only their lengths are used). ``h`` is the target element edge;
    ``region_h`` overrides it for named regions (``{"le_su": 0.02, "le_pr":
    0.02}`` refines the leading edge). Every region has at least ``minimum``
    nodes.
    """
    sizes = {str(k): float(v) for k, v in (region_h or {}).items()}
    if h <= 0.0 or any(v <= 0.0 for v in sizes.values()):
        raise ValueError("element size h must be positive")
    if len(templates) != len(samples):
        raise ValueError("one sample set per template")
    out: list[ShellTopologyTemplate] = []
    held: dict[str, int] = {}
    for tmpl, smp in zip(templates, samples):
        held = {name: _held_edges(_length(pts), held.get(name), sizes.get(name, h), minimum)
                for name, pts in smp.items()}
        least = _even_edges(0.0, h, minimum)
        loop = tuple(replace(b, n=held.get(b.name, least) + 1, sized=True) for b in tmpl.loop)
        webs = tuple(replace(w, n=held.get(w.name, least) + 1, sized=True) for w in tmpl.webs)
        out.append(replace(tmpl, loop=loop, webs=webs))
    return out


__all__ = ["region_counts"]
