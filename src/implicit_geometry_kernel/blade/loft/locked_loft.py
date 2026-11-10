"""Zone-locked loft of ``ShellGraph`` runs: held counts and shared partition vertices.

Wave 2 already refuses to lerp missing names. This layer stops adjacent independent
resamples from inventing a second copy of each partition edge.

Authored region identity, ``map_vertex``, and figure loft live in
:mod:`implicit_geometry_kernel.blade.topology`. ``shared_loft_names`` stays here
as the Wave 2 event helper; ``write_mid_surfaces`` no longer joins on it alone.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Sequence, Tuple

import numpy as np

from implicit_geometry_kernel.blade.shell_graph import ShellGraph, ShellRun

LOFTABLE_SURFACES = frozenset({"skin", "web", "core", "insert"})


def loftable_runs(graph: ShellGraph) -> Tuple[ShellRun, ...]:
    """Skin/web (and opt-in solid) runs in graph order; ribs and stringers omitted."""
    return tuple(r for r in graph.runs() if r.surface in LOFTABLE_SURFACES)


def loft_signature(graph: ShellGraph) -> Tuple[str, ...]:
    """Ordered loftable names. A change is a topology event for loft, not a lerp."""
    return tuple(r.name for r in loftable_runs(graph))


def shell_loft_zones(graphs: Sequence[ShellGraph]) -> list[slice]:
    """Maximal contiguous bands of constant loftable signature.

    Rib/stringer faces do not split a zone. A dropped web, cap runout, or bay
    rename does. Empty input yields no zones.
    """
    graphs = list(graphs)
    if not graphs:
        return []
    zones: list[slice] = []
    start = 0
    sig = loft_signature(graphs[0])
    for i in range(1, len(graphs)):
        nxt = loft_signature(graphs[i])
        if nxt != sig:
            zones.append(slice(start, i))
            start = i
            sig = nxt
    zones.append(slice(start, len(graphs)))
    return zones


def arc_resample_open(pts: np.ndarray, n: int) -> np.ndarray:
    """Open polyline resample that keeps both endpoints (``n >= 2``)."""
    pts = np.asarray(pts, float)
    n = max(int(n), 2)
    if len(pts) == 0:
        return np.zeros((n, 3))
    if pts.ndim != 2:
        raise ValueError("points must be (n, d)")
    if np.allclose(pts[0], pts[-1]) and len(pts) > 2:
        pts = pts[:-1]
    if len(pts) == 1:
        out = np.repeat(pts[:1], n, axis=0)
        if out.shape[1] == 2:
            out = np.column_stack([out, np.zeros(n)])
        return out
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 1e-15:
        out = np.repeat(pts[:1], n, axis=0)
        if out.shape[1] == 2:
            out = np.column_stack([out, np.zeros(n)])
        return out
    t = np.linspace(0.0, s[-1], n)
    cols = pts.shape[1]
    out = np.empty((n, cols))
    for j, tj in enumerate(t):
        k = int(np.searchsorted(s, tj, side="right") - 1)
        k = min(max(k, 0), len(seg) - 1)
        u = (tj - s[k]) / max(seg[k], 1e-15)
        out[j] = pts[k] + u * (pts[k + 1] - pts[k])
    if cols == 2:
        out = np.column_stack([out, np.zeros(n)])
    out[0] = pts[0] if pts.shape[1] == out.shape[1] else np.append(pts[0], 0.0)
    last = pts[-1] if pts.shape[1] == out.shape[1] else np.append(pts[-1], 0.0)
    out[-1] = last
    return out


def _held_counts(runs: Sequence[ShellRun], n_chord: int) -> Tuple[int, ...]:
    """Same held count on every run so a name can loft across a topology event."""
    n = max(int(n_chord), 2)
    return tuple(n for _ in runs)


@dataclass(frozen=True)
class LockedLoftPlan:
    """Held point counts for one constant-signature loft zone.

    Adjacent loop runs share partition-edge vertices at resample time; web
    midline ends are those same loop samples (shared nodes).
    """

    names: Tuple[str, ...]
    counts: Tuple[int, ...]
    loop_names: Tuple[str, ...]
    web_names: Tuple[str, ...]

    def count_of(self, name: str) -> int:
        return self.counts[self.names.index(name)]

    @classmethod
    def from_graph(cls, graph: ShellGraph, n_chord: int = 81) -> "LockedLoftPlan":
        runs = loftable_runs(graph)
        loop_names = tuple(r.name for r in graph.loop if r.surface in LOFTABLE_SURFACES)
        web_names = tuple(r.name for r in graph.webs if r.surface in LOFTABLE_SURFACES)
        solids = tuple(r.name for r in runs if r.name not in loop_names and r.name not in web_names)
        ordered = loop_names + web_names + solids
        by_name = {r.name: r for r in runs}
        held_runs = tuple(by_name[n] for n in ordered if n in by_name)
        return cls(tuple(r.name for r in held_runs), _held_counts(held_runs, n_chord),
                   loop_names, web_names)


def _closest_end(pts: np.ndarray, point: np.ndarray) -> int:
    d0 = float(np.linalg.norm(pts[0] - point))
    d1 = float(np.linalg.norm(pts[-1] - point))
    return 0 if d0 <= d1 else -1


def resample_locked(graph: ShellGraph, plan: LockedLoftPlan) -> Dict[str, np.ndarray]:
    """Resample ``plan`` names on ``graph`` with shared partition vertices.

    ``graph`` must carry every plan name (callers zone first). Raises ``KeyError``
    if a planned run is missing — no label lerp.
    """
    by_name = {r.name: r for r in graph.runs()}
    out: Dict[str, np.ndarray] = {}
    for name, n in zip(plan.names, plan.counts):
        if name not in by_name:
            raise KeyError(name)
        out[name] = arc_resample_open(by_name[name].points, n)
    for edge in graph.partition_edges():
        if edge.name_a not in out or edge.name_b not in out:
            continue
        point = np.asarray(edge.point, float)
        pa, pb = out[edge.name_a], out[edge.name_b]
        ia, ib = _closest_end(pa, point), _closest_end(pb, point)
        pa[ia] = point
        pb[ib] = point
    loop_pts = [out[n] for n in plan.loop_names if n in out]
    if loop_pts:
        loop_xy = np.vstack(loop_pts)
        for name in plan.web_names:
            if name not in out:
                continue
            pts = out[name]
            for i in (0, -1):
                j = int(np.argmin(np.linalg.norm(loop_xy - pts[i], axis=1)))
                pts[i] = loop_xy[j]
    return out


def locked_mid_surface_patches(
    graphs: Sequence[ShellGraph],
    n_chord: int = 81,
    members=None,
) -> list[tuple[str, np.ndarray, str]]:
    """Locked loft patches: ``(name, rings, surface)``.

    Each station is resampled with a plan for *that* graph. A name that survives
    a topology event (web drop, cap runout) is lofted through it. A name present
    on only one side of the event starts or ends there — no lerp. ``rings`` is
    ``(n_span_run, n_chord, 3)`` for one contiguous presence of ``name``.
    """
    want = None if members is None else set(members)
    graphs = list(graphs)
    samples: list[Dict[str, np.ndarray]] = []
    surface_of: Dict[str, str] = {}
    order: list[str] = []
    for graph in graphs:
        plan = LockedLoftPlan.from_graph(graph, n_chord=n_chord)
        samples.append(resample_locked(graph, plan))
        for run in loftable_runs(graph):
            surface_of[run.name] = run.surface
            if run.name not in order:
                order.append(run.name)
    patches: list[tuple[str, np.ndarray, str]] = []
    for name in order:
        if want is not None and name not in want:
            continue
        i = 0
        n = len(samples)
        while i < n:
            if name not in samples[i]:
                i += 1
                continue
            j = i + 1
            while j < n and name in samples[j]:
                j += 1
            if j - i >= 2:
                rings = np.stack([samples[k][name] for k in range(i, j)], axis=0)
                patches.append((name, rings, surface_of.get(name, "skin")))
            i = j
    return patches


__all__ = [
    "LOFTABLE_SURFACES",
    "LockedLoftPlan",
    "arc_resample_open",
    "loft_signature",
    "loftable_runs",
    "locked_mid_surface_patches",
    "resample_locked",
    "shell_loft_zones",
]
