"""Hard invariants of a shell mesh: the pass/fail gate under the quality report.

:func:`~implicit_geometry_kernel.mesh.quality.mesh_check` reports counts and
extrema for an engineer to judge. The checks here are different: each one
holds for every mesh the kernel means to build, so a single violation is a bug
in a mesher, never a matter of judgement.

* every per-element field has one entry per element;
* every element has at least three distinct node indices, all in range;
* nodes, thickness and eccentricity are finite, thickness is positive;
* every ``surface`` is one of :data:`~implicit_geometry_kernel.mesh.shell.SURFACES`;
* every node is used, and no two nodes coincide (nodes are shared, not duplicated);
* no two elements have the same node set, and none has zero area;
* on each surface kind, two elements that share an edge run it in opposite
  directions (consistent normals);
* an edge used by one element only lies on the root station, or on a web
  (a web that ends without a rib has a free end).

:func:`check_shell_mesh` returns the problems as messages;
:func:`assert_shell_mesh` raises :class:`ShellMeshInvariantError` on any.
``shell_mesh(..., strict=True)`` runs the gate on the mesh it returns.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from implicit_geometry_kernel.mesh.shell.types import ShellMesh

#: Surface kinds a shell mesh may carry (same as ``mesh.shell.SURFACES``).
_SURFACES = ("skin", "web", "rib", "tip_closure")

#: Per-element fields that must have one entry per element. ``constructor`` is
#: checked only when it is set.
_PER_ELEMENT = ("name", "surface", "partition", "thickness", "eccentricity", "layers", "layer_thicknesses")


class ShellMeshInvariantError(ValueError):
    """A shell mesh broke one or more hard invariants; ``problems`` lists them."""

    def __init__(self, problems: list[str]):
        self.problems = list(problems)
        head = f"shell mesh breaks {len(problems)} invariant(s)"
        super().__init__(head + ":\n  " + "\n  ".join(problems[:20]) + ("\n  ..." if len(problems) > 20 else ""))


def _scale(nodes: np.ndarray) -> float:
    if len(nodes) == 0:
        return 1.0
    span = float(np.linalg.norm(nodes.max(axis=0) - nodes.min(axis=0)))
    return span if span > 0.0 else 1.0


def _area(p: np.ndarray) -> float:
    """Area of a (possibly non-planar) polygon from its Newell normal."""
    c = p.mean(axis=0)
    q = p - c
    return 0.5 * float(np.linalg.norm(np.cross(q, np.roll(q, -1, axis=0)).sum(axis=0)))


def _close_pairs(nodes: np.ndarray, tol: float) -> list[tuple[int, int]]:
    """Node pairs closer than ``tol``: a sweep along x, NumPy only."""
    order = np.argsort(nodes[:, 0], kind="stable")
    xs = nodes[order, 0]
    hi = np.searchsorted(xs, xs + tol, side="right")
    out = []
    for a in (int(v) for v in np.flatnonzero(hi - np.arange(len(xs)) > 1)):
        i = int(order[a])
        for b in range(a + 1, int(hi[a])):
            j = int(order[b])
            if float(np.linalg.norm(nodes[i] - nodes[j])) < tol:
                out.append((min(i, j), max(i, j)))
    return sorted(out)


def check_shell_mesh(mesh: ShellMesh, *, tol: float | None = None, limit: int = 5) -> list[str]:
    """Problems with ``mesh`` as messages; empty when every invariant holds.

    ``tol`` is the node-coincidence distance in metres (default ``1e-9`` of
    the mesh's bounding-box diagonal); zero-area elements are those below
    ``tol**2``. Each check reports at most ``limit`` offending items, then a
    count.
    """
    out: list[str] = []
    nodes = np.asarray(mesh.nodes, dtype=float)
    n_nodes, n_cells = len(nodes), len(mesh.elements)
    scale = _scale(nodes) if nodes.ndim == 2 and len(nodes) else 1.0
    tol = 1e-9 * scale if tol is None else float(tol)

    def report(what: str, items: list) -> None:
        if items:
            shown = ", ".join(str(i) for i in items[:limit])
            more = f" (+{len(items) - limit} more)" if len(items) > limit else ""
            out.append(f"{what}: {shown}{more}")

    if nodes.ndim != 2 or (n_nodes and nodes.shape[1] != 3):
        return [f"nodes must be an (n, 3) array, got shape {nodes.shape}"]

    # Per-element fields.
    for name in _PER_ELEMENT:
        n = len(getattr(mesh, name))
        if n != n_cells:
            out.append(f"{name} has {n} entries for {n_cells} elements")
    if len(mesh.constructor) not in (0, n_cells):
        out.append(f"constructor has {len(mesh.constructor)} entries for {n_cells} elements")
    if out:
        return out  # the checks below index these fields by element

    # Connectivity.
    bad_index, short = [], []
    for i, e in enumerate(mesh.elements):
        if any(k < 0 or k >= n_nodes for k in e):
            bad_index.append(i)
        elif len(set(e)) < 3 or len(set(e)) != len(e):
            short.append(i)
    report("elements with a node index out of range", bad_index)
    report("elements with fewer than three distinct nodes or a repeated node", short)
    if bad_index:
        return out

    # Values.
    if not np.all(np.isfinite(nodes)):
        report("non-finite nodes", [int(i) for i in np.flatnonzero(~np.isfinite(nodes).all(axis=1))])
    thickness = np.asarray(mesh.thickness, dtype=float)
    report("elements with non-finite or non-positive thickness",
           [int(i) for i in np.flatnonzero(~(np.isfinite(thickness) & (thickness > 0.0)))])
    ecc = np.asarray(mesh.eccentricity, dtype=float).reshape(n_cells, -1)
    report("elements with non-finite eccentricity", [int(i) for i in np.flatnonzero(~np.isfinite(ecc).all(axis=1))])
    report("elements with an unknown surface", [i for i, s in enumerate(mesh.surface) if s not in _SURFACES])

    # Nodes: all used, none coincident.
    used = np.zeros(n_nodes, dtype=bool)
    for e in mesh.elements:
        used[list(e)] = True
    report("unused nodes", [int(i) for i in np.flatnonzero(~used)])
    if n_nodes > 1 and np.all(np.isfinite(nodes)):
        pairs = _close_pairs(nodes, tol)
        report(f"coincident node pairs (closer than {tol:.3g} m)", [tuple(int(v) for v in p) for p in pairs])

    # Elements: unique, with area.
    seen: dict[frozenset, int] = {}
    dup = []
    for i, e in enumerate(mesh.elements):
        cell = frozenset(e)
        if cell in seen:
            dup.append((seen[cell], i))
        else:
            seen[cell] = i
    report("duplicate elements", dup)
    if np.all(np.isfinite(nodes)):
        report("zero-area elements", [i for i, e in enumerate(mesh.elements)
                                      if len(set(e)) >= 3 and _area(nodes[list(e)]) <= tol * tol])

    # Edges: orientation per surface kind, and where free edges may lie.
    directed: dict[str, dict[tuple[int, int], list[bool]]] = defaultdict(lambda: defaultdict(list))
    owners: dict[tuple[int, int], list[int]] = defaultdict(list)
    for i, (e, surface) in enumerate(zip(mesh.elements, mesh.surface)):
        for a, b in zip(e, e[1:] + e[:1]):
            key = (a, b) if a < b else (b, a)
            directed[surface][key].append(a < b)
            owners[key].append(i)
    flipped = []
    for surface, edges in directed.items():
        for key, dirs in edges.items():
            if len(dirs) == 2 and dirs[0] == dirs[1]:
                flipped.append((surface, key))
    report("edges run the same way by two elements of one surface (inconsistent normals)", flipped)

    z_root = float(min(mesh.stations)) if len(mesh.stations) else (float(nodes[:, 2].min()) if n_nodes else 0.0)
    free = []
    for key, cells in owners.items():
        if len(cells) != 1:
            continue
        if mesh.surface[cells[0]] == "web":
            continue
        if abs(nodes[key[0], 2] - z_root) > tol or abs(nodes[key[1], 2] - z_root) > tol:
            free.append(key)
    report("free edges away from the root station", free)
    return out


def assert_shell_mesh(mesh: ShellMesh, *, tol: float | None = None) -> ShellMesh:
    """Return ``mesh`` if every invariant holds, else raise :class:`ShellMeshInvariantError`."""
    problems = check_shell_mesh(mesh, tol=tol)
    if problems:
        raise ShellMeshInvariantError(problems)
    return mesh


__all__ = ["ShellMeshInvariantError", "assert_shell_mesh", "check_shell_mesh"]
