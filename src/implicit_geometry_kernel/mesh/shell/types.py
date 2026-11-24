"""The shell mesh record types: :class:`ShellMesh` and :class:`MeshRecipe`.

Kept apart from the meshers so that each mesher (``shell_mesh``,
``partition_mesh``) imports the types without importing the other.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from implicit_geometry_kernel.mesh.quality import QUALITY_DTYPE, MeshRunQuality
from implicit_geometry_kernel.mesh.tables import CORNER_DTYPE, EDGE_DTYPE


@dataclass(frozen=True)
class MeshRecipe:
    """The controls that built a shell mesh, and the edge they resolved.

    ``tip_h_over_c`` is the fraction in force (the default or a caller's
    ``h_over_c``). It is ``None`` when the mesh was sized by ``h`` alone or by
    a station count. ``edge`` is that target length at each station, and empty
    in count mode. ``h_over_c`` records only a fraction the caller passed.
    """

    tip_h_over_c: float | None
    h: float | None
    h_over_c: float | None
    n_span: int | None
    region_h: tuple[tuple[str, float], ...]
    refine: tuple[tuple[float, float, float], ...]
    edge: tuple[float, ...]


@dataclass(frozen=True)
class ShellMesh:
    """Nodes and elements of the blade shell, with per-element shell data.

    ``elements[i]`` lists node indices (3 for a rib or tip triangle, 4 for a quad)
    in loop order. The per-element tuples share that index. Normals (right-hand
    rule on that order) point out of the section on the skin, towards the
    trailing edge on webs, and along +z on ribs and the tip closure; a cap
    stack's mid-plane lies ``eccentricity`` inside the skin, along -normal.
    """

    nodes: np.ndarray
    elements: tuple[tuple[int, ...], ...]
    name: tuple[str, ...]
    """Mesh region: a skin region (``cap_su``, ``bay_pr_0``, …), a web, a rib cell or the tip closure."""
    surface: tuple[str, ...]
    """``skin``, ``web``, ``rib`` or ``tip_closure``."""
    partition: tuple[str, ...]
    """The cut's partition the element lies on (``spar_cap_suction_surface``, a rib name, …)."""
    thickness: np.ndarray
    eccentricity: np.ndarray
    """Cap eccentricity of the partition (``ShellRun.eccentricity_m``, m): half the cap stack; 0 off the caps."""
    layers: tuple[tuple[str, ...], ...]
    layer_thicknesses: tuple[tuple[float, ...], ...]
    stations: tuple[float, ...] = field(default=())
    quality: tuple[MeshRunQuality, ...] = field(default=())
    """Per-partition rollup of :attr:`element_quality`."""
    recipe: MeshRecipe | None = None
    """The size controls and the resolved edge at each station."""
    element_quality: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=QUALITY_DTYPE))
    """One quality row per element: aspect, skew, angles, Jacobian, area, shell volume, volume ratio, warp."""
    n_source: int = 0
    """Skin and web nodes occupy ``nodes[:n_source]``. Interior rib and tip nodes follow."""
    bay: tuple[int, ...] = ()
    """Span bay index of each element. Ribs and the tip use the station index."""
    z0: tuple[float, ...] = ()
    z1: tuple[float, ...] = ()
    constructor: tuple[str, ...] = ()
    """How the element was built: ``ruled``, ``midspan_split``, ``step_patch``, ``transition``, ``fan``, ``interior``."""
    n0: tuple[int, ...] = ()
    """Node count of the region on the root side of the bay. Zero on a rib or the tip."""
    n1: tuple[int, ...] = ()
    flipped: tuple[int, ...] = ()
    """1 when orientation reversed the element's loop."""
    node_origin: tuple[str, ...] = ()
    """Per node: ``station``, ``span_mid``, ``centroid``, ``station_as_midpoint``, or ``interior``."""
    section_distance: np.ndarray = field(default_factory=lambda: np.empty(0))
    """Distance of each node from the section at its span. Station samples are zero."""
    valence: tuple[int, ...] = ()
    edges: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=EDGE_DTYPE))
    """One :data:`~implicit_geometry_kernel.mesh.tables.EDGE_DTYPE` row per undirected pair."""
    loops: tuple = ()
    counts: tuple = ()
    station_log: tuple = ()
    strips: tuple = ()
    near_miss: tuple = ()
    junctions: tuple = ()
    valence_defects: tuple[int, ...] = ()
    edge_sagitta: dict = field(default_factory=dict)
    corners: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=CORNER_DTYPE))
    """One :data:`~implicit_geometry_kernel.mesh.tables.CORNER_DTYPE` row per element corner. Not recomputed on export."""
    diagonal_min: tuple[float, ...] = ()
    """Shortest distance between non-adjacent nodes of each element, metres. A collapsed fan is a tiny diagonal."""

    @property
    def quads(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) == 4)

    @property
    def triangles(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) == 3)

    @property
    def polygons(self) -> tuple[int, ...]:
        return tuple(i for i, e in enumerate(self.elements) if len(e) > 4)

    def edge_use(self) -> dict[tuple[int, int], int]:
        """How many elements use each edge (2 on a smooth sheet, more at a junction)."""
        out: dict[tuple[int, int], int] = {}
        for e in self.elements:
            for a, b in zip(e, e[1:] + e[:1]):
                if a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                out[key] = out.get(key, 0) + 1
        return out

    def digest(self) -> str:
        """SHA-256 of the mesh content: nodes, connectivity, surfaces, partitions, thickness, eccentricity.

        Coordinates and properties are rounded to 12 decimals (metres) before
        hashing, so noise in the last bits changes the digest only when a value
        sits on a rounding boundary. Two meshes with the same digest are the
        same mesh for a solver.
        """
        import hashlib

        h = hashlib.sha256()

        def floats(a) -> None:
            arr = np.round(np.asarray(a, dtype=np.float64), 12) + 0.0  # +0.0 folds -0.0
            h.update(np.ascontiguousarray(arr, dtype="<f8").tobytes())

        def strings(values) -> None:
            h.update("\x1f".join(str(v) for v in values).encode("utf-8"))
            h.update(b"\x1e")

        floats(self.nodes)
        h.update(np.asarray([len(e) for e in self.elements], dtype="<i8").tobytes())
        h.update(np.asarray([i for e in self.elements for i in e], dtype="<i8").tobytes())
        strings(self.surface)
        strings(self.partition)
        floats(self.thickness)
        floats(self.eccentricity)
        return h.hexdigest()

    def provenance(self) -> dict:
        """What built this mesh: kernel version, mesher, size controls, stations and :meth:`digest`.

        A plain dict of strings, numbers and lists, ready for JSON or a VTK
        field. ``mesher`` is ``"partition"`` when the (t, z) partition mesher
        built the skin and web quads, else ``"ruled"``.
        """
        from implicit_geometry_kernel import __version__

        sheet = {c for c, s in zip(self.constructor, self.surface) if s in ("skin", "web")}
        recipe = self.recipe
        return {
            "kernel_version": __version__,
            "mesher": "partition" if sheet == {"partition"} else "ruled",
            "h": None if recipe is None else recipe.h,
            "h_over_c": None if recipe is None else recipe.h_over_c,
            "n_span": None if recipe is None else recipe.n_span,
            "stations": [float(z) for z in self.stations],
            "n_nodes": int(len(self.nodes)),
            "n_elements": int(len(self.elements)),
            "digest": self.digest(),
        }


__all__ = ["MeshRecipe", "ShellMesh"]
