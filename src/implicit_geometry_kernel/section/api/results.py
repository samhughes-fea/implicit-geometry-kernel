"""Section geometry worker result dataclasses."""

from __future__ import annotations

from collections.abc import Collection, Mapping
from dataclasses import dataclass
from typing import Any, Literal

from implicit_geometry_kernel.contract.cell_wall_loop import CellWallLoop
from implicit_geometry_kernel.contract.shell_midline_strip import ShellMidlineStrip
from implicit_geometry_kernel.integration.quadtree.tree import QuadTreeMesh
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
from implicit_geometry_kernel.section.layout.system_layout import SystemLayoutSpec
from implicit_geometry_kernel.section.midlines.host import SectionMidlineHost

LayoutInput = SystemLayoutSpec | SectionLayout

OutputGate = Literal["properties", "midlines", "cell_wall_loops"]
ALL_OUTPUT_GATES: frozenset[str] = frozenset(
    {"properties", "midlines", "cell_wall_loops"}
)


def normalize_output_gates(outputs: Collection[str]) -> frozenset[str]:
    """Validate ``outputs`` gate names.

    Unknown names raise ``ValueError``.
    """
    raw = {str(o) for o in outputs}
    unknown = raw - ALL_OUTPUT_GATES
    if unknown:
        raise ValueError(f"Unknown section geometry outputs: {sorted(unknown)}")
    return frozenset(raw)


@dataclass(frozen=True)
class DroppedWeb:
    """Web dropped while resolving spanwise anchors at a station."""

    x_m: float
    reason: str


@dataclass(frozen=True)
class SectionPropertiesMeta:
    """Frozen bag of property-backend metadata from :func:`compute_section_properties`."""

    backend: str
    grid: dict[str, int] | None = None
    quadtree_spec: dict[str, Any] | None = None
    constructive_labels: tuple[str, ...] | None = None
    wall_solid_labels: tuple[str, ...] | None = None
    residual_labels: tuple[str, ...] | None = None
    residual_maturity: dict[str, str] | None = None
    constructive_s: float | None = None
    polygon_count: int | None = None
    residual_quadtree_s: float | None = None
    quadtree_stats: dict[str, Any] | None = None
    constructive_polygons: dict[str, str] | None = None
    profiling: dict[str, Any] | None = None

    @classmethod
    def from_mapping(cls, meta: Mapping[str, Any]) -> SectionPropertiesMeta:
        def _str_tuple(key: str) -> tuple[str, ...] | None:
            raw = meta.get(key)
            if raw is None:
                return None
            return tuple(str(x) for x in raw)

        grid_raw = meta.get("grid")
        grid: dict[str, int] | None = None
        if isinstance(grid_raw, Mapping):
            grid = {str(k): int(v) for k, v in grid_raw.items()}

        maturity_raw = meta.get("residual_maturity")
        maturity: dict[str, str] | None = None
        if isinstance(maturity_raw, Mapping):
            maturity = {str(k): str(v) for k, v in maturity_raw.items()}

        polys_raw = meta.get("constructive_polygons")
        polys: dict[str, str] | None = None
        if isinstance(polys_raw, Mapping):
            polys = {str(k): str(v) for k, v in polys_raw.items()}

        qt_stats = meta.get("quadtree_stats")
        profiling = meta.get("profiling")
        return cls(
            backend=str(meta.get("backend", "")),
            grid=grid,
            quadtree_spec=(
                dict(meta["quadtree_spec"])
                if isinstance(meta.get("quadtree_spec"), Mapping)
                else None
            ),
            constructive_labels=_str_tuple("constructive_labels"),
            wall_solid_labels=_str_tuple("wall_solid_labels"),
            residual_labels=_str_tuple("residual_labels"),
            residual_maturity=maturity,
            constructive_s=(
                float(meta["constructive_s"])
                if meta.get("constructive_s") is not None
                else None
            ),
            polygon_count=(
                int(meta["polygon_count"])
                if meta.get("polygon_count") is not None
                else None
            ),
            residual_quadtree_s=(
                float(meta["residual_quadtree_s"])
                if meta.get("residual_quadtree_s") is not None
                else None
            ),
            quadtree_stats=dict(qt_stats) if isinstance(qt_stats, Mapping) else None,
            constructive_polygons=polys,
            profiling=dict(profiling) if isinstance(profiling, Mapping) else None,
        )


@dataclass(frozen=True)
class SectionGeometryResult:
    """Stable worker output; gated carriers are ``None`` when not requested.

    When midlines are gated, also fills shell-handoff scalars (``chord_m``,
    ``twist_rad``, ``leading_edge_twisted_pitch_axis_blade_frame``, ``trailing_edge_twisted_pitch_axis_blade_frame``) so FE can build
    ``ShellStripMeshingInputs`` without re-evaluating SDF.

    ``section`` implements :class:`SectionMidlineHost` for midline / mesh export.

    ``properties_map`` (when gated) carries area / first / second moments in the
    **solids' authoring frame as built** — with host twist that is LE-origin
    rotated-NACA, **not** frame **B** (`twisted_pitch_axis_blade_frame`). Do not
    silently re-twist the map to match emit midlines.

    ``quadtree_mesh`` is **G-PROPS debug only** (populated when ``retain_mesh``
    is true). It is not a shell-mesh product; density / Q4 / +x3 stay in
    ``section_mesh``.

    ``stage_timings`` is a small ``{stage: seconds}`` dict populated when
    ``SECTION_GEOMETRY_PROFILE=1`` (assemble / assemble_* / properties /
    properties_* / emit_* / midlines / junction_preparation /
    junction_preparation_* / cell_wall_loops).
    """

    section: SectionMidlineHost
    prep_key: tuple[Any, ...]
    properties_map: dict[str, dict[str, float]] | None = None
    properties_meta: SectionPropertiesMeta | None = None
    midline_strips: tuple[ShellMidlineStrip, ...] | None = None
    cell_wall_loops: tuple[CellWallLoop, ...] | None = None
    quadtree_mesh: QuadTreeMesh | None = None
    chord_m: float | None = None
    twist_rad: float | None = None
    leading_edge_twisted_pitch_axis_blade_frame: tuple[float, float] | None = None
    trailing_edge_twisted_pitch_axis_blade_frame: tuple[float, float] | None = None
    stage_timings: dict[str, float] | None = None


__all__ = [
    "ALL_OUTPUT_GATES",
    "DroppedWeb",
    "LayoutInput",
    "OutputGate",
    "SectionGeometryResult",
    "SectionPropertiesMeta",
    "normalize_output_gates",
]
