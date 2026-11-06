"""Computational choices for building a section. Does not change geometry names."""

from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from implicit_geometry_kernel.section.api.results import normalize_output_gates


def _env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    return default if raw is None else int(raw)


def _env_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    return default if raw is None else float(raw)


@dataclass(frozen=True)
class SectionOptions:
    """Outputs, props backend, quadtree settings and sample counts.

    ``SECTION_GEOMETRY_PROPS_BACKEND``, ``SECTION_GEOMETRY_QT_DEPTH``,
    ``SECTION_GEOMETRY_QT_MAX_CELLS`` and ``SECTION_GEOMETRY_QT_BAND_SCALE``
    are the defaults when the corresponding field is left unset. Profiling
    still reads ``SECTION_GEOMETRY_PROFILE``; that variable does not change
    results.
    """

    outputs: tuple[str, ...] = ("properties", "midlines", "cell_wall_loops")
    props_backend: Any | None = None
    n_web_samples: int = 20
    n_cap_samples: int = 80
    retain_mesh: bool = False
    quadtree_spec: Any | None = None
    airfoil_for_spec: Any | None = None
    qt_depth: int | None = None
    qt_max_cells: int | None = None
    qt_band_scale: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "outputs", tuple(normalize_output_gates(self.outputs)))
        object.__setattr__(self, "n_web_samples", int(self.n_web_samples))
        object.__setattr__(self, "n_cap_samples", int(self.n_cap_samples))
        object.__setattr__(self, "retain_mesh", bool(self.retain_mesh))
        if self.qt_depth is None:
            object.__setattr__(self, "qt_depth", _env_int("SECTION_GEOMETRY_QT_DEPTH", 7))
        if self.qt_max_cells is None:
            object.__setattr__(self, "qt_max_cells", _env_int("SECTION_GEOMETRY_QT_MAX_CELLS", 150000))
        if self.qt_band_scale is None:
            object.__setattr__(self, "qt_band_scale", _env_float("SECTION_GEOMETRY_QT_BAND_SCALE", 1.5))

    @classmethod
    def from_kwargs(
        cls,
        *,
        outputs: Iterable[str] = ("properties", "midlines", "cell_wall_loops"),
        props_backend: Any | None = None,
        n_web_samples: int = 20,
        n_cap_samples: int = 80,
        retain_mesh: bool = False,
        quadtree_spec: Any | None = None,
        airfoil_for_spec: Any | None = None,
        **_ignored: Any,
    ) -> SectionOptions:
        return cls(
            outputs=tuple(outputs),
            props_backend=props_backend,
            n_web_samples=n_web_samples,
            n_cap_samples=n_cap_samples,
            retain_mesh=retain_mesh,
            quadtree_spec=quadtree_spec,
            airfoil_for_spec=airfoil_for_spec,
        )
