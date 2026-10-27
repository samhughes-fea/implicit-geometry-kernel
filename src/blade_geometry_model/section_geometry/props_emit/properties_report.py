"""Section properties report container and tabular outputs."""

from __future__ import annotations

import csv
import json

import numpy as np

from blade_geometry_model.section_geometry.sdf.cache import SectionEvalCache
from blade_geometry_model.section_geometry.props_emit.props_backend import (
    PropsBackend,
    props_backend_from_env,
)
from blade_geometry_model.section_geometry.props_emit.props_compute import compute_section_properties


class SectionPropertiesReport:
    """Compute and collect section properties for all subcomponents.

    Parameters
    ----------
    section_geometry : MultiCellSection or SectionMidlineHost
    grid : SDFGrid, optional
        Required when ``backend='grid'``; optional for quadtree (used for spec bbox fallback).
    backend : {'grid', 'quadtree', 'hybrid'}
        Property integration backend (default from ``SECTION_GEOMETRY_PROPS_BACKEND`` env;
        default ``hybrid``).
    eval_cache : SectionEvalCache, optional
        Memoize grid evaluations when ``backend='grid'``.
    quadtree_spec : QuadtreeSpec, optional
        Override adaptive mesh settings for quadtree backend.
    chord, twist_rad, airfoil_for_spec
        Used to size the quadtree domain when ``backend='quadtree'``.

    Usage
    -----
        report = SectionPropertiesReport(bsg, grid, backend='quadtree')
        df     = report.to_dataframe()
        report.to_json("section_props.json")
    """

    def __init__(
        self,
        section_geometry,
        grid=None,
        *,
        backend: PropsBackend | None = None,
        eval_cache: SectionEvalCache | None = None,
        quadtree_spec=None,
        chord: float = 1.0,
        twist_rad: float = 0.0,
        airfoil_for_spec=None,
        _precomputed: tuple | None = None,
    ):
        self._bsg = section_geometry
        self._grid = grid
        self._eval_cache = eval_cache
        self._backend = backend or props_backend_from_env()
        self._meta: dict = {}
        self._data = {}
        self._mesh = None
        if _precomputed is not None:
            self._data, self._meta, self._mesh = _precomputed
        else:
            self._compute(
                quadtree_spec=quadtree_spec,
                chord=float(chord),
                twist_rad=float(twist_rad),
                airfoil_for_spec=airfoil_for_spec,
            )

    @classmethod
    def from_precomputed(cls, section_geometry, data, meta=None, mesh=None):
        """Wrap already-computed properties / mesh (skip rebuild in export)."""
        return cls(
            section_geometry,
            _precomputed=(dict(data), dict(meta or {}), mesh),
        )

    def _compute(self, *, quadtree_spec=None, chord, twist_rad, airfoil_for_spec):
        self._data, self._meta, self._mesh = compute_section_properties(
            self._bsg,
            backend=self._backend,
            grid=self._grid,
            eval_cache=self._eval_cache,
            quadtree_spec=quadtree_spec,
            chord=chord,
            twist_rad=twist_rad,
            airfoil_for_spec=airfoil_for_spec,
            with_mesh=True,
        )

    @property
    def backend_meta(self) -> dict:
        """Backend metadata (quadtree stats, grid size, profiling counters)."""
        return dict(self._meta)

    @property
    def mesh(self):
        """Shared quadtree mesh (``QuadTreeMesh``) for the quadtree backend, else None.

        Exposed so the plot renderer can reuse the already-built mesh
        (``mesh.phi_by_label``) instead of evaluating a second dense grid. Kept off
        ``backend_meta`` / ``to_json`` so it is never serialised.
        """
        return self._mesh

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------

    def __getitem__(self, label):
        return self._data[label]

    def __iter__(self):
        return iter(self._data)

    def summary(self):
        """Return a formatted string summary."""
        lines = ["=" * 64, f"{'Section Properties':^64}", "=" * 64]
        for label, props in self._data.items():
            lines.append(f"\n  [{label}]")
            lines.append(f"    Area   : {props['area']:.6g}")
            lines.append(f"    Centroid: ({props['cx']:.6g}, {props['cy']:.6g})")
            lines.append(f"    Ixx    : {props['Ixx']:.6g}")
            lines.append(f"    Iyy    : {props['Iyy']:.6g}")
            lines.append(f"    Ixy    : {props['Ixy']:.6g}")
        lines.append("\n" + "=" * 64)
        return "\n".join(lines)

    def to_dict(self):
        """Return all properties as a plain dict."""
        return dict(self._data)

    def to_json(self, filepath, *, job_meta=None):
        """Serialise to JSON file.

        Parameters
        ----------
        filepath : str or Path
            Output path.
        job_meta : dict, optional
            Merged under top-level key ``job_meta`` (orchestration provenance).
        """
        # Convert numpy types to Python natives
        def _convert(obj):
            if isinstance(obj, dict):
                return {str(k): _convert(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple)):
                return [_convert(v) for v in obj]
            if isinstance(obj, (np.floating, float)):
                return float(obj)
            if isinstance(obj, (np.integer, int)):
                return int(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj

        data = {
            label: {k: _convert(v) for k, v in props.items()}
            for label, props in self._data.items()
        }
        if job_meta is not None:
            data["job_meta"] = _convert(job_meta)
        with open(filepath, "w") as fh:
            json.dump(data, fh, indent=2)
        return filepath

    def to_csv(self, filepath):
        """Write properties to a flat CSV file."""
        fieldnames = ["label", "area", "cx", "cy", "Ixx", "Iyy", "Ixy",
                      "r_gyr_x", "r_gyr_y"]
        with open(filepath, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                    extrasaction="ignore")
            writer.writeheader()
            for props in self._data.values():
                row = {k: (float(props[k]) if k != "label" else props[k])
                       for k in fieldnames}
                writer.writerow(row)
        return filepath

    def to_dataframe(self):
        """Return a pandas DataFrame (requires pandas)."""
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("pandas is required for to_dataframe(). "
                              "Install with: pip install pandas")
        rows = []
        for label, props in self._data.items():
            row = {"label": label}
            row.update({k: float(v) for k, v in props.items()
                        if k != "label"})
            rows.append(row)
        return pd.DataFrame(rows).set_index("label")


__all__ = ["SectionPropertiesReport"]
