"""Station-local cache for grid-evaluated component SDF fields."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .ir import compile_expr


def _eikonal_redistance_2d(phi, dx, dy, *, n_iter=12, dt=0.15):
    """Few iterations of first-order reinitialization to bring |∇φ| → 1."""
    f = np.asarray(phi, dtype=np.float64, copy=True)
    s0 = np.sign(np.clip(f, a_min=-1e-15, a_max=1e-15))
    s0 = np.where(s0 == 0.0, 1.0, s0)
    for _ in range(int(n_iter)):
        gy = np.gradient(f, float(dy), axis=0)
        gx = np.gradient(f, float(dx), axis=1)
        gmag = np.sqrt(gx * gx + gy * gy) + 1e-20
        f = f - float(dt) * s0 * (gmag - 1.0)
    return f


@dataclass
class SectionEvalCache:
    """Memoize grid IR evaluation outputs by ``(grid_id, label)``."""

    _phi_by_key: dict[tuple[int, str], np.ndarray] = field(default_factory=dict)
    _expr_cache_by_grid: dict[int, dict[Any, np.ndarray]] = field(default_factory=dict)

    def get(self, grid: Any, label: str) -> np.ndarray | None:
        return self._phi_by_key.get((id(grid), str(label)))

    def set(self, grid: Any, label: str, phi: np.ndarray) -> np.ndarray:
        arr = np.asarray(phi, dtype=float)
        self._phi_by_key[(id(grid), str(label))] = arr
        return arr

    def get_or_eval(self, label: str, sdf_callable: Any, grid: Any) -> np.ndarray:
        return self.get_or_eval_with_owner(label, sdf_callable, grid, owner=None)

    def get_or_eval_with_owner(
        self, label: str, sdf_callable: Any, grid: Any, owner: Any | None
    ) -> np.ndarray:
        key = (id(grid), str(label))
        cached = self._phi_by_key.get(key)
        if cached is not None:
            return cached

        from blade_geometry_model.section_geometry.midlines.section_midline_host import (
            resolve_section_host,
        )

        phi = None
        expr = None
        host = resolve_section_host(owner) if owner is not None else None
        if host is not None:
            if str(label) in {"_skin_outer_boundary", "_skin_inner_boundary"}:
                attr = (
                    "_skin_outer_boundary_expr"
                    if str(label) == "_skin_outer_boundary"
                    else "_skin_inner_boundary_expr"
                )
                expr = getattr(host, attr, None)
            else:
                comps = getattr(host, "components", None) or {}
                comp = comps.get(str(label)) if isinstance(comps, dict) else None
                expr = getattr(comp, "expr", getattr(comp, "_expr", None)) if comp is not None else None
        else:
            expr = getattr(sdf_callable, "expr", getattr(sdf_callable, "_expr", None))

        if expr is not None:
            gid = id(grid)
            expr_cache = self._expr_cache_by_grid.setdefault(gid, {})
            phi = np.asarray(compile_expr(expr, grid, cache=expr_cache), dtype=float)

        if phi is None and owner is not None:
            twist = float(getattr(owner, "twist_rad", getattr(host, "twist_rad", 0.0) if host is not None else 0.0) or 0.0)
            if abs(twist) > 1e-10 and hasattr(grid, "rotated_eval"):
                base_components = getattr(owner, "components_unrotated", None)
                if isinstance(base_components, dict) and label in base_components:
                    phi = np.asarray(grid.rotated_eval(base_components[label], twist), dtype=float)
                elif str(label) == "_skin_outer_boundary":
                    base_sdf = getattr(owner, "_skin_outer_boundary_unrotated_sdf", None)
                    if base_sdf is not None:
                        phi = np.asarray(grid.rotated_eval(base_sdf, twist), dtype=float)
                elif str(label) == "_skin_inner_boundary":
                    base_sdf = getattr(owner, "_skin_inner_boundary_unrotated_sdf", None)
                    if base_sdf is not None:
                        phi = np.asarray(grid.rotated_eval(base_sdf, twist), dtype=float)
        if phi is None:
            phi = np.asarray(grid.eval(sdf_callable), dtype=float)

        import os

        reflag = os.getenv("SECTION_GEOMETRY_REDISTANCE", "0").lower() in {"1", "true", "yes"}
        if reflag and phi.ndim == 2 and hasattr(grid, "dx") and hasattr(grid, "dy"):
            try:
                phi = _eikonal_redistance_2d(phi, float(grid.dx), float(grid.dy))
            except Exception:
                pass
        self._phi_by_key[key] = phi
        return phi

    def clear(self) -> None:
        self._phi_by_key.clear()
        self._expr_cache_by_grid.clear()
