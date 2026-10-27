"""Section-level IR evaluator with shared subtree memoization."""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

from blade_geometry_model.section_geometry.midlines.section_midline_host import (
    resolve_section_host,
)
from blade_geometry_model.section_geometry.sdf.ir import CallableField, Expr, eval_expr_batch, union_expr

MultiSampler = Callable[[np.ndarray], dict[str, np.ndarray]]


class SectionFieldEvaluator:
    """Evaluate labelled section components with shared IR memo per point batch."""

    def __init__(
        self,
        expr_by_label: dict[str, Expr],
        *,
        twist_rad: float = 0.0,
    ):
        self._expr_by_label = dict(expr_by_label)
        self._twist = float(twist_rad)
        self._union_expr = (
            union_expr(list(expr_by_label.values())) if expr_by_label else None
        )

    @classmethod
    def from_section(cls, section: Any, twist_rad: float = 0.0) -> SectionFieldEvaluator:
        del twist_rad
        exprs: dict[str, Expr] = {}
        host = resolve_section_host(section)
        components = getattr(host, "components", None) or {}
        for lbl, comp in components.items():
            expr = getattr(comp, "expr", getattr(comp, "_expr", None))
            if expr is None:
                # Skin-only / opaque callables: fall through to CallableField path.
                continue
            key = str(getattr(lbl, "value", lbl))
            exprs[key] = expr
        if not exprs:
            # Fallback when no IR exprs were collected: labelled SDF callables
            # only. Prefer each callable's own IR expr; otherwise wrap as an
            # opaque CallableField so expr_bbox / sampling still work.
            raw_labels = list(getattr(section, "labels", ()) or ()) or list(section)
            for lbl in raw_labels:
                key = str(getattr(lbl, "value", lbl))
                callable_sdf = section[lbl]
                expr = getattr(callable_sdf, "expr", getattr(callable_sdf, "_expr", None))
                exprs[key] = expr if expr is not None else CallableField(key=key, fn=callable_sdf)
        twist = float(getattr(host, "twist_rad", 0.0) or 0.0)
        return cls(exprs, twist_rad=twist)

    @property
    def labels(self) -> list[str]:
        return list(self._expr_by_label.keys())

    def expr_for(self, label: str) -> Expr | None:
        """Return the IR expression for ``label`` (or ``None`` if unknown)."""
        return self._expr_by_label.get(str(label))

    def exprs_for(self, labels: list[str]) -> dict[str, Expr]:
        """Return the IR expressions for ``labels`` that are registered."""
        out: dict[str, Expr] = {}
        for lbl in labels:
            key = str(lbl)
            expr = self._expr_by_label.get(key)
            if expr is not None:
                out[key] = expr
        return out

    def eval_labels(
        self,
        labels: list[str],
        pts: np.ndarray,
        *,
        memo: dict[Expr, np.ndarray] | None = None,
    ) -> dict[str, np.ndarray]:
        """Evaluate multiple labels on the same point batch with one shared memo."""
        arr = np.asarray(pts, dtype=float).reshape(-1, 2)
        if arr.size == 0:
            return {str(lbl): np.empty(0, dtype=float) for lbl in labels}
        x = arr[:, 0]
        y = arr[:, 1]
        store: dict[Expr, np.ndarray] = {} if memo is None else memo
        out: dict[str, np.ndarray] = {}
        for lbl in labels:
            key = str(lbl)
            expr = self._expr_by_label.get(key)
            if expr is None:
                continue
            out[key] = np.asarray(
                eval_expr_batch(expr, x, y, store), dtype=float
            ).reshape(-1)
        return out

    def eval(self, label: str, x, y, memo: dict[Expr, np.ndarray] | None = None) -> np.ndarray:
        expr = self._expr_by_label[str(label)]
        return eval_expr_batch(expr, x, y, memo)

    def sample(self, label: str, pts: np.ndarray) -> np.ndarray:
        arr = np.asarray(pts, dtype=float).reshape(-1, 2)
        if arr.size == 0:
            return np.empty(0, dtype=float)
        memo: dict[Expr, np.ndarray] = {}
        out = self.eval(str(label), arr[:, 0], arr[:, 1], memo=memo).reshape(-1)
        from blade_geometry_model.section_geometry.sdf.profiling import get_stats

        stats = get_stats()
        if stats is not None:
            stats.phi_eval_count += int(out.size)
        return out

    def eval_union(self, x, y, memo: dict[Expr, np.ndarray] | None = None) -> np.ndarray:
        if self._union_expr is None:
            raise ValueError("No component expressions registered.")
        return eval_expr_batch(self._union_expr, x, y, memo)

    def make_multi_sampler(self, labels: list[str] | None = None) -> MultiSampler:
        """Return a callable that evaluates all labels on one point batch."""
        lbls = [str(l) for l in (labels if labels is not None else self.labels)]

        def multi_sampler(pts: np.ndarray) -> dict[str, np.ndarray]:
            arr = np.asarray(pts, dtype=float).reshape(-1, 2)
            if arr.size == 0:
                return {lbl: np.empty(0, dtype=float) for lbl in lbls}
            memo: dict[Expr, np.ndarray] = {}
            result = self.eval_labels(lbls, arr, memo=memo)
            from blade_geometry_model.section_geometry.sdf.profiling import get_stats

            stats = get_stats()
            if stats is not None:
                for vals in result.values():
                    stats.phi_eval_count += int(np.asarray(vals).size)
            return result

        return multi_sampler

    def make_sampler(self, label: str):
        lbl = str(label)

        def sampler(pts: np.ndarray) -> np.ndarray:
            return self.sample(lbl, pts)

        return sampler

    def make_samplers(self, labels: list[str]) -> dict[str, Any]:
        return {str(lbl): self.make_sampler(str(lbl)) for lbl in labels}
