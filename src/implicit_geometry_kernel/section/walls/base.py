"""Base class for IR-backed section subcomponents."""

from __future__ import annotations

from typing import Any


from implicit_geometry_kernel.sdf.nodes import Expr, eval_expr_batch


def expr_callable(expr: Expr):
    """Standalone callable wrapper for a single IR subtree."""

    def _fn(x, y, memo=None):
        return eval_expr_batch(expr, x, y, memo)

    _fn._expr = expr  # type: ignore[attr-defined]
    return _fn


class SectionComponent:
    """IR-backed SDF component; ``__call__`` evaluates ``self.expr``."""

    label: Any = None

    def __init__(self, expr: Expr):
        self.expr = expr
        self._expr = expr

    def __call__(self, x, y, memo=None):
        return eval_expr_batch(self.expr, x, y, memo)

    def with_expr(self, expr: Expr):
        """Same component type with a replaced IR subtree."""
        out = object.__new__(type(self))
        out.__dict__.update(vars(self))
        out.expr = expr
        out._expr = expr
        return out

    def twisted(self, twist_rad: float):
        """Rotate ``expr`` into the B-frame; preserve the concrete subclass."""
        from implicit_geometry_kernel.sdf.nodes import Rotate

        if abs(float(twist_rad)) <= 1e-10:
            return self
        return self.with_expr(Rotate(self.expr, angle=float(twist_rad), cx=0.0, cy=0.0))


def twist_component(comp: SectionComponent, twist_rad: float) -> SectionComponent:
    """Return a component whose expr is rotated into the B-frame."""
    return comp.twisted(twist_rad)
