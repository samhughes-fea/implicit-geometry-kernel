"""Base class for IR-backed section subcomponents."""

from __future__ import annotations

from typing import Any

import numpy as np

from blade_geometry_model.section_geometry.sdf.ir import Expr, eval_expr_batch


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

    def eval_batch(self, x, y, memo: dict | None = None):
        return eval_expr_batch(self.expr, x, y, memo)


def twist_component(comp: SectionComponent, twist_rad: float) -> SectionComponent:
    """Return a component whose expr is rotated into the B-frame."""
    from blade_geometry_model.section_geometry.sdf.ir import Rotate

    if abs(float(twist_rad)) <= 1e-10:
        return comp
    twisted = SectionComponent.__new__(SectionComponent)
    twisted.expr = Rotate(comp.expr, angle=float(twist_rad), cx=0.0, cy=0.0)
    twisted._expr = twisted.expr
    twisted.label = getattr(comp, "label", None)
    for attr in (
        "surface",
        "x_start",
        "x_end",
        "cap_height",
        "skin_thickness",
        "twist_angle",
        "_af",
        "_strip_clip",
        "_camber_xy",
        "x_top",
        "y_top",
        "x_bot",
        "y_bot",
        "x_top_ext",
        "y_top_ext",
        "x_bot_ext",
        "y_bot_ext",
        "thickness",
        "alignment",
        "inner_sdf",
        "outer_sdf",
        "y_min",
        "y_max",
        "bay_mode",
        "flapwise_nx",
        "flapwise_ny",
        "flapwise_c_lo",
        "flapwise_c_hi",
        "flapwise_clip_mode",
    ):
        if hasattr(comp, attr):
            setattr(twisted, attr, getattr(comp, attr))
    if hasattr(comp, "midline_polyline"):
        twisted.midline_polyline = comp.midline_polyline  # type: ignore[method-assign]
    return twisted
