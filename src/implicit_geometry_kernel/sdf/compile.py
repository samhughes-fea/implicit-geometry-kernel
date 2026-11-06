"""Former compile step from CSG shapes to evaluator nodes (kept for old imports).

The shapes in :mod:`implicit_geometry_kernel.sdf.nodes` are ``SDFNode`` subclasses and
evaluate themselves, so there is nothing to compile. These names remain until 0.6.0.
"""

from __future__ import annotations

from typing import Dict

import numpy as np

from implicit_geometry_kernel.sdf.nodes.csg import eval_expr_batch
from implicit_geometry_kernel.sdf.nodes.primitives import Expr


def compile_expr_to_node(expr: Expr) -> Expr:
    """Return ``expr`` itself: every shape is already an ``SDFNode``."""
    return expr


def eval_via_nodes(expr: Expr, x, y, memo: Dict[Expr, np.ndarray] | None = None) -> np.ndarray:
    """Same as :func:`implicit_geometry_kernel.sdf.nodes.eval_expr_batch`."""
    return eval_expr_batch(expr, x, y, memo=memo)


__all__ = ["compile_expr_to_node", "eval_via_nodes"]
