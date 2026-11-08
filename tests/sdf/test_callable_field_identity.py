"""Opaque callable leaves must not be confused by the compile cache."""

import numpy as np

from implicit_geometry_kernel.sdf.nodes import CallableField, eval_expr_batch


def test_same_token_different_callable_is_a_different_leaf():
    a = CallableField(key="clip_1", fn=lambda x, y: np.asarray(x) - 1.0)
    b = CallableField(key="clip_1", fn=lambda x, y: np.asarray(x) - 2.0)
    assert a != b
    x = np.array([0.0])
    y = np.array([0.0])
    assert eval_expr_batch(a, x, y)[0] == -1.0
    assert eval_expr_batch(b, x, y)[0] == -2.0


def test_same_token_same_callable_is_shared():
    def fn(x, y):
        return np.asarray(x)

    assert CallableField(key="k", fn=fn) == CallableField(key="k", fn=fn)
    assert hash(CallableField(key="k", fn=fn)) == hash(CallableField(key="k", fn=fn))
