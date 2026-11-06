"""
Rigid-body and affine transforms for SDF callables.

All transforms operate in the pull-back sense: instead of moving the
geometry, we transform the query coordinates into the geometry's local frame.
This preserves the SDF property exactly for rotations and translations.
"""

import numpy as np


def _ensure_f64(a):
    if isinstance(a, np.ndarray) and a.dtype == np.float64:
        return a
    return np.asarray(a, dtype=float)


def rotate_field(f, angle, cx=0.0, cy=0.0):
    """Wrap an SDF callable with a CCW rotation of `angle` radians about (cx, cy)."""
    cos_a = np.cos(-angle)
    sin_a = np.sin(-angle)

    def _rotated(x, y):
        x = _ensure_f64(x)
        y = _ensure_f64(y)
        xr = x - cx
        yr = y - cy
        xp = cos_a * xr - sin_a * yr + cx
        yp = sin_a * xr + cos_a * yr + cy
        return f(xp, yp)

    return _rotated


def translate_field(f, dx, dy):
    """Wrap an SDF callable with a translation by (dx, dy)."""
    def _translated(x, y):
        x = _ensure_f64(x)
        y = _ensure_f64(y)
        return f(x - dx, y - dy)
    return _translated


def scale_field(f, sx, sy=None, cx=0.0, cy=0.0):
    """Wrap an SDF callable with anisotropic scaling about (cx, cy).

    Anisotropic scaling (sx ≠ sy) distorts the SDF — use uniform scaling
    (sx == sy) if exact distances are required.
    """
    if sy is None:
        sy = sx
    sx = float(sx)
    sy = float(sy)
    if (not np.isfinite(sx)) or (not np.isfinite(sy)):
        raise ValueError(f"sx and sy must be finite; got sx={sx!r}, sy={sy!r}.")
    if sx <= 0.0 or sy <= 0.0:
        raise ValueError(
            f"sx and sy must be strictly positive; got sx={sx!r}, sy={sy!r}."
        )

    def _scaled(x, y):
        x = _ensure_f64(x)
        y = _ensure_f64(y)
        xp = (x - cx) / sx + cx
        yp = (y - cy) / sy + cy
        return f(xp, yp) * min(sx, sy)

    return _scaled


def mirror_field_x(f):
    """Mirror the SDF field about the x-axis (y → -y)."""
    def _mirrored(x, y):
        return f(_ensure_f64(x), -_ensure_f64(y))
    return _mirrored


def mirror_field_y(f):
    """Mirror the SDF field about the y-axis (x → -x)."""
    def _mirrored(x, y):
        return f(-_ensure_f64(x), _ensure_f64(y))
    return _mirrored


class SDFFrame:
    """Composable, ordered sequence of 2-D rigid-body transforms.

    Call ``apply(f)`` to wrap an SDF callable with the full chain.
    Transforms are applied to query coordinates in reverse order (pull-back).

    Example
    -------
        frame = SDFFrame().translate(0.1, 0.0).rotate(np.radians(15), cx=0.3)
        web_global = frame.apply(web_local_sdf)
    """

    def __init__(self):
        self._ops = []

    def rotate(self, angle, cx=0.0, cy=0.0):
        """Add a CCW rotation of `angle` radians about (cx, cy)."""
        self._ops.append(("rotate", {"angle": angle, "cx": cx, "cy": cy}))
        return self

    def translate(self, dx, dy):
        self._ops.append(("translate", {"dx": dx, "dy": dy}))
        return self

    def scale(self, sx, sy=None, cx=0.0, cy=0.0):
        self._ops.append(("scale", {"sx": sx, "sy": sy, "cx": cx, "cy": cy}))
        return self

    def mirror_x(self):
        self._ops.append(("mirror_x", {}))
        return self

    def mirror_y(self):
        self._ops.append(("mirror_y", {}))
        return self

    def apply(self, f):
        """Wrap callable f with this transform chain."""
        result = f
        for op_type, kwargs in self._ops:
            if op_type == "rotate":
                result = rotate_field(result, **kwargs)
            elif op_type == "translate":
                result = translate_field(result, **kwargs)
            elif op_type == "scale":
                result = scale_field(result, **kwargs)
            elif op_type == "mirror_x":
                result = mirror_field_x(result)
            elif op_type == "mirror_y":
                result = mirror_field_y(result)
        return result

    def inverse(self):
        """Return a new SDFFrame representing the inverse transform chain."""
        inv = SDFFrame()
        for op_type, kwargs in self._ops:
            if op_type == "rotate":
                inv._ops.insert(0, ("rotate", {
                    "angle": -kwargs["angle"], "cx": kwargs["cx"], "cy": kwargs["cy"],
                }))
            elif op_type == "translate":
                inv._ops.insert(0, ("translate", {
                    "dx": -kwargs["dx"], "dy": -kwargs["dy"],
                }))
            elif op_type == "scale":
                sx = kwargs["sx"]
                sy = kwargs.get("sy") or sx
                inv._ops.insert(0, ("scale", {
                    "sx": 1.0 / sx, "sy": 1.0 / sy,
                    "cx": kwargs.get("cx", 0.0), "cy": kwargs.get("cy", 0.0),
                }))
            elif op_type in ("mirror_x", "mirror_y"):
                inv._ops.insert(0, (op_type, {}))
        return inv

    def compose(self, other):
        """Return a new SDFFrame that applies self then other."""
        new_frame = SDFFrame()
        new_frame._ops = self._ops + other._ops
        return new_frame

    def __repr__(self):
        ops_str = ", ".join(
            f"{t}({', '.join(f'{k}={v:.4g}' if isinstance(v, float) else f'{k}={v}' for k, v in kw.items())})"
            for t, kw in self._ops
        )
        return f"SDFFrame([{ops_str}])"


def rotate_points(x, y, angle, cx=0.0, cy=0.0):
    """Rotate (x, y) point arrays CCW by `angle` radians about (cx, cy)."""
    x = _ensure_f64(x)
    y = _ensure_f64(y)
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    xr = cos_a * (x - cx) - sin_a * (y - cy) + cx
    yr = sin_a * (x - cx) + cos_a * (y - cy) + cy
    return xr, yr


def flapwise_aligned_web_angle(twist_angle):
    """Return the in-plane rotation to align a web with the flapwise axis.

    A chord-normal web is vertical at zero twist; to stay flapwise-aligned
    in a section rotated by ``twist_angle``, counter-rotate by ``-twist_angle``.
    """
    return -twist_angle
