"""Moved to :mod:`implicit_geometry_kernel.sdf.kernels.polygon_numba`; this path is an alias until 0.6.0."""

import sys

import implicit_geometry_kernel.sdf.kernels.polygon_numba as _moved

sys.modules[__name__] = _moved
