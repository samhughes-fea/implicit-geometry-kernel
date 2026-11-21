"""Deprecated flat path: use ``implicit_geometry_kernel.mesh.export.vtk``.

Kept for one 0.9.x release so callers get a warning before the path goes
(removal in 0.9.5, docs/MIGRATING.md). Importing this module does not warn;
reading a name from it does.
"""

from __future__ import annotations

from implicit_geometry_kernel.mesh._flat import forward

__getattr__ = forward(__name__, "implicit_geometry_kernel.mesh.export.vtk")
