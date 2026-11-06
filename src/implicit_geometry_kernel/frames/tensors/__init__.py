"""Pure S <-> B tensor and vector maps: second moments (inertia), moments, forces and
thickness. The package keeps every name the former ``frames/tensors.py`` module had."""

from implicit_geometry_kernel.frames.tensors.inertia import (
    _GEOM_I_FLOOR,
    structural_product_moment,
    transform_second_moments_B_to_S,
    transform_second_moments_S_to_B,
)
from implicit_geometry_kernel.frames.tensors.moments import (
    transform_moments_B_to_S,
    transform_moments_S_to_B,
)
from implicit_geometry_kernel.frames.tensors.thickness import (
    transform_thickness_B_to_Sprime,
    transform_thickness_Sprime_to_B,
)
from implicit_geometry_kernel.frames.tensors.forces import (
    transform_forces_B_to_S,
    transform_forces_S_to_B,
)

__all__ = [
    "_GEOM_I_FLOOR",
    "structural_product_moment",
    "transform_second_moments_B_to_S",
    "transform_second_moments_S_to_B",
    "transform_moments_B_to_S",
    "transform_moments_S_to_B",
    "transform_thickness_B_to_Sprime",
    "transform_thickness_Sprime_to_B",
    "transform_forces_B_to_S",
    "transform_forces_S_to_B",
]
