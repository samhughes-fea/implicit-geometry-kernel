"""Shell frames C / E / M + product load channels.

See ``README.md`` for mathematics and Exact I/O. Prefer the package barrel
``coordinate_frames`` for public imports.
"""

from .load_channels import ShellSolverChannels, blade_resultants_to_shell_channels
from .ply import (
    ply_strain_transform_3d_E_to_M,
    ply_strain_transform_3d_M_to_E,
    ply_strain_transform_E_to_M,
    ply_strain_transform_M_to_E,
    ply_stress_transform_3d_E_to_M,
    ply_stress_transform_3d_M_to_E,
    ply_stress_transform_E_to_M,
    ply_stress_transform_M_to_E,
    rotate_ply_strain_3d_E_to_M,
    rotate_ply_strain_3d_M_to_E,
    rotate_ply_strain_E_to_M,
    rotate_ply_strain_M_to_E,
    rotate_ply_stress_3d_E_to_M,
    rotate_ply_stress_3d_M_to_E,
    rotate_ply_stress_E_to_M,
    rotate_ply_stress_M_to_E,
    rotate_ply_vectors,
    rotate_transverse_shear_E_to_M,
)

__all__ = [
    "ShellSolverChannels",
    "blade_resultants_to_shell_channels",
    "ply_strain_transform_3d_E_to_M",
    "ply_strain_transform_3d_M_to_E",
    "ply_strain_transform_E_to_M",
    "ply_strain_transform_M_to_E",
    "ply_stress_transform_3d_E_to_M",
    "ply_stress_transform_3d_M_to_E",
    "ply_stress_transform_E_to_M",
    "ply_stress_transform_M_to_E",
    "rotate_ply_strain_3d_E_to_M",
    "rotate_ply_strain_3d_M_to_E",
    "rotate_ply_strain_E_to_M",
    "rotate_ply_strain_M_to_E",
    "rotate_ply_stress_3d_E_to_M",
    "rotate_ply_stress_3d_M_to_E",
    "rotate_ply_stress_E_to_M",
    "rotate_ply_stress_M_to_E",
    "rotate_ply_vectors",
    "rotate_transverse_shear_E_to_M",
]
