"""The geometry hand-off types that downstream mesh and section models consume."""

from implicit_geometry_kernel.contract.cell_wall_loop import CellWallLoop
from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
from implicit_geometry_kernel.contract.laminate_offset_region import LaminateOffsetRegion
from implicit_geometry_kernel.contract.shell_midline_strip import ShellMidlineStrip, StripBond

__all__ = [
    "CellWallLoop",
    "InfeasibleSectionGeometry",
    "LaminateOffsetRegion",
    "ShellMidlineStrip",
    "StripBond",
]
