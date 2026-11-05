"""Private SDF IR projection internals for wall recipes (G-RECIPE).

Former ``sections.components``. Not a public authoring surface — assemble
authors frozen wall recipes; these classes are ``to_sdf_component`` targets.
"""

from implicit_geometry_kernel.section.recipes.sdf_ir.base import SectionComponent, expr_callable, twist_component
from implicit_geometry_kernel.section.recipes.sdf_ir.inserts import LEInsert, TEInsert
from implicit_geometry_kernel.section.recipes.sdf_ir.outer_skin import OuterSkin
from implicit_geometry_kernel.section.recipes.sdf_ir.sandwich_core import SandwichCore
from implicit_geometry_kernel.section.recipes.sdf_ir.shear_web import ShearWeb
from implicit_geometry_kernel.section.recipes.sdf_ir.spar_cap import SparCap

__all__ = [
    "SectionComponent",
    "expr_callable",
    "twist_component",
    "OuterSkin",
    "SparCap",
    "ShearWeb",
    "SandwichCore",
    "TEInsert",
    "LEInsert",
]
