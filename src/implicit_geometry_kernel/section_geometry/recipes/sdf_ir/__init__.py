"""Private SDF IR projection internals for wall recipes (G-RECIPE).

Former ``sections.components``. Not a public authoring surface — assemble
authors frozen wall recipes; these classes are ``to_sdf_component`` targets.
"""

from .base import SectionComponent, expr_callable, twist_component
from .inserts import LEInsert, TEInsert
from .outer_skin import OuterSkin
from .sandwich_core import SandwichCore
from .shear_web import ShearWeb
from .spar_cap import SparCap

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
