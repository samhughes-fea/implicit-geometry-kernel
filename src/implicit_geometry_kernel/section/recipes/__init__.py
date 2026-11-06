"""Dual-emission wall recipes.

Authoring surface is ``SkinWall`` / ``WebWall`` / ``CapWall`` / ``CoreWall`` /
``InsertWall``. SDF IR classes (``OuterSkin``, ``SparCap``, …) live beside
``walls.py`` and are imported from those modules, not this barrel.
"""

from implicit_geometry_kernel.section.recipes.walls import CapWall, CoreWall, InsertWall, SkinWall, WebWall

__all__ = [
    "CapWall",
    "CoreWall",
    "InsertWall",
    "SkinWall",
    "WebWall",
]
