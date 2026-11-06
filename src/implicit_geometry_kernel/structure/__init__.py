"""Members of a thin-walled structure and how they combine (geometry only).

A structure is an outer mould line along a span with a wall and internal members.
Builders (:class:`MemberBuilder`) turn parameters into four kinds of contribution
(wall, branch, region, transverse), and :class:`StructureAssembler` combines them
into one :class:`Structure` of signed-distance fields. Thicknesses are span laws in
metres; plies, materials and meshes belong to downstream packages.

This layer sits between ``sdf`` and ``section``: it imports only ``sdf``.
"""

from implicit_geometry_kernel.structure.assembler import (
    BranchField,
    Structure,
    StructureAssembler,
    StructureError,
    WallThickness,
    assemble,
)
from implicit_geometry_kernel.structure.builders import BuildContext, BuilderRegistry, MemberBuilder
from implicit_geometry_kernel.structure.contributions import (
    BranchContribution,
    Contribution,
    RegionContribution,
    TransverseContribution,
    WallContribution,
)
from implicit_geometry_kernel.structure.laws import (
    ConstantLaw,
    ProductLaw,
    ScaledLaw,
    SpanLaw,
    SumLaw,
    TableLaw,
    ThicknessLaw,
    as_law,
)

__all__ = [
    "BranchContribution", "BranchField", "BuildContext", "BuilderRegistry", "ConstantLaw",
    "Contribution", "MemberBuilder", "ProductLaw", "RegionContribution", "ScaledLaw", "SpanLaw",
    "Structure", "StructureAssembler", "StructureError", "SumLaw", "TableLaw", "ThicknessLaw",
    "TransverseContribution", "WallContribution", "WallThickness", "as_law", "assemble",
]
