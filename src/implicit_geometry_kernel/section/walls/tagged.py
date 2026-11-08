"""Stamp ``FeatureTable`` ids onto section wall solids (recipes and IR)."""

from __future__ import annotations

from implicit_geometry_kernel.sdf.nodes.base import FeatureTable, SDFNode
from implicit_geometry_kernel.sdf.nodes.tagged import CallableSDF, ExprSDF, TagNode
from implicit_geometry_kernel.section.walls.base import SectionComponent

def kind_from_label(label: str) -> str:
    from implicit_geometry_kernel.section.assembly.members import wall_type_for_label

    wt = wall_type_for_label(label)
    if wt is not None:
        return wt.kind
    return "skin"


def tagged_component(comp: SectionComponent, table: FeatureTable | None = None) -> tuple[SDFNode, FeatureTable]:
    """Wrap a recipe component as a tagged ``SDFNode`` (``phi`` from its IR ``expr``)."""
    table = FeatureTable() if table is None else table
    label = str(getattr(comp, "label", None) or type(comp).__name__)
    spec = table.add(label, kind_from_label(label))
    expr = getattr(comp, "expr", None)
    node: SDFNode = ExprSDF(expr, spec.id) if expr is not None else CallableSDF(comp, spec.id)
    return TagNode(node, spec.id), table
