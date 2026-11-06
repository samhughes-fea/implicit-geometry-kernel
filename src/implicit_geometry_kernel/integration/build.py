"""Build 2-D area meshes from scalar field samplers."""

from __future__ import annotations

from typing import Mapping

from implicit_geometry_kernel.integration.quadtree.marching import (
    contour_from_quadtree_marching_squares,
    summarize_quadtree,
)
from implicit_geometry_kernel.integration.quadtree.soa import build_shared_quadtree
from implicit_geometry_kernel.integration.spec import QuadtreeSpec
from implicit_geometry_kernel.integration.types import AreaMesh2D, FieldMesh2D, PhiSampler


def build_area_mesh_2d(phi_sampler: PhiSampler, spec: QuadtreeSpec) -> FieldMesh2D:
    """Refine a quadtree on ``phi_sampler`` and optionally extract the φ=0 contour."""
    mesh = build_shared_quadtree(spec, {"phi": phi_sampler})
    stats = summarize_quadtree(mesh)
    contour = None
    try:
        contour = contour_from_quadtree_marching_squares(mesh.leaf_cells())
    except ValueError:
        pass
    return FieldMesh2D(mesh=mesh, contour=contour, stats=stats)


def build_area_meshes_2d(
    fields: Mapping[str, PhiSampler],
    spec: QuadtreeSpec,
) -> AreaMesh2D:
    """Build one quadtree mesh per named scalar field on the same domain spec."""
    built = {label: build_area_mesh_2d(sampler, spec) for label, sampler in fields.items()}
    return AreaMesh2D(fields=built)
