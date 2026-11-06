"""implicit_geometry_kernel - geometry of blade and wing cross-sections.

The kernel describes the geometry of blade and wing sections and lofts it along the
span. It holds nothing about materials, plies, stiffness, loads or solvers; the
laminate and section models build those on top of it.

Public API
----------
Airfoils   Airfoil (from coordinates, ``Airfoil.naca4``, ``Airfoil.ellipse``)
Layout     Station (span position, chord, twist, pitch axis, offset),
           SkinRegion (chordwise skin bands), SpanWeb (webs at chord fractions),
           SYSTEM_TYPE_KEYS and resolve_system_type (the 40 structural layouts)
Loft       BladeGeometry (a layout lofted between stations; optional span tables for
           chord, twist, pitch axis, offset), BladeGeometry.interpolate(z) (the
           interpolator), section_at(z) -> SectionCut, oml_polyline(z) (closed OML),
           oml_sdf() (exact polygon SDF at each z). Scalars interpolate linearly
           (PCHIP later, scalars only); airfoils stay 2-D linear blends.
Feature    FieldSample.meta ids (FeatureTable, TagNode, PriorityUnionNode),
           GeometrySolver (in-plane project / intersect), wall_meta / wall_mid_sdf,
           shell_graph / partition_edges / station_network (merged shell: caps are partitions).

Subpackages
-----------
blade         airfoils, stations, lofting and cut sections (the API above)
section       one 2-D section: the 40 layouts, placement, midlines, junctions,
              cell-wall loops, area and second moments
sdf           2-D signed-distance modelling: primitives, CSG, offsets, feature meta,
              lofted polygon SDF, in-plane Newton solver
frames        geometric coordinate frames and point transforms
integration   quadtree and shoelace area and moment integrals
contract      the hand-off types (ShellMidlineStrip, CellWallLoop)
extract       in-plane zero loops and feature curves (not a mesh)

Dependencies: numpy, shapely and numba; scipy (``[splines]``) for spline-resampled
airfoils. Matplotlib (``[plots]``) is for the sibling ``plotting`` package, not
the kernel. Nothing else in the suite.
"""

from implicit_geometry_kernel.blade import (
    Airfoil,
    BladeGeometry,
    CellLoop,
    CellPolygon,
    MidSurface,
    SectionCut,
    SectionWall,
    ShearWeb,
    PartitionEdge,
    ShellGraph,
    ShellRun,
    TopologyDiff,
    TopologyEvent,
    SkinRegion,
    SpanWeb,
    Station,
    StationNetwork,
    partition_edge_curve,
    section_at,
    shared_loft_names,
    shell_graph,
    shell_strips_from_graph,
    topology_diff,
    station_network,
    station_network_from_meta,
    wall_meta_from_cut,
)
from implicit_geometry_kernel.sdf import (
    FeatureTable,
    FieldSample,
    GeometrySolver,
    LoftedPolygonSDF,
    NO_FEATURE,
    PriorityUnionNode,
    TagNode,
)
from implicit_geometry_kernel.section.layout.system_layout import SYSTEM_TYPE_KEYS, resolve_system_type

__version__ = "0.4.0"

__all__ = [
    "Airfoil", "BladeGeometry", "CellLoop", "CellPolygon", "FeatureTable", "FieldSample", "GeometrySolver",
    "LoftedPolygonSDF", "MidSurface", "NO_FEATURE", "PriorityUnionNode", "SYSTEM_TYPE_KEYS",
    "PartitionEdge", "SectionCut", "SectionWall", "ShearWeb", "ShellGraph", "ShellRun", "SkinRegion",
    "SpanWeb", "Station", "StationNetwork", "TagNode", "TopologyDiff", "TopologyEvent",
    "partition_edge_curve", "resolve_system_type", "section_at", "shared_loft_names", "shell_graph",
    "shell_strips_from_graph",     "station_network", "station_network_from_meta", "topology_diff",
    "wall_meta_from_cut",
]
