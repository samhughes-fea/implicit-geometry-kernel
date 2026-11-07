"""Blades and wings: airfoils, spanwise stations, lofting and cut sections."""

from implicit_geometry_kernel.blade.airfoil import Airfoil
from implicit_geometry_kernel.blade.cut import (
    CellPolygon, MidSegment, MidSurface, PartitionEdge, SectionCut, SectionWall, ShellGraph, ShellRun,
    TopologyDiff, TopologyEvent, partition_edge_curve, shared_loft_names, shell_graph,
    shell_strips_from_graph, topology_diff,
)
from implicit_geometry_kernel.blade.loft import BladeGeometry, SpanTable, Thickness, section_at
from implicit_geometry_kernel.blade.windio import (
    cap_chord_windows_from_windio_structure,
    nd_arc_to_chord_frac,
    web_positions_from_windio_structure,
)
from implicit_geometry_kernel.blade.network import (
    StationNetwork, WallBranch, station_network, station_network_from_meta, wall_meta_from_cut,
)
from implicit_geometry_kernel.blade.io import read_planform_csv, write_windio
from implicit_geometry_kernel.blade.stations import PlyDrop, ShearWeb, SkinRegion, Station, Stringer, TransverseMember

__all__ = ["Airfoil", "BladeGeometry", "CellLoop", "CellPolygon", "MidSegment", "MidSurface", "PartitionEdge",
           "SectionCut", "SectionWall",
           "ShellGraph", "ShellRun",
           "PlyDrop", "ShearWeb", "SkinRegion", "SpanTable", "SpanWeb", "Station", "StationNetwork",
           "Stringer", "Thickness", "TransverseMember", "WallBranch",
           "read_planform_csv", "write_windio",
           "TopologyDiff", "TopologyEvent", "partition_edge_curve", "section_at", "shared_loft_names",
           "shell_graph", "shell_strips_from_graph", "station_network", "station_network_from_meta",
           "topology_diff",
           "cap_chord_windows_from_windio_structure", "nd_arc_to_chord_frac",
           "wall_meta_from_cut", "web_positions_from_windio_structure"]


def __getattr__(name: str):
    if name == "CellLoop":
        import warnings

        warnings.warn(
            "CellLoop is deprecated; use CellPolygon. Removed in 0.6.0.",
            DeprecationWarning,
            stacklevel=2,
        )
        return CellPolygon
    if name == "SpanWeb":
        import warnings

        warnings.warn(
            "SpanWeb is deprecated; use ShearWeb. Removed in 0.6.0.",
            DeprecationWarning,
            stacklevel=2,
        )
        return ShearWeb
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
