"""``blade.cut`` no longer re-exports shell-graph names."""

import implicit_geometry_kernel.blade.cut as cut

SHELL_GRAPH_NAMES = (
    "MidSegment",
    "MidSurface",
    "PartitionEdge",
    "ShellGraph",
    "ShellRun",
    "TopologyDiff",
    "TopologyEvent",
    "partition_edge_curve",
    "shared_loft_names",
    "shell_graph",
    "shell_strips_from_graph",
    "topology_diff",
    "LockedLoftPlan",
    "resample_locked",
    "shell_loft_zones",
)


def test_cut_does_not_reexport_shell_graph_names() -> None:
    for name in SHELL_GRAPH_NAMES:
        assert not hasattr(cut, name), name
