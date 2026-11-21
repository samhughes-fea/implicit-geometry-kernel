"""Shared pytest hooks."""

from __future__ import annotations

import pytest

# Tests that take 4 s or more on a laptop (``pytest --durations=60``). Pull requests
# run the full suite on one Python version and ``-m "not slow"`` on the others;
# pushes to main run everything everywhere. Update the list when timings move.
SLOW = {
    "tests/section/test_plot_section_with_quadtree.py::test_plot_section_with_quadtree_two_axes_and_rectangles",
    "tests/test_mapped_strip_quads.py::test_parametric_blade_quads_match_python_loops",
    "tests/test_mesh_health.py::test_write_mesh_tables_emits_failures_and_health_counts",
    "tests/test_reference_transition.py::test_iea10_cap_loft_key_stable_along_span",
    "tests/test_reference_transition.py::test_iea10_loft_keepset_matches_intercept",
    "tests/test_reference_transition.py::test_iea22_te_corners_do_not_chatter",
    "tests/test_reference_transition.py::test_windio_span_metrics",
    "tests/test_shell_mesh.py::test_a_dropped_web_ends_on_its_station_and_the_tip_uses_the_webs_that_remain",
    "tests/test_shell_mesh.py::test_elements_are_quads_and_triangles",
    "tests/test_shell_mesh.py::test_elements_carry_their_partition_data",
    "tests/test_shell_mesh.py::test_leading_edge_and_span_refinement",
    "tests/test_shell_mesh.py::test_nrel_te_cap_seam_shares_endpoints_on_probe_stations",
    "tests/test_shell_mesh.py::test_sized_quads_have_no_slivers",
    "tests/test_shell_mesh.py::test_the_default_edge_is_five_percent_of_chord_and_every_face_uses_the_sheet",
    "tests/test_shell_mesh.py::test_vtk_keeps_triangles",
    "tests/test_shell_topology.py::test_09_le_loft_tracks_web_indices",
    "tests/test_shell_topology.py::test_09_pa_4d_pressure_skin_covers_web_events",
    "tests/test_shell_topology.py::test_09_pa_4d_web3_onset_and_drop_are_skin_ruled",
    "tests/test_shell_topology.py::test_09_pa_loft_keeps_surviving_web_index",
    "tests/test_shell_topology.py::test_09_te_chatter_and_tip_span",
    "tests/test_shell_topology.py::test_09_te_wrap_held_without_flicker",
    "tests/test_shell_topology.py::test_catalogue_2df_te_trace_is_smooth",
}


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    # Tests that need the private reference-blades catalogue carry a skipif whose
    # reason names it. Mark them ``reference`` so ``pytest -m reference`` selects them.
    for item in items:
        if item.nodeid.split("[")[0] in SLOW:
            item.add_marker(pytest.mark.slow)
        for mark in item.iter_markers("skipif"):
            if "catalogue" in str(mark.kwargs.get("reason", "")):
                item.add_marker(pytest.mark.reference)
                break
