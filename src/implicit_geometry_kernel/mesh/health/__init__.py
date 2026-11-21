"""Mesh health registry: classify and flag for engineer judgement."""

from implicit_geometry_kernel.mesh.health.registry import (
    HEALTH_CLASSES, HEALTH_CODES, HEALTH_LAYERS, HEALTH_SCOPES, HEALTH_SEVERITIES,
    MISSING_ID, HealthClass, HealthHit, HealthThresholds,
    detect_edge_over_target, detect_mid_surface_fold, detect_ridge_mismatch,
    detect_sheet_intersection, diff_failure_elements, diff_health_counts,
    diff_study_folders, edge_health_codes, element_health_codes,
    format_health_summary, health_class, health_counts, health_counts_table,
    hits_to_rows, read_failures_csv, read_health_counts_csv, run_mesh_health,
    top_hits,
)

__all__ = [
    "HEALTH_CLASSES", "HEALTH_CODES", "HEALTH_LAYERS", "HEALTH_SCOPES", "HEALTH_SEVERITIES",
    "MISSING_ID", "HealthClass", "HealthHit", "HealthThresholds",
    "detect_edge_over_target", "detect_mid_surface_fold", "detect_ridge_mismatch",
    "detect_sheet_intersection", "diff_failure_elements", "diff_health_counts",
    "diff_study_folders", "edge_health_codes", "element_health_codes",
    "format_health_summary", "health_class", "health_counts", "health_counts_table",
    "hits_to_rows", "read_failures_csv", "read_health_counts_csv", "run_mesh_health",
    "top_hits",
]
