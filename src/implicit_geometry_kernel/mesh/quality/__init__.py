"""Per-element quality metrics and transparent mesh check rollups."""

from implicit_geometry_kernel.mesh.quality.metrics import (
    MeshCheck, MeshRunQuality, QUALITY_DTYPE, element_quality, mesh_check, mesh_quality,
    quality_rollup, volume_ratio_from_edges,
)

__all__ = [
    "MeshCheck", "MeshRunQuality", "QUALITY_DTYPE", "element_quality", "mesh_check",
    "mesh_quality", "quality_rollup", "volume_ratio_from_edges",
]


def __getattr__(name: str):
    """Expose metrics helpers (including private) for compatibility imports."""
    import implicit_geometry_kernel.mesh.quality.metrics as metrics
    try:
        return getattr(metrics, name)
    except AttributeError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
