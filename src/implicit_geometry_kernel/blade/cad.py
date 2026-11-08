"""Optional STEP export of mid-surface patches (extra ``[cad]``, cadquery)."""

from __future__ import annotations

from pathlib import Path


def cadquery_available() -> bool:
    try:
        import cadquery  # noqa: F401
    except ImportError:
        return False
    return True


def write_mid_surface_step(path, patches) -> str:
    """Write lofted mid-surface patches as a STEP file.

    ``patches`` is ``[(name, rings)]`` with ``rings`` shaped ``(n_span, n_chord, 3)``.
    """
    try:
        import cadquery as cq
    except ImportError as exc:
        raise ImportError(
            "STEP export needs extra [cad] (cadquery). Install implicit-geometry-kernel[cad]."
        ) from exc

    path = Path(path)
    solids = []
    for name, rings in patches:
        rings = list(rings)
        wires = []
        for ring in rings:
            pts = [cq.Vector(float(p[0]), float(p[1]), float(p[2])) for p in ring]
            pts.append(pts[0])
            wires.append(cq.Wire.makePolygon(pts))
        solids.append(cq.Solid.makeLoft(wires))
    if not solids:
        raise ValueError("no mid-surface patches to export")
    compound = solids[0]
    for extra in solids[1:]:
        compound = compound.fuse(extra)
    cq.exporters.export(compound, str(path))
    return str(path)


__all__ = ["cadquery_available", "write_mid_surface_step"]
