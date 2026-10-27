"""Stamp placement identity onto ``ShellMidlineStrip`` (WP3 D3.4).

Uses string construction defaults from :mod:`.anchoring` so this module stays
inside the ``section_geometry`` import boundary. Callers that already hold a
layup ``PlyPlacementPolicy`` may pass ``policy_id`` / ``digest`` / ``reference_face``
explicitly via :func:`stamp_strip_placement_fields`.
"""

from __future__ import annotations

from dataclasses import is_dataclass, replace
from typing import Any, Mapping

from blade_geometry_model.section_geometry.placement.anchoring import (
    default_construction_for_kind,
)


def placement_stamp_for_kind(wall_kind: str) -> dict[str, str]:
    """Default placement stamp fields for a strip ``kind`` (policy id + reference)."""
    reference_face, _direction, policy_id = default_construction_for_kind(wall_kind)
    return {
        "placement_policy_id": policy_id,
        "reference_face": reference_face,
    }


def require_placement_stamp(
    *,
    placement_policy_id: str | None,
    placement_digest: str | None,
    reference_face: str | None,
) -> None:
    """Refuse a partial placement claim (policy without digest/face)."""
    claimed = placement_policy_id is not None
    if not claimed:
        return
    missing = [
        name
        for name, value in (
            ("placement_digest", placement_digest),
            ("reference_face", reference_face),
        )
        if value is None or str(value).strip() == ""
    ]
    if missing:
        raise ValueError(
            f"placement_policy_id={placement_policy_id!r} requires {', '.join(missing)}"
        )


def stamp_strip_placement_fields(
    strip: Any,
    *,
    placement_policy_id: str,
    placement_digest: str,
    reference_face: str,
) -> Any:
    """Return ``strip`` with the three placement fields populated."""
    payload = {
        "placement_policy_id": str(placement_policy_id),
        "placement_digest": str(placement_digest),
        "reference_face": str(reference_face),
    }
    require_placement_stamp(**payload)
    if is_dataclass(strip) and hasattr(strip, "placement_policy_id"):
        return replace(strip, **{k: v for k, v in payload.items() if hasattr(strip, k)})
    return strip


def stamp_strip_placement(
    strip: Any,
    *,
    kind: str | None = None,
    placement_digest: str | None = None,
    placement_policy_id: str | None = None,
    reference_face: str | None = None,
) -> Any:
    """Stamp defaults for ``kind``, optionally overriding with explicit fields.

    Digest is required for a complete stamp. When only the wall kind is known,
    the policy id and reference face are filled from the geometry defaults and
    ``placement_digest`` must be supplied by the caller (typically the layup
    family's digest) — otherwise the strip is left unstamped rather than
    inventing a digest inside geometry.
    """
    wall = kind or str(getattr(strip, "kind", "") or "")
    if placement_policy_id is None or reference_face is None:
        if not wall:
            return strip
        ref, _direction, policy_id = default_construction_for_kind(wall)
        placement_policy_id = placement_policy_id or policy_id
        reference_face = reference_face or ref
    if placement_digest is None:
        return strip
    return stamp_strip_placement_fields(
        strip,
        placement_policy_id=str(placement_policy_id),
        placement_digest=str(placement_digest),
        reference_face=str(reference_face),
    )


def stamp_strips_placement(
    strips: list[Any],
    *,
    digest_by_kind: Mapping[str, str] | None = None,
) -> list[Any]:
    """Stamp strips when a per-kind placement digest is provided."""
    digests = dict(digest_by_kind or {})
    out: list[Any] = []
    for strip in strips:
        kind = str(getattr(strip, "kind", "") or "")
        out.append(
            stamp_strip_placement(
                strip,
                kind=kind or None,
                placement_digest=digests.get(kind),
            )
        )
    return out


__all__ = [
    "placement_stamp_for_kind",
    "require_placement_stamp",
    "stamp_strip_placement",
    "stamp_strip_placement_fields",
    "stamp_strips_placement",
]
