"""Map placement construction strings onto wall-normal offsets (WP3 D3.1 / D3.2).

Offsets are signed distances along the wall build direction, measured from the
fixed reference face into the wall stack. Positive values always point into the
material (or outward from a midsurface origin for webs).

This module deliberately speaks in *string* reference/build tokens so
``section_geometry`` does not import ``material_utilities`` (G-row boundary).
The layup kernel's ``PlyPlacementPolicy`` values match these tokens 1:1.
"""

from __future__ import annotations

from dataclasses import dataclass

# Accepted constructions — mirrors layup.types.placement._DROP_FACE_BY_CONSTRUCTION.
_ONE_SIDED = frozenset(
    {
        ("outer_mould_line", "inward"),
        ("inner_mould_line", "outward"),
        ("skin_bond_line", "into_cavity"),
    }
)
_MIDSURFACE_ORIGIN = frozenset(
    {
        ("midsurface", "outward_symmetric"),
        ("midsurface", "centred"),
    }
)
_DEFAULT_BY_KIND = {
    "skin": ("outer_mould_line", "inward", "skin_from_outer_mould_line"),
    "cap": ("skin_bond_line", "into_cavity", "spar_cap_from_skin_bond_line"),
    "spar": ("skin_bond_line", "into_cavity", "spar_cap_from_skin_bond_line"),
    "web": ("midsurface", "outward_symmetric", "web_from_midsurface"),
}


@dataclass(frozen=True)
class AnchoredWallGeometry:
    """Physical faces and midsurface of one wall relative to its reference face."""

    reference_face: str
    build_direction: str
    thickness_m: float
    reference_offset_m: float
    free_face_offset_m: float
    midsurface_offset_m: float
    second_free_face_offset_m: float | None = None
    placement_policy_id: str | None = None
    placement_digest: str | None = None

    @property
    def policy_key(self) -> tuple[str, str]:
        return (self.reference_face, self.build_direction)


def _normalise_construction(
    reference_face: str, build_direction: str
) -> tuple[str, str]:
    ref = str(reference_face).strip().lower()
    direction = str(build_direction).strip().lower()
    key = (ref, direction)
    if key not in _ONE_SIDED and key not in _MIDSURFACE_ORIGIN:
        raise ValueError(
            f"unsupported placement construction {ref!r}/{direction!r}; "
            f"accepted: {sorted(_ONE_SIDED | _MIDSURFACE_ORIGIN)}"
        )
    return key


def midsurface_offset_from_reference_m(
    reference_face: str, build_direction: str, thickness_m: float
) -> float:
    """Signed depth of the wall midsurface from the fixed reference face."""
    t = float(thickness_m)
    if not (t > 0.0):
        raise ValueError(f"thickness_m must be positive, got {thickness_m!r}")
    key = _normalise_construction(reference_face, build_direction)
    if key in _MIDSURFACE_ORIGIN:
        return 0.0
    return 0.5 * t


def physical_face_offsets_from_reference_m(
    reference_face: str, build_direction: str, thickness_m: float
) -> tuple[float, float, float | None]:
    """``(reference, free[, second_free])`` offsets from the fixed reference face."""
    t = float(thickness_m)
    if not (t > 0.0):
        raise ValueError(f"thickness_m must be positive, got {thickness_m!r}")
    key = _normalise_construction(reference_face, build_direction)
    if key in _MIDSURFACE_ORIGIN:
        half = 0.5 * t
        return (0.0, half, -half)
    return (0.0, t, None)


def anchored_wall_geometry(
    reference_face: str,
    build_direction: str,
    thickness_m: float,
    *,
    placement_policy_id: str | None = None,
    placement_digest: str | None = None,
) -> AnchoredWallGeometry:
    """Derive face and midsurface offsets for a placement construction."""
    ref, direction = _normalise_construction(reference_face, build_direction)
    free_ref, free, second = physical_face_offsets_from_reference_m(
        ref, direction, thickness_m
    )
    return AnchoredWallGeometry(
        reference_face=ref,
        build_direction=direction,
        thickness_m=float(thickness_m),
        reference_offset_m=free_ref,
        free_face_offset_m=free,
        midsurface_offset_m=midsurface_offset_from_reference_m(
            ref, direction, thickness_m
        ),
        second_free_face_offset_m=second,
        placement_policy_id=placement_policy_id,
        placement_digest=placement_digest,
    )


def anchored_wall_geometry_from_policy(policy: object, thickness_m: float) -> AnchoredWallGeometry:
    """Duck-typed entry for a layup ``PlyPlacementPolicy`` (or any object with the fields)."""
    return anchored_wall_geometry(
        str(getattr(getattr(policy, "reference_face", None), "value", getattr(policy, "reference_face", ""))),
        str(getattr(getattr(policy, "build_direction", None), "value", getattr(policy, "build_direction", ""))),
        thickness_m,
        placement_policy_id=getattr(policy, "name", None),
        placement_digest=getattr(policy, "digest", None),
    )


def default_construction_for_kind(wall_kind: str) -> tuple[str, str, str]:
    """``(reference_face, build_direction, policy_id)`` for a strip kind."""
    key = str(wall_kind or "").strip().lower()
    if key in ("spar_cap",):
        key = "cap"
    if key not in _DEFAULT_BY_KIND:
        raise KeyError(
            f"no default placement construction for wall kind {wall_kind!r}; "
            f"expected one of {sorted(_DEFAULT_BY_KIND)}"
        )
    return _DEFAULT_BY_KIND[key]


def airfoil_inward_depth_for_midsurface(
    reference_face: str,
    thickness_m: float,
    *,
    reference_iml_depth_m: float = 0.0,
) -> float:
    """Depth from the authoring OML to the wall midsurface (airfoil SDF convention)."""
    t = float(thickness_m)
    if not (t > 0.0):
        raise ValueError(f"thickness_m must be positive, got {thickness_m!r}")
    ref = str(reference_face).strip().lower()
    if ref == "outer_mould_line":
        return 0.5 * t
    if ref == "inner_mould_line":
        iml = float(reference_iml_depth_m) if reference_iml_depth_m > 0.0 else t
        return iml - 0.5 * t
    if ref == "midsurface":
        return 0.5 * t if float(reference_iml_depth_m) <= 0.0 else float(reference_iml_depth_m)
    if ref == "skin_bond_line":
        return float(reference_iml_depth_m) + 0.5 * t
    raise ValueError(f"no airfoil mid depth for reference {ref!r}")


__all__ = [
    "AnchoredWallGeometry",
    "airfoil_inward_depth_for_midsurface",
    "anchored_wall_geometry",
    "anchored_wall_geometry_from_policy",
    "default_construction_for_kind",
    "midsurface_offset_from_reference_m",
    "physical_face_offsets_from_reference_m",
]
