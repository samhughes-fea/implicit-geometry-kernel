"""Unit tests for coordinate_frames short-tag SoT."""

from __future__ import annotations

import pytest

from blade_geometry_model.coordinate_frames import (
    FRAME_B,
    FRAME_C,
    FRAME_E,
    FRAME_M,
    FRAME_R,
    FRAME_S,
    FRAME_SPRIME,
    FRAME_TAGS,
    require_frame_tag,
)


def test_named_constants_are_catalogue_members() -> None:
    named = (FRAME_R, FRAME_B, FRAME_S, FRAME_SPRIME, FRAME_E, FRAME_M, FRAME_C)
    assert named == FRAME_TAGS
    for tag in named:
        assert require_frame_tag(tag) == tag


@pytest.mark.parametrize("tag", list(FRAME_TAGS))
def test_require_frame_tag_accepts_exact_members(tag: str) -> None:
    assert require_frame_tag(tag) == tag


@pytest.mark.parametrize(
    "bad_tag",
    [None, "", "Sprime", "S′", "s", "BB", "geometry", "local"],
)
def test_require_frame_tag_rejects_aliases_and_unknown(bad_tag) -> None:
    with pytest.raises(ValueError, match="ERROR: frame does not exist"):
        require_frame_tag(bad_tag)
