"""Spar-cap construction for structural families A/B/C/D.

Thin re-export surface over :mod:`cap_feasibility` and :mod:`build_spar_caps`.
"""

from __future__ import annotations

from .cap_feasibility import (
    CAMBER_CLIP_DEV_BUDGET_M,
    MIN_CAP_VOID_ABS_M,
    MIN_CAP_VOID_FRAC,
    _cap_feasibility_clamp,
    cap_feasibility_clamp,
    default_camber_clip_n_samples,
)
from .build_spar_caps import build_spar_caps

__all__ = [
    "CAMBER_CLIP_DEV_BUDGET_M",
    "MIN_CAP_VOID_ABS_M",
    "MIN_CAP_VOID_FRAC",
    "default_camber_clip_n_samples",
    "cap_feasibility_clamp",
    "_cap_feasibility_clamp",
    "build_spar_caps",
]
