"""Thin ``CoordinateFrameTransforms`` facade over the span and tensor kernels."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union

import numpy as np

from .span import SpanFrameTransforms
from . import tensors as _tensors


@dataclass
class CoordinateFrameTransforms:
    """Stable API facade: R↔B span maps plus S↔B tensor transforms.

    Holds a :class:`SpanFrameTransforms` and forwards span methods via
    ``__getattr__``. S↔B second moments, moments, forces, thickness, are
    free functions in ``tensors``; this class delegates so existing
    ``CoordinateFrameTransforms().transform_*`` call sites keep working.
    """

    R: float = 0.0
    rR_grids: Optional[Union[np.ndarray, List[np.ndarray]]] = None
    _span: SpanFrameTransforms = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "_span", SpanFrameTransforms(R=self.R, rR_grids=self.rR_grids)
        )

    def __getattr__(self, name: str):
        return getattr(object.__getattribute__(self, "_span"), name)

    def transform_second_moments_S_to_B(
        self,
        Ixx_S: Union[float, np.ndarray],
        Iyy_S: Union[float, np.ndarray],
        Ixy_S: Union[float, np.ndarray] = 0.0,
        twist_deg: Union[float, np.ndarray] = 0.0,
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray], Union[float, np.ndarray]]:
        return _tensors.transform_second_moments_S_to_B(
            Ixx_S, Iyy_S, Ixy_S, twist_deg
        )

    def transform_second_moments_B_to_S(
        self,
        Ixx_B: Union[float, np.ndarray],
        Iyy_B: Union[float, np.ndarray],
        Ixy_B: Union[float, np.ndarray] = 0.0,
        twist_deg: Union[float, np.ndarray] = 0.0,
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray], Union[float, np.ndarray]]:
        return _tensors.transform_second_moments_B_to_S(
            Ixx_B, Iyy_B, Ixy_B, twist_deg
        )

    def transform_moments_B_to_S(
        self,
        M_YB: Union[float, np.ndarray],
        M_XB: Union[float, np.ndarray] = 0.0,
        twist_deg: Union[float, np.ndarray] = 0.0,
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        return _tensors.transform_moments_B_to_S(M_YB, M_XB, twist_deg)

    def transform_moments_S_to_B(
        self,
        M_YS: Union[float, np.ndarray],
        M_XS: Union[float, np.ndarray] = 0.0,
        twist_deg: Union[float, np.ndarray] = 0.0,
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        return _tensors.transform_moments_S_to_B(M_YS, M_XS, twist_deg)

    def transform_thickness_B_to_Sprime(
        self,
        thickness_B: Union[float, np.ndarray],
        twist_deg: Union[float, np.ndarray],
    ) -> Union[float, np.ndarray]:
        return _tensors.transform_thickness_B_to_Sprime(thickness_B, twist_deg)

    def transform_thickness_Sprime_to_B(
        self,
        thickness_Sprime: Union[float, np.ndarray],
        twist_deg: Union[float, np.ndarray],
    ) -> Union[float, np.ndarray]:
        return _tensors.transform_thickness_Sprime_to_B(thickness_Sprime, twist_deg)

    def transform_forces_B_to_S(
        self,
        V_XB: Union[float, np.ndarray],
        V_YB: Union[float, np.ndarray] = 0.0,
        twist_deg: Union[float, np.ndarray] = 0.0,
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        return _tensors.transform_forces_B_to_S(V_XB, V_YB, twist_deg)

    def transform_forces_S_to_B(
        self,
        V_XS: Union[float, np.ndarray],
        V_YS: Union[float, np.ndarray] = 0.0,
        twist_deg: Union[float, np.ndarray] = 0.0,
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        return _tensors.transform_forces_S_to_B(V_XS, V_YS, twist_deg)
