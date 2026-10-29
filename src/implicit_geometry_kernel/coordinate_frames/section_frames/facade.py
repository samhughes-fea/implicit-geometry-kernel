"""Thin ``CoordinateFrameTransforms`` facade over span / tensor / stiffness kernels."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Union

import numpy as np

from .span import SpanFrameTransforms
from . import stiffness as _stiffness
from . import tensors as _tensors


@dataclass
class CoordinateFrameTransforms(SpanFrameTransforms):
    """Stable API facade: R↔B span maps plus S↔B tensor / stiffness transforms.

    Span state and converters live on :class:`SpanFrameTransforms`. S↔B second
    moments, moments, forces, thickness, and B→S′ stiffness are free functions
    in ``tensors`` / ``stiffness``; this class delegates so existing
    ``CoordinateFrameTransforms().transform_*`` call sites keep working.
    """

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

    def transform_stiffness_B_to_Sprime(
        self,
        E: float,
        EI_flap_B: Union[float, np.ndarray],
        scaled_dimensions,
        *,
        EI_edge_B: Optional[Union[float, np.ndarray]] = None,
        station_indices: Optional[np.ndarray] = None,
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        return _stiffness.transform_stiffness_B_to_Sprime(
            E,
            EI_flap_B,
            scaled_dimensions,
            EI_edge_B=EI_edge_B,
            station_indices=station_indices,
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
