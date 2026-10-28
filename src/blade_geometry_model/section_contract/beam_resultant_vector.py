r"""
BeamResultantVector
===================
FE-agnostic DTO for a Saint-Venant six-vector of beam/section resultants,
carrying **canonical component order**, a **structural frame tag** (B / S / S′),
and the **moment reference point** in that frame.

Frame sense and S↔B polarity are owned by
``geometry_utilities.coordinate_frames`` — this module does not reimplement
twist matrices. See that package README and ``section_frames``
home.

**Boundary contract:**

* Storage order is always :data:`BEAM_IO_ORDER`
  ``(N, Vy, Vz, My, Mz, T)``. These ``My``/``Mz``/``Vy``/``Vz`` names are the
  **product shell solver channels** from
  :func:`~blade_geometry_model.coordinate_frames.blade_resultants_to_shell_channels`
  (``Mz``/``Vz`` ← flap, ``My``/``Vy`` ← edge). Consumers that need another
  conjugate order (e.g. shell K6 ``[N, My, Mz, T, Vz, Vy]``) must call
  :meth:`to_mode_conjugate` with an explicit named order — never silently
  reinterpret positional slots.
* ``about_xy`` is the section-plane moment/torque reference in the labelled
  frame, stored as **canonical** ``(X_*, Y_*)`` — same column order as
  ``ShellMidlineStrip.midline_twisted_pitch_axis_blade_frame``. Legacy product
  name ``about_yz`` remains as a deprecated read alias. Rigid transport uses
  :meth:`transported_to`.
* S↔B rotations of the in-plane shear/moment pairs and of ``about_*`` use
  :meth:`rotated_S_to_B` / :meth:`rotated_B_to_S` (delegating to
  ``coordinate_frames``). Generic :meth:`rotated` is the S→B map for
  geometric twist ``theta`` (radians) without retagging.
* This module has no FE / shell / strip dependency; it is a pure handoff contract.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from blade_geometry_model.coordinate_frames import (
    FRAME_B,
    FRAME_S,
    FRAME_SPRIME,
    CoordinateFrameTransforms,
    FrameTag,
    blade_resultants_to_shell_channels,
    require_frame_tag,
    rotate_B_to_S,
    rotate_S_to_B,
)
from blade_geometry_model.section_contract.load_orders import (
    BEAM_IO_ORDER,
    SHELL_K6_MODE_ORDER,
    STRIP_OPTIM_ORDER,
)

# Named component → index in BEAM_IO_ORDER
_IO_INDEX: dict[str, int] = {name: i for i, name in enumerate(BEAM_IO_ORDER)}

# Legacy alias — same value as FRAME_B (midlines / FE resultants live in B).
FRAME_SECTION_YZ = FRAME_B

# Default moment reference: pitch / load-stacking origin in the labelled
# structural frame (pivot chordline at (0, 0) after naca_to_S / chord_to_blade).
BEAM_LOAD_AXIS_XY: tuple[float, float] = (0.0, 0.0)
BEAM_LOAD_AXIS_YZ = BEAM_LOAD_AXIS_XY  # legacy alias

_FRAME_TRANSFORMS = CoordinateFrameTransforms()

def _rotate_about_xy(
    about_xy: tuple[float, float],
    twist_rad: float,
    *,
    s_to_b: bool,
) -> tuple[float, float]:
    """Rotate a canonical ``(X, Y)`` about-point with the geometry polarity."""
    pt = np.asarray([[float(about_xy[0]), float(about_xy[1])]], dtype=np.float64)
    out = rotate_S_to_B(pt, twist_rad) if s_to_b else rotate_B_to_S(pt, twist_rad)
    return float(out[0, 0]), float(out[0, 1])

def _rotate_channels_S_to_B(
    Vy: float,
    Vz: float,
    My: float,
    Mz: float,
    twist_rad: float,
) -> tuple[float, float, float, float]:
    """Product channels through coordinate_frames force/moment S→B maps.

    Channel map (see ``blade_resultants_to_shell_channels``)::

        Mz, Vz ← flap  (M_Y*, V_X*)
        My, Vy ← edge  (M_X*, V_Y*)
    """
    twist_deg = float(np.rad2deg(twist_rad))
    # Moments: (M_Y, M_X) = (Mz, My)
    M_YB, M_XB = _FRAME_TRANSFORMS.transform_moments_S_to_B(
        float(Mz), float(My), twist_deg
    )
    # Forces: (V_X, V_Y) = (Vz, Vy)
    V_XB, V_YB = _FRAME_TRANSFORMS.transform_forces_S_to_B(
        float(Vz), float(Vy), twist_deg
    )
    return float(V_YB), float(V_XB), float(M_XB), float(M_YB)

def _rotate_channels_B_to_S(
    Vy: float,
    Vz: float,
    My: float,
    Mz: float,
    twist_rad: float,
) -> tuple[float, float, float, float]:
    """Inverse of :func:`_rotate_channels_S_to_B`."""
    twist_deg = float(np.rad2deg(twist_rad))
    M_YS, M_XS = _FRAME_TRANSFORMS.transform_moments_B_to_S(
        float(Mz), float(My), twist_deg
    )
    V_XS, V_YS = _FRAME_TRANSFORMS.transform_forces_B_to_S(
        float(Vz), float(Vy), twist_deg
    )
    return float(V_YS), float(V_XS), float(M_XS), float(M_YS)

@dataclass(frozen=True)
class BeamResultantVector:
    """Saint-Venant six-vector with frame and moment reference metadata.

    Parameters
    ----------
    R
        Length-6 array in :data:`BEAM_IO_ORDER` ``(N, Vy, Vz, My, Mz, T)``.
    frame
        Structural frame tag for ``R`` and ``about_xy`` (``\"B\"``, ``\"S\"``,
        or ``\"S'\"``). Prefer :data:`FRAME_B` on the geometry→shell handoff.
    about_xy
        Canonical ``(X_*, Y_*)`` moment/torque reference in ``frame``.
    """

    R: NDArray[np.float64]
    frame: FrameTag
    about_xy: tuple[float, float]

    def __post_init__(self) -> None:
        r = np.asarray(self.R, dtype=np.float64).ravel()
        if r.size != 6:
            raise ValueError(f"R must have length 6 (BEAM_IO_ORDER), got {r.size}")
        if not np.all(np.isfinite(r)):
            raise ValueError("R components must be finite")
        object.__setattr__(self, "R", r.copy())
        x, y = float(self.about_xy[0]), float(self.about_xy[1])
        if not (np.isfinite(x) and np.isfinite(y)):
            raise ValueError("about_xy must be finite")
        object.__setattr__(self, "about_xy", (x, y))
        object.__setattr__(self, "frame", require_frame_tag(self.frame))

    # ------------------------------------------------------------------
    # Named component access
    # ------------------------------------------------------------------

    @property
    def about_yz(self) -> tuple[float, float]:
        """Deprecated alias of :attr:`about_xy` (legacy product ``yz`` name)."""
        warnings.warn(
            "BeamResultantVector.about_yz is deprecated; use about_xy "
            "(canonical (X_*, Y_*) moment/torque reference)",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.about_xy

    def component(self, name: str) -> float:
        """Return one component by BEAM_IO_ORDER name (e.g. ``\"My\"``)."""
        try:
            return float(self.R[_IO_INDEX[name]])
        except KeyError as exc:
            raise KeyError(
                f"unknown component {name!r}; valid names: {BEAM_IO_ORDER}"
            ) from exc

    @property
    def N(self) -> float:
        return float(self.R[0])

    @property
    def Vy(self) -> float:
        return float(self.R[1])

    @property
    def Vz(self) -> float:
        return float(self.R[2])

    @property
    def My(self) -> float:
        return float(self.R[3])

    @property
    def Mz(self) -> float:
        return float(self.R[4])

    @property
    def T(self) -> float:
        return float(self.R[5])

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_named(
        cls,
        *,
        N: float = 0.0,
        Vy: float = 0.0,
        Vz: float = 0.0,
        My: float = 0.0,
        Mz: float = 0.0,
        T: float = 0.0,
        frame: FrameTag = FRAME_B,
        about_xy: tuple[float, float] = BEAM_LOAD_AXIS_XY,
    ) -> BeamResultantVector:
        """Build from named Saint-Venant / product-channel components."""
        r = np.array([N, Vy, Vz, My, Mz, T], dtype=np.float64)
        return cls(R=r, frame=frame, about_xy=about_xy)

    @classmethod
    def from_flap_edge(
        cls,
        *,
        M_flapwise: float = 0.0,
        M_edgewise: float = 0.0,
        V_flapwise: float = 0.0,
        V_edgewise: float = 0.0,
        N: float = 0.0,
        T: float = 0.0,
        frame: FrameTag = FRAME_B,
        about_xy: tuple[float, float] = BEAM_LOAD_AXIS_XY,
    ) -> BeamResultantVector:
        """Build via :func:`blade_resultants_to_shell_channels` (flap/edge → My/Mz)."""
        ch = blade_resultants_to_shell_channels(
            M_flapwise=M_flapwise,
            M_edgewise=M_edgewise,
            V_flapwise=V_flapwise,
            V_edgewise=V_edgewise,
            N=N,
            T=T,
        )
        return cls.from_named(
            N=ch.N,
            Vy=ch.Vy,
            Vz=ch.Vz,
            My=ch.My,
            Mz=ch.Mz,
            T=ch.T,
            frame=frame,
            about_xy=about_xy,
        )

    @classmethod
    def from_row(
        cls,
        row: Sequence[float] | NDArray[np.float64],
        order: Sequence[str],
        *,
        frame: FrameTag = FRAME_B,
        about_xy: tuple[float, float] = BEAM_LOAD_AXIS_XY,
    ) -> BeamResultantVector:
        """Build from a positional row given an **explicit** source ``order``.

        Never pass a bare array without ``order`` — that is the slot-swap bug
        class this contract eliminates.
        """
        names = tuple(str(n) for n in order)
        if len(names) != 6:
            raise ValueError(f"order must have length 6, got {len(names)}")
        if set(names) != set(BEAM_IO_ORDER):
            raise ValueError(
                f"order must be a permutation of {BEAM_IO_ORDER}, got {names}"
            )
        arr = np.asarray(row, dtype=np.float64).ravel()
        if arr.size != 6:
            raise ValueError(f"row must have length 6, got {arr.size}")
        r = np.empty(6, dtype=np.float64)
        for i, name in enumerate(names):
            r[_IO_INDEX[name]] = arr[i]
        return cls(R=r, frame=frame, about_xy=about_xy)

    @classmethod
    def from_extreme_loads_station(
        cls,
        extreme_loads: Any,
        station_index: int,
        *,
        frame: FrameTag = FRAME_B,
        about_xy: tuple[float, float] = BEAM_LOAD_AXIS_XY,
    ) -> BeamResultantVector:
        """Adapt one station from :class:`~blade_precompute.beam_resultants.core.types.ExtremeLoads`."""
        i = int(station_index)
        return cls.from_named(
            N=float(extreme_loads.N[i]),
            Vy=float(extreme_loads.Vy[i]),
            Vz=float(extreme_loads.Vz[i]),
            My=float(extreme_loads.My[i]),
            Mz=float(extreme_loads.Mz[i]),
            T=float(extreme_loads.T[i]),
            frame=frame,
            about_xy=about_xy,
        )

    @classmethod
    def from_integrated_resultants_station(
        cls,
        integrated: Any,
        station_index: int,
        *,
        frame: FrameTag = FRAME_B,
        about_xy: tuple[float, float] = BEAM_LOAD_AXIS_XY,
    ) -> BeamResultantVector:
        """Adapt one station from :class:`~blade_utilities.loads.IntegratedResultants`."""
        i = int(station_index)
        return cls.from_named(
            N=float(integrated.N[i]),
            Vy=float(integrated.Vy[i]),
            Vz=float(integrated.Vz[i]),
            My=float(integrated.My[i]),
            Mz=float(integrated.Mz[i]),
            T=float(integrated.T[i]),
            frame=frame,
            about_xy=about_xy,
        )

    @classmethod
    def from_mapping(
        cls,
        values: Mapping[str, float],
        *,
        frame: FrameTag = FRAME_B,
        about_xy: tuple[float, float] = BEAM_LOAD_AXIS_XY,
    ) -> BeamResultantVector:
        """Build from a mapping of component name → value (missing → 0)."""
        return cls.from_named(
            N=float(values.get("N", 0.0)),
            Vy=float(values.get("Vy", 0.0)),
            Vz=float(values.get("Vz", 0.0)),
            My=float(values.get("My", 0.0)),
            Mz=float(values.get("Mz", 0.0)),
            T=float(values.get("T", 0.0)),
            frame=frame,
            about_xy=about_xy,
        )

    # ------------------------------------------------------------------
    # Transforms
    # ------------------------------------------------------------------

    def transported_to(self, x: float, y: float) -> BeamResultantVector:
        r"""Rigid equipollent transport of moments about a new reference point.

        Arguments are canonical section-plane ``(X_*, Y_*)`` (legacy call sites
        may still name them ``y, z`` — same slots).

        With \(\Delta X = X_{\mathrm{old}}-X\), \(\Delta Y = Y_{\mathrm{old}}-Y\)
        and product channels ``Vz``/``Vy`` = flap/edge = ``F_{X}``/``F_{Y}``:

        * \(T' = T + \Delta X\, V_y - \Delta Y\, V_z\)
          (geometric \(T = X F_Y - Y F_X\))
        * \(M_y' = M_y + \Delta Y\, N\)
        * \(M_z' = M_z - \Delta X\, N\)

        Forces are unchanged. Frame is preserved.
        """
        x_new, y_new = float(x), float(y)
        x_old, y_old = self.about_xy
        dx = x_old - x_new
        dy = y_old - y_new
        if dx == 0.0 and dy == 0.0:
            return self
        N, Vy, Vz, My, Mz, T = (float(v) for v in self.R)
        T_p = T + dx * Vy - dy * Vz
        My_p = My + dy * N
        Mz_p = Mz - dx * N
        return BeamResultantVector(
            R=np.array([N, Vy, Vz, My_p, Mz_p, T_p], dtype=np.float64),
            frame=self.frame,
            about_xy=(x_new, y_new),
        )

    def rotated(self, theta: float) -> BeamResultantVector:
        """Apply geometric twist ``theta`` (radians) as an S→B map (no retag).

        Product channels and ``about_xy`` use
        ``coordinate_frames`` polarity (``transform_*_S_to_B`` /
        ``rotate_S_to_B``). Prefer :meth:`rotated_S_to_B` /
        :meth:`rotated_B_to_S` when crossing named frames.
        """
        return self._rotated_impl(float(theta), s_to_b=True, frame=self.frame)

    def rotated_S_to_B(self, twist_rad: float) -> BeamResultantVector:
        """Map components and ``about_xy`` from S → B; retag ``frame`` to :data:`FRAME_B`."""
        return self._rotated_impl(float(twist_rad), s_to_b=True, frame=FRAME_B)

    def rotated_B_to_S(self, twist_rad: float) -> BeamResultantVector:
        """Map components and ``about_xy`` from B → S; retag ``frame`` to :data:`FRAME_S`."""
        return self._rotated_impl(float(twist_rad), s_to_b=False, frame=FRAME_S)

    def _rotated_impl(
        self,
        twist_rad: float,
        *,
        s_to_b: bool,
        frame: FrameTag,
    ) -> BeamResultantVector:
        t = float(twist_rad)
        if abs(t) < 1e-15 and frame == self.frame:
            return self
        N, Vy, Vz, My, Mz, T = (float(v) for v in self.R)
        if abs(t) < 1e-15:
            Vy_p, Vz_p, My_p, Mz_p = Vy, Vz, My, Mz
            about = self.about_xy
        elif s_to_b:
            Vy_p, Vz_p, My_p, Mz_p = _rotate_channels_S_to_B(Vy, Vz, My, Mz, t)
            about = _rotate_about_xy(self.about_xy, t, s_to_b=True)
        else:
            Vy_p, Vz_p, My_p, Mz_p = _rotate_channels_B_to_S(Vy, Vz, My, Mz, t)
            about = _rotate_about_xy(self.about_xy, t, s_to_b=False)
        return BeamResultantVector(
            R=np.array([N, Vy_p, Vz_p, My_p, Mz_p, T], dtype=np.float64),
            frame=frame,
            about_xy=about,
        )

    def with_frame(self, frame: FrameTag) -> BeamResultantVector:
        """Return a copy with a new frame tag (components / about unchanged)."""
        if frame == self.frame:
            return self
        return BeamResultantVector(R=self.R, frame=frame, about_xy=self.about_xy)

    def with_about_xy(self, x: float, y: float) -> BeamResultantVector:
        """Retag ``about_xy`` without transporting (rare; prefer :meth:`transported_to`)."""
        return BeamResultantVector(
            R=self.R, frame=self.frame, about_xy=(float(x), float(y))
        )

    def with_about_yz(self, x: float, y: float) -> BeamResultantVector:
        """Deprecated alias of :meth:`with_about_xy`."""
        warnings.warn(
            "BeamResultantVector.with_about_yz is deprecated; use with_about_xy",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.with_about_xy(x, y)

    def to_mode_conjugate(self, mode_order: Sequence[str]) -> NDArray[np.float64]:
        """Project components into a named conjugate order (length-6 array).

        Example (shell periodic K6)::

            f_beam = vec.to_mode_conjugate(SHELL_K6_MODE_ORDER)
            # → [N, My, Mz, T, Vz, Vy]
        """
        names = tuple(str(n) for n in mode_order)
        if len(names) != 6:
            raise ValueError(f"mode_order must have length 6, got {len(names)}")
        if set(names) != set(BEAM_IO_ORDER):
            raise ValueError(
                f"mode_order must be a permutation of {BEAM_IO_ORDER}, got {names}"
            )
        return np.array([self.component(n) for n in names], dtype=np.float64)

    def to_shell_k6_conjugate(self) -> NDArray[np.float64]:
        """Shortcut: shell periodic mode order ``[N, My, Mz, T, Vz, Vy]``."""
        return self.to_mode_conjugate(SHELL_K6_MODE_ORDER)

    def to_beam_io_array(self) -> NDArray[np.float64]:
        """Return a copy of ``R`` in :data:`BEAM_IO_ORDER`."""
        return self.R.copy()

    def to_named_dict(self) -> dict[str, float]:
        """Named components in BEAM_IO_ORDER."""
        return {name: float(self.R[i]) for i, name in enumerate(BEAM_IO_ORDER)}

__all__ = [
    "BEAM_IO_ORDER",
    "SHELL_K6_MODE_ORDER",
    "STRIP_OPTIM_ORDER",
    "FRAME_B",
    "FRAME_S",
    "FRAME_SPRIME",
    "FRAME_SECTION_YZ",
    "BEAM_LOAD_AXIS_XY",
    "BEAM_LOAD_AXIS_YZ",
    "BeamResultantVector",
]
