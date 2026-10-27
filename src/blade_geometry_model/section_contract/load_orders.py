"""Named beam six-vector conjugate orders (load / FE handoff — not geometry policy).

``BeamResultantVector`` storage remains :data:`BEAM_IO_ORDER`. Product channel
sense (``Mz``/``Vz`` ← flap, ``My``/``Vy`` ← edge) is owned by
``coordinate_frames.blade_resultants_to_shell_channels``. Consumers that need
another conjugate order call ``to_mode_conjugate`` with these named tuples.
"""

from __future__ import annotations

# Beam six-vector I/O index convention (Saint-Venant storage order).
BEAM_IO_ORDER: tuple[str, ...] = ("N", "Vy", "Vz", "My", "Mz", "T")

# Shell periodic K6 mode conjugate order (FE / load channel).
SHELL_K6_MODE_ORDER: tuple[str, ...] = ("N", "My", "Mz", "T", "Vz", "Vy")

# Strip / optimiser R_sec row order historically used in DesignEvaluator.
STRIP_OPTIM_ORDER: tuple[str, ...] = ("N", "My", "Mz", "T", "Vy", "Vz")

__all__ = [
    "BEAM_IO_ORDER",
    "SHELL_K6_MODE_ORDER",
    "STRIP_OPTIM_ORDER",
]
