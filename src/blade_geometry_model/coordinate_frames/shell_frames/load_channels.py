"""Product shell FE load-channel adapter (B resultants → solver DOFs)."""

from typing import NamedTuple


class ShellSolverChannels(NamedTuple):
    """Product shell FE solver DOF channels (adapter names — not element \(M_{11}\)).

    Attributes
    ----------
    My, Mz, Vy, Vz, N, T
        Saint-Venant components in the product map: ``Mz``/``Vz`` = flap,
        ``My``/``Vy`` = edge.
    """

    N: float
    Vy: float
    Vz: float
    My: float
    Mz: float
    T: float


def blade_resultants_to_shell_channels(
    *,
    M_flapwise: float = 0.0,
    M_edgewise: float = 0.0,
    V_flapwise: float = 0.0,
    V_edgewise: float = 0.0,
    N: float = 0.0,
    T: float = 0.0,
) -> ShellSolverChannels:
    """Map named B-frame / physical beam resultants to product shell solver DOFs.

    Production shell sizing and BSF real-loads use::

        Mz, Vz ← flap   (m_YB, v_XB)
        My, Vy ← edge   (m_XB, v_YB)

    Prefer calling this helper instead of hard-coding ``My``/``Mz`` at call sites
    so the flap/edge map cannot drift.
    """
    return ShellSolverChannels(
        N=float(N),
        Vy=float(V_edgewise),
        Vz=float(V_flapwise),
        My=float(M_edgewise),
        Mz=float(M_flapwise),
        T=float(T),
    )
