"""E ↔ M ply frame transforms (laminate/element ↔ ply material)."""

import numpy as np


# **E** (element) axes ≡ laminate axes ≡ \((x_1,x_2,x_3)\): \(x_1\) strip/laminate
# longitudinal, \(x_2\) contour tangent, \(x_3\) through-thickness. "Laminate frame"
# in ``material_utilities`` and "element frame" in the shell model name the same frame.
#
# **M** (ply material) is E rotated in-plane about \(x_3\) by the ply angle θ; axes
# stay \((1,2,3)\) with 1 = fibre.
#
# Voigt stress \([\sigma_{11},\sigma_{22},\tau_{12}]\) and engineering strain
# \([\varepsilon_{11},\varepsilon_{22},\gamma_{12}]\) **do not share a transformation
# matrix** — they differ by factors of 2 on the shear terms. Using one for the other
# silently corrupts the shear component. See README.md §7.


def ply_stress_transform_E_to_M(theta_deg) -> np.ndarray:
    r"""Voigt **stress** transform \(T_\sigma\): ``sigma_M = T @ sigma_E``.

    Parameters
    ----------
    theta_deg : float or np.ndarray
        Ply angle(s) [degrees], rotation about \(x_3\) from E to M.

    Returns
    -------
    np.ndarray
        Shape ``theta_deg.shape + (3, 3)``. Scalar θ gives ``(3, 3)``.

    Notes
    -----
    Do not use this on strain — see :func:`ply_strain_transform_E_to_M`.
    """
    t = np.deg2rad(np.asarray(theta_deg, dtype=np.float64))
    c, s = np.cos(t), np.sin(t)
    cc, ss, sc = c * c, s * s, s * c
    out = np.empty(t.shape + (3, 3), dtype=np.float64)
    out[..., 0, 0], out[..., 0, 1], out[..., 0, 2] = cc, ss, 2.0 * sc
    out[..., 1, 0], out[..., 1, 1], out[..., 1, 2] = ss, cc, -2.0 * sc
    out[..., 2, 0], out[..., 2, 1], out[..., 2, 2] = -sc, sc, cc - ss
    return out


def ply_strain_transform_E_to_M(theta_deg) -> np.ndarray:
    r"""Engineering **strain** transform \(T_\varepsilon\): ``eps_M = T @ eps_E``.

    Differs from :func:`ply_stress_transform_E_to_M` by the factor-2 shear terms
    (``sc`` vs ``2sc`` off-diagonal, ``-2sc`` vs ``-sc`` on the shear row) — a
    consequence of engineering shear strain \(\\gamma_{12} = 2\\varepsilon_{12}\).

    Parameters
    ----------
    theta_deg : float or np.ndarray
        Ply angle(s) [degrees], rotation about \(x_3\) from E to M.

    Returns
    -------
    np.ndarray
        Shape ``theta_deg.shape + (3, 3)``. Scalar θ gives ``(3, 3)``.
    """
    t = np.deg2rad(np.asarray(theta_deg, dtype=np.float64))
    c, s = np.cos(t), np.sin(t)
    cc, ss, sc = c * c, s * s, s * c
    out = np.empty(t.shape + (3, 3), dtype=np.float64)
    out[..., 0, 0], out[..., 0, 1], out[..., 0, 2] = cc, ss, sc
    out[..., 1, 0], out[..., 1, 1], out[..., 1, 2] = ss, cc, -sc
    out[..., 2, 0], out[..., 2, 1], out[..., 2, 2] = -2.0 * sc, 2.0 * sc, cc - ss
    return out


def ply_stress_transform_M_to_E(theta_deg) -> np.ndarray:
    """Inverse of :func:`ply_stress_transform_E_to_M` (rotation by ``-θ``)."""
    return ply_stress_transform_E_to_M(-np.asarray(theta_deg, dtype=np.float64))


def ply_strain_transform_M_to_E(theta_deg) -> np.ndarray:
    """Inverse of :func:`ply_strain_transform_E_to_M` (rotation by ``-θ``)."""
    return ply_strain_transform_E_to_M(-np.asarray(theta_deg, dtype=np.float64))


def rotate_ply_vectors(vec: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Apply per-ply ``(3, 3)`` transforms to per-ply 3-vectors.

    Low-level applier shared by the convenience rotators and by
    ``material_utilities.laminate.clpt_recovery.rotate_plies_to_material``.

    Parameters
    ----------
    vec : np.ndarray
        ``(..., n_comp, n_ply, 3)`` or ``(..., n_comp, n_ply, n_face, 3)``.
    T : np.ndarray
        ``(..., n_comp, n_ply, 3, 3)`` matching ``vec``'s component/ply axes.

    Returns
    -------
    np.ndarray
        Same shape as ``vec``.

    Notes
    -----
    The face axis is detected as ``vec.ndim >= 4``; a bare ``(n_ply, n_face, 3)``
    array is therefore treated as ``(n_comp, n_ply, 3)``. Callers with a face axis
    and no component axis must add a leading axis first (as
    ``ply_strains_stresses_material_frame`` does).
    """
    vec = np.asarray(vec, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    if vec.ndim >= 4:
        return np.einsum("...cpij,...cpfj->...cpfi", T, vec)
    return np.einsum("...cpij,...cpj->...cpi", T, vec)


def _matvec(T: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """Broadcasting ``(..., 3, 3) @ (..., 3)`` matrix-vector product."""
    return np.einsum("...ij,...j->...i", np.asarray(T, dtype=np.float64),
                     np.asarray(vec, dtype=np.float64))


def rotate_ply_stress_E_to_M(sigma_E: np.ndarray, theta_deg) -> np.ndarray:
    """Rotate Voigt stress from laminate/element axes into ply material axes.

    ``sigma_E`` is ``(..., 3)`` = ``[σ11, σ22, τ12]``; ``theta_deg`` must broadcast
    against ``sigma_E.shape[:-1]`` (e.g. pass ``(n_ply, 1)`` for a ``(n_ply, n_face, 3)``
    array).
    """
    return _matvec(ply_stress_transform_E_to_M(theta_deg), sigma_E)


def rotate_ply_stress_M_to_E(sigma_M: np.ndarray, theta_deg) -> np.ndarray:
    """Inverse of :func:`rotate_ply_stress_E_to_M`."""
    return _matvec(ply_stress_transform_M_to_E(theta_deg), sigma_M)


def rotate_ply_strain_E_to_M(eps_E: np.ndarray, theta_deg) -> np.ndarray:
    """Rotate engineering strain from laminate/element axes into ply material axes.

    ``eps_E`` is ``(..., 3)`` = ``[ε11, ε22, γ12]``; ``theta_deg`` must broadcast
    against ``eps_E.shape[:-1]``.
    """
    return _matvec(ply_strain_transform_E_to_M(theta_deg), eps_E)


def rotate_ply_strain_M_to_E(eps_M: np.ndarray, theta_deg) -> np.ndarray:
    """Inverse of :func:`rotate_ply_strain_E_to_M`."""
    return _matvec(ply_strain_transform_M_to_E(theta_deg), eps_M)


# Canonical 3D Voigt order: [11, 22, 33, 23, 13, 12]
_CANONICAL_3D = ("11", "22", "33", "23", "13", "12")


def ply_stress_transform_3d_E_to_M(theta_deg) -> np.ndarray:
    r"""Six-component Voigt stress transform about \(x_3\): ``sigma_M = T @ sigma_E``.

    Ordering is ``[σ11, σ22, σ33, τ23, τ13, τ12]``.
    """
    t = np.deg2rad(np.asarray(theta_deg, dtype=np.float64))
    c, s = np.cos(t), np.sin(t)
    cc, ss, sc = c * c, s * s, s * c
    out = np.zeros(t.shape + (6, 6), dtype=np.float64)
    out[..., 0, 0], out[..., 0, 1], out[..., 0, 5] = cc, ss, 2.0 * sc
    out[..., 1, 0], out[..., 1, 1], out[..., 1, 5] = ss, cc, -2.0 * sc
    out[..., 2, 2] = 1.0
    out[..., 3, 3], out[..., 3, 4] = c, -s
    out[..., 4, 3], out[..., 4, 4] = s, c
    out[..., 5, 0], out[..., 5, 1], out[..., 5, 5] = -sc, sc, cc - ss
    return out


def ply_strain_transform_3d_E_to_M(theta_deg) -> np.ndarray:
    r"""Six-component engineering-strain transform about \(x_3\).

    Ordering is ``[ε11, ε22, ε33, γ23, γ13, γ12]``.
    """
    t = np.deg2rad(np.asarray(theta_deg, dtype=np.float64))
    c, s = np.cos(t), np.sin(t)
    cc, ss, sc = c * c, s * s, s * c
    out = np.zeros(t.shape + (6, 6), dtype=np.float64)
    out[..., 0, 0], out[..., 0, 1], out[..., 0, 5] = cc, ss, sc
    out[..., 1, 0], out[..., 1, 1], out[..., 1, 5] = ss, cc, -sc
    out[..., 2, 2] = 1.0
    out[..., 3, 3], out[..., 3, 4] = c, -s
    out[..., 4, 3], out[..., 4, 4] = s, c
    out[..., 5, 0], out[..., 5, 1], out[..., 5, 5] = -2.0 * sc, 2.0 * sc, cc - ss
    return out


def ply_stress_transform_3d_M_to_E(theta_deg) -> np.ndarray:
    """Inverse of :func:`ply_stress_transform_3d_E_to_M`."""
    return ply_stress_transform_3d_E_to_M(-np.asarray(theta_deg, dtype=np.float64))


def ply_strain_transform_3d_M_to_E(theta_deg) -> np.ndarray:
    """Inverse of :func:`ply_strain_transform_3d_E_to_M`."""
    return ply_strain_transform_3d_E_to_M(-np.asarray(theta_deg, dtype=np.float64))


def rotate_ply_stress_3d_E_to_M(sigma_E: np.ndarray, theta_deg) -> np.ndarray:
    """Rotate canonical 3D Voigt stress from E to M."""
    return _matvec(ply_stress_transform_3d_E_to_M(theta_deg), sigma_E)


def rotate_ply_stress_3d_M_to_E(sigma_M: np.ndarray, theta_deg) -> np.ndarray:
    """Rotate canonical 3D Voigt stress from M to E."""
    return _matvec(ply_stress_transform_3d_M_to_E(theta_deg), sigma_M)


def rotate_ply_strain_3d_E_to_M(eps_E: np.ndarray, theta_deg) -> np.ndarray:
    """Rotate canonical 3D engineering strain from E to M."""
    return _matvec(ply_strain_transform_3d_E_to_M(theta_deg), eps_E)


def rotate_ply_strain_3d_M_to_E(eps_M: np.ndarray, theta_deg) -> np.ndarray:
    """Rotate canonical 3D engineering strain from M to E."""
    return _matvec(ply_strain_transform_3d_M_to_E(theta_deg), eps_M)


def rotate_transverse_shear_E_to_M(
    tau13_E: np.ndarray,
    tau23_E: np.ndarray,
    theta_deg,
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate ``(τ13, τ23)`` from laminate E into ply material M via the 3D map."""
    tau13 = np.asarray(tau13_E, dtype=np.float64)
    tau23 = np.asarray(tau23_E, dtype=np.float64)
    zeros = np.zeros_like(tau13)
    vec = np.stack([zeros, zeros, zeros, tau23, tau13, zeros], axis=-1)
    rotated = rotate_ply_stress_3d_E_to_M(vec, theta_deg)
    return rotated[..., 4], rotated[..., 3]

