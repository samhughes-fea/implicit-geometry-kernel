"""Paint-stable loop roles for spanwise template compatibility."""

from __future__ import annotations

from collections.abc import Sequence


def role_atom(name: str) -> tuple:
    """Paint-stable role for compatibility (cap flicker shares a cap atom)."""
    if name == "te" or name.startswith("te_"):
        return ("TE",)
    if name.startswith("le_"):
        return ("LE", name[-2:])
    if name.startswith("cap_"):
        return ("cap", name[4:6])
    if name.startswith("bay_"):
        parts = name.split("_")
        side = parts[1] if len(parts) > 1 else ""
        return ("bay", side)
    if name.startswith("web"):
        return ("web", name)
    return ("panel", name)


def signatures_compatible(a, b) -> bool:
    """True if role atoms match after TE wrap / cap-paint flicker.

    Web add/drop is **not** compatible (new zone) even if loop atoms match.
    """
    if a == b:
        return True
    if not isinstance(a, tuple) or not isinstance(b, tuple):
        return False
    if len(a) == 2 and len(b) == 2 and isinstance(a[0], tuple):
        ra, wa = a
        rb, wb = b
        if wa != wb:
            return False
        return _loop_atoms(ra) == _loop_atoms(rb)
    return _loop_atoms(a) == _loop_atoms(b)


def _loop_atoms(names: Sequence[str]) -> frozenset:
    atoms = []
    for n in names:
        atom = role_atom(n)
        if atom[0] == "bay":
            atoms.append(("panel", atom[1]))
        elif atom[0] == "cap":
            atoms.append(("panel", atom[1]))
        else:
            atoms.append(atom)
    return frozenset(atoms)


__all__ = ["role_atom", "signatures_compatible"]
