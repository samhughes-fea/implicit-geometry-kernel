"""Wall-solid / core / residual label classification for hybrid properties."""

from __future__ import annotations

from typing import Any, Iterable


def _label_key(label: Any) -> str:
    val = getattr(label, "value", label)
    return str(val)


def is_core_label(label: Any) -> bool:
    """True for sandwich-core bay labels (constructive polygon path)."""
    key = _label_key(label)
    return key.startswith("sandwich_core_") or key in {"sandwich_core"}


# LE/TE inserts use constructive polygons by default (parity locked).
INSERTS_CONSTRUCTIVE_BY_DEFAULT = True


def is_insert_label(label: Any) -> bool:
    """True for LE/TE insert labels."""
    key = _label_key(label)
    return key in {
        "leading_edge_insert",
        "trailing_edge_insert",
    }


def is_residual_label(label: Any) -> bool:
    """True for LE/TE inserts on the residual φ path. Cores are constructive."""
    if not is_insert_label(label):
        return False
    return not INSERTS_CONSTRUCTIVE_BY_DEFAULT


def is_skin_label(label: Any) -> bool:
    """True for outer-skin wall labels (full or suction/pressure)."""
    key = _label_key(label)
    return (
        key == "outer_skin"
        or key.startswith("skin_suction_surface")
        or key.startswith("skin_pressure_surface")
    )


def is_wall_solid_label(label: Any) -> bool:
    """True for skin / web / spar-cap / constructive-insert wall solids."""
    key = _label_key(label)
    if is_insert_label(key) and INSERTS_CONSTRUCTIVE_BY_DEFAULT:
        return True
    if is_residual_label(key) or is_core_label(key):
        return False
    if is_skin_label(label):
        return True
    if key.startswith("shear_web_"):
        return True
    if key.startswith("spar_cap"):
        return True
    return False




def is_constructive_label(label: Any) -> bool:
    """True for any label integrated via constructive polygons (walls + cores)."""
    return is_wall_solid_label(label) or is_core_label(label)


def classify_property_labels(labels: Iterable[Any]) -> tuple[list[str], list[str]]:
    """Split section labels into (constructive_labels, residual_labels).

    ``constructive_labels`` includes wall solids and sandwich cores.
    ``residual_labels`` is LE/TE inserts on residual φ (and any unknown labels).
    """
    constructive: list[str] = []
    residual: list[str] = []
    for raw in labels:
        key = _label_key(raw)
        if is_constructive_label(key):
            constructive.append(key)
        elif is_residual_label(key):
            residual.append(key)
        else:
            residual.append(key)
    return constructive, residual
