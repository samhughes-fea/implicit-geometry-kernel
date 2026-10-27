"""SDF subcomponent panel + taxonomy caption helpers (geometry-owned)."""

from __future__ import annotations

from typing import Any

import numpy as np


def taxonomy_plot_caption(layout_key: str, layout: Any) -> str:
    """Human-readable X/Y/Z caption for a registry layout (``layout_key`` for display)."""
    del layout_key  # kept for call-site symmetry with figure titles
    y = str(layout.structural_family)
    y_desc = {
        "A": "no spar caps",
        "B": "single cap band (pitch / max-t)",
        "C": "discrete caps per web",
        "D": "continuous box (first web to last web)",
    }.get(y, y)
    z = (
        "flapwise (B-frame vertical webs)"
        if layout.web_orientation == "flapwise"
        else "chord-normal (S-frame vertical webs)"
    )
    parts = [
        f"X = {layout.n_webs} web(s)",
        f"Y = {y} ({y_desc})",
        f"Z = {z}",
    ]
    return " | ".join(parts)


def _section_query_from_blade(
    pts_b: np.ndarray,
    *,
    twist_rad: float,
    chord_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Map canonical B ``(X_B, Y_B)`` → coordinates expected by twisted section SDFs.

    Delegates to :func:`plot_section.section_query_from_blade`.
    """
    from blade_geometry_model.section_geometry.visualisation.plot_section import (
        section_query_from_blade,
    )

    return section_query_from_blade(
        pts_b, twist_rad=float(twist_rad), chord_m=float(chord_m)
    )


def draw_section_sdf_panel(ax_geom: Any, section: Any, airfoil: Any, twist_rad: float) -> None:
    """Fill/contour each section label's SDF in **canonical B**, for midline overlay.

    Contours are drawn with horizontal = ``Y_B``, vertical = ``X_B`` (pair with
    :func:`apply_structural_section_axes` for LE-on-RHS).
    """
    from matplotlib.patches import Patch

    from blade_geometry_model.coordinate_frames import (
        naca_to_twisted_pitch_axis_blade_frame,
    )

    kind_colour = {"skin": "#7eb8da", "web": "#e74c3c", "cap": "#f1c40f"}
    kind_alpha = {"skin": 0.30, "web": 0.55, "cap": 0.55}

    chord_m = float(getattr(airfoil, "chord", 0.0) or 0.0)
    if not np.isfinite(chord_m) or chord_m <= 0.0:
        raise ValueError("draw_section_sdf_panel requires airfoil.chord > 0")

    verts_naca = np.asarray(airfoil.vertices, dtype=float)
    verts_b = naca_to_twisted_pitch_axis_blade_frame(
        verts_naca, float(twist_rad), chord_m
    )
    pad = 0.05 * float(max(np.ptp(verts_b[:, 0]), np.ptp(verts_b[:, 1]), 1e-9))
    xb_min = float(verts_b[:, 0].min()) - pad
    xb_max = float(verts_b[:, 0].max()) + pad
    yb_min = float(verts_b[:, 1].min()) - pad
    yb_max = float(verts_b[:, 1].max()) + pad

    xb = np.linspace(xb_min, xb_max, 360)
    yb = np.linspace(yb_min, yb_max, 600)
    XB, YB = np.meshgrid(xb, yb)
    pts_b = np.column_stack([XB.ravel(), YB.ravel()])
    qx, qy = _section_query_from_blade(
        pts_b, twist_rad=float(twist_rad), chord_m=chord_m
    )
    QX = qx.reshape(XB.shape)
    QY = qy.reshape(XB.shape)

    label_to_kind: dict[str, str] = {}
    for label in section.labels:
        s = str(label).lower()
        if "skin" in s:
            label_to_kind[label] = "skin"
        elif "web" in s:
            label_to_kind[label] = "web"
        elif "cap" in s:
            label_to_kind[label] = "cap"
    n_sub = len(list(getattr(section, "labels", []) or []))
    for label, kind in label_to_kind.items():
        try:
            phi = np.asarray(section[label](QX, QY), dtype=float)
        except Exception:
            continue
        col = kind_colour[kind]
        # Display: horizontal = Y_B, vertical = X_B.
        ax_geom.contourf(
            YB, XB, phi, levels=[-1e9, 0.0], colors=[col], alpha=kind_alpha[kind], zorder=1,
        )
        ax_geom.contour(
            YB, XB, phi, levels=[0.0], colors=[col], linewidths=0.8, alpha=0.9, zorder=2,
        )

    ax_geom.legend(
        handles=[
            Patch(facecolor=kind_colour["skin"], alpha=kind_alpha["skin"], label="skin SDF"),
            Patch(facecolor=kind_colour["cap"], alpha=kind_alpha["cap"], label="cap SDF"),
            Patch(facecolor=kind_colour["web"], alpha=kind_alpha["web"], label="web SDF"),
        ],
        loc="lower right",
        fontsize=8,
        facecolor="#1a1a2e",
        edgecolor="#2a2a3a",
        labelcolor="#e0e0e0",
        ncol=1,
    )
    ax_geom.set_title(
        f"Subcomponents (SDF)\n{n_sub} subcomponents · coloured by kind",
        color="#e0e0e0",
        fontsize=10,
    )


__all__ = [
    "draw_section_sdf_panel",
    "taxonomy_plot_caption",
]
