"""Planform CSV in and a windIO-shaped dict out so a geometry round-trips."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from implicit_geometry_kernel.blade.airfoil import Airfoil
from implicit_geometry_kernel.blade.loft import BladeGeometry
from implicit_geometry_kernel.blade.stations import Station


def read_planform_csv(path, airfoil: Airfoil | None = None, **blade_kw) -> BladeGeometry:
    """Read ``z,chord,twist,pitch_axis,offset_x,offset_y[,airfoil]`` rows.

    ``airfoil`` is used when the file has no airfoil column (NACA 4-digit codes
    or ellipse). Layout and thicknesses come from ``blade_kw``.
    ``pivot_fraction`` is required on ``blade_kw`` unless a pitch-axis column
    is present (then that column is passed through as ``pivot_fraction``).
    """
    text = Path(path).read_text(encoding="utf-8")
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.replace(";", ",").split(",")]
        if parts[0].lower() in {"z", "span"}:
            continue
        rows.append(parts)
    if len(rows) < 2:
        raise ValueError("planform CSV needs at least two stations")
    has_pitch_col = all(len(parts) > 3 for parts in rows)
    stations = []
    for parts in rows:
        z, chord = float(parts[0]), float(parts[1])
        twist = float(parts[2]) if len(parts) > 2 else 0.0
        pitch = float(parts[3]) if len(parts) > 3 else 0.25
        ox = float(parts[4]) if len(parts) > 4 else 0.0
        oy = float(parts[5]) if len(parts) > 5 else 0.0
        foil = airfoil
        if len(parts) > 6 and parts[6]:
            token = parts[6]
            if token.lower().startswith("ellipse"):
                foil = Airfoil.ellipse(float(token.split(":")[-1]) if ":" in token else 1.0)
            else:
                code = token.upper().replace("NACA", "").strip()
                foil = Airfoil.naca5(code) if len(code) == 5 else Airfoil.naca4(code)
        if foil is None:
            raise ValueError("give an airfoil column or the airfoil= argument")
        stations.append(Station(z, chord, foil, twist=twist, pitch_axis=pitch, offset=(ox, oy)))
    layout = blade_kw.pop("layout", "2D-F-LC-0R")
    skin = blade_kw.pop("skin_thickness", 0.004)
    if "pivot_fraction" not in blade_kw and "pitch_axis" not in blade_kw:
        if not has_pitch_col:
            raise TypeError("pivot_fraction is required (or a pitch-axis CSV column)")
        blade_kw["pivot_fraction"] = (
            [s.z for s in stations],
            [s.pitch_axis for s in stations],
        )
    return BladeGeometry(stations, layout, skin, **blade_kw)


def write_windio(geo: BladeGeometry) -> dict[str, Any]:
    """Minimal windIO 2.0-shaped outer_shape_bem from the loft stations."""
    zs = np.array([s.z for s in geo.stations], float)
    span = float(zs[-1] - zs[0]) or 1.0
    grid = ((zs - zs[0]) / span).tolist()
    chords, twists, pitches, ox, oy = [], [], [], [], []
    for z in zs:
        c, t, p, off, _ = geo.interpolate(float(z))
        chords.append(c)
        twists.append(t)
        pitches.append(p)
        ox.append(off[0])
        oy.append(off[1])
    return {
        "outer_shape_bem": {
            "reference_axis": {
                "x": {"grid": grid, "values": ox},
                "y": {"grid": grid, "values": oy},
                "z": {"grid": grid, "values": zs.tolist()},
            },
            "chord": {"grid": grid, "values": chords},
            "twist": {"grid": grid, "values": twists},
            "pitch_axis": {"grid": grid, "values": pitches},
            "airfoil_position": {
                "grid": grid,
                "labels": [s.airfoil.name for s in geo.stations],
            },
        }
    }
