"""Gallery hosts go through ``SectionGeometry``, not ``MultiCellSection`` constructors.

The worker still returns a ``SectionMidlineHost`` (usually ``MultiCellSection``)
on ``result.section``. Catalogue keys plus ``web_offsets`` are the public door.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from implicit_geometry_kernel.catalogue import resolve_system_type
from implicit_geometry_kernel.section import SectionGeometry


def layout_token(
    n_webs: int,
    family: str,
    *,
    alignment: str = "chord_normal",
) -> str:
    """WX-Y-Z-NR token for a station gallery case."""
    fam = str(family).upper()
    n = int(n_webs)
    if n == 0:
        return f"0{fam}-LC-0R"
    y = "F" if str(alignment) == "flapwise" else "CN"
    return f"{n}{fam}-{y}-LC-0R"


def _cap_m(cap_height: Any, spar_cap_thickness: float | None) -> float:
    if spar_cap_thickness is not None:
        return float(spar_cap_thickness)
    if cap_height is None:
        return 0.012
    if isinstance(cap_height, (int, float)):
        return float(cap_height)
    return float(cap_height[0])


def assemble_host(
    airfoil: Any,
    *,
    layout: Any | None = None,
    n_webs: int | None = None,
    family: str = "D",
    alignment: str = "chord_normal",
    web_offsets: Sequence[float] | None = None,
    web_fracs: Sequence[float] | None = None,
    skin_thickness: float = 0.003,
    spar_cap_thickness: float | None = None,
    cap_height: Any = None,
    shear_web_thickness: Any = 0.004,
    twist_rad: float = 0.0,
    outputs: Sequence[str] = (),
):
    """Assemble a plot/property host through ``SectionGeometry.execute()``."""
    xs: tuple[float, ...] | None = None
    if web_offsets is not None:
        xs = tuple(float(x) for x in web_offsets)
    elif web_fracs is not None:
        chord = float(airfoil.chord)
        xs = tuple(float(f) * chord for f in web_fracs)
    if layout is None:
        n = n_webs
        if n is None:
            n = 0 if xs is None else len(xs)
        layout = resolve_system_type(layout_token(int(n), family, alignment=alignment))
    cap = _cap_m(cap_height, spar_cap_thickness)
    return (
        SectionGeometry(
            layout=layout,
            airfoil=airfoil,
            skin_thickness=skin_thickness,
            spar_cap_thickness=cap,
            shear_web_thickness=shear_web_thickness,
            twist_rad=float(twist_rad),
            web_offsets=xs,
            outputs=tuple(outputs),
        )
        .execute()
        .get_results()
        .section
    )
