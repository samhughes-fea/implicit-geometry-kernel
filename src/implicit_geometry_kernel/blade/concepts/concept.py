"""``BladeConcept`` / ``SpanwiseScheme``: full-blade authoring compiled to loft kwargs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Union

import numpy as np

from implicit_geometry_kernel.blade.stations import PlyDrop, Station, Stringer, TransverseMember
from implicit_geometry_kernel.catalogue import SystemLayoutSpec, resolve_system_type
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout

LayoutInput = Union[str, SectionLayout, SystemLayoutSpec]


def catalogue_ribs(stations: Sequence[Station], n_ribs: int) -> tuple[TransverseMember, ...]:
    """Evenly spaced rib faces from letter R, including a tip closure."""
    n = int(n_ribs)
    if n <= 0:
        return ()
    z0, z1 = float(stations[0].z), float(stations[-1].z)
    if abs(z1 - z0) < 1e-12:
        return (TransverseMember("tip_closure", z1, 0.012),)
    zs = np.linspace(z0, z1, n + 1)[1:]
    out = []
    for i, z in enumerate(zs):
        name = "tip_closure" if i == n - 1 else f"rib_{i + 1}"
        out.append(TransverseMember(name, float(z), 0.012))
    return tuple(out)


def _n_ribs(layout: LayoutInput) -> int:
    n = getattr(layout, "n_ribs", None)
    if n is not None:
        return int(n)
    if isinstance(layout, str):
        return int(resolve_system_type(layout).n_ribs)
    return 0


@dataclass(frozen=True)
class SpanwiseScheme:
    """Typed list of spanwise extras already modeled on ``BladeGeometry``."""

    transverse: tuple[TransverseMember, ...] = ()
    stringers: tuple[Stringer, ...] = ()
    ply_drops: tuple[PlyDrop, ...] = ()


@dataclass(frozen=True)
class BladeConcept:
    """A catalogue token or parametric layout plus a spanwise scheme.

    W/X/Y are the station slice (``SystemLayoutSpec``). Z is web freeze. R is
    transverse faces, filled from the token when ``scheme.transverse`` is empty
    and ``to_geometry_kwargs(stations=...)`` is given.
    """

    layout: LayoutInput
    scheme: SpanwiseScheme = SpanwiseScheme()

    @classmethod
    def from_key(cls, key: str, scheme: SpanwiseScheme | None = None) -> "BladeConcept":
        return cls(resolve_system_type(key), scheme or SpanwiseScheme())

    def to_geometry_kwargs(self, stations: Sequence[Station] | None = None) -> dict:
        transverse = self.scheme.transverse
        if stations is not None and not transverse:
            n = _n_ribs(self.layout)
            if n:
                transverse = catalogue_ribs(stations, n)
        return {
            "layout": self.layout,
            "transverse": transverse,
            "stringers": self.scheme.stringers,
            "ply_drops": self.scheme.ply_drops,
        }

    def region_graph(self, **kwargs):
        """Compile this concept's layout to a keyed :class:`WallRegionGraph`."""
        from implicit_geometry_kernel.blade.topology import WallRegionGraph

        return WallRegionGraph.from_section_layout(self.layout, **kwargs)


__all__ = ["BladeConcept", "SpanwiseScheme", "catalogue_ribs"]
