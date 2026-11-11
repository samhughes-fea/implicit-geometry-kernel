"""``BladeConcept`` / ``SpanwiseScheme``: full-blade authoring compiled to loft kwargs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Union

import numpy as np

from implicit_geometry_kernel.blade.stations import PlyDrop, Station, Stringer, TipClosure, TransverseMember
from implicit_geometry_kernel.catalogue import SystemLayoutSpec, resolve_system_type
from implicit_geometry_kernel.section.layout.section_layout import SectionLayout

LayoutInput = Union[str, SectionLayout, SystemLayoutSpec]


def catalogue_ribs(stations: Sequence[Station], n_ribs: int) -> tuple[TransverseMember, ...]:
    """``n_ribs`` interior rib faces, evenly spaced between the root and the tip.

    The tip is closed by the blade's :class:`TipClosure`, which is not a rib, so
    a ``-10R`` layout has 10 ribs here plus the tip closure.
    """
    n = int(n_ribs)
    if n <= 0:
        return ()
    z0, z1 = float(stations[0].z), float(stations[-1].z)
    zs = np.linspace(z0, z1, n + 2)[1:-1]
    return tuple(TransverseMember(f"rib_{i + 1}", float(z), 0.012) for i, z in enumerate(zs))


def _n_ribs(layout: LayoutInput) -> int:
    n = getattr(layout, "n_ribs", None)
    if n is not None:
        return int(n)
    if isinstance(layout, str):
        return int(resolve_system_type(layout).n_ribs)
    return 0


@dataclass(frozen=True)
class SpanwiseScheme:
    """Typed list of spanwise extras already modeled on ``BladeGeometry``.

    ``tip_closure`` of ``None`` gives the blade the default :class:`TipClosure`.
    """

    transverse: tuple[TransverseMember, ...] = ()
    stringers: tuple[Stringer, ...] = ()
    ply_drops: tuple[PlyDrop, ...] = ()
    tip_closure: TipClosure | None = None


@dataclass(frozen=True)
class BladeConcept:
    """A catalogue token or parametric layout plus a spanwise scheme.

    W/X/Y are the station slice (``SystemLayoutSpec``). Z is web freeze. R is
    the number of interior rib faces, filled from the token when
    ``scheme.transverse`` is empty and ``to_geometry_kwargs(stations=...)`` is
    given. The tip closure is not counted in R.
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
        out: dict[str, object] = {
            "layout": self.layout,
            "transverse": transverse,
            "stringers": self.scheme.stringers,
            "ply_drops": self.scheme.ply_drops,
        }
        if self.scheme.tip_closure is not None:
            out["tip_closure"] = self.scheme.tip_closure
        return out

    def region_graph(self, **kwargs):
        """Compile this concept's layout to a keyed :class:`WallRegionGraph`."""
        from implicit_geometry_kernel.blade.topology import WallRegionGraph

        return WallRegionGraph.from_section_layout(self.layout, **kwargs)


__all__ = ["BladeConcept", "SpanwiseScheme", "catalogue_ribs"]
