"""Web-anchor and station-layout helpers for BladeGeometry."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from implicit_geometry_kernel.blade.loft.loft_tables import _with_web_fractions
from implicit_geometry_kernel.blade.stations import Stringer, TransverseMember


class LoftLayoutMixin:
    def _web_chord_frac_law(self, index: int):
        """Chord-fraction law for loft identity: LC catalogue, or PA/LE metres."""
        if self._web_x_from_pitch is not None:
            x_pa = float(self._web_x_from_pitch[int(index)])

            def frac(z, _x=x_pa):
                zz = float(z)
                chord = max(self._span_scalar(self._chord, [s.chord for s in self.stations], zz), 1e-12)
                return float(self._pitch(zz)) + _x / chord

            return frac
        if self.web_offsets is not None:
            x_le = float(self.web_offsets[int(index)])

            def frac(z, _x=x_le):
                zz = float(z)
                chord = max(self._span_scalar(self._chord, [s.chord for s in self.stations], zz), 1e-12)
                return _x / chord

            return frac
        fracs = tuple(float(f) for f in getattr(self.spec, "web_chord_fracs", ()) or ())
        f = float(fracs[int(index)])
        return lambda z, _f=f: _f

    def _with_anchored_web_laws(self, graph):
        """PA/LE loft anchors follow the metre freeze, not catalogue LC fractions."""
        if self._web_x_from_pitch is None and self.web_offsets is None:
            return graph
        webs = tuple(replace(w, frac=self._web_chord_frac_law(w.index)) for w in graph.webs)
        pad = float(getattr(graph, "cap_pad", None) or 0.0)
        anchors = []
        for a in graph.anchors:
            if a.name.startswith("web_"):
                try:
                    i = int(a.name.split("_", 1)[1])
                except ValueError:
                    anchors.append(a)
                    continue
                anchors.append(replace(a, frac=self._web_chord_frac_law(i)))
            elif pad > 0.0:
                parts = str(a.name).split("_")
                if (len(parts) >= 3 and parts[0] == "cap" and parts[1].isdigit()
                        and parts[2] in ("lo", "hi")):
                    i = int(parts[1])
                    law = self._web_chord_frac_law(i)
                    sign = -1.0 if parts[2] == "lo" else 1.0
                    anchors.append(replace(
                        a,
                        frac=lambda z, _law=law, _p=pad, _s=sign: float(
                            np.clip(_law(z) + _s * _p, 0.02, 0.98)
                        ),
                    ))
                    continue
                anchors.append(a)
            else:
                anchors.append(a)
        return replace(graph, webs=webs, anchors=tuple(anchors))

    def _transverse_at(self, z: float) -> tuple[TransverseMember, ...]:
        key = round(float(z), 9)
        return tuple(m for m in self.transverse if round(float(m.z), 9) == key)

    def _stringers_at(self, z: float) -> tuple[Stringer, ...]:
        key = round(float(z), 9)
        out = []
        for stringer in self.stringers:
            if stringer.extent is None:
                out.append(stringer)
                continue
            z0, z1 = round(float(stringer.extent[0]), 9), round(float(stringer.extent[1]), 9)
            if z0 <= key <= z1:
                out.append(stringer)
        return tuple(out)

    def _layout_at(self, z: float):
        from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
        from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
        from implicit_geometry_kernel.catalogue import SystemLayoutSpec

        spec = self.spec
        if self._web_positions_fn is None and self._cap_windows_fn is None:
            return spec
        if isinstance(spec, SystemLayoutSpec):
            layout = SectionLayout.from_system_layout_spec(spec)
        else:
            layout = spec
        fracs = None
        if self._web_positions_fn is not None:
            raw = self._web_positions_fn(float(z))
            fracs = tuple(float(f) for f in np.atleast_1d(raw).reshape(-1))
        if fracs is not None:
            layout = _with_web_fractions(layout, fracs)
        if self._cap_windows_fn is not None:
            windows = self._cap_windows_fn(float(z))
            if windows is not None:
                raw = tuple(windows) if not isinstance(windows, tuple) else windows
                if raw and not isinstance(raw[0], (tuple, list)):
                    raw = (tuple(float(x) for x in raw),)
                layout = SectionLayout(
                    web_chord_fracs=layout.web_chord_fracs,
                    orientation=layout.orientation,
                    cap_placement=CapPlacementRule.from_chord_windows(raw),
                    catalogue_key=layout.catalogue_key,
                    web_stations_m=layout.web_stations_m,
                    web_axes_m=layout.web_axes_m,
                    structural_family=layout.structural_family,
                    flags=layout.flags,
                )
        return layout

