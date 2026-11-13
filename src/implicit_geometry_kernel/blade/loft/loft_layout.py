"""Web-anchor and station-layout helpers for BladeGeometry."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace

import numpy as np

from implicit_geometry_kernel.blade.concepts.web_anchoring import chord_fraction_from_pitch
from implicit_geometry_kernel.blade.loft.loft_tables import _metre_rows, _with_web_fractions
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

    def _pitch_metres_at(self, z: float):
        fn = getattr(self, "_web_x_from_pitch_fn", None)
        if fn is None:
            return None
        return _metre_rows(fn(float(z)), "web_x_from_pitch")

    def _cap_metres_at(self, z: float, chord: float, pitch: float):
        fn = getattr(self, "_cap_x_from_pitch_fn", None)
        if fn is None:
            return None
        raw = _metre_rows(fn(float(z)), "cap_x_from_pitch")
        if raw is None:
            return None
        if raw and not isinstance(raw[0], tuple):
            raise ValueError("cap_x_from_pitch needs (x0, x1) pairs")
        return tuple(
            (
                chord_fraction_from_pitch(a, chord, pitch),
                chord_fraction_from_pitch(b, chord, pitch),
            )
            for a, b in raw
        )

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

    def _tip_closure_at(self, z: float):
        """The tip closure, with its thickness resolved, on the tip station only."""
        from dataclasses import replace as _replace

        if round(float(z), 9) != round(float(self.tip_closure_z), 9):
            return None
        return _replace(self.tip_closure, thickness=self.tip_closure_thickness())

    def _stringers_at(self, z: float) -> tuple[Stringer, ...]:
        key = round(float(z), 9)
        out = []
        for stringer in self.stringers:
            if stringer.extent is None:
                out.append(self._stringer_on_chord(stringer, z))
                continue
            z0, z1 = round(float(stringer.extent[0]), 9), round(float(stringer.extent[1]), 9)
            if z0 <= key <= z1:
                out.append(self._stringer_on_chord(stringer, z))
        return tuple(out)

    def _stringer_on_chord(self, stringer: Stringer, z: float) -> Stringer:
        if getattr(stringer, "position_m", None) is None:
            return stringer
        chord = max(self._span_scalar(self._chord, [s.chord for s in self.stations], float(z)), 1e-12)
        pitch = float(self._pitch(float(z)))
        frac = chord_fraction_from_pitch(stringer.position_m, chord, pitch)
        return replace(stringer, position=float(frac))

    def _web_feet_at(self, z: float):
        """Two-foot webs at z: ``((suction, pressure), ...)`` chord fractions, LE to TE.

        ``None`` when ``web_positions`` gives one fraction per web (a straight line
        in the layout's orientation). A web given by two feet is the line through the
        outer-mould-line points at those chord fractions on each side.
        """
        if self._web_positions_fn is None:
            return None
        raw = self._web_positions_fn(float(z))
        if raw is None or not len(raw):
            return None
        rows = [r for r in raw]
        if not all(isinstance(r, (tuple, list, np.ndarray)) and np.size(r) == 2 for r in rows):
            if any(isinstance(r, (tuple, list, np.ndarray)) for r in rows):
                raise ValueError("web_positions entries must all be fractions or all "
                                 "(suction, pressure) pairs")
            return None
        feet = sorted(((float(r[0]), float(r[1])) for r in rows), key=lambda p: p[0] + p[1])
        for a, b in feet:
            if not (0.0 < a < 1.0 and 0.0 < b < 1.0):
                raise ValueError(f"web feet must be chord fractions inside (0, 1); got {(a, b)}")
        return tuple(feet)

    def _outer_fraction_to_midplane(self, z: float, frac: float, side: str) -> float:
        """Chord fraction on the cap mid-plane whose outer foot is ``frac`` on this side."""
        chord, _, _, _, (x, upper, lower) = self.interpolate(float(z))
        y_of = upper if side == "suction" else lower
        cap = self._cap(float(z))
        t_cap = float(cap[0] if isinstance(cap, (tuple, list)) else cap)
        inset = float(self._skin(float(z))) + 0.5 * t_cap
        eps = 1e-3
        f0 = float(np.clip(frac - eps, 0.0, 1.0))
        f1 = float(np.clip(frac + eps, 0.0, 1.0))
        p0 = np.array([f0, float(np.interp(f0, x, y_of))]) * chord
        p1 = np.array([f1, float(np.interp(f1, x, y_of))]) * chord
        p = np.array([float(frac), float(np.interp(frac, x, y_of))]) * chord
        tangent = p1 - p0
        norm = float(np.linalg.norm(tangent))
        if norm <= 1e-15:
            return float(frac)
        tangent /= norm
        normal = np.array([-tangent[1], tangent[0]])
        if float(np.dot(normal, np.array([0.5 * chord, 0.0]) - p)) < 0.0:
            normal = -normal
        landed = p + normal * inset
        return float(landed[0] / chord)

    def _cap_windows_from_outer_surface(self, z: float, windows):
        def pair(window, side):
            lo, hi = float(window[0]), float(window[1])
            return (
                self._outer_fraction_to_midplane(z, lo, side),
                self._outer_fraction_to_midplane(z, hi, side),
            )

        if isinstance(windows, Mapping):
            return {
                "suction": pair(windows["suction"], "suction"),
                "pressure": pair(windows["pressure"], "pressure"),
            }
        raw = tuple(windows)
        if raw and not isinstance(raw[0], (tuple, list)):
            raw = (tuple(float(v) for v in raw),)
        suction = tuple(pair(w, "suction") for w in raw)
        pressure = tuple(pair(w, "pressure") for w in raw)
        if len(suction) == 1:
            return {"suction": suction[0], "pressure": pressure[0]}
        return suction

    def _layout_at(self, z: float):
        from implicit_geometry_kernel.section.layout.cap_placement import CapPlacementRule
        from implicit_geometry_kernel.section.layout.section_layout import SectionLayout
        from implicit_geometry_kernel.catalogue import SystemLayoutSpec

        spec = self.spec
        if (self._web_positions_fn is None and self._cap_windows_fn is None
                and getattr(self, "_cap_x_from_pitch_fn", None) is None):
            return spec
        if isinstance(spec, SystemLayoutSpec):
            layout = SectionLayout.from_system_layout_spec(spec)
        else:
            layout = spec
        fracs = None
        if self._web_positions_fn is not None:
            feet = self._web_feet_at(z)
            if feet is not None:
                fracs = tuple(0.5 * (a + b) for a, b in feet)
            else:
                raw = self._web_positions_fn(float(z))
                fracs = tuple(float(f) for f in np.atleast_1d(raw).reshape(-1))
        if fracs is not None:
            layout = _with_web_fractions(layout, fracs)
            if not fracs and layout.cap_placement.kind in ("continuous_box", "per_web_pads"):
                # no webs live here (a file whose webs start outboard of the root): a
                # web-pad or box rule has nothing to place caps by
                layout = replace(layout, cap_placement=CapPlacementRule.none(),
                                 structural_family="A", catalogue_key=None,
                                 flags=tuple(f for f in layout.flags if f != "has_spar_caps"))
        chord = max(self._span_scalar(self._chord, [s.chord for s in self.stations], float(z)), 1e-12)
        pitch = float(self._pitch(float(z)))
        windows = self._cap_metres_at(z, chord, pitch)
        if windows is None and self._cap_windows_fn is not None:
            windows = self._cap_windows_fn(float(z))
        if windows is not None and getattr(self, "cap_edge_reference", "midplane") == "outer":
            windows = self._cap_windows_from_outer_surface(z, windows)
        if windows is not None:
            if isinstance(windows, Mapping):
                rule = CapPlacementRule.from_side_windows(windows["suction"],
                                                          windows["pressure"])
            else:
                raw = tuple(windows) if not isinstance(windows, tuple) else windows
                if raw and not isinstance(raw[0], (tuple, list)):
                    raw = (tuple(float(x) for x in raw),)
                rule = CapPlacementRule.from_chord_windows(raw)
            layout = SectionLayout(
                web_chord_fracs=layout.web_chord_fracs,
                orientation=layout.orientation,
                cap_placement=rule,
                catalogue_key=layout.catalogue_key,
                web_stations_m=layout.web_stations_m,
                web_axes_m=layout.web_axes_m,
                structural_family=layout.structural_family,
                flags=layout.flags,
            )
        return layout

