"""Interpolate, OML polyline and section_at for BladeGeometry."""
from __future__ import annotations

import re

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.loft import BladeGeometry

from implicit_geometry_kernel.blade.cut import CellPolygon, SectionCut, SectionWall, _twist_offset
from implicit_geometry_kernel.blade.loft.loft_tables import _FORK_TWIST_SIGN, _match, _morph_surfaces
from implicit_geometry_kernel.contract.errors import InfeasibleSectionGeometry
from implicit_geometry_kernel.section.layout.thickness_limits import (
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M as _MIN_SKIN,
    MIN_REALISTIC_SPAR_WALL_THICKNESS_M as _MIN_CAP,
    MIN_REALISTIC_WEB_WALL_THICKNESS_M as _MIN_WEB,
)


def _gauge_values(value) -> tuple[float, ...]:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return (float(value),)
    return tuple(float(v) for v in value)


def _scale_gauge(value, k: float):
    vals = tuple(v * k for v in _gauge_values(value))
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return vals[0]
    return vals


def _min_thickness_scale(t_skin, t_cap, t_web) -> float:
    """Scale-up so manufacturing minima are met; the cut is placed back by ``1/k``."""
    k = max(
        1.0,
        *(
            lo / t
            for lo, values in ((_MIN_SKIN, t_skin), (_MIN_WEB, t_web), (_MIN_CAP, t_cap))
            for t in _gauge_values(values)
            if t > 0.0
        ),
    )
    if k > 1.0:
        k *= 1.0 + 1e-9
    return k


#: how far inside its outer-mould-line foot a two-foot web's axis point is placed, in skin
#: thicknesses along the web's line: past the skin mid-surface (half a skin thickness in
#: along the wall normal) for webs within about 75 degrees of the normal, and within the
#: section builder's search for the skin mid-surface from it (four skin thicknesses)
_FOOT_AXIS_DEPTH = 2.0


def _with_web_axes(layout, feet, x: np.ndarray, upper: np.ndarray, lower: np.ndarray,
                   chord: float, skin: float):
    """``layout`` with each web the line through the OML points at its two feet.

    ``feet`` are ``(suction, pressure)`` chord fractions; ``x``, ``upper``, ``lower`` the
    unit-chord surfaces of the section, ``chord`` its (scaled) chord and ``skin`` its
    (scaled) skin thickness. The section builder extends a web from its axis points to
    the skin mid-surface along the line, so each axis point is taken on the line just
    inside its foot (``_FOOT_AXIS_DEPTH`` skin thicknesses).
    """
    import dataclasses

    axes = []
    for f_su, f_pr in feet:
        top = np.array([f_su, float(np.interp(f_su, x, upper))]) * chord
        bot = np.array([f_pr, float(np.interp(f_pr, x, lower))]) * chord
        span = float(np.linalg.norm(top - bot))
        if span <= 2.0 * _FOOT_AXIS_DEPTH * skin:
            raise InfeasibleSectionGeometry(
                f"web feet {(f_su, f_pr)} are {span:.4g} m apart, too close for a "
                f"skin of {skin:.4g} m")
        u = (bot - top) / span
        a = top + _FOOT_AXIS_DEPTH * skin * u
        b = bot - _FOOT_AXIS_DEPTH * skin * u
        axes.append(((float(a[0]), float(a[1])), (float(b[0]), float(b[1]))))
    return dataclasses.replace(layout, web_axes_m=tuple(axes))


def _surface_crossing(frac: float, theta: float, x: np.ndarray, y: np.ndarray, upper: bool) -> float:
    """Chord fraction where the line through ``(frac, 0)`` at ``theta`` from +y meets the surface."""
    residual = (x - frac) * np.cos(theta) - y * np.sin(theta)
    best = None
    for i in range(len(x) - 1):
        a, b = float(residual[i]), float(residual[i + 1])
        if a == 0.0 or a * b < 0.0:
            den = a - b
            t = 0.0 if den == 0.0 else a / den
            t = float(np.clip(t, 0.0, 1.0))
            xx = float(x[i] + t * (x[i + 1] - x[i]))
            yy = float(y[i] + t * (y[i + 1] - y[i]))
            if best is None or (yy > best[1] if upper else yy < best[1]):
                best = (xx, yy)
    if best is None:
        return float(frac)
    return float(best[0])


def _tilted_feet(fracs, tilt_deg, x, upper, lower):
    """``(suction, pressure)`` chord fractions for webs tilted ``tilt_deg`` from the thickness direction."""
    feet = []
    for frac, deg in zip(fracs, tilt_deg):
        # Positive tilt walks the suction foot toward the trailing edge.
        theta = -float(np.radians(deg))
        if abs(theta) < 1e-12:
            feet.append((float(frac), float(frac)))
            continue
        su = _surface_crossing(float(frac), theta, x, upper, True)
        pr = _surface_crossing(float(frac), theta, x, lower, False)
        feet.append((su, pr))
    return tuple(feet)


class LoftCutMixin:
    def interpolate(self: BladeGeometry, z: float):
        """(chord, twist, pitch_axis, offset, (x, upper, lower)) interpolated at z.

        Circle-like parents morph with a smoothstep blend; otherwise surfaces follow
        ``Airfoil.interpolate_by_thickness``. Chord, twist and offset come from span
        tables when given, otherwise from the stations. Pitch axis is always
        ``pivot_fraction``.

        A blade built with :meth:`BladeGeometry.from_section` already carries the
        consumer's blend in its samples. Between those stations the surfaces and the
        scalars are linear, and ``_morph_surfaces`` is not called.
        """
        zs = np.array([s.z for s in self.stations])
        if z < zs[0] - 1e-12 * max(1.0, abs(zs[0])) or z > zs[-1] + 1e-12 * max(1.0, abs(zs[-1])):
            raise ValueError(f"z = {z} lies outside the stations ({zs[0]} to {zs[-1]})")
        k = int(np.clip(np.searchsorted(zs, z, side="right") - 1, 0, len(zs) - 2))
        a, b = self.stations[k], self.stations[k + 1]
        w = float(np.clip((z - a.z) / (b.z - a.z), 0.0, 1.0))
        samples = getattr(self, "_section_samples", None)
        if samples is not None:
            # The callable was already the consumer's blend. Lerp the samples.
            # Do not run _morph_surfaces on this path.
            from implicit_geometry_kernel.blade.loft.section_callable import lerp_sample

            blended = lerp_sample(samples[k], samples[k + 1], w, float(z))
            return (
                blended.chord,
                blended.twist,
                blended.pivot_fraction,
                blended.offset,
                (np.asarray(blended.x, float), np.asarray(blended.upper, float), np.asarray(blended.lower, float)),
            )
        x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, self.n_airfoil_points)))
        upper, lower = _morph_surfaces(a.airfoil, b.airfoil, w, x, getattr(self, "airfoil_morph", "auto"))
        zz = float(z)
        chord = self._span_scalar(self._chord, [s.chord for s in self.stations], zz)
        twist = self._span_scalar(self._twist, [s.twist for s in self.stations], zz)
        pitch = float(self._pitch(zz))
        if self._offset is not None:
            off = self._offset(zz)
            offset = (float(off[0]), float(off[1]))
        else:
            ox = self._span_scalar(None, [s.offset[0] for s in self.stations], zz)
            oy = self._span_scalar(None, [s.offset[1] for s in self.stations], zz)
            offset = (ox, oy)
        return (chord, twist, pitch, offset, (x, upper, lower))

    def _placement_twist(self: BladeGeometry, twist_deg: float) -> float:
        """Degrees passed to the nose-up placement. ``nose_down`` negates them."""
        twist = float(twist_deg)
        if getattr(self, "twist_sense", "nose_up") == "nose_down":
            return -twist
        return twist

    def _span_scalar(self: BladeGeometry, table, values, z: float) -> float:
        if table is not None:
            return float(table(z))
        zs = np.asarray([s.z for s in self.stations], dtype=np.float64)
        return float(np.interp(z, zs, np.asarray(values, dtype=np.float64)))

    def oml_polyline(self: BladeGeometry, z: float) -> np.ndarray:
        """Closed outer mould line in the blade frame at ``z``, shape ``(N, 3)``.

        Same interpolator and placement as ``section_at`` (pitch axis, twist, offset).
        """
        z = float(z)
        chord, twist, pitch, offset, (x, upper, lower) = self.interpolate(z)
        naca = np.vstack([np.column_stack([x[::-1], upper[::-1]]), np.column_stack([x[1:], lower[1:]])]) * chord
        if np.allclose(naca[0], naca[-1]):
            naca = naca[:-1]
        local = np.column_stack([naca[:, 0] - pitch * chord, naca[:, 1]])
        xy = _twist_offset(local, self._placement_twist(twist), offset)
        pts = np.column_stack([xy, np.full(len(xy), z)])
        return np.vstack([pts, pts[0]])

    # ------------------------------------------------------------ cutting
    def section_at(self: BladeGeometry, z: float) -> SectionCut:
        """The wall midlines, thicknesses, labels and cell loops of the section at z.

        Cuts are cached per ``round(z, 9)`` (LRU 256). ``CutWallMetaSDF`` and
        ``StackThickness`` use the same key; this is the G7 cost fix, not a BVH.
        """
        from implicit_geometry_kernel.frames import twisted_pitch_axis_blade_frame_to_naca
        from implicit_geometry_kernel.sdf.shapes.airfoil import AirfoilSDF
        from implicit_geometry_kernel.section.spec import SectionGeometry, SectionOptions, SectionSpec

        z = round(float(z), 9)
        hit = self._section_cache.get(z)
        if hit is not None:
            self._section_cache.move_to_end(z)
            return hit
        chord, twist, pitch, offset, (x, upper, lower) = self.interpolate(z)
        layout = self._layout_at(z)
        t_skin, t_cap, t_web = self._skin(z), (self._cap(z) if layout.has_spar_caps else 0.0), self._web(z)
        # section_geometry floors wall thicknesses at wind-turbine manufacturing minima
        # (2.5 mm skin, 4 mm web, 5 mm cap). The geometry is similar under scaling, so a
        # thinner-walled section is built k times larger and scaled back: exact, and the
        # frozen geometry code is used as it is.
        k = _min_thickness_scale(t_skin, t_cap, t_web)
        signed_twist = self._placement_twist(twist)
        twist_rad = _FORK_TWIST_SIGN * np.radians(signed_twist)
        pitch_m = self._pitch_metres_at(z)
        feet = self._web_feet_at(z)
        if pitch_m is not None and pitch_m and isinstance(pitch_m[0], tuple):
            feet = tuple(
                (float(pitch) + float(a) / chord, float(pitch) + float(b) / chord)
                for a, b in pitch_m
            )
        tilt = self._web_tilts_at(z)
        if tilt is not None and any(abs(float(v)) > 1e-12 for v in tilt):
            if feet is not None:
                raise ValueError(
                    "web_tilt_deg tilts a single chord fraction; two-foot webs already set the line"
                )
            fracs = tuple(float(f) for f in getattr(layout, "web_chord_fracs", ()) or ())
            if len(tilt) == 1 and len(fracs) != 1:
                tilt = tuple([tilt[0]] * len(fracs))
            if len(tilt) != len(fracs):
                raise ValueError(f"layout has {len(fracs)} webs; got {len(tilt)} tilt angles")
            feet = _tilted_feet(fracs, tilt, x, upper, lower)
        if feet is not None:
            from implicit_geometry_kernel.catalogue import SystemLayoutSpec
            from implicit_geometry_kernel.section.layout.section_layout import SectionLayout

            if isinstance(layout, SystemLayoutSpec):
                layout = SectionLayout.from_system_layout_spec(layout)
            layout = _with_web_axes(layout, feet, x, upper, lower, chord * k, t_skin * k)
        web_offsets = None
        if pitch_m is not None and not (pitch_m and isinstance(pitch_m[0], tuple)):
            web_offsets = tuple((float(xm) + pitch * chord) * k for xm in pitch_m)
        elif self.web_offsets is not None:
            web_offsets = tuple(v * k for v in self.web_offsets)
        elif self._web_x_from_pitch is not None:
            web_offsets = tuple((xm + pitch * chord) * k for xm in self._web_x_from_pitch)
        spec = SectionSpec.from_kwargs(
            layout,
            AirfoilSDF.from_surface_graphs(x, upper, lower, chord * k),
            skin_thickness=_scale_gauge(t_skin, k),
            spar_cap_thickness=_scale_gauge(t_cap, k),
            shear_web_thickness=_scale_gauge(t_web, k),
            twist_rad=twist_rad,
            pivot_fraction=pitch,
            web_offsets=web_offsets,
        )
        runner = SectionGeometry(
            spec, SectionOptions(outputs=("midlines", "cell_wall_loops"))
        ).execute()
        result = runner.get_results()
        resolution = runner.get_layout_resolution()

        def place(pts_B):
            pts = np.asarray(pts_B, float) / k
            naca = twisted_pitch_axis_blade_frame_to_naca(
                pts, _FORK_TWIST_SIGN * np.radians(signed_twist), chord, pitch,
            )
            local = np.column_stack([naca[:, 0] - pitch * chord, naca[:, 1]])
            xy = _twist_offset(local, signed_twist, offset)
            return np.column_stack([xy, np.full(len(xy), z)])

        from implicit_geometry_kernel.section.junctions.cap_skin_bonds import (
            bond_caps_to_skin,
        )
        from implicit_geometry_kernel.section.midlines.export import build_solid_cap_strips
        from implicit_geometry_kernel.section.midlines.skin_split import (
            split_skin_at_cap_footprints,
        )

        mesh_strips = tuple(result.midline_strips or ())
        cap_strips = build_solid_cap_strips(
            result.section,
            twist_rad=float(twist_rad),
            chord_m=float(chord * k),
            pivot_fraction=pitch,
        )
        # Solids path: cap recipe midlines + skin/web mesh strips, bonded for
        # graph partitions. Split host skins at footprint ends so cap colouring
        # cannot emit two loop runs with one wall name.
        if cap_strips:
            strips = list(
                split_skin_at_cap_footprints(
                    bond_caps_to_skin(mesh_strips + cap_strips)
                )
            )
        else:
            strips = list(mesh_strips)
        by_label = {s.label: s for s in strips}
        name = _lone_pad_names(by_label)
        walls = []
        for s in strips:
            pts = place(s.midline_twisted_pitch_axis_blade_frame)
            side = s.skin_side or s.surface
            side = None if side is None else side.replace("_surface", "")
            bonds = {}
            for bond in s.bonds:
                host = by_label.get(bond.host_label)
                if host is None:
                    continue
                bonds[name.get(bond.host_label, bond.host_label)] = (
                    _match(s.midline_twisted_pitch_axis_blade_frame, bond.points),
                    _match(host.midline_twisted_pitch_axis_blade_frame, bond.host_points),
                )
            label = name.get(s.label, s.label)
            walls.append(SectionWall(label, s.kind, pts, float(s.thickness_m) / k, side, bonds))
        cells = tuple(
            CellPolygon(int(c.id), place(c.polyline_twisted_pitch_axis_blade_frame), float(c.area_m2) / k**2,
                     tuple(name.get(w, w) for w in c.wall_labels))
            for c in (result.cell_wall_loops or ())
        )
        dropped = tuple(f"web at {d.x_m / k:.6g} m from the LE: {d.reason}" for d in resolution.dropped)
        dropped_idx = tuple(
            int(idx) for d in resolution.dropped if (idx := getattr(d, "index", None)) is not None
        )
        recipe = resolution.layout
        n_live = int(getattr(recipe, "n_webs", 0) or len(resolution.web_indices))
        fam = str(getattr(recipe, "structural_family", None) or "")
        cut = SectionCut(
            z, chord, signed_twist, pitch, (float(offset[0]), float(offset[1])), str(self.layout), tuple(walls), cells,
            dropped, self._transverse_at(z), self._stringers_at(z), result.section, k,
            dropped_web_indices=dropped_idx,
            station_webs=n_live,
            station_family=fam,
            tip_closure=self._tip_closure_at(z),
        )
        object.__setattr__(cut, "_wall_reference", getattr(self, "wall_reference", "mid"))
        if self.merge_cap_and_skin:
            from implicit_geometry_kernel.blade.shell_graph import shell_graph, shell_strips_from_graph
            from dataclasses import replace as _replace

            # Keep the graph built from the stacked walls. Rebuilding it from the
            # merged runs drops the caps: those runs have no skin bonds, and the
            # mid-surface only paints a cap where a bond says so (BG-040).
            partition_graph = shell_graph(cut)
            merged = []
            for strip in shell_strips_from_graph(cut):
                xy = np.asarray(strip.midline_twisted_pitch_axis_blade_frame, float)
                if xy.shape[1] == 2:
                    pts = np.column_stack([xy, np.full(len(xy), z)])
                else:
                    pts = xy
                stack_t = (
                    float(sum(strip.layer_thicknesses_m))
                    if strip.layer_thicknesses_m
                    else float(strip.thickness_m)
                )
                merged.append(SectionWall(strip.label, strip.kind, pts, stack_t, None, {}))
            cut = _replace(cut, walls=tuple(merged))
            object.__setattr__(cut, "_partition_graph", partition_graph)
        self._section_cache[z] = cut
        while len(self._section_cache) > self._section_cache_size:
            self._section_cache.popitem(last=False)
        return cut


_PAD = {side: re.compile(rf"^spar_cap_{side}_surface_(\d+)$") for side in ("suction", "pressure")}


def _lone_pad_names(by_label) -> dict[str, str]:
    """Publish a lone indexed cap pad under its surface's canonical label.

    Family C stores pads as ``spar_cap_<side>_surface_<i>``. When a surface has
    one pad and no canonical cap, that pad is the surface's cap, so it keeps the
    name a family-D box uses at neighbouring stations and loft keys stay stable.
    """
    out: dict[str, str] = {}
    for side, pattern in _PAD.items():
        canon = f"spar_cap_{side}_surface"
        pads = [str(lb) for lb in by_label if pattern.match(str(lb))]
        if len(pads) == 1 and canon not in by_label:
            out[pads[0]] = canon
    return out
