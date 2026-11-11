"""Interpolate, OML polyline and section_at for BladeGeometry."""
from __future__ import annotations

import re

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from implicit_geometry_kernel.blade.loft import BladeGeometry

from implicit_geometry_kernel.blade.cut import CellPolygon, SectionCut, SectionWall, _twist_offset
from implicit_geometry_kernel.blade.loft.loft_tables import _FORK_TWIST_SIGN, _match, _morph_surfaces
from implicit_geometry_kernel.section.layout.thickness_limits import (
    MIN_REALISTIC_SKIN_WALL_THICKNESS_M as _MIN_SKIN,
    MIN_REALISTIC_SPAR_WALL_THICKNESS_M as _MIN_CAP,
    MIN_REALISTIC_WEB_WALL_THICKNESS_M as _MIN_WEB,
)


def _min_thickness_scale(t_skin: float, t_cap: float, t_web: float) -> float:
    """Scale-up so manufacturing minima are met; the cut is placed back by ``1/k``."""
    k = max(
        1.0,
        *(
            lo / t
            for lo, t in ((_MIN_SKIN, t_skin), (_MIN_WEB, t_web), (_MIN_CAP, t_cap))
            if t > 0.0
        ),
    )
    if k > 1.0:
        k *= 1.0 + 1e-9
    return k


class LoftCutMixin:
    def interpolate(self: BladeGeometry, z: float):
        """(chord, twist, pitch_axis, offset, (x, upper, lower)) interpolated at z.

        Circle-like parents morph with a smoothstep blend; otherwise surfaces follow
        ``Airfoil.interpolate_by_thickness``. Chord, twist and offset come from span
        tables when given, otherwise from the stations. Pitch axis is always
        ``pivot_fraction``.
        """
        zs = np.array([s.z for s in self.stations])
        if z < zs[0] - 1e-12 * max(1.0, abs(zs[0])) or z > zs[-1] + 1e-12 * max(1.0, abs(zs[-1])):
            raise ValueError(f"z = {z} lies outside the stations ({zs[0]} to {zs[-1]})")
        k = int(np.clip(np.searchsorted(zs, z, side="right") - 1, 0, len(zs) - 2))
        a, b = self.stations[k], self.stations[k + 1]
        w = float(np.clip((z - a.z) / (b.z - a.z), 0.0, 1.0))
        x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, self.n_airfoil_points)))
        upper, lower = _morph_surfaces(a.airfoil, b.airfoil, w, x)
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
        xy = _twist_offset(local, twist, offset)
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
        twist_rad = _FORK_TWIST_SIGN * np.radians(twist)
        web_offsets = None
        if self.web_offsets is not None:
            web_offsets = tuple(v * k for v in self.web_offsets)
        elif self._web_x_from_pitch is not None:
            web_offsets = tuple((x + pitch * chord) * k for x in self._web_x_from_pitch)
        spec = SectionSpec.from_kwargs(
            layout,
            AirfoilSDF.from_surface_graphs(x, upper, lower, chord * k),
            skin_thickness=t_skin * k,
            spar_cap_thickness=t_cap * k,
            shear_web_thickness=t_web * k,
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
                pts, _FORK_TWIST_SIGN * np.radians(twist), chord, pitch,
            )
            local = np.column_stack([naca[:, 0] - pitch * chord, naca[:, 1]])
            xy = _twist_offset(local, twist, offset)
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
            z, chord, twist, pitch, (float(offset[0]), float(offset[1])), str(self.layout), tuple(walls), cells,
            dropped, self._transverse_at(z), self._stringers_at(z), result.section, k,
            dropped_web_indices=dropped_idx,
            station_webs=n_live,
            station_family=fam,
        )
        if self.merge_cap_and_skin:
            from implicit_geometry_kernel.blade.shell_graph import shell_strips_from_graph
            from dataclasses import replace as _replace

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
