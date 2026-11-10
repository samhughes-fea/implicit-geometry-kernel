"""Split host skins at cap-footprint ends (loop-run uniqueness)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from implicit_geometry_kernel.contract.shell_midline_strip import ShellMidlineStrip
from implicit_geometry_kernel.section.midlines.ring_geometry import open_polyline_prefix_lengths


def split_skin_at_cap_footprints(
    strips: tuple[ShellMidlineStrip, ...],
    *,
    merge_s_tol_m: float | None = None,
) -> tuple[ShellMidlineStrip, ...]:
    """Split host skins at cap-footprint ends so loop run names stay unique.

    Web splits already isolate LE / bay / TE. A cap sitting *inside* a panel
    (one-web root/tip) would otherwise colour ``skin → cap → skin`` with the
    same wall label twice. Extra cuts are merged when they coincide with web
    feet (two-web 2D-F is a no-op).
    """
    from implicit_geometry_kernel.section.assembly.members import wall_cell_face
    from implicit_geometry_kernel.section.junctions.cap_skin_bonds import _chain
    from implicit_geometry_kernel.section.midlines.cap_skin_ray_resample import (
        _surface_side,
    )
    from implicit_geometry_kernel.section.midlines.skin_split import (
        _MERGE_S_TOL_FRAC,
        _PRESSURE,
        _SUCTION,
        _classify_side_arcs,
        _label_key,
        _mean_dist_to_polyline,
        _normalise_side,
        _s_on_polyline,
    )

    if not any(wall_cell_face(s.label) == "cap" and s.bonds for s in strips):
        return strips

    side_hosts: dict[str, list[tuple[int, ShellMidlineStrip]]] = {
        _SUCTION: [],
        _PRESSURE: [],
    }
    for i, s in enumerate(strips):
        if wall_cell_face(s.label) != "skin":
            continue
        side = _normalise_side(s.skin_side, label=_label_key(s.label))
        if side in side_hosts:
            side_hosts[side].append((i, s))

    new_skins: list[ShellMidlineStrip] = []
    for side in (_SUCTION, _PRESSURE):
        hosts = side_hosts[side]
        if not hosts:
            continue
        try:
            chain = _chain(hosts)
        except ValueError:
            new_skins.extend(s for _i, s in hosts)
            continue
        chunks = [np.asarray(arc.oriented, dtype=float) for arc in chain]
        pts = np.vstack([chunks[0]] + [p[1:] for p in chunks[1:]])
        L = float(open_polyline_prefix_lengths(pts)[-1]) if pts.shape[0] > 1 else 0.0
        tol = (
            float(merge_s_tol_m)
            if merge_s_tol_m is not None
            else _MERGE_S_TOL_FRAC * max(L, 1e-12)
        )
        joints: list[float] = []
        acc = 0.0
        for piece in chunks[:-1]:
            acc += float(open_polyline_prefix_lengths(piece)[-1])
            joints.append(acc)
        panels = [0.0] + joints + [L]
        extra: list[float] = []
        for st in strips:
            if wall_cell_face(st.label) != "cap":
                continue
            if _surface_side(st.surface) != _surface_side(side):
                continue
            for bond in st.bonds:
                hp = np.asarray(bond.host_points, dtype=float)
                if hp.shape[0] < 2:
                    continue
                s_lo, s_hi = sorted(
                    (_s_on_polyline(hp[0], pts), _s_on_polyline(hp[-1], pts))
                )
                min_remnant = max(tol, 0.01 * max(L, 1e-12), 1e-3)
                if (s_hi - s_lo) < min_remnant:
                    continue
                overlap = []
                for k in range(len(panels) - 1):
                    s0, s1 = float(panels[k]), float(panels[k + 1])
                    lo, hi = max(s_lo, s0), min(s_hi, s1)
                    if hi - lo > 0.0:
                        overlap.append((hi - lo, k, s0, s1))
                if not overlap:
                    continue
                _span, _k, s0, s1 = max(overlap, key=lambda item: item[0])
                if s_lo > s0 + min_remnant:
                    extra.append(s_lo)
                if s_hi < s1 - min_remnant:
                    extra.append(s_hi)
        if not extra:
            new_skins.extend(s for _i, s in hosts)
            continue
        classified = _classify_side_arcs(
            side,
            pts,
            float(hosts[0][1].thickness_m),
            joints + extra,
            merge_s_tol_m=merge_s_tol_m,
        )
        new_skins.extend(classified if classified else [s for _i, s in hosts])

    if len(new_skins) < 2:
        return strips

    remapped: list[ShellMidlineStrip] = []
    for st in strips:
        if wall_cell_face(st.label) == "skin":
            continue
        if wall_cell_face(st.label) != "cap" or not st.bonds:
            remapped.append(st)
            continue
        new_bonds = []
        for bond in st.bonds:
            hp = np.asarray(bond.host_points, dtype=float)
            best = min(
                new_skins,
                key=lambda sk: _mean_dist_to_polyline(
                    hp, np.asarray(sk.midline_twisted_pitch_axis_blade_frame, float)
                ),
            )
            new_bonds.append(replace(bond, host_label=best.label))
        remapped.append(replace(st, bonds=tuple(new_bonds)))
    return tuple(new_skins + remapped)


__all__ = ["split_skin_at_cap_footprints"]
