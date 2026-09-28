"""Split outer skin midlines into hierarchical suction/pressure × LE/bay/TE arcs."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from numpy.typing import NDArray

from blade_structure_utilities.geometry_utilities.section_contract.shell_midline_strip import (
    ShellMidlineStrip,
)
from blade_structure_utilities.geometry_utilities.section_geometry.recipes.sdf_ir.outer_skin import (
    split_closed_ring_ss_ps,
)

from .ring_geometry import (
    merge_sorted_s_values,
    open_polyline_prefix_lengths,
    point_to_open_polyline_nearest,
)

#: Split positions closer than this fraction of a skin half-ring's arc length
#: are treated as one cut. Chosen to reproduce the historic absolute 1e-4 m
#: exactly at the 1 m reference chord, where a half-ring is ~1.0 m, while
#: scaling with the section: everything else in the junction path is
#: chord-relative, and an absolute tolerance loses headroom against the minimum
#: web spacing (``min_web_gap_frac`` = 2 % of chord) as the chord shrinks —
#: 200x at 1.0 m, 48x at the 0.241 m tip chord, 10x at 0.05 m. Callers that
#: need the old behaviour still pass ``merge_s_tol_m`` explicitly.
_MERGE_S_TOL_FRAC: float = 1e-4

_SUCTION = "suction_surface"
_PRESSURE = "pressure_surface"
_LEADING = "leading_edge"
_TRAILING = "trailing_edge"


def _label_key(label: object) -> str:
    return str(getattr(label, "value", label))


def _normalise_side(side: str | None, *, label: str = "", pts: NDArray | None = None) -> str | None:
    if side is not None:
        s = str(side).strip().lower()
        if s in (_SUCTION, "ss", "upper"):
            return _SUCTION
        if s in (_PRESSURE, "ps", "lower"):
            return _PRESSURE
    lab = label.lower()
    if "suction" in lab or lab.endswith("_ss") or ":ss" in lab:
        return _SUCTION
    if "pressure" in lab or lab.endswith("_ps") or ":ps" in lab:
        return _PRESSURE
    if pts is not None and pts.size:
        # Canonical B: X_B is flapwise and the suction surface is on +X_B. This
        # read Y_B (chordwise) before SG-015, i.e. it classified by LE/TE half.
        return _SUCTION if float(np.mean(pts[:, 0])) >= 0.0 else _PRESSURE
    return None


def _open_arc_between(
    pts: NDArray[np.float64],
    prefix: NDArray[np.float64],
    s_lo: float,
    s_hi: float,
) -> NDArray[np.float64]:
    """Extract open polyline samples for ``s_lo <= s <= s_hi`` along ``pts``."""
    if s_hi <= s_lo + 1e-15:
        return np.zeros((0, 2), dtype=float)
    n = int(pts.shape[0])
    out: list[NDArray[np.float64]] = []

    def _at(s: float) -> NDArray[np.float64]:
        if s <= float(prefix[0]) + 1e-15:
            return pts[0].copy()
        if s >= float(prefix[-1]) - 1e-15:
            return pts[-1].copy()
        j = int(np.searchsorted(prefix, s, side="right") - 1)
        j = int(np.clip(j, 0, n - 2))
        seg = float(prefix[j + 1] - prefix[j])
        t = 0.0 if seg < 1e-30 else (s - float(prefix[j])) / seg
        return (1.0 - t) * pts[j] + t * pts[j + 1]

    out.append(_at(s_lo))
    for j in range(1, n - 1):
        sj = float(prefix[j])
        if s_lo + 1e-9 < sj < s_hi - 1e-9:
            out.append(pts[j].copy())
    out.append(_at(s_hi))
    return np.stack(out, axis=0)


def _ensure_ss_ps_skins(
    skins: list[ShellMidlineStrip],
) -> list[ShellMidlineStrip]:
    """Normalise to two open suction/pressure half-rings."""
    if not skins:
        return []
    if len(skins) == 1 and skins[0].closed:
        ring = np.asarray(skins[0].midline_twisted_pitch_axis_blade_frame, dtype=float)
        # split_closed_ring_ss_ps reasons chord-first (LE = min x, suction = +y);
        # these strips are canonical B (X_B flap, Y_B chord), so transpose around
        # the call (SG-015). OuterSkin calls it in its chord-first authoring frame.
        ss_chord_first, ps_chord_first = split_closed_ring_ss_ps(ring[:, ::-1])
        ss, ps = ss_chord_first[:, ::-1].copy(), ps_chord_first[:, ::-1].copy()
        t = float(skins[0].thickness_m)
        return [
            ShellMidlineStrip(
                label="skin_suction_surface",
                kind="skin",
                midline_twisted_pitch_axis_blade_frame=ss,
                thickness_m=t,
                closed=False,
                skin_side=_SUCTION,
            ),
            ShellMidlineStrip(
                label="skin_pressure_surface",
                kind="skin",
                midline_twisted_pitch_axis_blade_frame=ps,
                thickness_m=t,
                closed=False,
                skin_side=_PRESSURE,
            ),
        ]
    out: list[ShellMidlineStrip] = []
    for s in skins:
        lab = _label_key(s.label)
        side = _normalise_side(
            s.skin_side,
            label=lab,
            pts=np.asarray(s.midline_twisted_pitch_axis_blade_frame, dtype=float),
        )
        if side is None:
            side = _SUCTION
        long_label = (
            "skin_suction_surface" if side == _SUCTION else "skin_pressure_surface"
        )
        out.append(
            replace(
                s,
                closed=False,
                skin_side=side,
                label=long_label,
            )
        )
    return out


def _classify_side_arcs(
    side: str,
    pts: NDArray[np.float64],
    thickness_m: float,
    split_s: list[float],
    *,
    merge_s_tol_m: float | None,
) -> list[ShellMidlineStrip]:
    """Split one open half-ring at web-foot arc parameters; name LE/bay/TE."""
    prefix = open_polyline_prefix_lengths(pts)
    L = float(prefix[-1]) if prefix.size else 0.0
    if L < 1e-12 or pts.shape[0] < 2:
        return []
    if merge_s_tol_m is None:
        merge_s_tol_m = _MERGE_S_TOL_FRAC * L

    cuts = merge_sorted_s_values([0.0, L] + list(split_s), merge_s_tol_m)
    cuts = [0.0] + [c for c in cuts if 1e-9 < c < L - 1e-9] + [L]
    cuts = merge_sorted_s_values(cuts, merge_s_tol_m)
    if len(cuts) < 2:
        return []

    segs: list[tuple[float, float, NDArray[np.float64]]] = []
    for i in range(len(cuts) - 1):
        s0, s1 = float(cuts[i]), float(cuts[i + 1])
        arc = _open_arc_between(pts, prefix, s0, s1)
        if arc.shape[0] >= 2:
            segs.append((s0, s1, arc))
    if not segs:
        return []

    root = "skin_suction_surface" if side == _SUCTION else "skin_pressure_surface"
    n = len(segs)
    out: list[ShellMidlineStrip] = []
    bay_i = 0
    for i, (_s0, _s1, arc) in enumerate(segs):
        # Zero-web: a single half-ring is a full side arc, not an LE tip panel.
        if n == 1:
            label = root
            out.append(
                ShellMidlineStrip(
                    label=label,
                    kind="skin",
                    midline_twisted_pitch_axis_blade_frame=arc,
                    thickness_m=float(thickness_m),
                    closed=False,
                    skin_side=side,
                    skin_arc=None,
                    bay_index=None,
                )
            )
            continue
        # Suction is TE→LE (index 0 = TE); pressure is LE→TE (index 0 = LE).
        if i == 0:
            skin_arc, bay_index = (
                (_TRAILING, None) if side == _SUCTION else (_LEADING, None)
            )
        elif i == n - 1:
            skin_arc, bay_index = (
                (_LEADING, None) if side == _SUCTION else (_TRAILING, None)
            )
        else:
            skin_arc, bay_index = "bay", bay_i
            bay_i += 1

        if skin_arc == "bay":
            label = f"{root}:bay_{int(bay_index)}"
        else:
            label = f"{root}:{skin_arc}"
        out.append(
            ShellMidlineStrip(
                label=label,
                kind="skin",
                midline_twisted_pitch_axis_blade_frame=arc,
                thickness_m=float(thickness_m),
                closed=False,
                skin_side=side,
                skin_arc=skin_arc,
                bay_index=bay_index,
            )
        )
    return out


def split_skin_at_junctions(
    strips: tuple[ShellMidlineStrip, ...],
    *,
    merge_s_tol_m: float | None = None,
) -> tuple[ShellMidlineStrip, ...]:
    """Snap web/cap feet and panelise skin into suction/pressure × LE/bay/TE arcs.

    Hierarchical labels: ``skin_{suction|pressure}_surface:{leading_edge|bay_i|trailing_edge}``.
    """
    skin_list = [s for s in strips if s.kind == "skin"]
    if not skin_list:
        return strips

    side_skins = _ensure_ss_ps_skins(skin_list)
    by_side: dict[str, ShellMidlineStrip] = {}
    for s in side_skins:
        side = _normalise_side(s.skin_side, label=_label_key(s.label))
        if side is not None:
            by_side[side] = s
    if _SUCTION not in by_side or _PRESSURE not in by_side:
        others = [s for s in strips if s.kind != "skin"]
        return tuple(side_skins + others)

    out_non_skin: list[ShellMidlineStrip] = []
    web_feet_by_side: dict[str, list[float]] = {_SUCTION: [], _PRESSURE: []}
    web_polys: list[NDArray[np.float64]] = []

    for st in strips:
        if st.kind == "skin":
            continue
        if st.kind == "web":
            arr = np.asarray(st.midline_twisted_pitch_axis_blade_frame, dtype=float).copy()
            if arr.shape[0] < 2:
                out_non_skin.append(st)
                continue
            # The foot is already the junction: ShearWeb.midline_polyline emits
            # the web mid-surface out to the skin mid-surface, so there is
            # nothing here to relocate. This pass only *indexes* — it finds the
            # arc length at which the skin must be split. Moving the foot is
            # what used to require guessing a direction, and guessing it
            # perpendicular to the skin is what kinked the web.
            for row in (0, -1):
                p = arr[row]
                best_side = None
                best_d = 1e300
                best_s: float | None = None
                for side, sk in by_side.items():
                    pts = np.asarray(sk.midline_twisted_pitch_axis_blade_frame, dtype=float)
                    d, q = point_to_open_polyline_nearest(p, pts)
                    if d >= best_d:
                        continue
                    pref = open_polyline_prefix_lengths(pts)
                    a = pts[:-1]
                    b = pts[1:]
                    ab = b - a
                    l2 = np.sum(ab * ab, axis=1) + 1e-30
                    t = np.clip(np.sum((q.reshape(1, 2) - a) * ab, axis=1) / l2, 0.0, 1.0)
                    proj = a + t[:, None] * ab
                    j = int(np.argmin(np.sum((q.reshape(1, 2) - proj) ** 2, axis=1)))
                    best_d = float(d)
                    best_side = side
                    best_s = float(pref[j] + float(t[j]) * float(np.linalg.norm(ab[j])))
                if best_side is not None and best_s is not None:
                    web_feet_by_side[best_side].append(float(best_s))
            out_non_skin.append(
                replace(st, midline_twisted_pitch_axis_blade_frame=arr)
            )
            web_polys.append(arr)
            continue

        # Caps are passed through untouched. A cap end is either at a web
        # junction — where it already sits exactly, measured 0.000 mm from the
        # web mid-surface for every family-D layout — or it is a free end
        # terminating mid-bay, 40-149 mm from the nearest web on families B and
        # C. Neither belongs on the skin: a cap mid-surface runs a full
        # (skin + cap) / 2 inboard of the skin mid-surface and never meets it.
        # This block used to snap any end that failed an "is it on a web"
        # tolerance onto the nearest skin point, displacing free B/C cap ends by
        # ~10 mm and inflating their arc length by ~18 %.
        out_non_skin.append(st)

    new_skins: list[ShellMidlineStrip] = []
    for side in (_SUCTION, _PRESSURE):
        sk = by_side[side]
        pts = np.asarray(sk.midline_twisted_pitch_axis_blade_frame, dtype=float)
        new_skins.extend(
            _classify_side_arcs(
                side,
                pts,
                float(sk.thickness_m),
                web_feet_by_side.get(side, []),
                merge_s_tol_m=merge_s_tol_m,
            )
        )

    if len(new_skins) < 2:
        new_skins = list(by_side.values())

    return tuple(new_skins + out_non_skin)


__all__ = [
    "split_skin_at_junctions",
]
