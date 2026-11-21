"""Loft mapped strips between sampled stations."""
from __future__ import annotations

import math
from typing import Mapping

import numpy as np

from implicit_geometry_kernel.blade.loft.locked_loft import arc_resample_open
from implicit_geometry_kernel.blade.topology.sizing import resample_open_h
from implicit_geometry_kernel.blade.topology.template import (
    LoopBranch,
    ShellTopologyTemplate,
    WebBranch,
)

class LoftFace:
    """One lofted face. Iterates as ``(name, face, collapsing)`` so existing callers keep working."""

    __slots__ = ("name", "face", "collapsing", "constructor", "origins")

    def __init__(self, name: str, face: np.ndarray, collapsing: bool, constructor: str,
                 origins: tuple[str, ...] | None = None):
        self.name = name
        self.face = np.asarray(face, float)
        self.collapsing = bool(collapsing)
        self.constructor = constructor
        self.origins = tuple(origins) if origins is not None else ("station",) * len(self.face)

    def __iter__(self):
        yield self.name
        yield self.face
        yield self.collapsing

    def __getitem__(self, index):
        return (self.name, self.face, self.collapsing)[index]


def _resample_polyline(pts: np.ndarray, count: int) -> np.ndarray:
    """Arc-length resample that keeps the section endpoints."""
    pts = np.asarray(pts, float)
    count = int(count)
    if len(pts) == count or count < 2 or len(pts) < 2:
        return pts
    out = arc_resample_open(pts, count)
    out[0] = pts[0]
    out[-1] = pts[-1]
    return out


def align_event_polylines(samples: list, templates: list) -> list:
    """Give both ends of a topology bay the station chain's node count.

    When a region is born or dies, the neighbour's station chain picks up that
    region's nodes on one side only. Resampling the neighbour's own polyline to
    the chain count, from root to tip, makes the bay ruled. Later bays see the
    updated station.

    The return is one band of pairs per bay. A resample replaces the tip
    polyline in the pair already built; the next bay reads the updated station
    on its own first pass.
    """
    bands = []
    for index in range(len(samples) - 1):
        pairs = _band_pairs(samples[index], samples[index + 1], templates[index], templates[index + 1])
        updated = []
        for name, pa, pb in pairs:
            if (
                name in samples[index + 1]
                and len(pa) >= 2 and len(pb) >= 2 and len(pa) != len(pb)
            ):
                samples[index + 1][name] = _resample_polyline(samples[index + 1][name], len(pa))
                pb = np.asarray(samples[index + 1][name], float)
            updated.append((name, pa, pb))
        bands.append(updated)
    return bands


def _polyline_covered(pts: np.ndarray, pairs, own: str) -> bool:
    """True when every point of ``pts`` is already a node of another strip."""
    pts = np.asarray(pts, float)
    if len(pts) < 2:
        return False
    for name, pa, pb in pairs:
        if name == own:
            continue
        for host in (pa, pb):
            host = np.asarray(host, float)
            if len(host) < 2:
                continue
            if all(float(np.min(np.linalg.norm(host - point, axis=1))) <= 1e-8 for point in pts):
                return True
    return False


def loft_mapped_strips(s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
                       t0: ShellTopologyTemplate, t1: ShellTopologyTemplate,
                       *, on_section=None, strip_log: list | None = None,
                       pairs: list | None = None) -> list:
    """All-quad strips between two sampled stations. Omit one-sided webs.

    Each t1 skin branch is ruled from the t0 skin walk between its mapped
    endpoints (web onset/drop splits a panel along the mould, not a chord).
    Both sides of every strip are the stations' own nodes, so neighbouring
    bands share every station node.

    A region present on both stations at the same count is ruled quads. Where a
    region is born, dies, or the trailing face closes, the band gets a mid-span
    row: equal strips split into two quads through one shared seam node, and
    the collapsed face reuses that node. Its centroid stays on the triangle.

    Returns ``(name, face, collapsing)``; ``face`` is a quad ``(4, 3)``.
    ``collapsing`` marks the faces of a transition (and the triangle halves of
    a folded quad, which are ``(3, 3)``).
    """
    if pairs is None:
        pairs = _band_pairs(s0, s1, t0, t1)
    event = any(len(pa) != len(pb) for _name, pa, pb in pairs)
    trigger = ""
    if event:
        for name, pa, pb in pairs:
            if len(pa) != len(pb):
                trigger = name
                break
    # A ruled bay never reads the seams, and building them asks for a section
    # at the mid-span of every bay.
    seams = _band_seams(pairs, on_section) if event else {}
    faces: list[LoftFace] = []
    for name, pa, pb in pairs:
        stats: dict = {}
        if not event:
            made = _ruled_quads(name, pa, pb)
            kind = "ruled"
        elif len(pa) == len(pb):
            made = _split_quads(name, pa, pb, on_section=on_section, seams=seams)
            kind = "midspan_split"
        elif min(len(pa), len(pb)) == 1 and _polyline_covered(
                pa if len(pa) > len(pb) else pb, pairs, name):
            made = []
            kind = "omitted"
            note = "absorbed"
        elif min(len(pa), len(pb)) == 1:
            made = _transition_quads(name, pa, pb, on_section=None, stats=stats, seams=seams)
            kind = "fan"
        else:
            made = _transition_quads(name, pa, pb, on_section=on_section, stats=stats, seams=seams)
            kind = "transition"
        note = ""
        if kind == "fan":
            note = "apex"
        elif kind == "omitted":
            note = "absorbed"
        elif any(face.constructor == "bowtie" for face in made):
            note = "bowtie"
        if strip_log is not None:
            strip_log.append({
                "region": name, "n0": int(len(pa)), "n1": int(len(pb)),
                "constructor": kind, "trigger": trigger,
                "adv_a": int(stats.get("adv_a", 0)), "adv_b": int(stats.get("adv_b", 0)),
                "ties": int(stats.get("ties", 0)), "note": note,
                "site0": int(stats.get("site0", -1)), "site1": int(stats.get("site1", -1)),
            })
        faces.extend(made)
    if strip_log is not None:
        names1 = {w.name for w in t1.webs}
        names0 = {w.name for w in t0.webs}
        for web in t0.webs:
            if web.name not in s1 or web.name not in names1:
                strip_log.append({
                    "region": web.name, "n0": int(len(s0.get(web.name, ()))), "n1": 0,
                    "constructor": "omitted", "trigger": trigger,
                    "adv_a": 0, "adv_b": 0, "ties": 0, "note": "web omitted",
                    "site0": -1, "site1": -1,
                })
        for web in t1.webs:
            if web.name not in s0 or web.name not in names0:
                strip_log.append({
                    "region": web.name, "n0": 0, "n1": int(len(s1.get(web.name, ()))),
                    "constructor": "omitted", "trigger": trigger,
                    "adv_a": 0, "adv_b": 0, "ties": 0, "note": "web omitted",
                    "site0": -1, "site1": -1,
                })
    return faces


def _band_pairs(s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
                t0: ShellTopologyTemplate, t1: ShellTopologyTemplate
                ) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """``(name, t0 nodes, t1 nodes)`` for every strip of the band."""
    raw: list[tuple[str, np.ndarray, np.ndarray]] = []
    t0_wraps = any(b.wrap for b in t0.loop)
    groups = _met_key_groups(s1, t1)
    for branch in t0.loop:
        apex = _dies_at(branch, s0, s1, t0, t1, groups)
        if apex is not None:
            raw.append((branch.name, np.asarray(s0[branch.name], float), apex.reshape(1, -1)))
    for branch in t1.loop:
        if branch.wrap and not t0_wraps and branch.name in s1:
            ridge0 = _sample_te_ridge(s0)
            if ridge0 is not None:
                raw.append((branch.name, ridge0.reshape(1, -1), np.asarray(s1[branch.name], float)))
            continue
        if branch.wrap or branch.name not in s1:
            continue
        pb = np.asarray(s1[branch.name], float)
        pa = _skin_parent(
            _outer_ends(branch, s0, t0, groups), s0, s1, t0, t1, pb, groups,
        )
        if pa is not None:
            raw.append((branch.name, pa, pb))
    for branch in t0.loop:
        if not branch.wrap or branch.name not in s0:
            continue
        dest = _strip_dest(branch, s1, t0, t1)
        pa = np.asarray(s0[branch.name], float)
        if dest is not None and dest in s1:
            pb = np.asarray(s1[dest], float)
        else:
            ridge = _sample_te_ridge(s1)
            if ridge is None:
                continue
            pb, dest = ridge.reshape(1, -1), branch.name
        raw.append((dest, pa, pb))
    _absorb_rename_gaps(raw, s0, s1, t0, t1, groups)
    names1 = {w.name for w in t1.webs}
    for web in t0.webs:
        if web.name not in s0 or web.name not in s1 or web.name not in names1:
            continue
        raw.append((web.name, np.asarray(s0[web.name], float), np.asarray(s1[web.name], float)))
    out: list[tuple[str, np.ndarray, np.ndarray]] = []
    for name, pa, pb in raw:
        if len(pa) == 0 or len(pb) == 0 or not (np.all(np.isfinite(pa)) and np.all(np.isfinite(pb))):
            continue
        pa, pb = _dedupe_open(np.asarray(pa, float)), _dedupe_open(np.asarray(pb, float))
        if len(pa) == 1 and len(pb) == 1:
            continue
        out.append((name, pa, pb))
    return out


def _met_key_groups(s1: Mapping[str, np.ndarray], t1: ShellTopologyTemplate) -> list[set]:
    """Sets of t1 keys that meet at one point: the ends of each region with no length.

    The sampler gives a region with no length no nodes (its two keys meet, for
    example where a cap edge passes a web foot); the keys of each set are one
    point on the t1 station.
    """
    groups: list[set] = []
    for b in t1.loop:
        if b.wrap or b.name in s1 or b.start == b.end:
            continue
        q0, q1 = _xyz_of_vid(s1, t1, b.start), _xyz_of_vid(s1, t1, b.end)
        if q0 is None or q1 is None or float(np.linalg.norm(q0 - q1)) > 0.0:
            continue
        merged = {(b.start, b.side), (b.end, b.side)}
        for g in [g for g in groups if g & merged]:
            merged |= g
            groups.remove(g)
        groups.append(merged)
    return groups


def _group_of(groups: list[set], vid, side: str) -> set | None:
    for g in groups:
        if (vid, side) in g:
            return g
    return None


def _le_order(s0: Mapping[str, np.ndarray], t0: ShellTopologyTemplate, side: str, vid) -> float | None:
    """Arc distance of a t0 key from the LE end of its side (orders the keys of a group)."""
    chain = _concat_side_samples(s0, t0, side)
    p = _te_key_xyz(s0, t0, vid, side)
    if chain is None or p is None:
        return None
    seg = np.linalg.norm(np.diff(chain, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    i = int(np.argmin(np.linalg.norm(chain - p, axis=1)))
    le = _xyz_of_vid(s0, t0, ("anchor", "le"))
    i_le = 0 if le is None else int(np.argmin(np.linalg.norm(chain - le, axis=1)))
    return float(abs(arc[i] - arc[i_le]))


def _outer_ends(branch: LoopBranch, s0: Mapping[str, np.ndarray], t0: ShellTopologyTemplate,
                groups: list[set]) -> LoopBranch:
    """A t1 branch whose end keys meet others on t1 takes the outer t0 key of each set.

    Branches run LE to TE on their side, so a branch ending in a set ends at the
    set's most LE-ward key on t0 and one starting in a set starts at its most
    TE-ward key. The t0 regions between the keys of a set are left to die at
    the meeting point (:func:`_dies_at`), so no t0 node is covered twice.
    """
    from dataclasses import replace

    def pick(vid, last: bool):
        g = _group_of(groups, vid, branch.side)
        if g is None:
            return vid
        ranked: list[tuple[float, object]] = []
        for v, side in g:
            d = _le_order(s0, t0, side, v) if side == branch.side else None
            if d is not None:
                ranked.append((d, v))
        if not ranked:
            return vid
        ranked.sort(key=lambda dv: dv[0])
        return ranked[-1][1] if last else ranked[0][1]

    start, end = pick(branch.start, True), pick(branch.end, False)
    if start == branch.start and end == branch.end:
        return branch
    return replace(branch, start=start, end=end)


def _dies_at(branch: LoopBranch, s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
             t0: ShellTopologyTemplate, t1: ShellTopologyTemplate, groups: list[set]
             ) -> np.ndarray | None:
    """The t1 point a t0 region shrinks to, when both its keys meet on t1."""
    if branch.wrap or branch.name not in s0 or not len(s0[branch.name]):
        return None
    g = _group_of(groups, branch.start, branch.side)
    if g is None or (branch.end, branch.side) not in g:
        return None
    for vid, _side in g:
        p = _xyz_of_vid(s1, t1, vid)
        if p is not None:
            return p
    return None


def _ruled_quad_array(pa: np.ndarray, pb: np.ndarray) -> np.ndarray:
    """One ``(n, 4, 3)`` array of ruled quads, corner order ``a, b, c, d``."""
    return np.stack([pa[:-1], pa[1:], pb[1:], pb[:-1]], axis=1)


def _ruled_quad_masks(quads: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Which ruled quads are real faces, and which of those are bowties.

    ``quads`` is ``(n, 4, 3)`` in corner order ``a, b, c, d``. A face needs three
    distinct corners, matching :func:`_finite_loft_face`. A bowtie is the same
    fold :func:`_quad_bowtie` splits into two triangles.
    """
    quads = np.asarray(quads, float)
    if len(quads) == 0:
        empty = np.zeros(0, dtype=bool)
        return empty, empty
    a, b, c, d = quads[:, 0], quads[:, 1], quads[:, 2], quads[:, 3]

    def apart(u: np.ndarray, v: np.ndarray) -> np.ndarray:
        return np.linalg.norm(u - v, axis=1) > 1e-12

    keep_b = apart(b, a)
    last = np.where(keep_b[:, None], b, a)
    keep_c = apart(c, last)
    last = np.where(keep_c[:, None], c, last)
    keep_d = apart(d, last) & apart(d, a)
    keep = (1 + keep_b.astype(np.int64) + keep_c.astype(np.int64) + keep_d.astype(np.int64)) >= 3
    n1 = np.cross(b - a, c - a)
    n2 = np.cross(c - a, d - a)
    ln1 = np.linalg.norm(n1, axis=1)
    ln2 = np.linalg.norm(n2, axis=1)
    bow = keep & (ln1 >= 1e-18) & (ln2 >= 1e-18) & (np.sum(n1 * n2, axis=1) < -0.2 * ln1 * ln2)
    return keep, bow


def _ruled_quads(name: str, pa: np.ndarray, pb: np.ndarray) -> list[LoftFace]:
    pa = np.asarray(pa, float)
    pb = np.asarray(pb, float)
    out: list[LoftFace] = []
    if len(pa) < 2:
        return out
    for quad in _ruled_quad_array(pa, pb):
        if _quad_bowtie(quad):
            a, b, c, d = quad
            _append_face(out, name, np.stack([a, b, c]), True, "bowtie")
            _append_face(out, name, np.stack([a, c, d]), True, "bowtie")
        else:
            _append_face(out, name, np.array(quad, dtype=float, copy=True), False, "ruled")
    return out


def _mid(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Midpoint of a spanwise or interior edge, ``p`` on t0 first (same bits on both sides)."""
    return 0.5 * (p + q)


def _split_quad_arrays(pa: np.ndarray, pb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Station-side and far-side quads through the mid-span row, each ``(n, 4, 3)``."""
    mids = 0.5 * (pa + pb)
    near = np.stack([pa[:-1], pa[1:], mids[1:], mids[:-1]], axis=1)
    far = np.stack([mids[:-1], mids[1:], pb[1:], pb[:-1]], axis=1)
    return near, far


def _xyz_key(p: np.ndarray) -> tuple[float, float, float]:
    return (float(p[0]), float(p[1]), float(p[2]))


def _seam_key(a: np.ndarray, b: np.ndarray) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Identity of one spanwise edge, independent of which face asks for it."""
    ka, kb = _xyz_key(a), _xyz_key(b)
    return (ka, kb) if ka <= kb else (kb, ka)


def _band_seams(pairs: list[tuple[str, np.ndarray, np.ndarray]], on_section) -> dict:
    """One mid-span node per spanwise edge a full strip already owns.

    The first strip that has both endpoints writes the node. A collapsed face
    and the neighbouring split read it back, so the seam is one point.
    """
    # A three-node skin collapse is one coarse edge. Snapping its seam to the
    # nearest key lands on a corner and the fan quad folds flat, so those
    # seams stay on the midline. The trailing edge keeps the key.
    short: list[np.ndarray] = []
    for name, pa, pb in pairs:
        if name == "te" or min(len(pa), len(pb)) != 1:
            continue
        long = pa if len(pa) >= len(pb) else pb
        if 2 <= len(long) <= 3:
            short.append(np.asarray(long, float))

    def _on_short(point: np.ndarray) -> bool:
        point = np.asarray(point, float)
        for curve in short:
            if float(np.min(np.linalg.norm(curve - point, axis=1))) <= 1e-8:
                return True
        return False

    seams: dict = {}
    for name, pa, pb in pairs:
        if len(pa) < 2 or len(pb) < 2:
            continue
        for a, b in ((pa[0], pb[0]), (pa[-1], pb[-1])):
            key = _seam_key(a, b)
            if key not in seams:
                on_key = not (_on_short(a) or _on_short(b))
                seams[key] = _project(name, _mid(a, b), on_section, on_key)
    return seams


def _project(name: str, xyz: np.ndarray, on_section, at_end: bool) -> np.ndarray:
    if on_section is None:
        return np.asarray(xyz, float)
    return np.asarray(on_section(name, np.asarray(xyz, float), at_end), float)


def _span_mids(name: str, pa: np.ndarray, pb: np.ndarray, on_section) -> np.ndarray:
    """Mid-span node of each spanwise edge, on the section when ``on_section`` is given.

    The two chordwise ends keep the station-key lookup. The edges between them
    share one projection onto the station loop.
    """
    mids = 0.5 * (np.asarray(pa, float) + np.asarray(pb, float))
    if on_section is None or len(mids) == 0:
        return np.asarray(mids, float)
    out = np.array(mids, dtype=float, copy=True)
    out[0] = on_section(name, out[0], True)
    if len(out) > 1:
        out[-1] = on_section(name, out[-1], True)
    if len(out) > 2:
        out[1:-1] = np.asarray(on_section(name, out[1:-1], False), float)
    return out


def _split_quads(name: str, pa: np.ndarray, pb: np.ndarray, *, on_section=None,
                 seams: dict | None = None) -> list[LoftFace]:
    """Two quads per ruled quad, through the mid-span row of an event band."""
    pa = np.asarray(pa, float)
    pb = np.asarray(pb, float)
    out: list[LoftFace] = []
    if len(pa) < 2 or len(pa) != len(pb):
        return out
    mids = _span_mids(name, pa, pb, on_section)
    if seams:
        for index in (0, len(mids) - 1):
            hit = seams.get(_seam_key(pa[index], pb[index]))
            if hit is not None:
                mids[index] = hit
    near = np.stack([pa[:-1], pa[1:], mids[1:], mids[:-1]], axis=1)
    far = np.stack([mids[:-1], mids[1:], pb[1:], pb[:-1]], axis=1)
    near_origin = ("station", "station", "span_mid", "span_mid")
    far_origin = ("span_mid", "span_mid", "station", "station")
    for i in range(len(near)):
        _append_face(out, name, np.array(near[i], dtype=float, copy=True), False, "midspan_split", near_origin)
        _append_face(out, name, np.array(far[i], dtype=float, copy=True), False, "midspan_split", far_origin)
    return out


def _quarter_is_shallow(edge: float, rise: float, skew: float = 0.0) -> bool:
    """True when the quarter-bay template would leave a corner under 20°.

    ``rise`` is the part of the end ruling square to the coarse edge, and
    ``skew`` is how far that ruling slides along the edge, in edge lengths.
    A square ruling is under 20° when the rise is shorter than about three
    quarters of the edge.     A slide of more than half an edge keeps the sliver
    under 20° out to a rise about one and a half edges. A taller, squarer bay
    keeps the quarter-bay template.
    """
    if edge <= 1e-12:
        return False
    ratio = rise / edge
    if ratio < math.tan(math.radians(20.0)) / 0.5:
        return True
    if ratio < 1.55 and skew >= 0.55:
        return True
    # A mild slide at this height puts more than twice the median into the
    # coarse quad of the quarter-bay split. The full-bay sample stays legal
    # there. A harder slide at the same height does not.
    return 0.8 <= ratio < 1.05 and 0.30 <= skew < 0.45


def _ruling_frame(across: np.ndarray, ruling: np.ndarray) -> tuple[float, float]:
    """Perpendicular rise and along-edge slide of one end ruling, each over the edge."""
    across = np.asarray(across, float)
    ruling = np.asarray(ruling, float)
    edge = float(np.linalg.norm(across))
    if edge <= 1e-12:
        return 1.0, 0.0
    along = float(np.dot(ruling, across)) / edge
    perp = float(np.linalg.norm(ruling - across * (along / edge)))
    return perp / edge, abs(along) / edge


def _site_fractions(rise_over_edge: float, skew: float) -> tuple[float, float]:
    """Fractions of the coarse edge and the end ruling for one interior sample.

    The sample is the point that keeps a short bay's site quads at 20° or more
    and within twice the median area of their bay. A square ruling uses the
    knots. A ruling that slides along the edge moves the sample toward the fine
    station, which is what keeps that skewed site off the coarse edge.
    """
    r = float(np.clip(rise_over_edge, 0.24, 1.20))
    knots = np.array([0.24, 0.30, 0.40, 0.50, 0.70, 1.20])
    alpha = float(np.interp(r, knots, [0.60, 0.68, 0.68, 0.68, 0.68, 0.50]))
    beta = float(np.interp(r, knots, [0.63, 0.44, 0.37, 0.34, 0.34, 0.55]))
    alpha_slide = float(np.interp(r, knots, [0.82, 0.82, 0.82, 0.82, 0.82, 0.50]))
    beta_slide = float(np.interp(r, knots, [0.90, 0.90, 0.90, 0.90, 0.90, 0.70]))
    blend = float(np.clip((skew - 0.25) / 0.45, 0.0, 1.0))
    return (1.0 - blend) * alpha + blend * alpha_slide, (1.0 - blend) * beta + blend * beta_slide


def _site_sample(origin: np.ndarray, across: np.ndarray, ruling: np.ndarray) -> np.ndarray:
    """Interior sample in the plane of the coarse edge and the end ruling."""
    rise, skew = _ruling_frame(across, ruling)
    alpha, beta = _site_fractions(rise, skew)
    return np.asarray(origin, float) + alpha * np.asarray(across, float) + beta * np.asarray(ruling, float)


def _step_quads(name: str, pa: np.ndarray, pb: np.ndarray, *, on_section=None,
                stats: dict | None = None) -> list[LoftFace]:
    """A count change of two edges: 1:1 quads, plus two adjacent 2:1 sites.

    The two sites sit symmetric about the middle of the arc and share the
    spanwise edge between them. Each coarse edge faces two fine edges. No
    coarse edge faces three. A region with only two coarse edges uses both.
    The coarse edge stays one quad edge, so the next bay shares it, and the
    fine middle node stays a corner. On a tall bay the extra nodes gather
    above the mid-span line. On a short bay that quarter-height sliver would
    be under 20°, so the width change is a trapezoid of the full bay and the
    end mid-span nodes stay where the neighbouring strip expects them.
    """
    pa = np.asarray(pa, float)
    pb = np.asarray(pb, float)
    if abs(len(pa) - len(pb)) != 2 or min(len(pa), len(pb)) < 3:
        return []
    if len(pb) > len(pa):
        coarse, fine = pa, pb
    else:
        coarse, fine = pb, pa
    nc = len(coarse)
    n_edges = nc - 1
    if n_edges < 2:
        return []
    start = (n_edges - 2) // 2
    if stats is not None:
        stats["site0"] = start
        stats["site1"] = start + 1
    out: list[LoftFace] = []
    coarse_is_pa = bool(np.shares_memory(coarse, pa))

    def seg(c0, c1, f0, f1, end0: bool, end1: bool) -> None:
        """One 1:1 spanwise quad, split through the mid-span row.

        ``end0`` and ``end1`` mark the region's two boundary edges, so the
        mid-span node there is the same projection the neighbouring strip uses.
        """
        if coarse_is_pa:
            a0, a1, b0, b1 = c0, c1, f0, f1
        else:
            a0, a1, b0, b1 = f0, f1, c0, c1
        m0 = _project(name, _mid(a0, b0), on_section, end0)
        m1 = _project(name, _mid(a1, b1), on_section, end1)
        _append_face(out, name, np.stack([a0, a1, m1, m0]), False, "midspan_split",
                     ("station", "station", "span_mid", "span_mid"))
        _append_face(out, name, np.stack([m0, m1, b1, b0]), False, "midspan_split",
                     ("span_mid", "span_mid", "station", "station"))

    def lifted(samples: np.ndarray, c0, c2, index: int) -> np.ndarray:
        """Section sample toward one fine node, off the mid-span line."""
        arc = _arc_fracs(samples)
        t = float(arc[index])
        on_coarse = (1.0 - t) * c0 + t * c2
        rise, _skew = _ruling_frame(c2 - c0, samples[0] - c0)
        span = 0.72 if rise > 1.55 else 0.75
        raw = (1.0 - span) * on_coarse + span * samples[index]
        return _project(name, raw, on_section, False)

    def emit(verts, origins) -> None:
        if not coarse_is_pa:
            verts = verts[::-1]
            origins = origins[::-1]
        _append_face(out, name, np.stack(verts), True, "step_patch", origins)

    e = 0
    fi = 0
    while e < n_edges:
        if e == start:
            c0, c1, c2 = coarse[e], coarse[e + 1], coarse[e + 2]
            samples = fine[fi:fi + 5]
            f0, f1, f2, f3, f4 = samples
            ml = _project(name, _mid(c0, f0), on_section, e == 0)
            mr = _project(name, _mid(c2, f4), on_section, e + 2 == n_edges)
            left_edge = float(np.linalg.norm(c1 - c0))
            right_edge = float(np.linalg.norm(c2 - c1))
            left_rise, left_skew = _ruling_frame(c1 - c0, f0 - c0)
            right_rise, right_skew = _ruling_frame(c1 - c2, f4 - c2)
            if (_quarter_is_shallow(left_edge, left_rise * left_edge, left_skew)
                    or _quarter_is_shallow(right_edge, right_rise * right_edge, right_skew)):
                gl = _project(name, _site_sample(c0, c1 - c0, f0 - c0), on_section, False)
                gr = _project(name, _site_sample(c2, c1 - c2, f4 - c2), on_section, False)
                station = ("station", "station", "span_mid", "span_mid")
                emit((c0, c1, gl, ml), station)
                emit((ml, gl, f1, f0), ("span_mid", "span_mid", "station", "station"))
                emit((gl, c1, f2, f1), ("span_mid", "station", "station", "station"))
                emit((c1, c2, mr, gr), station)
                emit((f4, f3, gr, mr), station)
                emit((gr, f3, f2, c1), ("span_mid", "station", "station", "station"))
            else:
                mm = _project(name, _mid(c1, f2), on_section, False)
                p = lifted(samples, c0, c2, 1)
                s = lifted(samples, c0, c2, 2)
                q = lifted(samples, c0, c2, 3)
                station = ("station", "station", "span_mid", "span_mid")
                interior = ("span_mid",) * 4
                emit((f0, f1, p, ml), station)
                emit((f1, f2, s, p), station)
                emit((ml, p, s, mm), interior)
                emit((c0, c1, mm, ml), station)
                emit((f2, f3, q, s), station)
                emit((f3, f4, mr, q), station)
                emit((mm, s, q, mr), interior)
                emit((c1, c2, mr, mm), station)
            e += 2
            fi += 4
        else:
            seg(coarse[e], coarse[e + 1], fine[fi], fine[fi + 1], e == 0, e + 1 == n_edges)
            e += 1
            fi += 1
    return out


def _transition_quads(name: str, pa: np.ndarray, pb: np.ndarray, *,
                      on_section=None, stats: dict | None = None,
                      seams: dict | None = None) -> list:
    """All-quad strip between station polylines with different node counts.

    Both sides have an even number of edges (each region has an odd node
    count, and onset keys land on even nodes), so every other node is a
    coarse node. The coarse sides are zipped into triangles and each triangle
    splits into three quads through its centroid and edge midpoints. On a
    station line the midpoint of a coarse edge is the station node between
    its ends, and on the band's spanwise sides it is the mid-span row node,
    so the strip shares every node with its neighbours.     A side of one node
    (a panel born from a point) is the apex of a fan. That fan stays on the
    triangle rulings: the section has no face left to project the centroid
    and the spanwise mids onto.
    """
    if (len(pa) - 1) % 2 or (len(pb) - 1) % 2:
        raise ValueError(
            f"transition strip {name!r} has {len(pa) - 1} and {len(pb) - 1} edges; "
            "an all-quad transition needs an even number on each station"
        )
    ca, cb = pa[::2], pb[::2]

    def point(v: tuple[str, int]) -> np.ndarray:
        return ca[v[1]] if v[0] == "a" else cb[v[1]]

    def mid(u: tuple[str, int], v: tuple[str, int]) -> tuple[np.ndarray, str]:
        if u[0] == v[0] and abs(u[1] - v[1]) == 1:
            fine = pa if u[0] == "a" else pb
            return fine[2 * min(u[1], v[1]) + 1], "station"
        if u[0] == "b" and v[0] == "a":
            u, v = v, u
        at_end = u[1] in (0, len(ca) - 1) and v[1] in (0, len(cb) - 1)
        if seams:
            hit = seams.get(_seam_key(point(u), point(v)))
            if hit is not None:
                return np.asarray(hit, float), "span_mid"
        return _project(name, _mid(point(u), point(v)), section, at_end), "span_mid"

    adv_a = adv_b = ties = 0
    tris: list[tuple[tuple[str, int], tuple[str, int], tuple[str, int]]] = []
    fan = len(ca) == 1 or len(cb) == 1
    section = None if fan else on_section
    if len(ca) == 1:
        tris = [(("a", 0), ("b", j + 1), ("b", j)) for j in range(len(cb) - 1)]
    elif len(cb) == 1:
        tris = [(("a", i), ("a", i + 1), ("b", 0)) for i in range(len(ca) - 1)]
    else:
        # advance the side whose next coarse edge is centred earlier along the
        # strip, so the zip is the same whichever station is the finer one
        fa, fb = _arc_fracs(ca), _arc_fracs(cb)
        i = j = 0
        while i < len(ca) - 1 or j < len(cb) - 1:
            left = fa[i] + fa[i + 1] if i < len(ca) - 1 else math.inf
            right = fb[j] + fb[j + 1] if j < len(cb) - 1 else math.inf
            if j == len(cb) - 1 or (i < len(ca) - 1 and left <= right):
                if i < len(ca) - 1 and j < len(cb) - 1 and abs(left - right) <= 1e-12:
                    ties += 1
                tris.append((("a", i), ("a", i + 1), ("b", j)))
                i += 1
                adv_a += 1
            else:
                tris.append((("a", i), ("b", j + 1), ("b", j)))
                j += 1
                adv_b += 1
    if stats is not None:
        stats["adv_a"] = adv_a
        stats["adv_b"] = adv_b
        stats["ties"] = ties
    kind = "fan" if fan else "transition"
    out: list[LoftFace] = []
    for v0, v1, v2 in tris:
        p0, p1, p2 = point(v0), point(v1), point(v2)
        raw_c = (p0 + p1 + p2) / 3.0
        c = _project(name, raw_c, section, False)
        c_origin = "centroid" if on_section is None else "span_mid"
        m01, o01 = mid(v0, v1)
        m12, o12 = mid(v1, v2)
        m20, o20 = mid(v2, v0)
        _append_face(out, name, np.stack([p0, m01, c, m20]), True, kind, ("station", o01, c_origin, o20))
        _append_face(out, name, np.stack([p1, m12, c, m01]), True, kind, ("station", o12, c_origin, o01))
        _append_face(out, name, np.stack([p2, m20, c, m12]), True, kind, ("station", o20, c_origin, o12))
    return out


def _dedupe_open(pts: np.ndarray) -> np.ndarray:
    keep = [pts[0]]
    for q in pts[1:]:
        if not np.allclose(keep[-1], q, atol=1e-12, rtol=0.0):
            keep.append(q)
    return np.asarray(keep, float)


def _arc_fracs(pts: np.ndarray) -> np.ndarray:
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    return s / s[-1] if s[-1] > 0.0 else np.linspace(0.0, 1.0, len(pts))


def _match_open_counts(pa: np.ndarray, pb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Arc-resample both open polylines to the same vertex count (keep endpoints)."""
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    n = max(len(pa), len(pb), 2)
    if len(pa) != n:
        pa = arc_resample_open(pa, n)
    if len(pb) != n:
        pb = arc_resample_open(pb, n)
    return pa, pb


def _finite_loft_face(face: np.ndarray) -> bool:
    """Drop collapsed or non-finite strips (they poison 3-D autoscale)."""
    p = np.asarray(face, float)
    if p.ndim != 2 or len(p) < 3 or p.shape[1] < 3:
        return False
    if not np.all(np.isfinite(p)):
        return False
    uniq = [p[0]]
    for q in p[1:]:
        if not np.allclose(uniq[-1], q, atol=1e-12, rtol=0.0):
            uniq.append(q)
    if len(uniq) >= 2 and np.allclose(uniq[0], uniq[-1], atol=1e-12, rtol=0.0):
        uniq = uniq[:-1]
    return len(uniq) >= 3


def _append_face(faces: list, name: str, face: np.ndarray, collapsing: bool,
                 constructor: str = "ruled", origins: tuple[str, ...] | None = None) -> None:
    if _finite_loft_face(face):
        faces.append(LoftFace(name, face, collapsing, constructor, origins))


def _xyz_of_vid(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate, vid) -> np.ndarray | None:
    for b in tmpl.loop:
        if b.name not in samples:
            continue
        pts = np.asarray(samples[b.name], float)
        if len(pts) == 0:
            continue
        if b.start == vid:
            return pts[0].copy()
        if b.end == vid:
            return pts[-1].copy()
    return None


def _concat_side_samples(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate,
                         side: str) -> np.ndarray | None:
    """LE→TE skin polyline on one side from sampled loop branches."""
    parts: list[np.ndarray] = []
    for b in tmpl.loop:
        if b.wrap or b.side != side or b.name not in samples:
            continue
        _append_ring(parts, np.asarray(samples[b.name], float))
    if not parts:
        return None
    return np.vstack(parts)


def _resample_keys(p0: np.ndarray, p1: np.ndarray, arc: np.ndarray, n: int,
                   *, bunch_end: bool = False, size_field: bool = False,
                   chord: float = 1.0, t: float = 0.0) -> np.ndarray:
    """Arc-length resample between two keys that already lie on ``arc`` (TE→LE)."""
    n = max(int(n), 2)
    if len(arc) < 2:
        out = np.repeat(np.asarray(p0, float).reshape(1, -1), n, axis=0)
        out[-1] = p1
        return out
    i0 = int(np.argmin(np.linalg.norm(arc - p0, axis=1)))
    i1 = int(np.argmin(np.linalg.norm(arc - p1, axis=1)))
    if i0 == i1:
        sub = np.vstack([p0, p1])
    elif i0 < i1:
        sub = np.vstack([p0, arc[i0 + 1 : i1], p1])
    else:
        sub = np.vstack([p0, arc[i0 - 1 : i1 : -1], p1])
    if size_field:
        sampled = resample_open_h(sub, chord, t)
    elif bunch_end and n > 2:
        sampled = _arc_resample_bunched(sub, n)
    else:
        sampled = arc_resample_open(sub, n)
    sampled[0] = np.asarray(p0, float)
    sampled[-1] = np.asarray(p1, float)
    return sampled


def _arc_resample_bunched(pts: np.ndarray, n: int) -> np.ndarray:
    """Arc-length resample clustered toward the last vertex (TE)."""
    pts = np.asarray(pts, float)
    n = max(int(n), 2)
    if len(pts) < 2:
        return arc_resample_open(pts, n)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 1e-15:
        return arc_resample_open(pts, n)
    from implicit_geometry_kernel.blade.loft.span_place import resample_open_many

    return resample_open_many(pts[None, :, :], n, bunch=True)[0]


def _on_chain(chain: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Nearest even node of the t0 side chain to a key that t0 lacks (XY).

    Region boundaries sit on even chain nodes (every region has an even number
    of edges), so a key placed on an even node keeps both strips beside it at
    an even edge count, which the all-quad transition needs.
    """
    pts = np.asarray(chain, float)[::2]
    d = np.linalg.norm(pts[:, :2] - np.asarray(p, float)[:2], axis=1)
    return pts[int(np.argmin(d))].copy()


def _absorb_rename_gaps(raw: list, s0: Mapping[str, np.ndarray],
                        s1: Mapping[str, np.ndarray], t0: ShellTopologyTemplate,
                        t1: ShellTopologyTemplate, groups: list) -> None:
    """Prepend each rename-gap t0 region onto exactly one surviving strip.

    A region absent on t1 with no meet-set apex is a rename gap (cap→bay,
    bay→te). Prefer the strip whose ``pa`` starts at the gap's TE end; else the
    strip whose ``pa`` starts on an interior node of the gap. Skip gaps already
    covered (``_outer_ends`` kept the long parent).
    """
    names1 = {b.name for b in t1.loop}
    dying = [
        b for b in t0.loop
        if not b.wrap and b.name not in names1 and b.name in s0
        and _dies_at(b, s0, s1, t0, t1, groups) is None
    ]
    # TE-ward first so a chain of dying bays prepends onto the same survivor.
    for branch in reversed(dying):
        pts = np.asarray(s0[branch.name], float)
        if len(pts) < 2 or _polyline_covered(pts, raw, branch.name):
            continue
        # Contiguous TE-ward survivor: pa starts where the gap ends.
        hit_index = None
        prefix = None
        for index, (name, pa, pb) in enumerate(raw):
            pa = np.asarray(pa, float)
            if len(pa) < 1:
                continue
            if float(np.linalg.norm(pa[0] - pts[-1])) <= 1e-8:
                hit_index, prefix = index, pts[:-1]
                break
        if hit_index is None:
            # Survivor landed mid-gap (cap LE key inside a dying bay).
            for index, (name, pa, pb) in enumerate(raw):
                pa = np.asarray(pa, float)
                if len(pa) < 1:
                    continue
                distances = np.linalg.norm(pts - pa[0], axis=1)
                hit = int(np.argmin(distances))
                if float(distances[hit]) > 1e-8 or hit == 0:
                    continue
                hit_index, prefix = index, pts[:hit]
                break
        if hit_index is None or prefix is None or len(prefix) == 0:
            continue
        name, pa, pb = raw[hit_index]
        merged = np.vstack([prefix, np.asarray(pa, float)])
        raw[hit_index] = (name, merged, pb)
        # A LE neighbour may already have kept a long parent across this gap;
        # clip it at the new join so the two strips share one station edge.
        join = np.asarray(merged[0], float)
        for index, (other, oa, ob) in enumerate(raw):
            if index == hit_index:
                continue
            oa = np.asarray(oa, float)
            if len(oa) < 3:
                continue
            distances = np.linalg.norm(oa - join, axis=1)
            at = int(np.argmin(distances))
            if float(distances[at]) > 1e-8 or at == 0 or at >= len(oa) - 1:
                continue
            raw[index] = (other, oa[: at + 1], ob)


def _skin_parent(branch: LoopBranch, s0: Mapping[str, np.ndarray],
                 s1: Mapping[str, np.ndarray], t0: ShellTopologyTemplate,
                 t1: ShellTopologyTemplate, pb: np.ndarray, groups: list
                 ) -> np.ndarray | None:
    """t0 skin arc between the mapped endpoints of a t1 branch (not a chord).

    A web that only exists on t1 is placed on the t0 mould by nearest XY, so the
    TE panel splits at the new foot instead of collapsing to the TE point.
    """
    pb = np.asarray(pb, float)
    n = max(len(pb), 2)
    chain = _concat_side_samples(s0, t0, branch.side)
    p0 = _te_key_xyz(s0, t0, branch.start, branch.side)
    p1 = _te_key_xyz(s0, t0, branch.end, branch.side)
    if chain is None:
        return _collapsed_parent(branch, s0, t0, t1, n)
    if p0 is None:
        p0 = _on_chain(chain, pb[0])
    if p1 is None:
        p1 = _on_chain(chain, pb[-1])
    own = _chain_nodes_between(chain, p0, p1)
    if own is None:
        return _resample_keys(p0, p1, chain, n)
    # Keys that meet on the next station expand this chain across a collapsing
    # meet-set fan; keep the short sample then. Any other longer parent (rename
    # gap, same-name panel that slid) must stay long or station edges open.
    sample = s0.get(branch.name)
    if sample is not None and len(sample) >= 2 and len(own) > len(sample):
        if _extras_only_meet_set(own, sample, branch.side, s0, s1, t0, t1, groups):
            return np.asarray(sample, float).copy()
    return own


def _extras_only_meet_set(own: np.ndarray, sample: np.ndarray, side: str,
                          s0: Mapping[str, np.ndarray], s1: Mapping[str, np.ndarray],
                          t0: ShellTopologyTemplate, t1: ShellTopologyTemplate,
                          groups: list) -> bool:
    """True when every node of ``own`` beyond ``sample`` lies on a meet-set fan."""
    sample_keys = {tuple(np.round(p, 9)) for p in np.asarray(sample, float)}
    extras = [
        p for p in np.asarray(own, float)
        if tuple(np.round(p, 9)) not in sample_keys
    ]
    if not extras:
        return True
    for point in extras:
        on_meet = False
        for branch in t0.loop:
            if branch.wrap or branch.side != side or branch.name not in s0:
                continue
            pts = np.asarray(s0[branch.name], float)
            if not any(float(np.linalg.norm(p - point)) <= 1e-8 for p in pts):
                continue
            if _dies_at(branch, s0, s1, t0, t1, groups) is None:
                return False
            on_meet = True
            break
        if not on_meet:
            return False
    return True


def _chain_nodes_between(chain: np.ndarray, p0: np.ndarray, p1: np.ndarray) -> np.ndarray | None:
    """The t0 station's own nodes from ``p0`` to ``p1``.

    A strip's lower edge must be the t0 station's nodes, or the strip below and
    the strip above meet at different nodes along the station line. Where the
    two stations carry different counts (a topology event), the band takes an
    all-quad transition instead of resampling either side.
    """
    chain = np.asarray(chain, float)
    d0 = np.linalg.norm(chain - np.asarray(p0, float), axis=1)
    d1 = np.linalg.norm(chain - np.asarray(p1, float), axis=1)
    i0, i1 = int(np.argmin(d0)), int(np.argmin(d1))
    if d0[i0] > 1e-12 or d1[i1] > 1e-12:
        return None
    step = 1 if i1 >= i0 else -1
    return chain[i0 : i1 + step if i1 + step >= 0 else None : step].copy()


def _te_key_xyz(samples: Mapping[str, np.ndarray], tmpl: ShellTopologyTemplate, vid, side: str) -> np.ndarray | None:
    """A t1 key on the t0 station, with the TE ridge and TE corner as one key.

    Where the trailing edge goes from blunt to sharp, a sharp branch ends at the
    ridge ``("TE",)`` and the blunt station has a corner ``("TE", side)`` per
    side instead (and the other way round). The strip must end at that corner,
    or it stops one node short of it and leaves a sliver beside the TE strip.
    """
    if vid in tmpl.vertices:
        got = _xyz_of_vid(samples, tmpl, vid)
        if got is not None:
            return got
    if vid == ("TE",):
        return _xyz_of_vid(samples, tmpl, ("TE", side))
    if isinstance(vid, tuple) and len(vid) == 2 and vid[0] == "TE":
        return _xyz_of_vid(samples, tmpl, ("TE",))
    return None


def _collapsed_parent(branch: LoopBranch, s0: Mapping[str, np.ndarray],
                      t0: ShellTopologyTemplate, t1: ShellTopologyTemplate, n: int) -> np.ndarray | None:
    try:
        start0 = t1.map_vertex(branch.start, t0)
        end0 = t1.map_vertex(branch.end, t0)
    except KeyError:
        return None
    p0 = _xyz_of_vid(s0, t0, start0)
    p1 = _xyz_of_vid(s0, t0, end0)
    if p0 is None:
        p0 = p1
    if p1 is None:
        p1 = p0
    if p0 is None or p1 is None:
        return None
    n = max(int(n), 2)
    if np.allclose(p0, p1):
        return np.repeat(np.asarray(p0, float).reshape(1, -1), n, axis=0)
    return np.linspace(p0, p1, n)


def _same_web_set(t0: ShellTopologyTemplate, t1: ShellTopologyTemplate) -> bool:
    s0, s1 = t0.signature, t1.signature
    if isinstance(s0, tuple) and isinstance(s1, tuple) and len(s0) == 2 and len(s1) == 2:
        return s0[1] == s1[1]
    return s0 == s1


def _sample_te_ridge(samples: Mapping[str, np.ndarray]) -> np.ndarray | None:
    for name in ("te_su", "te_pr", "te"):
        pts = np.asarray(samples.get(name, []), float)
        if len(pts):
            return np.asarray(pts[-1], float)
    return None


def _strip_dest(branch: LoopBranch, s1: Mapping[str, np.ndarray],
                t0: ShellTopologyTemplate, t1: ShellTopologyTemplate) -> str | None:
    if _same_web_set(t0, t1) and branch.name in s1:
        return branch.name
    try:
        start_m = t0.map_vertex(branch.start, t1)
        end_m = t0.map_vertex(branch.end, t1)
    except KeyError:
        return branch.name if branch.name in s1 else None
    if start_m == end_m and not branch.wrap:
        return None
    for other in t1.loop:
        if other.name not in s1 or other.side != branch.side:
            continue
        if other.start == start_m and other.end == end_m:
            return other.name
    return None


def _quad_bowtie(quad: np.ndarray) -> bool:
    """True if the two triangles of a loft quad have opposing normals (a fold)."""
    a, b, c, d = np.asarray(quad, float)
    n1 = np.cross(b - a, c - a)
    n2 = np.cross(c - a, d - a)
    if np.linalg.norm(n1) < 1e-18 or np.linalg.norm(n2) < 1e-18:
        return False
    return float(np.dot(n1, n2)) < -0.2 * float(np.linalg.norm(n1) * np.linalg.norm(n2))


def _append_ring(parts: list[np.ndarray], chain: np.ndarray | None) -> None:
    if chain is None:
        return
    c = np.asarray(chain, float)
    if len(c) == 0:
        return
    if parts and np.allclose(parts[-1][-1], c[0], atol=1e-12, rtol=0.0):
        c = c[1:]
    if len(c):
        parts.append(c)


def _web_pts_from_su(samples: Mapping[str, np.ndarray], web: WebBranch | None,
                     su_end: np.ndarray) -> np.ndarray | None:
    if web is None or web.name not in samples:
        return None
    pts = np.asarray(samples[web.name], float)
    if len(pts) < 2:
        return None
    if np.linalg.norm(pts[-1] - su_end) < np.linalg.norm(pts[0] - su_end):
        pts = pts[::-1].copy()
    return pts

