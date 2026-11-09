"""SystemType-style layout keys → frozen layout specs (SDF-first; no Shapely).

Catalogue tokens are ``[W][X]-[Y]-[Z]-[R]`` with required web anchoring and a
required rib count, e.g. ``2D-F-LC-0R``, ``2D-F-LE-10R``, ``2D-F-PA-10R``.
Orientation cores (``2D-F``, ``0A``) live in the registry.
:func:`resolve_system_type` requires ``LC``, ``LE`` or ``PA`` and ``nR``
(including ``0R``). Bare ``2D-F-LC`` and bare ``-R`` are invalid.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from typing import Literal

from implicit_geometry_kernel.contract.errors import (
    InfeasibleSectionGeometry,
)
from implicit_geometry_kernel.catalogue.structural import StructuralFamily

WebOrientation = Literal["chord_normal", "flapwise"]
WebPositioning = Literal["uniform", "structural", "fixed", "none"]
WebAnchoring = Literal["local_chord", "leading_edge", "pitch_axis"]

_CORE = r"(0[AB]|[1-5][ABCD]-(?:CN|F))"
_ANCHOR_TAG = {"local_chord": "LC", "leading_edge": "LE", "pitch_axis": "PA"}
_TAG_ANCHOR = {v: k for k, v in _ANCHOR_TAG.items()}
_TOKEN = re.compile(rf"^{_CORE}-(LC|LE|PA)-(\d+)R$")
_MISSING_RIBS = re.compile(rf"^{_CORE}-(LC|LE|PA)$")


@dataclass(frozen=True)
class SystemLayoutSpec:
    """Catalogue layout intent for one ``[W][X]-[Y]-[Z]-[R]`` token.

    ``web_chord_fracs`` are fractions in ``[0, 1]`` measured along chord from the
    airfoil leading edge in section *x* (same convention as :class:`AirfoilSDF`).

    When ``n_webs == 0`` (cores ``0A`` / ``0B``), the station worker still
    builds a multi-cell host with empty web positions (family A or B).
    ``geometry_mode="zero_web_catalogue"`` remains catalogue metadata for those
    keys; it does not select a separate host class.

    ``key`` is the full token (``2D-F-LC-0R``, ``2D-F-PA-10R``). ``layout_key``
    is the orientation core (``2D-F``). ``web_anchoring`` is local-chord, leading
    edge, or pitch-axis. ``n_ribs`` is the required ``-NR`` count (``0R`` when
    there are no ribs).
    """

    key: str
    n_cells: int
    n_webs: int
    has_spar_caps: bool
    web_orientation: WebOrientation
    web_positioning: WebPositioning
    web_chord_fracs: tuple[float, ...]
    spar_cap_positioning: Literal["max_thickness", "fixed", "box_spar", "none"]
    geometry_mode: Literal["multicell", "zero_web_catalogue"]
    structural_family: StructuralFamily = "D"
    description: str = ""
    layout_key: str = ""
    web_anchoring: WebAnchoring = "local_chord"
    n_ribs: int = 0


def _uniform_fracs(n_webs: int) -> tuple[float, ...]:
    if n_webs <= 0:
        return ()
    edges = [i / (n_webs + 1) for i in range(1, n_webs + 1)]
    return tuple(float(x) for x in edges)


def _default_web_fracs(n_webs: int) -> tuple[float, ...]:
    # Keep legacy fixed stations for X=2; use centerline for X=1; uniform otherwise.
    if n_webs == 1:
        return (0.50,)
    if n_webs == 2:
        return (0.15, 0.50)
    return _uniform_fracs(n_webs)


def _default_web_positioning(n_webs: int) -> WebPositioning:
    return "fixed" if n_webs <= 2 else "uniform"


def _family_desc(y: str) -> str:
    return {
        "A": "no spar caps",
        "B": "single cap band (pitch / max-thickness)",
        "C": "discrete spar caps per web",
        "D": "continuous box spar caps",
    }[y]


_LAYOUT_REGISTRY: dict[str, SystemLayoutSpec] = {
    # 1 cell, 0 webs — outer airfoil SDF only (matches “skin-only” medial use).
    "0A": SystemLayoutSpec(
        key="0A",
        n_cells=1,
        n_webs=0,
        has_spar_caps=False,
        web_orientation="chord_normal",
        web_positioning="none",
        web_chord_fracs=(),
        spar_cap_positioning="none",
        geometry_mode="zero_web_catalogue",
        structural_family="A",
        description="Type 0A analogue: one cell, no webs, no spar caps (SDF = outer skin).",
    ),
    "0B": SystemLayoutSpec(
        key="0B",
        n_cells=1,
        n_webs=0,
        has_spar_caps=True,
        web_orientation="chord_normal",
        web_positioning="none",
        web_chord_fracs=(),
        spar_cap_positioning="max_thickness",
        geometry_mode="zero_web_catalogue",
        structural_family="B",
        description="Type 0B: spar caps without webs (MultiCellSection family B, max-thickness anchor).",
    ),
}

for _x in range(1, 6):
    _fracs = _default_web_fracs(_x)
    _pos = _default_web_positioning(_x)
    for _y in ("A", "B", "C", "D"):
        # D = continuous box between first and last web; with one web it degenerates
        # into zero-span caps and is not a valid physical configuration.
        if _x == 1 and _y == "D":
            continue
        for _z, _orient in (("CN", "chord_normal"), ("F", "flapwise")):
            _key = f"{_x}{_y}-{_z}"
            _has_caps = _y != "A"
            _cap_pos: Literal["max_thickness", "fixed", "box_spar", "none"] = (
                "none" if _y == "A" else "max_thickness"
            )
            _LAYOUT_REGISTRY[_key] = SystemLayoutSpec(
                key=_key,
                n_cells=_x + 1,
                n_webs=_x,
                has_spar_caps=_has_caps,
                web_orientation=_orient,  # type: ignore
                web_positioning=_pos,
                web_chord_fracs=_fracs,
                spar_cap_positioning=_cap_pos,
                geometry_mode="multicell",
                structural_family=_y,  # type: ignore[arg-type]
                description=(
                    f"SystemType{_x}{_y}-{_z}: {_x + 1} cells, {_x} "
                    f"{'chord-normal' if _z == 'CN' else 'flapwise'} webs, {_family_desc(_y)}."
                ),
            )


SYSTEM_TYPE_LAYOUT_KEYS: tuple[str, ...] = tuple(
    f"{core}-LC-0R" for core in sorted(_LAYOUT_REGISTRY)
)
SYSTEM_TYPE_KEYS: tuple[str, ...] = tuple(
    f"{core}-{z}-0R"
    for core in sorted(_LAYOUT_REGISTRY)
    for z in ("LC", "LE", "PA")
)


def _compose_key(core: str, anchoring: str, n_ribs: int = 0) -> str:
    tag = _ANCHOR_TAG.get(anchoring, "LC")
    return f"{core}-{tag}-{max(int(n_ribs), 0)}R"


def _bind_spec(core_spec: SystemLayoutSpec, *, anchoring: str, n_ribs: int = 0) -> SystemLayoutSpec:
    layout_key = core_spec.key if core_spec.key in _LAYOUT_REGISTRY else (
        core_spec.layout_key or core_spec.key
    )
    n = max(int(n_ribs), 0)
    return replace(
        core_spec,
        key=_compose_key(layout_key, anchoring, n),
        layout_key=layout_key,
        web_anchoring=anchoring,  # type: ignore[arg-type]
        n_ribs=n,
    )


def resolve_system_type(key: str) -> SystemLayoutSpec:
    k = (key or "").strip()
    if not k:
        raise ValueError("system_type key must be a non-empty string.")
    if k in _LAYOUT_REGISTRY:
        raise KeyError(
            f"Unknown system_type {k!r}. Anchoring and rib count are required; "
            f"try {k}-LC-0R or {k}-PA-10R."
        )
    if k.endswith("-R") and not re.search(r"-\d+R$", k):
        raise KeyError(
            f"Unknown system_type {k!r}. Rib count is required, e.g. {k[:-1]}0R."
        )
    if _MISSING_RIBS.match(k):
        raise KeyError(
            f"Unknown system_type {k!r}. Rib count is required; try {k}-0R or {k}-10R."
        )
    m = _TOKEN.match(k)
    if m is None:
        sample = ", ".join(SYSTEM_TYPE_LAYOUT_KEYS[:8])
        raise KeyError(
            f"Unknown system_type {k!r}. Expected [W][X]-[Y]-[Z]-[R] with Z in "
            f"LC|LE|PA and R a rib count including 0R, e.g. {sample}, …"
        )
    core, tag, rib = m.group(1), m.group(2), m.group(3)
    if core not in _LAYOUT_REGISTRY:
        raise KeyError(f"Unknown system_type {k!r}. No layout core {core!r}.")
    n_ribs = int(rib)
    return _bind_spec(_LAYOUT_REGISTRY[core], anchoring=_TAG_ANCHOR[tag], n_ribs=n_ribs)


# Families in descending structural elaboration. A section that cannot sustain
# its family at a given web count falls to the next family that can, never the
# other way, so the cascade is monotone.
_FAMILY_ORDER: tuple[str, ...] = ("D", "C", "B", "A")


def family_min_webs(family: str) -> int:
    """Fewest webs the catalogue defines for ``family``.

    ``D`` floors at 2: it is a continuous box between the first and last web, so
    with one web the caps degenerate to zero span. ``C`` and ``B`` floor at 1,
    ``A`` at 0.
    """
    fam = str(family).upper()
    counts = [
        spec.n_webs for spec in _LAYOUT_REGISTRY.values()
        if str(spec.structural_family).upper() == fam
    ]
    if not counts:
        raise KeyError(f"Unknown structural family {family!r}.")
    return min(counts)


def effective_layout_for_web_count(
    layout: SystemLayoutSpec, n_webs: int
) -> SystemLayoutSpec:
    """Assemble recipe for a station that can carry only ``n_webs``.

    Not a new blade concept: span freeze (Z) and rib count (R) stay on the
    parent token. Web positions are anchored spanwise, so as the section shrinks
    it stops being able to carry every web. This maps the surviving count onto
    a catalogue W/X slice, stepping **down** only.

    The rule is family-general: keep the parent family and orientation where the
    catalogue defines that web count, otherwise fall to the next less-elaborate
    family that does, and never cross the cap-bearing boundary — a section with
    spar caps does not silently lose them by stepping down. The resulting
    cascades are

    * ``A``: 5A → 4A → 3A → 2A → **1A**
    * ``B``: 5B → 4B → 3B → 2B → **1B**
    * ``C``: 5C → 4C → 3C → 2C → **1C**
    * ``D``: 5D → 4D → 3D → 2D → **1C**

    ``D`` crosses to ``C`` at one web because ``1D`` does not exist — a
    continuous box needs a fore and an aft web. ``1C`` keeps a shear web and a
    sized cap, which measurement showed is unthrottled where ``2D`` is not (see
    ``examples/output/spanwise_web_anchors/README.md``).

    **The cascade floors at one web.** ``0A`` / ``0B`` are valid layouts when
    *selected*, but they are not degradation targets: a section with no web has
    no shear path at all, so losing the last one is a different kind of event
    from losing one of several. A station that cannot anchor a single web is
    over-constrained and raises, rather than being quietly handed a web-less
    section.

    **The floor is reachable, and how much warning you get depends on where you
    start.** On the reference rotor it never fires — the survivor count bottoms
    out at one — but that is a property of that rotor, not of the rule. On a
    3.2 m blade sampled to ``r/R = 1.0`` it fires at ``r/R ~ 0.96``, where the
    chord is small enough that a minimum-gauge skin on both surfaces leaves no
    cavity to anchor into. Note this is absolute section size, not ``t/c``:
    feasibility thickening makes the section fatter and still does not save it.

    A **single-web layout has no ladder to spend**. ``step_down_cascade("1B-CN")``
    is the one-entry tuple ``('1B-CN',)``, so a run pinned there converts the
    first station that cannot anchor a web straight into a raise, with no
    degradation in between. Multi-web families absorb the same geometry by
    stepping (``5D-CN`` -> ... -> ``1C-CN``). If outboard stations must survive,
    that is an argument for starting above the floor, not for lowering it.
    """
    n = int(n_webs)
    if n < 0:
        raise ValueError(f"n_webs must be >= 0; got {n_webs!r}.")
    if n > int(layout.n_webs):
        raise ValueError(
            f"step-down only: {layout.key} has {layout.n_webs} webs, asked for {n}."
        )
    if n == int(layout.n_webs):
        if layout.key in _LAYOUT_REGISTRY:
            return _bind_spec(layout, anchoring=layout.web_anchoring, n_ribs=layout.n_ribs)
        return layout

    want_caps = bool(layout.has_spar_caps)
    if n == 0:
        # Reachable only from a layout that had webs to lose: a 0-web parent is
        # returned unchanged by the equality branch above.
        raise InfeasibleSectionGeometry(
            f"{layout.key} cannot anchor a single web at this station. The "
            "step-down cascade floors at one web, because a section with no web "
            "has no shear path; 0A/0B are selectable layouts, not degradation "
            "targets. Report this station as over-constrained."
        )

    orient = "F" if str(layout.web_orientation) == "flapwise" else "CN"
    parent = str(layout.structural_family).upper()
    start = _FAMILY_ORDER.index(parent) if parent in _FAMILY_ORDER else 0
    for fam in _FAMILY_ORDER[start:]:
        candidate = _LAYOUT_REGISTRY.get(f"{n}{fam}-{orient}")
        if candidate is not None and bool(candidate.has_spar_caps) == want_caps:
            return _bind_spec(
                candidate, anchoring=layout.web_anchoring, n_ribs=layout.n_ribs,
            )
    raise KeyError(
        f"No catalogue layout for {n} web(s) stepping down from {layout.key} "
        f"(orientation {orient}, spar caps {want_caps})."
    )


def layout_without_spar_caps(layout: SystemLayoutSpec) -> SystemLayoutSpec:
    """The cap-less sibling of ``layout`` — same web count and orientation, family ``A``.

    **This is a different event from the web step-down cascade, and deliberately
    not part of it.** :func:`effective_layout_for_web_count` never crosses the
    cap-bearing boundary, because losing caps is not a consequence of losing
    webs. This function exists for the other case: the local cavity cannot
    contain a minimum-gauge spar cap, so the cap *runs out* — which is what real
    blades do outboard, and what the alternative (sizing a cap thinner than the
    material can be laid up) is not (see ``docs/section/bug_log.md``).

    The cascade's floor-at-one-web rule does not apply here and its rationale is
    worth separating: a section with no web has no shear path, so losing the last
    one is an over-constraint. A cap-less section keeps every web it had — family
    ``A`` is webbed — so the shear path is untouched and the section remains
    buildable. Only the flapwise bending contribution moves to the skin.

    Zero-web layouts are returned by their registry entry (``0B`` -> ``0A``);
    a layout that already has no caps is returned unchanged.
    """
    if not bool(layout.has_spar_caps):
        return layout
    n = int(layout.n_webs)
    if n == 0:
        capless = _LAYOUT_REGISTRY.get("0A")
        if capless is None:  # pragma: no cover - registry is static
            raise KeyError("No cap-less zero-web layout (0A) in the registry.")
        return _bind_spec(capless, anchoring=layout.web_anchoring, n_ribs=layout.n_ribs)
    orient = "F" if str(layout.web_orientation) == "flapwise" else "CN"
    key = f"{n}A-{orient}"
    capless = _LAYOUT_REGISTRY.get(key)
    if capless is None:
        raise KeyError(
            f"No cap-less layout {key!r} to run {layout.key} out into."
        )
    return _bind_spec(capless, anchoring=layout.web_anchoring, n_ribs=layout.n_ribs)


def step_down_cascade(layout: SystemLayoutSpec) -> tuple[str, ...]:
    """Full ladder of **assemble-recipe** keys from ``layout`` down to its floor.

    These tokens name station W/X slices (with parent Z and R copied), not new
    blade concepts. Floors at one web — see :func:`effective_layout_for_web_count`.
    A layout that already has no webs is its own terminal rung.
    """
    floor = 1 if int(layout.n_webs) >= 1 else 0
    return tuple(
        effective_layout_for_web_count(layout, n).key
        for n in range(int(layout.n_webs), floor - 1, -1)
    )


__all__ = [
    "SYSTEM_TYPE_KEYS",
    "SYSTEM_TYPE_LAYOUT_KEYS",
    "StructuralFamily",
    "SystemLayoutSpec",
    "WebAnchoring",
    "WebOrientation",
    "WebPositioning",
    "effective_layout_for_web_count",
    "family_min_webs",
    "layout_without_spar_caps",
    "resolve_system_type",
    "step_down_cascade",
]
