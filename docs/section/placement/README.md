# Placement geometry (WP3+)

**Audience: internal.** Maintainer notes; not part of the library manual.

Moved to `section/layout/` (`anchoring.py`, `collision.py`, `handoff.py`).
The `section.placement` import path was removed in 0.6.0.

Anchored wall geometry for surface-aware ply drop.

## Ownership

| Concern | Package |
|---------|---------|
| Ply order / drop face / policy digest | downstream layup library |
| Physical face offsets / midsurface depth | this package (`section.layout`) |
| Stamp onto `ShellMidlineStrip` | `section.layout.handoff` |

This package must not import a layup library. Construction tokens here are
plain strings.

## Defaults

| Wall | Reference | Build | Midline |
|------|-----------|-------|---------|
| skin | outer mould line | inward | `t/2` from OML |
| spar cap | skin bond line | into cavity | `t_skin + h/2` |
| shear web | midsurface | outward symmetric | web axis |

Named non-default skins: `inner_mould_line` / `outward`, and `midsurface` / `centred`.

## Invalidation

A placement digest change with the same zone set is **geometry-only** (rebuild
geometry → mesh → bind → FE → evidence). Zone identity / connectivity / presence
changes bump topology generation (`placement_lifecycle.classify_placement_change`).
