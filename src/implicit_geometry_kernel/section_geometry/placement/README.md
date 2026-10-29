# Placement geometry (WP3+)

Anchored wall geometry for surface-aware ply drop.

## Ownership

| Concern | Package |
|---------|---------|
| Ply order / drop face / policy digest | `material_utilities.layup.types.placement` |
| Physical face offsets / midsurface depth | this package (`section_geometry.placement`) |
| Stamp onto `ShellMidlineStrip` | `stamp_strip_identities` (gauge) + `placement.handoff` |
| Cache / lifecycle | `section_selection` (`placement_lifecycle`, `shell_prep_cache_key`) |

`section_geometry` must not import `material_utilities`. Construction tokens here are
plain strings that match the layup enums 1:1.

## Defaults (COTIDE Sheffield 2026 / BL-028)

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
