# Wall recipes

One recipe description, three projections (G-RECIPE SoT):

| Projection | Status |
|---|---|
| Boundary polygon | Delegates to constructive emitters (`emit_*_polygon`) |
| Midline | Wired via `to_midline()` / component `midline_polyline()` for skin / web / cap |
| SDF | `to_sdf_expr` / `to_sdf_component` → private IR under [`sdf_ir/`](sdf_ir/) |

Recipes hold parameters only. Assemble authors frozen recipes first; SDF IR
classes in `sdf_ir/` are private projection internals (not a second authoring
surface). Emit / hybrid props consume recipes stored on the host.

| Recipe | Boundary | Midline | SDF IR |
|---|---|---|---|
| `SkinWall` | constructive | yes | `OuterSkin` |
| `WebWall` | constructive | yes | `ShearWeb` |
| `CapWall` | constructive | yes | `SparCap` |
| `CoreWall` | constructive | n/a (bay fill) | `SandwichCore` |
| `InsertWall` | constructive (default) | n/a | `LEInsert` / `TEInsert` |

`INSERTS_CONSTRUCTIVE_BY_DEFAULT` is **True**: LE/TE inserts use the constructive
polygon path in hybrid properties. Set `emission_maturity='residual_phi'` only
for debug / residual-φ comparison.

`from_component` helpers on wall recipes are **test/migration only** — production
assemble must not call them.
