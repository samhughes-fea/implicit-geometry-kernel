# Wall recipes

One recipe description, three projections (G-RECIPE SoT):

| Projection | Status |
|---|---|
| Boundary polygon | Delegates to constructive emitters (`emit_*_polygon`) |
| Midline | Wired via `to_midline()` / component `midline_polyline()` for skin / web / cap |
| SDF | `to_sdf_expr` / `to_sdf_component` → private IR classes next to `walls.py` |

Recipes hold parameters only. Assemble authors frozen recipes first; SDF IR
classes (`OuterSkin`, `SparCap`, `ShearWebComponent`, `SandwichCore`, `LEInsert`,
`TEInsert`) are private projection internals (not a second authoring
surface). Import them from `recipes.outer_skin` and siblings, not the recipes
barrel. The barrel `__all__` is wall recipes only.

Emit / hybrid props consume recipes stored on the host.

| Recipe | Boundary | Midline | SDF IR |
|---|---|---|---|
| `SkinWall` | constructive | yes | `OuterSkin` |
| `WebWall` | constructive | yes | `ShearWebComponent` |
| `CapWall` | constructive | yes | `SparCap` |
| `CoreWall` | constructive | n/a (bay fill) | `SandwichCore` |
| `InsertWall` | constructive (default) | n/a | `LEInsert` / `TEInsert` |

`INSERTS_CONSTRUCTIVE_BY_DEFAULT` is **True**: LE/TE inserts use the constructive
polygon path in hybrid properties. Set `emission_maturity='residual_phi'` only
for debug / residual-φ comparison.

`from_component` helpers on wall recipes are **test/migration only** — production
assemble must not call them.

## Outer and inner mould lines

- **OML:** aerofoil contour `{φ_af = 0}`.
- **IML:** inward offset by `skin_thickness` (cavity boundary after the skin slab).

Skin is an inward slab from the OML. Spar caps sit on the IML. Cores fill the
cavity minus wall solids (skin / webs / caps).
