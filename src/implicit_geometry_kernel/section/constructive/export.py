"""Former path of the wall and core polygon export, now ``section.boundary.walls`` and
``section.boundary.cores``; kept as an alias until 0.6.0."""

from implicit_geometry_kernel.section.boundary.walls import (  # noqa: F401
    _label_key,
    _resolve_host,
    _component_map,
    _union_discrete_caps,
    _cap_from_recipes,
    _shared_deeper_by_surface,
    _web_union,
    _spar_cap_label_for_surface,
    _host_recipes,
    _host_cap_recipes,
    _recipe_for,
    _skin_via_recipe,
    _web_poly_via_recipe,
    _cap_poly_via_recipe,
    _core_poly_via_recipe,
    _insert_poly_via_recipe,
    _is_skin_key,
    _is_web_key,
    emit_wall_solid_polygons,
)
from implicit_geometry_kernel.section.boundary.cores import (  # noqa: F401
    emit_core_polygons,
    _keep_poly,
    emit_constructive_polygons,
)
