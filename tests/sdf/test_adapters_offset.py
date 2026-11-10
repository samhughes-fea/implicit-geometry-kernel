"""SDFCallable/SDFNode plane adapter, tagged IR, and spanwise offsets."""

import numpy as np

from implicit_geometry_kernel.sdf import sdf_circle
from implicit_geometry_kernel.sdf.nodes.tagged import FieldOffsetNode, MidSurfaceNode, as_node
from implicit_geometry_kernel.sdf.nodes import Circle, Offset
from implicit_geometry_kernel.section.walls.base import SectionComponent
from implicit_geometry_kernel.section.walls.tagged import tagged_component


def test_as_xy_matches_primitive():
    node = Circle(0.0, 0.0, 1.0)
    f = node.as_xy(z=4.0)
    x, y = np.linspace(-1.5, 1.5, 7), np.zeros(7)
    np.testing.assert_allclose(f(x, y), sdf_circle(x, y, 0.0, 0.0, 1.0))


def test_as_node_round_trip():
    node = as_node(lambda x, y: sdf_circle(x, y, 0.0, 0.0, 1.0), meta=3)
    s = node.evaluate(np.array([[0.0, 0.0, 2.0]]))
    assert s.phi[0] < 0.0 and s.meta[0] == 3


def test_tagged_component_stamps_skin():
    comp = SectionComponent(Offset(Circle(0.0, 0.0, 1.0), -0.05))
    comp.label = "skin_suction_surface"
    node, table = tagged_component(comp)
    assert table.by_name("skin_suction_surface").kind == "skin"
    inside = node.evaluate(np.array([[0.97, 0.0, 0.0]]))
    assert inside.meta[0] == 0


def test_field_offset_and_mid_surface():
    oml = Circle(0.0, 0.0, 1.0)
    cavity = FieldOffsetNode(oml, lambda p: np.full(len(p), 0.1), scale=1.0)
    mid = MidSurfaceNode(oml, cavity)
    on_mid = mid.evaluate(np.array([[0.95, 0.0, 0.0]]))
    assert abs(on_mid.phi[0]) < 1e-12
    np.testing.assert_allclose(mid.thickness(np.array([[0.95, 0.0, 0.0]])), 0.1)
