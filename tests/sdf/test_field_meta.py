"""Feature ids on SDF samples: TagNode and PriorityUnionNode."""

import numpy as np

from implicit_geometry_kernel.sdf import FeatureTable, PriorityUnionNode, TagNode
from implicit_geometry_kernel.sdf.nodes import Circle


def _xy(xs, ys, z=0.0):
    pts = np.column_stack([np.asarray(xs, float), np.asarray(ys, float), np.full(len(xs), z)])
    return pts


def test_feature_table_reuses_names():
    table = FeatureTable()
    a = table.add("skin_suction_surface", "skin")
    b = table.add("skin_suction_surface", "skin")
    c = table.add("shear_web_1", "web")
    assert a.id == b.id == 0
    assert c.id == 1
    assert c.priority > a.priority


def test_tag_node_stamps_meta():
    tagged = TagNode(Circle(0.0, 0.0, 1.0), 7)
    s = tagged.evaluate(_xy([0.0], [0.0]))
    assert s.meta[0] == 7
    assert s.phi[0] < 0.0


def test_priority_union_inside_prefers_cap():
    skin = TagNode(Circle(0.0, 0.0, 1.0), 1)
    cap = TagNode(Circle(0.0, 0.0, 0.4), 2)
    field = PriorityUnionNode([skin, cap], [1, 3])
    inside_cap = field.evaluate(_xy([0.0], [0.0]))
    in_skin_only = field.evaluate(_xy([0.7], [0.0]))
    outside = field.evaluate(_xy([3.0], [0.0]))
    assert inside_cap.meta[0] == 2 and inside_cap.phi[0] < 0.0
    assert in_skin_only.meta[0] == 1 and in_skin_only.phi[0] < 0.0
    assert outside.phi[0] > 0.0
    assert outside.meta[0] == 1  # nearest child
