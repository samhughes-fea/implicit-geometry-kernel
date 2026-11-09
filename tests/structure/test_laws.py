import pickle

import numpy as np
import pytest

from implicit_geometry_kernel.structure import ConstantLaw, ProductLaw, ScaledLaw, SumLaw, TableLaw, as_law


def test_constant_and_combinators():
    a, b = ConstantLaw(2.0), TableLaw([0.0, 1.0], [1.0, 3.0])
    z = np.array([0.0, 0.5, 1.0])
    assert np.allclose(SumLaw(a, b)(z), [3.0, 4.0, 5.0])
    assert np.allclose(ProductLaw(a, b)(z), [2.0, 4.0, 6.0])
    assert np.allclose(ScaledLaw(b, 2.0, 1.0)(z), [3.0, 5.0, 7.0])
    assert np.allclose((b * 0.5)(z), [0.5, 1.0, 1.5])
    assert as_law(a) is a and as_law(1.5).at(9.0) == 1.5


def test_single_plane_batch_matches_pointwise():
    law = TableLaw([0.0, 2.0, 5.0], [0.03, 0.02, 0.005])
    z = np.full(7, 1.3)
    assert np.allclose(law(z), law.at(1.3))


def test_table_law_hits_knots_is_linear_and_clamped():
    zk = [0.0, 1.0, 2.0, 4.0, 6.0]
    vk = [0.030, 0.030, 0.020, 0.008, 0.008]
    law = TableLaw(zk, vk)
    assert np.allclose(law(np.array(zk)), vk, rtol=0, atol=1e-15)
    z = np.linspace(-1.0, 7.0, 4001)
    v = law(z)
    assert np.all(np.diff(v) <= 1e-15)
    assert v.min() >= min(vk) - 1e-15 and v.max() <= max(vk) + 1e-15
    assert law.at(-1.0) == vk[0] and law.at(7.0) == vk[-1]
    assert law.at(1.5) == pytest.approx(0.025)
    assert law.at(3.0) == pytest.approx(0.014)


def test_table_law_ply_drop_plateaus():
    law = TableLaw([0.0, 2.0, 3.0, 6.0], [0.030, 0.030, 0.008, 0.008])
    z = np.linspace(0.0, 6.0, 61)
    v = law(z)
    assert np.allclose(v[z <= 2.0], 0.030)
    assert np.allclose(v[z >= 3.0], 0.008)
    mid = v[(z > 2.0) & (z < 3.0)]
    assert np.all(np.diff(mid) < 0.0)


@pytest.mark.parametrize("args", [([], []), ([0.0, 1.0], [1.0]), ([0.0, 0.0], [1.0, 2.0]), ([0.0, np.nan], [1.0, 2.0])])
def test_bad_tables_are_rejected(args):
    with pytest.raises(ValueError):
        TableLaw(*args)


def test_laws_pickle():
    law = SumLaw(TableLaw([0.0, 1.0], [1.0, 2.0]), ConstantLaw(0.5))
    again = pickle.loads(pickle.dumps(law))
    assert np.allclose(again(np.linspace(0, 1, 5)), law(np.linspace(0, 1, 5)))
