import pytest

oms = pytest.importorskip("pyopenms")

from openms_python.py_featuremap import Py_FeatureMap


def build_feature_map(count: int = 4) -> Py_FeatureMap:
    feature_map = oms.FeatureMap()
    for idx in range(count):
        feature = oms.Feature()
        feature.setUniqueId(idx)
        feature.setRT(100.0 + idx)
        feature.setMZ(400.0 + idx)
        feature_map.push_back(feature)
    return Py_FeatureMap(feature_map)


def test_py_featuremap_len_and_indexing():
    fmap = build_feature_map()

    assert len(fmap) == 4
    assert fmap[0].getUniqueId() == 0
    assert fmap[-1].getUniqueId() == 3

    sliced = fmap[1:4:2]
    assert isinstance(sliced, Py_FeatureMap)
    assert len(sliced) == 2
    assert [feat.getUniqueId() for feat in sliced] == [1, 3]

    reversed_slice = fmap[::-1]
    assert [feat.getUniqueId() for feat in reversed_slice] == [3, 2, 1, 0]

    assert len(fmap[10:15]) == 0

    with pytest.raises(IndexError):
        _ = fmap[10]

    with pytest.raises(TypeError):
        _ = fmap[None]
