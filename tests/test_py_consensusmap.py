import pytest

oms = pytest.importorskip("pyopenms")

from openms_python.py_consensusmap import Py_ConsensusMap


def build_consensus_map(count: int = 3) -> Py_ConsensusMap:
    cmap = oms.ConsensusMap()
    for idx in range(count):
        feature = oms.ConsensusFeature()
        feature.setUniqueId(idx)
        feature.setIntensity(100.0 + idx)
        feature.setMZ(400.0 + idx)
        cmap.push_back(feature)
    return Py_ConsensusMap(cmap)


def test_py_consensusmap_len_and_indexing():
    cmap = build_consensus_map()

    assert len(cmap) == 3
    assert cmap[0].getUniqueId() == 0
    assert cmap[-1].getUniqueId() == 2

    reversed_slice = cmap[::-1]
    assert [feat.getUniqueId() for feat in reversed_slice] == [2, 1, 0]

    sliced = cmap[0:2]
    assert isinstance(sliced, Py_ConsensusMap)
    assert len(sliced) == 2

    with pytest.raises(IndexError):
        _ = cmap[3]


def test_py_consensusmap_load_store_roundtrip(tmp_path):
    cmap = build_consensus_map()
    out_path = tmp_path / "consensus.consensusXML"

    cmap.store(out_path)
    assert out_path.exists()

    loaded = Py_ConsensusMap().load(out_path)
    assert len(loaded) == len(cmap)
    assert [feat.getUniqueId() for feat in loaded] == [feat.getUniqueId() for feat in cmap]


def test_py_consensusmap_rejects_unknown_extension(tmp_path):
    cmap = build_consensus_map()
    bad_path = tmp_path / "consensus.featureXML"

    with pytest.raises(ValueError):
        cmap.store(bad_path)

    with pytest.raises(ValueError):
        Py_ConsensusMap().load(bad_path)


def test_py_consensusmap_append_and_extend():
    cmap = Py_ConsensusMap()
    feature = oms.ConsensusFeature()
    feature.setUniqueId(42)
    cmap.append(feature)

    assert len(cmap) == 1
    assert cmap[0].getUniqueId() == 42

    other_features = []
    for uid in (43, 44):
        entry = oms.ConsensusFeature()
        entry.setUniqueId(uid)
        other_features.append(entry)

    cmap.extend(other_features)
    assert [feat.getUniqueId() for feat in cmap] == [42, 43, 44]


def test_py_consensusmap_remove_and_delete():
    cmap = build_consensus_map(5)

    cmap.remove(0)
    assert [feat.getUniqueId() for feat in cmap] == [1, 2, 3, 4]

    del cmap[-2]
    assert [feat.getUniqueId() for feat in cmap] == [1, 2, 4]

    extras = []
    for uid in (10, 11):
        feat = oms.ConsensusFeature()
        feat.setUniqueId(uid)
        extras.append(feat)
    cmap.extend(extras)

    del cmap[::2]
    assert [feat.getUniqueId() for feat in cmap] == [2, 10]

    cmap.append(oms.ConsensusFeature())
    cmap[-1].setUniqueId(12)

    del cmap[::-1]
    assert len(cmap) == 0

    with pytest.raises(IndexError):
        cmap.remove(0)

    with pytest.raises(TypeError):
        del cmap[object()]


def test_py_consensusmap_is_iterable():
    cmap = build_consensus_map(4)

    uids = [feature.getUniqueId() for feature in cmap]
    assert uids == [0, 1, 2, 3]
