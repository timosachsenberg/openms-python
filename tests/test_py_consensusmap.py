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
