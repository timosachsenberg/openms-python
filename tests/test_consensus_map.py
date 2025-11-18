import pyopenms as oms
import pytest

from openms_python import Py_FeatureMap, Py_ConsensusMap


def _simple_feature_map(rt: float) -> oms.FeatureMap:
    fmap = oms.FeatureMap()
    feature = oms.Feature()
    feature.setRT(rt)
    feature.setMZ(500.0)
    feature.setIntensity(100.0)
    fmap.push_back(feature)
    return fmap


def test_align_and_link_identity_creates_consensus():
    fmap_a = Py_FeatureMap(_simple_feature_map(10.0))
    fmap_b = Py_FeatureMap(_simple_feature_map(10.0))

    consensus = Py_ConsensusMap.align_and_link(
        [fmap_a, fmap_b],
        alignment_method="identity",
    )

    assert isinstance(consensus, Py_ConsensusMap)
    assert len(consensus) == 1
    headers = consensus.native.getColumnHeaders()
    assert len(headers) == 2


def test_align_and_link_invalid_method_raises():
    fmap = Py_FeatureMap(_simple_feature_map(10.0))
    with pytest.raises(ValueError):
        Py_ConsensusMap.align_and_link([fmap], alignment_method="unknown")
