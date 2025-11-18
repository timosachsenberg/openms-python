import pandas as pd
import pytest

oms = pytest.importorskip("pyopenms")

from openms_python.py_feature import Py_Feature
from openms_python.py_featuremap import Py_FeatureMap
from openms_python.py_identifications import Identifications


def build_feature_map(count: int = 4) -> Py_FeatureMap:
    feature_map = oms.FeatureMap()
    for idx in range(count):
        feature = oms.Feature()
        feature.setUniqueId(idx)
        feature.setRT(100.0 + idx)
        feature.setMZ(400.0 + idx)
        feature_map.push_back(feature)
    return Py_FeatureMap(feature_map)


def _protein(identifier: str, accession: str) -> oms.ProteinIdentification:
    entry = oms.ProteinIdentification()
    entry.setIdentifier(identifier)
    hit = oms.ProteinHit()
    hit.setAccession(accession)
    entry.setHits([hit])
    return entry


def _peptide(identifier: str, sequence: str, score: float, accession: str) -> oms.PeptideIdentification:
    entry = oms.PeptideIdentification()
    entry.setIdentifier(identifier)
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(sequence))
    hit.setScore(score)
    evidence = oms.PeptideEvidence()
    evidence.setProteinAccession(accession)
    hit.setPeptideEvidences([evidence])
    entry.setHits([hit])
    return entry


def test_py_featuremap_len_and_indexing():
    fmap = build_feature_map()

    assert len(fmap) == 4
    assert isinstance(fmap[0], Py_Feature)
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


def test_py_featuremap_iteration_returns_py_features():
    fmap = build_feature_map(3)

    collected = list(fmap)

    assert all(isinstance(feature, Py_Feature) for feature in collected)
    assert [feature.getUniqueId() for feature in collected] == [0, 1, 2]


def test_py_featuremap_load_store_roundtrip(tmp_path):
    fmap = build_feature_map()
    out_path = tmp_path / "features.featureXML"

    fmap.store(out_path)
    assert out_path.exists()

    loaded = Py_FeatureMap().load(out_path)
    assert len(loaded) == len(fmap)
    assert [feat.getUniqueId() for feat in loaded] == [feat.getUniqueId() for feat in fmap]


def test_py_featuremap_rejects_unknown_extension(tmp_path):
    fmap = build_feature_map()
    bad_path = tmp_path / "features.txt"

    with pytest.raises(ValueError):
        fmap.store(bad_path)

    with pytest.raises(ValueError):
        Py_FeatureMap().load(bad_path)


def test_py_featuremap_append_and_extend():
    fmap = Py_FeatureMap()
    feature = oms.Feature()
    feature.setUniqueId(10)
    fmap.append(feature)

    assert len(fmap) == 1
    assert fmap[0].getUniqueId() == 10

    second = oms.Feature()
    second.setUniqueId(11)
    third = oms.Feature()
    third.setUniqueId(12)

    fmap.extend([second, third])
    assert [feat.getUniqueId() for feat in fmap] == [10, 11, 12]


def test_py_featuremap_remove_and_delete():
    fmap = build_feature_map(5)

    fmap.remove(0)
    assert [feat.getUniqueId() for feat in fmap] == [1, 2, 3, 4]

    del fmap[-1]
    assert [feat.getUniqueId() for feat in fmap] == [1, 2, 3]

    new_feature = oms.Feature()
    new_feature.setUniqueId(10)
    other_feature = oms.Feature()
    other_feature.setUniqueId(11)
    fmap.extend([new_feature, other_feature])

    del fmap[1:4:2]
    assert [feat.getUniqueId() for feat in fmap] == [1, 3, 11]

    extra = []
    for uid in (12, 13):
        feat = oms.Feature()
        feat.setUniqueId(uid)
        extra.append(feat)
    fmap.extend(extra)

    del fmap[::-2]
    assert [feat.getUniqueId() for feat in fmap] == [3, 12]

    with pytest.raises(IndexError):
        fmap.remove(10)

    with pytest.raises(TypeError):
        del fmap[None]


def test_py_featuremap_dataframe_round_trip():
    df = pd.DataFrame(
        {
            "unique_id": [101, 102],
            "rt": [100.0, 101.0],
            "mz": [400.0, 401.0],
            "intensity": [1000.0, 2000.0],
            "charge": [2, 3],
            "width": [0.5, 0.6],
            "overall_quality": [0.9, 0.8],
            "annotation": ["feature_a", "feature_b"],
        }
    )

    fmap = Py_FeatureMap.from_dataframe(df)
    assert len(fmap) == 2
    assert fmap[0].getUniqueId() == 101
    assert fmap[0].getMetaValue("annotation") == "feature_a"

    recreated = fmap.to_dataframe()
    assert recreated.shape[0] == 2
    for column in ["unique_id", "rt", "mz", "intensity", "charge", "width", "overall_quality"]:
        assert column in recreated.columns
    assert recreated["rt"].tolist() == [100.0, 101.0]
    assert recreated["annotation"].tolist() == ["feature_a", "feature_b"]

    alias = Py_FeatureMap.from_df(df)
    assert len(alias) == 2


def test_py_featuremap_from_dataframe_missing_columns():
    df = pd.DataFrame({"rt": [1.0], "mz": [2.0]})

    with pytest.raises(ValueError):
        Py_FeatureMap.from_dataframe(df)


def test_py_featuremap_infer_proteins_collects_assigned_ids():
    protein = _protein("run", "P10")
    pep_a = _peptide("run", "PEPA", 5.0, "P10")
    pep_b = _peptide("run", "PEPB", 25.0, "P10")

    feature = oms.Feature()
    feature.setPeptideIdentifications([pep_a, pep_b])

    fmap = oms.FeatureMap()
    fmap.setProteinIdentifications([protein])
    fmap.push_back(feature)

    wrapper = Py_FeatureMap(fmap)
    inferred = wrapper.infer_proteins(algorithm="basic")

    assert isinstance(inferred, Identifications)
    assert inferred.protein_identifications[0].getHits()[0].getScore() == pytest.approx(25.0)
