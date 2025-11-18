import pandas as pd
import pytest

oms = pytest.importorskip("pyopenms")

from openms_python.py_consensusmap import Py_ConsensusMap
from openms_python.py_identifications import Identifications


def build_consensus_map(count: int = 3) -> Py_ConsensusMap:
    cmap = oms.ConsensusMap()
    for idx in range(count):
        feature = oms.ConsensusFeature()
        feature.setUniqueId(idx)
        feature.setIntensity(100.0 + idx)
        feature.setMZ(400.0 + idx)
        cmap.push_back(feature)
    return Py_ConsensusMap(cmap)


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


def test_py_consensusmap_is_iterable():
    cmap = build_consensus_map(4)

    collected = list(cmap)

    assert len(collected) == 4
    assert all(isinstance(feature, oms.ConsensusFeature) for feature in collected)
    assert [feature.getUniqueId() for feature in collected] == [0, 1, 2, 3]


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


def test_py_consensusmap_dataframe_roundtrip_preserves_meta():
    cmap = Py_ConsensusMap()
    for idx in range(2):
        feature = oms.ConsensusFeature()
        feature.setUniqueId(idx)
        feature.setRT(100.0 + idx)
        feature.setMZ(500.0 + idx)
        feature.setIntensity(1000.0 + (idx * 10))
        feature.setCharge(2)
        feature.setWidth(0.5)
        feature.setQuality(0.9)
        feature.setMetaValue("note", f"c{idx}")
        cmap.append(feature)

    df = cmap.to_dataframe()
    assert set(["rt", "mz", "intensity", "note"]).issubset(df.columns)

    df["mz"] += 0.1
    rebuilt = Py_ConsensusMap.from_dataframe(df)

    assert len(rebuilt) == len(cmap)
    assert rebuilt[0].getMetaValue("note") == "c0"
    assert rebuilt[1].getMZ() == pytest.approx(cmap[1].getMZ() + 0.1)


def test_py_consensusmap_from_dataframe_requires_columns():
    df = pd.DataFrame({"rt": [1.0], "mz": [2.0]})

    with pytest.raises(ValueError):
        Py_ConsensusMap.from_dataframe(df)


def test_py_consensusmap_infer_proteins_uses_consensus_map():
    protein = _protein("run", "P10")
    pep_a = _peptide("run", "PEPA", 5.0, "P10")
    pep_b = _peptide("run", "PEPB", 30.0, "P10")

    feature = oms.ConsensusFeature()
    feature.setPeptideIdentifications([pep_a, pep_b])

    cmap = oms.ConsensusMap()
    cmap.setProteinIdentifications([protein])
    cmap.push_back(feature)

    wrapper = Py_ConsensusMap(cmap)
    inferred = wrapper.infer_proteins(algorithm="basic", include_unassigned=True)

    assert isinstance(inferred, Identifications)
    assert inferred.protein_identifications[0].getHits()[0].getScore() == pytest.approx(30.0)


def test_py_consensusmap_infer_protein_quantities(monkeypatch):
    captured = {}

    class DummyInference:
        def __init__(self):
            captured["instantiated"] = True

        def infer(self, consensus_map, reference_map):
            captured["consensus_map"] = consensus_map
            captured["reference_map"] = reference_map

    monkeypatch.setattr(oms, "ProteinInference", lambda: DummyInference())

    wrapper = Py_ConsensusMap()
    result = wrapper.infer_protein_quantities(reference_map=3)

    assert result is wrapper
    assert captured["instantiated"] is True
    assert captured["consensus_map"] is wrapper.native
    assert captured["reference_map"] == 3
