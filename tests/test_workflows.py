import math
from pathlib import Path

import pandas as pd
import pyopenms as oms

from openms_python import (
    Py_MSExperiment,
    Py_FeatureMap,
    Identifications,
    ProteinStream,
    map_identifications_to_features,
    align_feature_maps,
    link_features,
    export_quant_table,
)


def _make_simple_experiment() -> Py_MSExperiment:
    exp = oms.MSExperiment()
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(1)
    spectrum.setRT(10.0)
    spectrum.set_peaks(([100.0, 101.0, 102.0], [1000.0, 1500.0, 1300.0]))
    exp.addSpectrum(spectrum)
    exp.updateRanges()
    return Py_MSExperiment(exp)


def test_detect_features_returns_feature_map():
    experiment = _make_simple_experiment()
    feature_map = experiment.detect_features()
    assert isinstance(feature_map, Py_FeatureMap)
    assert len(feature_map) >= 0


def test_protein_stream_pipeline_generates_spectra(tmp_path):
    fasta_content = ">test\nMPEPTIDEK\n"
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_content)

    stream = ProteinStream.from_fasta(fasta_path)
    peptides = stream.digest(min_length=3, max_length=15)
    spectra = peptides.theoretical_spectra(min_charge=1, max_charge=1)
    record = next(iter(spectra))
    assert record.peptide.size() >= 3
    assert record.spectrum.size() > 0


def _make_identifications() -> Identifications:
    target = oms.PeptideIdentification()
    target.setHigherScoreBetter(False)
    target_hit = oms.PeptideHit()
    target_hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    target_hit.setScore(0.01)
    target_hit.setMetaValue("target_decoy", "target")
    target.setHits([target_hit])

    decoy = oms.PeptideIdentification()
    decoy.setHigherScoreBetter(False)
    decoy_hit = oms.PeptideHit()
    decoy_hit.setSequence(oms.AASequence.fromString("PEPTIDD"))
    decoy_hit.setScore(0.02)
    decoy_hit.setMetaValue("target_decoy", "decoy")
    decoy.setHits([decoy_hit])

    return Identifications([], [target, decoy])


def test_identifications_filter_by_fdr_filters_decoys():
    ids = _make_identifications()
    filtered = ids.filter_by_fdr(threshold=0.5)
    assert len(filtered.peptides) == 1
    top_hit = filtered.peptides[0].getHits()[0]
    assert str(top_hit.getMetaValue("target_decoy")) == "target"


def test_identifications_filter_by_fdr_handles_protein_level():
    target = oms.ProteinIdentification()
    target.setHigherScoreBetter(False)
    target_hit = oms.ProteinHit()
    target_hit.setScore(0.01)
    target_hit.setAccession("P01234")
    target_hit.setMetaValue("target_decoy", "target")
    target.setHits([target_hit])

    decoy = oms.ProteinIdentification()
    decoy.setHigherScoreBetter(False)
    decoy_hit = oms.ProteinHit()
    decoy_hit.setScore(0.02)
    decoy_hit.setAccession("DECOY_P01234")
    decoy_hit.setMetaValue("target_decoy", "decoy")
    decoy.setHits([decoy_hit])

    ids = Identifications([target, decoy], [])
    filtered = ids.filter_by_fdr(level="protein", threshold=0.5)
    assert len(filtered.proteins) == 1
    assert filtered.proteins[0].getHits()[0].getAccession() == "P01234"


def test_map_identifications_to_features_sets_annotation():
    feature = oms.Feature()
    feature.setRT(10.0)
    feature.setMZ(500.0)
    feature.setIntensity(1000.0)
    feature_map = oms.FeatureMap()
    feature_map.push_back(feature)

    ids = _make_identifications()
    for pep in ids.peptides:
        pep.setRT(10.0)
        pep.setMZ(500.0)

    annotated = map_identifications_to_features(Py_FeatureMap(feature_map), ids)
    assert len(annotated) == 1
    assert annotated[0].native.getPeptideIdentifications()


def test_align_feature_maps_shifts_second_map():
    fmap1 = oms.FeatureMap()
    base = oms.Feature()
    base.setRT(10.0)
    base.setMZ(500.0)
    base.setIntensity(100.0)
    fmap1.push_back(base)

    fmap2 = oms.FeatureMap()
    shifted = oms.Feature()
    shifted.setRT(15.0)
    shifted.setMZ(500.0)
    shifted.setIntensity(80.0)
    fmap2.push_back(shifted)

    aligned = align_feature_maps([Py_FeatureMap(fmap1), Py_FeatureMap(fmap2)])
    assert math.isclose(aligned[1][0].native.getRT(), aligned[0][0].native.getRT(), rel_tol=1e-6, abs_tol=1e-6)


def test_link_features_and_export_quant_table():
    fmap_a = oms.FeatureMap()
    feat_a = oms.Feature()
    feat_a.setRT(10.0)
    feat_a.setMZ(500.0)
    feat_a.setIntensity(100.0)
    fmap_a.push_back(feat_a)

    fmap_b = oms.FeatureMap()
    feat_b = oms.Feature()
    feat_b.setRT(10.2)
    feat_b.setMZ(500.1)
    feat_b.setIntensity(110.0)
    fmap_b.push_back(feat_b)

    consensus = link_features([Py_FeatureMap(fmap_a), Py_FeatureMap(fmap_b)])
    df = export_quant_table(consensus)
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] >= 1
    # Expect one column per input map
    assert {col for col in df.columns if col.startswith("map_")}
