"""Tests for the high-level workflow helpers."""
from __future__ import annotations

import numpy as np
import pyopenms as oms

from openms_python import FastaWorkflow, Identifications, Py_FeatureMap, Py_MSExperiment


def _make_experiment() -> Py_MSExperiment:
    exp = oms.MSExperiment()
    for rt, peaks in [(10.0, [(100.0, 200.0), (101.0, 50.0)]), (15.0, [(200.0, 500.0)])]:
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(rt)
        mz, intensity = zip(*peaks)
        spec.set_peaks((np.array(mz), np.array(intensity)))
        exp.addSpectrum(spec)
    ms2 = oms.MSSpectrum()
    ms2.setMSLevel(2)
    exp.addSpectrum(ms2)
    return Py_MSExperiment(exp)


def _make_feature(rt: float, mz: float, intensity: float) -> oms.Feature:
    feature = oms.Feature()
    feature.setRT(rt)
    feature.setMZ(mz)
    feature.setIntensity(intensity)
    return feature


def test_detect_features_filters_intensity():
    exp = _make_experiment()
    fmap = exp.detect_features(min_intensity=60.0, mz_tolerance=0.5)
    assert isinstance(fmap, Py_FeatureMap)
    assert len(fmap) == 2
    feature = fmap[0]
    assert feature.metaValueExists("n_peaks")
    assert feature.getIntensity() > 0


def test_stream_theoretical_spectra(tmp_path):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">P1 Some protein\nPEPTIDEK\n", encoding="utf-8")
    workflow = FastaWorkflow.from_file(fasta)
    generator = workflow.theoretical_spectra(charges=(1, 2), min_length=4)
    first = next(generator)
    assert first["protein_id"] == "P1"
    assert first["peptide_sequence"].startswith("PEP")
    assert isinstance(first["spectrum"], oms.MSSpectrum)


def _make_identifications() -> Identifications:
    pep_good = oms.PeptideIdentification()
    pep_good.setHigherScoreBetter(False)
    pep_good.setScoreType("Posterior Error Probability")
    pep_good.setRT(10.0)
    pep_good.setMZ(500.0)
    hit = oms.PeptideHit()
    hit.setScore(0.001)
    hit.setMetaValue("target_decoy", "target")
    pep_good.setHits([hit])

    pep_bad = oms.PeptideIdentification()
    pep_bad.setHigherScoreBetter(False)
    pep_bad.setScoreType("Posterior Error Probability")
    pep_bad.setRT(11.0)
    pep_bad.setMZ(600.0)
    hit_bad = oms.PeptideHit()
    hit_bad.setScore(0.5)
    hit_bad.setMetaValue("target_decoy", "decoy")
    pep_bad.setHits([hit_bad])

    return Identifications(peptides=[pep_good, pep_bad])


def test_identifications_filter_by_fdr():
    ids = _make_identifications()
    filtered = ids.filter_by_fdr(0.05)
    assert len(filtered.peptide_identifications) == 1
    assert np.isclose(filtered.peptide_identifications[0].getRT(), 10.0)


def test_map_identifications_to_features_assigns_hits():
    ids = _make_identifications()
    fmap = oms.FeatureMap()
    fmap.push_back(_make_feature(10.0, 500.0, 1_000.0))
    wrapper = Py_FeatureMap(fmap)
    annotated = wrapper.annotate_with_identifications(ids, mz_tolerance=0.1, rt_tolerance=0.5, ppm=False)
    assert len(annotated[0].getPeptideIdentifications()) == 1


def test_align_feature_maps_normalises_rt():
    fmap_a = oms.FeatureMap()
    fmap_a.push_back(_make_feature(10.0, 500.0, 1_000.0))
    fmap_b = oms.FeatureMap()
    fmap_b.push_back(_make_feature(20.0, 500.0, 2_000.0))
    aligned = Py_FeatureMap.align([Py_FeatureMap(fmap_a), Py_FeatureMap(fmap_b)], reference_index=0)
    assert np.isclose(aligned[0][0].getRT(), 10.0)
    assert np.isclose(aligned[1][0].getRT(), 10.0)


def test_link_features_and_export_quant_table():
    fmap_a = oms.FeatureMap()
    fmap_a.push_back(_make_feature(10.0, 500.0, 1_000.0))
    fmap_b = oms.FeatureMap()
    fmap_b.push_back(_make_feature(10.2, 500.05, 2_000.0))
    consensus = Py_FeatureMap.link([Py_FeatureMap(fmap_a), Py_FeatureMap(fmap_b)], mz_tolerance=0.2, rt_tolerance=0.5)
    assert len(consensus) == 1
    table = consensus.to_quant_table()
    assert set(["map_0", "map_1"]).issubset(table.columns)
    assert np.isclose(table["map_0"].iloc[0], 1_000.0)
    assert np.isclose(table["map_1"].iloc[0], 2_000.0)
