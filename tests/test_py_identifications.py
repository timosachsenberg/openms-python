"""Tests for the identification wrappers."""
from __future__ import annotations

from pathlib import Path

import pyopenms as oms

from openms_python.py_identifications import (
    Identifications,
    PeptideIdentifications,
    ProteinIdentifications,
)


def _protein(identifier: str, accession: str) -> oms.ProteinIdentification:
    prot = oms.ProteinIdentification()
    prot.setIdentifier(identifier)
    hit = oms.ProteinHit()
    hit.setAccession(accession)
    prot.setHits([hit])
    return prot


def _peptide(
    identifier: str,
    sequence: str,
    score: float,
    *,
    higher_is_better: bool = True,
    accessions: tuple[str, ...] = (),
) -> oms.PeptideIdentification:
    pep = oms.PeptideIdentification()
    pep.setIdentifier(identifier)
    pep.setHigherScoreBetter(higher_is_better)
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(sequence))
    hit.setScore(score)
    if accessions:
        evidences = []
        for accession in accessions:
            evidence = oms.PeptideEvidence()
            evidence.setProteinAccession(accession)
            evidences.append(evidence)
        hit.setPeptideEvidences(evidences)
    pep.setHits([hit])
    return pep


def test_protein_identifications_behave_like_sequence() -> None:
    proteins = ProteinIdentifications([_protein("run1", "P1")])
    proteins.append(_protein("run2", "P2"))

    assert len(proteins) == 2
    assert proteins[0].getHits()[0].getAccession() == "P1"
    assert isinstance(proteins[:1], ProteinIdentifications)
    assert proteins.find_by_identifier("run2").getHits()[0].getAccession() == "P2"
    assert [p.getIdentifier() for p in proteins] == ["run1", "run2"]


def test_peptide_identifications_filter_and_find() -> None:
    pep_ids = PeptideIdentifications(
        [
            _peptide("run1", "PEPTIDE", 50.0, accessions=("P1",)),
            _peptide("run1", "PEPTIDER", 10.0, higher_is_better=False, accessions=("P2",)),
        ]
    )

    assert pep_ids.find_by_sequence("peptide").getHits()[0].getScore() == 50.0

    filtered = pep_ids.filter_by_score(5.0)
    assert len(filtered) == 1
    assert filtered[0].getHits()[0].getSequence().toString() == "PEPTIDE"
    assert [pep.getIdentifier() for pep in pep_ids] == ["run1", "run1"]


def test_identifications_roundtrip_and_helpers(tmp_path: Path) -> None:
    proteins = ProteinIdentifications([
        _protein("run", "P10"),
        _protein("run_alt", "P99"),
    ])
    peptides = PeptideIdentifications([
        _peptide("run", "PEPTIDE", 60.0, accessions=("P10",)),
        _peptide("run_alt", "OTHER", 5.0, higher_is_better=False, accessions=("P99",)),
    ])
    ids = Identifications(proteins, peptides)

    output = tmp_path / "test.idXML"
    ids.store(output)

    loaded = Identifications.from_idxml(output)

    assert loaded.summary() == {
        "proteins": 2,
        "peptides": 2,
        "protein_hits": 2,
        "peptide_hits": 2,
    }

    assert loaded.find_protein_by_accession("P10").getHits()[0].getAccession() == "P10"
    assert (
        loaded.find_peptide_by_sequence("peptide").getHits()[0].getSequence().toString() == "PEPTIDE"
    )

    matched = loaded.peptides_for_protein("P10")
    assert len(matched) == 1
    assert matched[0].getHits()[0].getSequence().toString() == "PEPTIDE"

    filtered_ids = loaded.filter_peptides_by_score(4.0)
    assert len(filtered_ids.peptide_identifications) == 1


def test_identifications_is_iterable() -> None:
    ids = Identifications(
        ProteinIdentifications([_protein("run", "P10")]),
        PeptideIdentifications([
            _peptide("run", "PEPTIDE", 60.0, accessions=("P10",)),
            _peptide("run", "OTHER", 5.0, higher_is_better=False, accessions=("P10",)),
        ]),
    )

    assert [pep.getHits()[0].getSequence().toString() for pep in ids] == [
        "PEPTIDE",
        "OTHER",
    ]
