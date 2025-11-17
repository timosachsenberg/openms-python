"""Wrappers for Protein- and PeptideIdentifications collections."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple, Union

import pyopenms as oms

from ._io_utils import ensure_allowed_suffix, IDENTIFICATION_EXTENSIONS


class ProteinIdentifications:
    """Sequence-like container for :class:`pyopenms.ProteinIdentification`."""

    def __init__(self, proteins: Optional[Iterable[oms.ProteinIdentification]] = None) -> None:
        self._proteins: List[oms.ProteinIdentification] = []
        if proteins:
            for protein in proteins:
                self.append(protein)

    @property
    def native(self) -> List[oms.ProteinIdentification]:
        """Return the mutable list of native protein identifications."""

        return self._proteins

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._proteins)

    def __iter__(self) -> Iterator[oms.ProteinIdentification]:  # pragma: no cover - trivial
        return iter(self._proteins)

    def __getitem__(self, key: Union[int, slice]) -> Union[oms.ProteinIdentification, "ProteinIdentifications"]:
        if isinstance(key, slice):
            return ProteinIdentifications(self._proteins[key])
        return self._proteins[key]

    def append(self, protein: oms.ProteinIdentification) -> "ProteinIdentifications":
        """Append a :class:`pyopenms.ProteinIdentification`."""

        self._proteins.append(oms.ProteinIdentification(protein))
        return self

    def extend(self, proteins: Iterable[oms.ProteinIdentification]) -> "ProteinIdentifications":
        for protein in proteins:
            self.append(protein)
        return self

    def find_by_identifier(self, identifier: str) -> Optional[oms.ProteinIdentification]:
        """Return the first protein identification with *identifier*."""

        for protein in self._proteins:
            if protein.getIdentifier() == identifier:
                return protein
        return None

    def find_by_accession(self, accession: str) -> Optional[oms.ProteinIdentification]:
        """Return the first protein identification containing *accession*."""

        for protein in self._proteins:
            for hit in protein.getHits():
                if hit.getAccession() == accession:
                    return protein
        return None

    def summary(self) -> List[tuple[str, int]]:
        """Return (identifier, #hits) tuples summarising the entries."""

        return [(protein.getIdentifier(), len(protein.getHits())) for protein in self._proteins]

    def copy(self) -> "ProteinIdentifications":
        """Return a shallow copy of this container."""

        return ProteinIdentifications(self._proteins)


class PeptideIdentifications:
    """Sequence-like container for :class:`pyopenms.PeptideIdentification`."""

    def __init__(self, peptides: Optional[Iterable[oms.PeptideIdentification]] = None) -> None:
        self._peptides: List[oms.PeptideIdentification] = []
        if peptides:
            for peptide in peptides:
                self.append(peptide)

    @property
    def native(self) -> List[oms.PeptideIdentification]:
        return self._peptides

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._peptides)

    def __iter__(self) -> Iterator[oms.PeptideIdentification]:  # pragma: no cover - trivial
        return iter(self._peptides)

    def __getitem__(self, key: Union[int, slice]) -> Union[oms.PeptideIdentification, "PeptideIdentifications"]:
        if isinstance(key, slice):
            return PeptideIdentifications(self._peptides[key])
        return self._peptides[key]

    def append(self, peptide: oms.PeptideIdentification) -> "PeptideIdentifications":
        self._peptides.append(oms.PeptideIdentification(peptide))
        return self

    def extend(self, peptides: Iterable[oms.PeptideIdentification]) -> "PeptideIdentifications":
        for peptide in peptides:
            self.append(peptide)
        return self

    def find_by_identifier(self, identifier: str) -> Optional[oms.PeptideIdentification]:
        for peptide in self._peptides:
            if peptide.getIdentifier() == identifier:
                return peptide
        return None

    def find_by_sequence(self, sequence: str) -> Optional[oms.PeptideIdentification]:
        target = sequence.upper()
        for peptide in self._peptides:
            for hit in peptide.getHits():
                if hit.getSequence().toString().upper() == target:
                    return peptide
        return None

    def filter_by_score(self, threshold: float) -> "PeptideIdentifications":
        """Return peptide IDs whose top hit passes *threshold*."""

        kept: List[oms.PeptideIdentification] = []
        for peptide in self._peptides:
            hits = peptide.getHits()
            if not hits:
                continue
            best = hits[0]
            if peptide.isHigherScoreBetter():
                condition = best.getScore() >= threshold
            else:
                condition = best.getScore() <= threshold
            if condition:
                kept.append(peptide)
        return PeptideIdentifications(kept)

    def copy(self) -> "PeptideIdentifications":
        return PeptideIdentifications(self._peptides)


class Identifications:
    """Convenience wrapper combining protein and peptide identifications."""

    def __init__(
        self,
        proteins: Optional[Union[ProteinIdentifications, Iterable[oms.ProteinIdentification]]] = None,
        peptides: Optional[Union[PeptideIdentifications, Iterable[oms.PeptideIdentification]]] = None,
    ) -> None:
        self.proteins = proteins if isinstance(proteins, ProteinIdentifications) else ProteinIdentifications(proteins)
        self.peptides = peptides if isinstance(peptides, PeptideIdentifications) else PeptideIdentifications(peptides)

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self.peptides)

    def __iter__(self) -> Iterator[oms.PeptideIdentification]:  # pragma: no cover - trivial
        return iter(self.peptides)

    @property
    def protein_identifications(self) -> ProteinIdentifications:
        return self.proteins

    @property
    def peptide_identifications(self) -> PeptideIdentifications:
        return self.peptides

    @classmethod
    def from_idxml(cls, filepath: Union[str, Path]) -> "Identifications":
        """Load protein and peptide identifications from an idXML file."""

        ensure_allowed_suffix(filepath, IDENTIFICATION_EXTENSIONS, "Identifications")
        proteins: List[oms.ProteinIdentification] = []
        peptides: List[oms.PeptideIdentification] = []
        oms.IdXMLFile().load(str(filepath), proteins, peptides)
        return cls(proteins, peptides)

    def store(self, filepath: Union[str, Path]) -> "Identifications":
        """Persist the identifications to an idXML file."""

        ensure_allowed_suffix(filepath, IDENTIFICATION_EXTENSIONS, "Identifications")
        oms.IdXMLFile().store(str(filepath), self.proteins.native, self.peptides.native)
        return self

    def to_idxml(self, filepath: Union[str, Path]) -> "Identifications":  # pragma: no cover - trivial
        return self.store(filepath)

    def find_protein(self, identifier: str) -> Optional[oms.ProteinIdentification]:
        return self.proteins.find_by_identifier(identifier)

    def find_protein_by_accession(self, accession: str) -> Optional[oms.ProteinIdentification]:
        return self.proteins.find_by_accession(accession)

    def find_peptide(self, identifier: str) -> Optional[oms.PeptideIdentification]:
        return self.peptides.find_by_identifier(identifier)

    def find_peptide_by_sequence(self, sequence: str) -> Optional[oms.PeptideIdentification]:
        return self.peptides.find_by_sequence(sequence)

    def peptides_for_protein(self, accession: str) -> PeptideIdentifications:
        """Return peptides that reference the provided protein accession."""

        matches: List[oms.PeptideIdentification] = []
        for peptide in self.peptides:
            for hit in peptide.getHits():
                for evidence in hit.getPeptideEvidences():
                    if evidence.getProteinAccession() == accession:
                        matches.append(peptide)
                        break
                else:
                    continue
                break
        return PeptideIdentifications(matches)

    def filter_peptides_by_score(self, threshold: float) -> "Identifications":
        """Return a copy with peptides filtered by score threshold."""

        filtered = self.peptides.filter_by_score(threshold)
        return Identifications(self.proteins.copy(), filtered)

    def filter_by_fdr(
        self,
        *,
        threshold: float = 0.01,
        level: str = "peptide",
        qvalue_keys: Sequence[str] = ("q-value", "FDR", "pep:score"),
    ) -> "Identifications":
        """Return a copy filtered by False Discovery Rate at the given level."""

        level = level.lower()
        if level not in {"peptide", "protein"}:
            raise ValueError("level must be either 'peptide' or 'protein'")

        if level == "peptide":
            entries = [oms.PeptideIdentification(entry) for entry in self.peptides.native]
            _run_fdr(entries)
            _ensure_q_values(entries, is_peptide=True)
            kept = [entry for entry in entries if _extract_q_value(entry, qvalue_keys, True) <= threshold]
            return Identifications(self.proteins.copy(), kept)

        proteins = [oms.ProteinIdentification(entry) for entry in self.proteins.native]
        _run_fdr(proteins)
        _ensure_q_values(proteins, is_peptide=False)
        kept_proteins = [entry for entry in proteins if _extract_q_value(entry, qvalue_keys, False) <= threshold]
        return Identifications(kept_proteins, self.peptides.copy())

    def summary(self) -> dict:
        """Return basic counts about the contained identifications."""

        protein_hits = sum(len(protein.getHits()) for protein in self.proteins)
        peptide_hits = sum(len(peptide.getHits()) for peptide in self.peptides)
        return {
            "proteins": len(self.proteins),
            "peptides": len(self.peptides),
            "protein_hits": protein_hits,
            "peptide_hits": peptide_hits,
        }


def _run_fdr(container: List[Union[oms.PeptideIdentification, oms.ProteinIdentification]]) -> None:
    if not container:
        return
    fdr = oms.FalseDiscoveryRate()
    try:
        fdr.apply(container)
    except RuntimeError:
        # Fallback to manual estimation below
        pass


def _ensure_q_values(container: List, is_peptide: bool) -> None:
    needs_assignment = not any(_has_q_value(entry, is_peptide) for entry in container)
    if not needs_assignment:
        return

    entries: List[Tuple[Union[oms.PeptideIdentification, oms.ProteinIdentification], float, bool, bool]] = []
    for entry in container:
        hits = entry.getHits()
        if not hits:
            continue
        label = _target_decoy_label(entry, hits[0], is_peptide)
        if label is None:
            continue
        entries.append((entry, hits[0].getScore(), label == "decoy", entry.isHigherScoreBetter()))

    if not entries:
        return

    higher_better = entries[0][3]
    sorted_entries = sorted(entries, key=lambda item: item[1], reverse=higher_better)
    ratios: List[float] = []
    target = 0
    decoy = 0
    for entry, _, is_decoy, _ in sorted_entries:
        if is_decoy:
            decoy += 1
        else:
            target += 1
        ratios.append(decoy / max(1, target))

    best = 1.0
    for idx in reversed(range(len(sorted_entries))):
        best = min(best, ratios[idx])
        entry = sorted_entries[idx][0]
        _set_q_value(entry, best, is_peptide)


def _has_q_value(entry, is_peptide: bool) -> bool:
    if entry.metaValueExists("q-value"):
        return True
    hits = entry.getHits()
    if not hits:
        return False
    return hits[0].metaValueExists("q-value")


def _extract_q_value(entry, keys: Sequence[str], is_peptide: bool) -> float:
    if entry.metaValueExists("q-value"):
        return float(entry.getMetaValue("q-value"))
    hits = entry.getHits()
    if hits:
        hit = hits[0]
        if hit.metaValueExists("q-value"):
            return float(hit.getMetaValue("q-value"))
        for key in keys:
            if hit.metaValueExists(key):
                return float(hit.getMetaValue(key))
    for key in keys:
        if entry.metaValueExists(key):
            return float(entry.getMetaValue(key))
    return float(hits[0].getScore()) if hits else 1.0


def _set_q_value(entry, value: float, is_peptide: bool) -> None:
    entry.setMetaValue("q-value", float(value))
    for hit in entry.getHits():
        hit.setMetaValue("q-value", float(value))


def _target_decoy_label(entry, hit, is_peptide: bool) -> Optional[str]:
    if hit.metaValueExists("target_decoy"):
        return str(hit.getMetaValue("target_decoy")).lower()
    if entry.metaValueExists("target_decoy"):
        return str(entry.getMetaValue("target_decoy")).lower()
    return None
