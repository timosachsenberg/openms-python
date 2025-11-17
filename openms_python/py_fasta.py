"""Chaining-friendly helpers for FASTA iteration and digestion."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Union

import pyopenms as oms


class FastaWorkflow:
    """Provide fluent helpers for protein iteration and digestion."""

    def __init__(self, entries: Iterable[oms.FASTAEntry]):
        self._entries: List[oms.FASTAEntry] = [oms.FASTAEntry(entry) for entry in entries]

    @classmethod
    def from_file(cls, fasta: Union[str, Path]) -> "FastaWorkflow":
        entries: List[oms.FASTAEntry] = []
        oms.FASTAFile().load(str(fasta), entries)
        return cls(entries)

    @classmethod
    def from_entries(cls, entries: Iterable[oms.FASTAEntry]) -> "FastaWorkflow":
        return cls(entries)

    def proteins(self) -> Iterator[oms.FASTAEntry]:
        """Yield copies of all stored FASTA entries."""

        for entry in self._entries:
            yield oms.FASTAEntry(entry)

    def digest(
        self,
        *,
        enzyme: str = "Trypsin",
        missed_cleavages: int = 0,
        min_length: int = 6,
        max_length: int = 60,
    ) -> Iterator[Dict[str, Union[str, oms.AASequence]]]:
        """Yield dictionaries describing peptides produced by in-silico digestion."""

        digestion = oms.ProteaseDigestion()
        digestion.setEnzyme(enzyme)
        digestion.setMissedCleavages(missed_cleavages)

        for entry in self.proteins():
            aa_sequence = oms.AASequence.fromString(entry.sequence)
            peptides: List[oms.AASequence] = []
            digestion.digest(aa_sequence, peptides, min_length, max_length)
            for peptide in peptides:
                yield {
                    "protein_id": entry.identifier,
                    "protein_description": entry.description,
                    "peptide_sequence": peptide,
                }

    def theoretical_spectra(
        self,
        *,
        enzyme: str = "Trypsin",
        missed_cleavages: int = 0,
        min_length: int = 6,
        max_length: int = 60,
        charges: Sequence[int] = (1,),
    ) -> Iterator[Dict[str, Union[str, int, oms.MSSpectrum]]]:
        """Yield theoretical spectra for every digested peptide."""

        generator = oms.TheoreticalSpectrumGenerator()

        for peptide_record in self.digest(
            enzyme=enzyme,
            missed_cleavages=missed_cleavages,
            min_length=min_length,
            max_length=max_length,
        ):
            peptide: oms.AASequence = peptide_record["peptide_sequence"]  # type: ignore[assignment]
            for charge in charges:
                spectrum = oms.MSSpectrum()
                generator.getSpectrum(spectrum, peptide, charge, charge)
                yield {
                    "protein_id": peptide_record["protein_id"],
                    "protein_description": peptide_record["protein_description"],
                    "peptide_sequence": peptide.toString(),
                    "charge": int(charge),
                    "spectrum": spectrum,
                }
