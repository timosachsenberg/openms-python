"""High level utilities wrapping the building blocks in :mod:`openms_python`."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

import pandas as pd
import pyopenms as oms

from .py_msexperiment import Py_MSExperiment
from .py_featuremap import Py_FeatureMap
from .py_consensusmap import Py_ConsensusMap
from .py_identifications import Identifications


_MSExperimentLike = Union[Py_MSExperiment, oms.MSExperiment]
_FeatureMapLike = Union[Py_FeatureMap, oms.FeatureMap]
_ConsensusMapLike = Union[Py_ConsensusMap, oms.ConsensusMap]


@dataclass
class DigestedPeptide:
    """Container storing peptide digestion results."""

    protein: oms.FASTAEntry
    peptide: oms.AASequence


@dataclass
class TheoreticalSpectrumRecord:
    """Result entry for theoretical spectrum generation."""

    protein: oms.FASTAEntry
    peptide: oms.AASequence
    spectrum: oms.MSSpectrum


class _StreamWrapper:
    """Small helper around generator factories to allow chaining."""

    def __init__(self, iterator_factory: Callable[[], Iterator]):
        self._iterator_factory = iterator_factory

    def __iter__(self) -> Iterator:
        return self._iterator_factory()


class ProteinStream(_StreamWrapper):
    """Expose proteins from a FASTA source and allow chaining operations."""

    @classmethod
    def from_fasta(cls, source: Union[str, Path, Iterable[oms.FASTAEntry]]) -> "ProteinStream":
        entries = _collect_fasta_entries(source)
        return cls(lambda: (entry for entry in entries))

    def digest(
        self,
        *,
        enzyme: str = "Trypsin",
        missed_cleavages: int = 0,
        min_length: int = 6,
        max_length: int = 60,
    ) -> "PeptideStream":
        """Digest every protein into peptides lazily."""

        def iterator() -> Iterator[DigestedPeptide]:
            digestion = oms.ProteaseDigestion()
            digestion.setEnzyme(enzyme)
            digestion.setMissedCleavages(int(missed_cleavages))
            for protein in self:
                peptides: List[oms.AASequence] = []
                sequence = oms.AASequence.fromString(protein.sequence)
                digestion.digest(sequence, peptides, int(min_length), int(max_length))
                for peptide in peptides:
                    yield DigestedPeptide(protein=protein, peptide=peptide)

        return PeptideStream(iterator)


class PeptideStream(_StreamWrapper):
    """Yield :class:`DigestedPeptide` entries and build theoretical spectra."""

    def theoretical_spectra(
        self,
        *,
        min_charge: int = 1,
        max_charge: int = 2,
        generator_params: Optional[Dict[str, Union[int, float, str]]] = None,
    ) -> "SpectrumStream":
        """Generate theoretical spectra on-the-fly for each peptide."""

        def iterator() -> Iterator[TheoreticalSpectrumRecord]:
            generator = oms.TheoreticalSpectrumGenerator()
            params = generator.getDefaults()
            if generator_params:
                for key, value in generator_params.items():
                    params.setValue(key, value)
            generator.setParameters(params)
            for entry in self:
                spectrum = oms.MSSpectrum()
                generator.getSpectrum(spectrum, entry.peptide, int(min_charge), int(max_charge))
                yield TheoreticalSpectrumRecord(entry.protein, entry.peptide, spectrum)

        return SpectrumStream(iterator)


class SpectrumStream(_StreamWrapper):
    """Terminal stage of the FASTA -> peptide -> spectrum pipeline."""

    pass


def stream_theoretical_spectra_from_fasta(
    fasta: Union[str, Path, Iterable[oms.FASTAEntry]],
    *,
    digestion: Optional[Dict[str, Union[int, str]]] = None,
    generator: Optional[Dict[str, Union[int, float, str]]] = None,
) -> SpectrumStream:
    """Convenience helper returning the streaming spectra pipeline."""

    digestion = digestion or {}
    generator = generator or {}
    return ProteinStream.from_fasta(fasta).digest(**digestion).theoretical_spectra(**generator)


def map_identifications_to_features(
    feature_map: _FeatureMapLike,
    identifications: Identifications,
    *,
    rt_tolerance: float = 5.0,
    mz_tolerance: float = 20.0,
    mz_unit: str = "ppm",
    use_centroid_rt: bool = False,
    use_centroid_mz: bool = False,
) -> Py_FeatureMap:
    """Annotate a feature map with identifications using :class:`pyopenms.IDMapper`."""

    mapper = oms.IDMapper()
    params = mapper.getDefaults()
    params.setValue("rt_tolerance", float(rt_tolerance))
    params.setValue("mz_tolerance", float(mz_tolerance))
    params.setValue("mz_measure", mz_unit)
    mapper.setParameters(params)

    native_feature_map = _as_feature_map(feature_map)
    mapper.annotate(
        native_feature_map,
        identifications.peptides.native,
        identifications.proteins.native,
        use_centroid_rt,
        use_centroid_mz,
        oms.MSExperiment(),
    )
    return Py_FeatureMap(native_feature_map)


def align_feature_maps(
    feature_maps: Sequence[_FeatureMapLike],
    *,
    reference_index: int = 0,
) -> List[Py_FeatureMap]:
    """Align retention times by matching the median RT of each map to the reference."""

    if not feature_maps:
        return []

    if reference_index < 0 or reference_index >= len(feature_maps):
        raise IndexError("reference_index is out of range")

    native_maps = [_as_feature_map(feature_map) for feature_map in feature_maps]
    reference_rts = _collect_rts(native_maps[reference_index])
    reference_median = median(reference_rts) if reference_rts else 0.0

    aligned: List[Py_FeatureMap] = []
    for idx, fmap in enumerate(native_maps):
        shift = 0.0
        if idx != reference_index:
            rts = _collect_rts(fmap)
            if rts:
                shift = median(rts) - reference_median
        adjusted = oms.FeatureMap()
        for feature_idx in range(fmap.size()):
            feature = oms.Feature(fmap[feature_idx])
            if shift:
                feature.setRT(feature.getRT() - shift)
            adjusted.push_back(feature)
        aligned.append(Py_FeatureMap(adjusted))
    return aligned


def link_features(
    feature_maps: Sequence[_FeatureMapLike],
    *,
    params: Optional[Dict[str, Union[int, float, str]]] = None,
) -> Py_ConsensusMap:
    """Group features across runs into a consensus map."""

    grouping = oms.FeatureGroupingAlgorithmQT()
    param_obj = grouping.getDefaults()
    if params:
        for key, value in params.items():
            param_obj.setValue(key, value)
    grouping.setParameters(param_obj)

    native_maps = [_as_feature_map(feature_map) for feature_map in feature_maps]
    consensus_map = oms.ConsensusMap()
    grouping.group(native_maps, consensus_map)
    if not consensus_map.getColumnHeaders():
        headers: Dict[int, oms.ColumnHeader] = {}
        for index in range(len(native_maps)):
            header = oms.ColumnHeader()
            header.label = f"map_{index}"
            headers[index] = header
        consensus_map.setColumnHeaders(headers)
    return Py_ConsensusMap(consensus_map)


def export_quant_table(
    consensus_map: _ConsensusMapLike,
    *,
    include_coordinates: bool = True,
) -> pd.DataFrame:
    """Convert a consensus map to a tidy quantification table."""

    cmap = consensus_map.native if isinstance(consensus_map, Py_ConsensusMap) else consensus_map
    headers = cmap.getColumnHeaders()
    label_lookup = {
        index: (header.label or header.filename or f"map_{index}")
        for index, header in headers.items()
    }
    sample_columns = [label_lookup[idx] for idx in sorted(label_lookup)]

    rows: List[Dict[str, Union[float, int]]] = []
    for feature in cmap:
        row: Dict[str, Union[float, int]] = {label: float("nan") for label in sample_columns}
        if include_coordinates:
            row.update(
                {
                    "rt": feature.getRT(),
                    "mz": feature.getMZ(),
                    "intensity": feature.getIntensity(),
                    "charge": feature.getCharge(),
                }
            )
        for handle in feature.getFeatureList():
            label = label_lookup.get(handle.getMapIndex(), f"map_{handle.getMapIndex()}")
            row[label] = handle.getIntensity()
        rows.append(row)

    columns: List[str] = sample_columns
    if include_coordinates:
        columns = ["rt", "mz", "intensity", "charge"] + columns
    return pd.DataFrame(rows, columns=columns if rows else columns)


def _collect_fasta_entries(source: Union[str, Path, Iterable[oms.FASTAEntry]]) -> Tuple[oms.FASTAEntry, ...]:
    if isinstance(source, (str, Path)):
        entries: List[oms.FASTAEntry] = []
        oms.FASTAFile().load(str(source), entries)
        return tuple(entries)
    return tuple(oms.FASTAEntry(entry) for entry in source)


def _as_feature_map(feature_map: Optional[_FeatureMapLike]) -> oms.FeatureMap:
    if feature_map is None:
        return oms.FeatureMap()
    if isinstance(feature_map, Py_FeatureMap):
        return oms.FeatureMap(feature_map.native)
    return oms.FeatureMap(feature_map)


def _collect_rts(feature_map: oms.FeatureMap) -> List[float]:
    rts: List[float] = []
    for idx in range(feature_map.size()):
        rts.append(feature_map[idx].getRT())
    return rts
