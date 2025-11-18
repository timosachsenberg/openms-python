from __future__ import annotations

import math
from pathlib import Path
from time import perf_counter
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd
import pyopenms as oms
from pyopenms import Constants

from openms_python import Py_MSExperiment, ProteinStream


GENERATOR_PARAMS = {
    "add_y_ions": "true",
    "add_b_ions": "true",
    "add_first_prefix_ion": "true",
    "add_precursor_peaks": "false",
    "add_losses": "false",
    "add_metainfo": "true",
}


class SimpleSearchEngine:
    def __init__(
        self,
        targets: Sequence[oms.AASequence],
        decoys: Sequence[oms.AASequence],
        *,
        mass_tolerance: float = 0.01,
    ) -> None:
        self.targets = list(targets)
        self.decoys = list(decoys)
        self.mass_tolerance = float(mass_tolerance)
        self.target_masses = np.array([pep.getMonoWeight() for pep in self.targets])
        self.decoy_masses = np.array([pep.getMonoWeight() for pep in self.decoys])
        self.generator = _make_theoretical_generator()
        self.aligner = _make_aligner()

    def search(self, experiment: Py_MSExperiment) -> pd.DataFrame:
        results: List[dict] = []
        for spectrum in experiment.ms2_spectra():
            precursor_mass = spectrum.precursor_mass
            target_candidates = self._select_candidates(
                precursor_mass, self.target_masses, self.targets
            )
            decoy_candidates = self._select_candidates(
                precursor_mass, self.decoy_masses, self.decoys
            )

            target_score, target_seq = self._score_candidates(spectrum, target_candidates)
            decoy_score, decoy_seq = self._score_candidates(spectrum, decoy_candidates)

            if target_score >= decoy_score:
                score = target_score
                label = 1
                sequence = target_seq
            else:
                score = decoy_score
                label = 0
                sequence = decoy_seq

            results.append(
                {
                    "Score": float(score),
                    "Label": label,
                    "Sequence": sequence.toString() if sequence is not None else None,
                }
            )
        return pd.DataFrame(results)

    def _score_candidates(
        self, spectrum: oms.MSSpectrum | "Py_MSSpectrum", candidates: Sequence[oms.AASequence]
    ) -> tuple[float, oms.AASequence | None]:
        best_score = 0.0
        best_sequence: oms.AASequence | None = None
        for candidate in candidates:
            score = self._hyperscore(spectrum, candidate)
            if score > best_score:
                best_score = score
                best_sequence = candidate
        return best_score, best_sequence

    def _hyperscore(
        self, spectrum: oms.MSSpectrum | "Py_MSSpectrum", peptide: oms.AASequence
    ) -> float:
        if isinstance(spectrum, oms.MSSpectrum):
            observed = spectrum
            precursor_charge = _infer_charge_from_native(observed)
        else:
            observed = spectrum.native  # type: ignore[attr-defined]
            precursor_charge = spectrum.precursor_charge or _infer_charge_from_native(observed)  # type: ignore[attr-defined]

        fragment_charge = max(1, min(int(precursor_charge) - 1, 2)) if precursor_charge else 1
        theoretical = oms.MSSpectrum()
        self.generator.getSpectrum(theoretical, peptide, 1, fragment_charge)
        theoretical.sortByPosition()

        alignment: List[tuple[int, int]] = []
        self.aligner.getSpectrumAlignment(alignment, theoretical, observed)
        if not alignment:
            return 0.0

        exp_indices = np.array([pair[1] for pair in alignment], dtype=int)
        theo_indices = np.array([pair[0] for pair in alignment], dtype=int)
        intensities = np.take(_get_intensities(spectrum), exp_indices)
        y_count, b_count = _count_fragment_types(theoretical, theo_indices)
        return math.log1p(float(np.sum(intensities))) + math.lgamma(y_count + 1.0) + math.lgamma(
            b_count + 1.0
        )

    def _select_candidates(
        self,
        precursor_mass: float | None,
        mass_array: np.ndarray,
        candidates: Sequence[oms.AASequence],
    ) -> List[oms.AASequence]:
        if precursor_mass is None:
            return []
        matches = np.where(np.isclose(mass_array, precursor_mass, atol=self.mass_tolerance))[0]
        return [candidates[int(idx)] for idx in matches]


def _make_theoretical_generator() -> oms.TheoreticalSpectrumGenerator:
    generator = oms.TheoreticalSpectrumGenerator()
    params = generator.getDefaults()
    for key, value in GENERATOR_PARAMS.items():
        params.setValue(key, value)
    generator.setParameters(params)
    return generator


def _make_aligner() -> oms.SpectrumAlignment:
    aligner = oms.SpectrumAlignment()
    params = aligner.getDefaults()
    params.setValue("tolerance", 0.05)
    params.setValue("is_relative_tolerance", "false")
    aligner.setParameters(params)
    return aligner


def _get_intensities(spectrum: oms.MSSpectrum | "Py_MSSpectrum") -> np.ndarray:
    if isinstance(spectrum, oms.MSSpectrum):
        _, intensity = spectrum.get_peaks()
        return np.asarray(intensity, dtype=float)
    return np.asarray(spectrum.intensity, dtype=float)  # type: ignore[attr-defined]


def _infer_charge_from_native(spectrum: oms.MSSpectrum) -> int:
    precursors = spectrum.getPrecursors()
    if precursors:
        return max(int(precursors[0].getCharge()), 1)
    return 1


def _count_fragment_types(theoretical: oms.MSSpectrum, indices: Iterable[int]) -> tuple[int, int]:
    arrays = theoretical.getStringDataArrays()
    if not arrays:
        return (0, 0)
    annotations = arrays[0]
    y_count = 0
    b_count = 0
    for idx in indices:
        label = annotations[int(idx)]
        if isinstance(label, bytes):
            label = label.decode()
        else:
            label = str(label)
        if label.startswith("y"):
            y_count += 1
        elif label.startswith("b"):
            b_count += 1
    return y_count, b_count


def _reverse_sequence(peptide: oms.AASequence) -> oms.AASequence:
    return oms.AASequence.fromString(str(peptide)[::-1])


def _build_fasta(tmp_path: Path) -> Path:
    fasta_path = tmp_path / "targets.fasta"
    fasta_path.write_text(
        ">pep1\nPEPTIDERK\n"
        ">pep2\nGILQTYK\n"
    )
    return fasta_path


def _experimental_from_theoretical(
    theoretical: oms.MSSpectrum, peptide: oms.AASequence, rt: float
) -> oms.MSSpectrum:
    mz, intensity = theoretical.get_peaks()
    mz_array = np.asarray(mz, dtype=float)
    intensity_array = np.asarray(intensity, dtype=float)
    if intensity_array.size == 0:
        intensity_array = np.array([100.0])
        mz_array = np.array([peptide.getMonoWeight() / 2.0])
    scaled = (intensity_array / intensity_array.max()) * 100.0
    noise_mz = np.array([mz_array.min() - 5.0, mz_array.max() + 5.0])
    noise_intensity = np.array([5.0, 3.0])
    combined_mz = np.concatenate([mz_array, noise_mz])
    combined_intensity = np.concatenate([scaled, noise_intensity])

    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(2)
    spectrum.setRT(rt)
    spectrum.set_peaks((combined_mz.tolist(), combined_intensity.tolist()))
    spectrum.sortByPosition()

    precursor = oms.Precursor()
    precursor.setCharge(2)
    precursor.setMZ((peptide.getMonoWeight() + 2 * Constants.PROTON_MASS_U) / 2.0)
    spectrum.setPrecursors([precursor])
    return spectrum


def _sparsify_spectrum(spectrum: oms.MSSpectrum, keep_every: int = 2) -> oms.MSSpectrum:
    if keep_every <= 1:
        return spectrum
    mz, intensity = spectrum.get_peaks()
    mz_array = np.asarray(mz, dtype=float)
    intensity_array = np.asarray(intensity, dtype=float)
    if mz_array.size <= keep_every:
        return spectrum
    mask = np.zeros_like(mz_array, dtype=bool)
    mask[::keep_every] = True
    filtered_mz = mz_array[mask]
    filtered_intensity = intensity_array[mask]
    spectrum.set_peaks((filtered_mz.tolist(), filtered_intensity.tolist()))
    spectrum.sortByPosition()
    return spectrum


def _build_experiment(
    records: Sequence, decoy_sequence: oms.AASequence
) -> Py_MSExperiment:
    experiment = Py_MSExperiment()
    for idx, record in enumerate(records):
        experimental = _experimental_from_theoretical(record.spectrum, record.peptide, 50.0 + idx)
        experiment.append(experimental)
    generator = _make_theoretical_generator()
    decoy_spectrum = oms.MSSpectrum()
    generator.getSpectrum(decoy_spectrum, decoy_sequence, 1, 2)
    decoy_spectrum.sortByPosition()
    decoy_exp = _experimental_from_theoretical(decoy_spectrum, decoy_sequence, 200.0)
    experiment.append(_sparsify_spectrum(decoy_exp, keep_every=3))
    return experiment


def _compute_q_values(psms: pd.DataFrame) -> pd.DataFrame:
    ordered = psms.sort_values("Score", ascending=False).reset_index(drop=True)
    is_target = (ordered["Label"] == 1).astype(int)
    is_decoy = (ordered["Label"] == 0).astype(int)
    ordered["targets"] = is_target.cumsum()
    ordered["decoys"] = is_decoy.cumsum()
    ordered["fdr"] = ordered["decoys"] / ordered["targets"].replace(0, np.nan)
    ordered["fdr"] = ordered["fdr"].fillna(0.0)
    ordered["q_value"] = ordered["fdr"][::-1].cummin()[::-1]
    return ordered.drop(columns=["targets", "decoys", "fdr"])


def _time_search(engine: SimpleSearchEngine, experiment: Py_MSExperiment) -> tuple[pd.DataFrame, float]:
    start = perf_counter()
    psms = engine.search(experiment)
    elapsed = perf_counter() - start
    return psms, elapsed


def test_simple_search_engine_identifies_targets(tmp_path):
    fasta_path = _build_fasta(tmp_path)
    targets = [
        entry.peptide
        for entry in ProteinStream.from_fasta(fasta_path).digest(min_length=6, max_length=10)
    ]
    assert len(targets) == 2
    decoys = [_reverse_sequence(peptide) for peptide in targets]

    processed = _build_experiment(
        list(
            ProteinStream.from_fasta(fasta_path)
            .digest(min_length=6, max_length=10)
            .theoretical_spectra(generator_params=GENERATOR_PARAMS)
        ),
        decoys[0],
    ).filter_top_n_peaks(200).normalize_to_tic()

    assert all(math.isclose(spec.total_ion_current, 1.0, rel_tol=1e-6) for spec in processed.ms2_spectra())
    assert all(spec.precursor_mass is not None for spec in processed.ms2_spectra())

    engine = SimpleSearchEngine(targets, decoys)
    psms = engine.search(processed)
    assert len(psms) == 3

    confident = _compute_q_values(psms).query("q_value <= 0.01")
    assert len(confident) == 2
    assert set(confident["Sequence"]) == {seq.toString() for seq in targets}

    # measure runtime when operating on an mzML file
    mzml_path = tmp_path / "experiment.mzML"
    processed.to_mzml(str(mzml_path))
    reloaded = Py_MSExperiment.from_mzml(str(mzml_path)).filter_top_n_peaks(200).normalize_to_tic()
    file_psms, elapsed = _time_search(engine, reloaded)
    pd.testing.assert_frame_equal(psms.reset_index(drop=True), file_psms.reset_index(drop=True))
    assert elapsed >= 0.0
