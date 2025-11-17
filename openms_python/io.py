"""I/O helpers for openms_python."""

from __future__ import annotations

from typing import Iterator, Optional

import pyopenms as oms

from .py_msexperiment import Py_MSExperiment
from .py_msspectrum import Py_MSSpectrum


def read_mzml(filepath: str) -> Py_MSExperiment:
    """Read an mzML file into a :class:`Py_MSExperiment`."""
    return Py_MSExperiment.from_mzml(filepath)


def write_mzml(experiment: Py_MSExperiment, filepath: str) -> None:
    """Write a :class:`Py_MSExperiment` to an mzML file."""
    experiment.to_mzml(filepath)


def stream_mzml(filepath: str, ms_level: Optional[int] = None) -> Iterator[Py_MSSpectrum]:
    """Yield spectra from an mzML file without loading it entirely into memory.

    This is a light wrapper around :func:`pyopenms.stream_mzML` that yields
    :class:`Py_MSSpectrum` objects so downstream code can keep using the
    Pythonic helper API while working in a streaming fashion.

    Args:
        filepath: Path to the mzML file that should be streamed.
        ms_level: If provided, only spectra whose MS level matches this value
            are yielded.

    Yields:
        ``Py_MSSpectrum`` objects, one spectrum at a time.
    """

    def _generator() -> Iterator[Py_MSSpectrum]:
        with oms.stream_mzML(filepath) as spectra:
            for spectrum in spectra:
                if ms_level is not None and spectrum.getMSLevel() != ms_level:
                    continue
                yield Py_MSSpectrum(spectrum)

    return _generator()
