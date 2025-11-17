"""Utilities for reading, writing, and streaming mzML files with pyOpenMS."""
from __future__ import annotations

from contextlib import AbstractContextManager
from pathlib import Path
from queue import SimpleQueue
from threading import Thread
from typing import Iterator, Union

try:
    import pyopenms as oms
except ImportError:  # pragma: no cover - handled gracefully at runtime
    oms = None  # type: ignore

from .py_msspectrum import Py_MSSpectrum


class _PyOpenMSNotAvailable(RuntimeError):
    """Raised when pyOpenMS is required but not installed."""


def _ensure_pyopenms_available() -> None:
    """Ensure pyOpenMS is importable before performing IO."""

    if oms is None:
        raise _PyOpenMSNotAvailable(
            "pyOpenMS is required for reading/writing mzML files. "
            "Install the 'pyopenms' package to use these helpers."
        )


def read_mzml(path: Union[str, Path], *, as_wrapper: bool = True):
    """Load an mzML file and optionally wrap it in :class:`Py_MSExperiment`."""

    from .py_msexperiment import Py_MSExperiment  # Local import to avoid cycles

    _ensure_pyopenms_available()
    mzml = oms.MzMLFile()
    experiment = oms.MSExperiment()
    mzml.load(str(path), experiment)
    return Py_MSExperiment(experiment) if as_wrapper else experiment


def write_mzml(experiment, path: Union[str, Path]) -> None:
    """Persist an MSExperiment or Py_MSExperiment to disk."""

    _ensure_pyopenms_available()
    native = getattr(experiment, "_experiment", experiment)
    if not isinstance(native, oms.MSExperiment):
        raise TypeError(
            "write_mzml expects a pyopenms.MSExperiment or Py_MSExperiment instance"
        )

    mzml = oms.MzMLFile()
    mzml.store(str(path), native)


class _QueueingSpectrumConsumer:
    """Minimal consumer that forwards streamed spectra to a queue."""

    def __init__(self, queue: SimpleQueue, as_wrapper: bool) -> None:
        self._queue = queue
        self._as_wrapper = as_wrapper

    def consumeSpectrum(self, spec: oms.MSSpectrum) -> None:  # type: ignore[override]
        payload: Union[oms.MSSpectrum, Py_MSSpectrum]
        copied = oms.MSSpectrum(spec)
        payload = Py_MSSpectrum(copied) if self._as_wrapper else copied
        self._queue.put(payload)

    def consumeChromatogram(self, chrom) -> None:  # pragma: no cover - passthrough
        # Chromatogram streaming is currently unsupported.
        return None

    def setExperimentalSettings(self, settings) -> None:  # pragma: no cover
        return None

    def setExpectedSize(self, spectra: int, chromatograms: int) -> None:  # pragma: no cover
        return None


class _MzMLStream(AbstractContextManager[Iterator[Union[Py_MSSpectrum, oms.MSSpectrum]]]):
    """Context manager/iterator that yields spectra as they are parsed."""

    def __init__(self, path: Union[str, Path], *, as_wrapper: bool = True) -> None:
        _ensure_pyopenms_available()
        self._path = str(path)
        self._as_wrapper = as_wrapper
        self._queue: SimpleQueue = SimpleQueue()
        self._sentinel = object()
        self._thread: Thread | None = None
        self._error: Exception | None = None
        self._closed = False
        self._sentinel_seen = False

    def _start(self) -> None:
        if self._thread is not None:
            return

        consumer = _QueueingSpectrumConsumer(self._queue, self._as_wrapper)

        def _run() -> None:
            try:
                mzml = oms.MzMLFile()
                mzml.transform(self._path, consumer)
            except Exception as exc:  # pragma: no cover - exercised indirectly
                self._error = exc
            finally:
                self._queue.put(self._sentinel)

        self._thread = Thread(target=_run, name="MzMLStream")
        self._thread.daemon = True
        self._thread.start()

    # Context manager protocol -------------------------------------------------
    def __enter__(self) -> "_MzMLStream":
        self._start()
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        self.close()
        return False

    # Iterator protocol --------------------------------------------------------
    def __iter__(self) -> "_MzMLStream":
        self._start()
        return self

    def __next__(self):
        if self._closed:
            raise StopIteration

        item = self._queue.get()
        if item is self._sentinel:
            self._sentinel_seen = True
            self.close()
            if self._error:
                raise self._error
            raise StopIteration

        return item

    def close(self) -> None:
        """Drain the stream and wait for the parser thread to exit."""

        if self._closed:
            return

        self._closed = True

        if self._thread is None:
            return

        if not self._sentinel_seen:
            while True:
                item = self._queue.get()
                if item is self._sentinel:
                    break
            self._sentinel_seen = True

        self._thread.join()
        if self._error:
            raise self._error


def stream_mzml(path: Union[str, Path], *, as_wrapper: bool = True) -> _MzMLStream:
    """Return a streaming iterator over spectra in an mzML file.

    Example:
        >>> with stream_mzml("big.mzML") as spectra:
        ...     for spectrum in spectra:
        ...         process(spectrum)

    Args:
        path: Path to the mzML file on disk.
        as_wrapper: When True (default) yield :class:`Py_MSSpectrum` objects;
            otherwise return native ``pyopenms.MSSpectrum`` instances.
    """

    return _MzMLStream(path, as_wrapper=as_wrapper)
