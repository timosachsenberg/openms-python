"""
Pythonic wrapper for pyOpenMS MSExperiment class.
"""

from pathlib import Path
from typing import Iterator, Optional, List, Union, Sequence, Dict, Any, Iterable, Callable, Set
import pandas as pd
import numpy as np
import pyopenms as oms
from .py_msspectrum import Py_MSSpectrum
from ._io_utils import ensure_allowed_suffix, MS_EXPERIMENT_EXTENSIONS


PEAK_PICKER_REGISTRY: Dict[str, Any] = {
    "hires": oms.PeakPickerHiRes,
    "cwt": getattr(oms, "PeakPickerCWT", oms.PeakPickerHiRes),
}


class _Py_MSExperimentSlicing:
    """Helper class for advanced slicing operations on Py_MSExperiment."""
    
    def __init__(self, experiment: 'Py_MSExperiment'):
        self.experiment = experiment
    
    def __getitem__(self, key) -> 'Py_MSExperiment':
        """Handle slice syntax for m/z ranges."""
        if isinstance(key, slice):
            if key.start is not None and key.stop is not None:
                # m/z slice: experiment.mz_filter[400.2:500.9]
                return self.experiment.filter_by_mz_range(key.start, key.stop)
            else:
                raise ValueError("Both start and stop values required for m/z slice")
        else:
            raise TypeError(f"Invalid slice type: {type(key)}")


class _Py_MSExperimentRTSlicing:
    """Helper class for retention time slicing operations on MSExperiment."""
    
    def __init__(self, experiment: 'Py_MSExperiment'):
        self.experiment = experiment
    
    def __getitem__(self, key) -> 'Py_MSExperiment':
        """Handle slice syntax for retention time ranges."""
        if isinstance(key, slice):
            if key.start is not None and key.stop is not None:
                # Retention time slice: experiment.rt_filter[2010.0:2010.5]
                return self.experiment.filter_by_rt_range(key.start, key.stop)
            else:
                raise ValueError("Both start and stop values required for retention time slice")
        else:
            raise TypeError(f"Invalid slice type: {type(key)}")
    
    def __call__(self, rt_min: float, rt_max: float) -> 'Py_MSExperiment':
        """Handle function call syntax: exp.rt_filter(rt_min, rt_max)"""
        return self.experiment.filter_by_rt_range(rt_min, rt_max)
    
class _Py_MSExperimentRTSlicing:
    """Helper class for retention time slicing operations on MSExperiment."""
    
    def __init__(self, experiment: 'Py_MSExperiment'):
        self.experiment = experiment
    
    def __getitem__(self, key) -> 'Py_MSExperiment':
        """Handle slice syntax for retention time ranges."""
        if isinstance(key, slice):
            if key.start is not None and key.stop is not None:
                # Retention time slice: experiment.rt[2010.0:2010.5]
                return self.experiment.filter_by_rt_range(key.start, key.stop)
            else:
                raise ValueError("Both start and stop values required for retention time slice")
        else:
            raise TypeError(f"Invalid slice type: {type(key)}")
    
    def __call__(self, rt_min: float, rt_max: float) -> 'Py_MSExperiment':
        """Handle function call syntax: exp.rt(rt_min, rt_max)"""
        return self.experiment.filter_by_rt_range(rt_min, rt_max)


class Py_MSExperiment:
    """
    A Pythonic wrapper around pyOpenMS MSExperiment.
    
    This class provides intuitive iteration, properties, and DataFrame
    conversion for mass spectrometry experiments.
    
    Example:
        >>> exp = MSExperiment.from_file('data.mzML')
        >>> print(f"Loaded {len(exp)} spectra")
        >>> 
        >>> # Iterate over MS1 spectra
        >>> for spec in exp.ms1_spectra():
        ...     print(f"RT: {spec.retention_time:.2f}")
        >>> 
        >>> # Convert to DataFrame
        >>> df = exp.to_dataframe()
        >>> 
        >>> # Method chaining
        >>> filtered = exp.filter_by_ms_level(1).filter_by_rt(100, 200)
    """
    
    def __init__(self, native_experiment: Optional[oms.MSExperiment] = None):
        """
        Initialize MSExperiment wrapper.
        
        Args:
            native_experiment: pyOpenMS MSExperiment object (creates new if None)
        """
        self._experiment = native_experiment if native_experiment is not None else oms.MSExperiment()
    
    # ==================== Class Methods for Creation ====================
    
    @classmethod
    def from_mzml(cls, filepath: str) -> 'Py_MSExperiment':
        """
        Load an MSExperiment from an mzML file.
        
        Args:
            filepath: Path to mzML file
            
        Returns:
            MSExperiment object
            
        Example:
            >>> exp = MSExperiment.from_file('data.mzML')
        """
        exp = oms.MSExperiment()
        oms.MzMLFile().load(filepath, exp)
        return cls(exp)
    
    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        mz: str = 'mz',
        intensity: str = 'intensity',
        rt: str = 'retention_time',
        ms_level: str = 'ms_level',
        native_id: str = 'native_id',
        spec_idx: str = 'spectrum_id',
    ) -> 'Py_MSExperiment':
        """Create an experiment from a peak-level :class:`pandas.DataFrame`."""

        required = {spec_idx, mz, intensity}
        missing = required - set(df.columns)
        if missing:
            missing_str = ", ".join(sorted(missing))
            raise ValueError(f"DataFrame is missing required columns: {missing_str}")

        exp = oms.MSExperiment()

        grouped = df.groupby(spec_idx, sort=False)
        for _, spectrum_df in grouped:
            spec = oms.MSSpectrum()
            spec.set_peaks((
                spectrum_df[mz].astype(float).tolist(),
                spectrum_df[intensity].astype(float).tolist(),
            ))

            if rt in spectrum_df.columns:
                spec.setRT(float(spectrum_df[rt].iloc[0]))

            if ms_level in spectrum_df.columns:
                spec.setMSLevel(int(spectrum_df[ms_level].iloc[0]))

            if native_id in spectrum_df.columns:
                spec.setNativeID(str(spectrum_df[native_id].iloc[0]))

            exp.addSpectrum(spec)

        return cls(exp)

    @classmethod
    def from_df(cls, df: pd.DataFrame, **kwargs) -> 'Py_MSExperiment':
        """Alias for :meth:`from_dataframe` matching :meth:`get_df`."""

        return cls.from_dataframe(df, **kwargs)
    
    def to_mzml(self, filepath: str):
        """
        Save MSExperiment to an mzML file.

        Args:
            filepath: Output path for mzML file

        Example:
            >>> exp.to_file('output.mzML')
        """
        oms.MzMLFile().store(filepath, self._experiment)

    def load(self, filepath: Union[str, Path]) -> 'Py_MSExperiment':
        """Load data from an MS file using its extension for detection."""

        ensure_allowed_suffix(filepath, MS_EXPERIMENT_EXTENSIONS, "MSExperiment")
        oms.FileHandler().loadExperiment(str(filepath), self._experiment)
        return self

    def store(self, filepath: Union[str, Path]) -> 'Py_MSExperiment':
        """Store the experiment to disk based on the output extension."""

        ensure_allowed_suffix(filepath, MS_EXPERIMENT_EXTENSIONS, "MSExperiment")
        oms.FileHandler().storeExperiment(str(filepath), self._experiment)
        return self
    
    # ==================== Pythonic Properties ====================
    
    @property
    def nr_spectra(self) -> int:
        """Get number of spectra in the experiment."""
        return self._experiment.getNrSpectra()
    
    @property
    def rt_range(self) -> tuple:
        """
        Get retention time range.
        
        Returns:
            Tuple of (min_rt, max_rt)
        """
        if len(self) == 0:
            return (0.0, 0.0)
        rts = [spec.retention_time for spec in self]
        return (min(rts), max(rts))
    
    @property
    def ms_levels(self) -> set:
        """Get set of MS levels present in the experiment."""
        return {spec.ms_level for spec in self}
    
    @property
    def mz_filter(self) -> _Py_MSExperimentSlicing:
        """
        Access m/z-based slicing operations.
        
        This property provides slice-like syntax for m/z range filtering.
        
        Example:
            >>> # Filter by m/z range 400-500
            >>> filtered = exp.mz_filter[400.2:500.9]
            >>>
            >>> # Chain with spectrum slicing
            >>> subset = exp[1:5].mz_filter[400.2:500.9]
        """
        return _Py_MSExperimentSlicing(self)
    
    @property
    def rt_filter(self) -> _Py_MSExperimentRTSlicing:
        """
        Access retention time-based slicing operations.
        
        This property provides slice-like syntax for retention time range filtering.
        
        Example:
            >>> # Filter by retention time range 2010-2010.5 seconds
            >>> filtered = exp.rt_filter[2010.0:2010.5]
            >>>
            >>> # Alternative function call syntax
            >>> filtered = exp.rt_filter(2010.0, 2010.5)
            >>>
            >>> # Chain with spectrum slicing
            >>> subset = exp[1:5].rt_filter[2010.0:2010.5]
        """
        return _Py_MSExperimentRTSlicing(self)
    
    # ==================== Magic Methods ====================
    
    def __len__(self) -> int:
        """Return number of spectra."""
        return self.nr_spectra
    
    def __iter__(self) -> Iterator[Py_MSSpectrum]:
        """Iterate over all spectra."""
        for i in range(len(self)):
            yield self[i]
    
    def __repr__(self) -> str:
        """Return string representation."""
        ms_levels_str = ', '.join(f"MS{level}" for level in sorted(self.ms_levels))
        rt_min, rt_max = self.rt_range
        return (
            f"MSExperiment(spectra={len(self)}, "
            f"levels=[{ms_levels_str}], "
            f"rt_range=[{rt_min:.2f}, {rt_max:.2f}]s)"
        )
    
    # ==================== Smart Iteration ====================
    
    def ms1_spectra(self) -> Iterator[Py_MSSpectrum]:
        """
        Iterate over MS1 spectra only.
        
        Example:
            >>> for spec in exp.ms1_spectra():
            ...     print(f"MS1 at RT={spec.retention_time:.2f}")
        """
        for spec in self:
            if spec.is_ms1:
                yield spec
    
    def ms2_spectra(self) -> Iterator[Py_MSSpectrum]:
        """
        Iterate over MS2 spectra only.
        
        Example:
            >>> for spec in exp.ms2_spectra():
            ...     print(f"MS2 precursor: {spec.precursor_mz:.4f}")
        """
        for spec in self:
            if spec.is_ms2:
                yield spec
    
    def spectra_by_level(self, ms_level: int) -> Iterator[Py_MSSpectrum]:
        """
        Iterate over spectra of a specific MS level.
        
        Args:
            ms_level: MS level to filter (1, 2, etc.)
            
        Example:
            >>> for spec in exp.spectra_by_level(2):
            ...     print(spec)
        """
        for spec in self:
            if spec.ms_level == ms_level:
                yield spec
    
    def spectra_in_rt_range(self, min_rt: float, max_rt: float) -> Iterator[Py_MSSpectrum]:
        """
        Iterate over spectra within retention time range.
        
        Args:
            min_rt: Minimum retention time (seconds)
            max_rt: Maximum retention time (seconds)
            
        Example:
            >>> for spec in exp.spectra_in_rt_range(100, 200):
            ...     print(f"RT: {spec.retention_time:.2f}")
        """
        for spec in self:
            if min_rt <= spec.retention_time <= max_rt:
                yield spec
    
    # ==================== DataFrame Conversion ====================
    
    def to_dataframe(self, include_peaks: bool = True, ms_level: Optional[int] = None) -> pd.DataFrame:
        """
        Convert experiment to pandas DataFrame.
        
        Args:
            include_peaks: If True, include individual peaks; if False, spectrum-level only
            ms_level: Optional MS level filter (None = all levels)
            
        Returns:
            DataFrame with spectrum and peak information
            
        Example:
            >>> # Spectrum-level DataFrame
            >>> df_spectra = exp.to_dataframe(include_peaks=False)
            >>> 
            >>> # Peak-level DataFrame for MS1 only
            >>> df_peaks = exp.to_dataframe(include_peaks=True, ms_level=1)
        """
        if include_peaks:
            return self._peaks_dataframe(ms_level)
        else:
            return self._spectra_dataframe(ms_level)

    def get_df(self, include_peaks: bool = True, ms_level: Optional[int] = None) -> pd.DataFrame:
        """Alias for :meth:`to_dataframe` for backwards compatibility."""

        return self.to_dataframe(include_peaks=include_peaks, ms_level=ms_level)
    
    def _spectra_dataframe(self, ms_level: Optional[int] = None) -> pd.DataFrame:
        """Create spectrum-level DataFrame."""
        data = []
        
        spectra = self.spectra_by_level(ms_level) if ms_level else self
        
        for idx, spec in enumerate(spectra):
            row = {
                'spectrum_index': idx,
                'native_id': spec.native_id,
                'retention_time': spec.retention_time,
                'ms_level': spec.ms_level,
                'n_peaks': len(spec),
                'total_ion_current': spec.total_ion_current,
                'base_peak_mz': spec.base_peak_mz,
                'base_peak_intensity': spec.base_peak_intensity,
            }
            
            # Add precursor info for MS2+
            if spec.ms_level >= 2:
                row['precursor_mz'] = spec.precursor_mz
                row['precursor_charge'] = spec.precursor_charge
            
            data.append(row)
        
        return pd.DataFrame(data)
    
    def _peaks_dataframe(self, ms_level: Optional[int] = None) -> pd.DataFrame:
        """Create peak-level DataFrame."""
        data = []
        
        spectra = self.spectra_by_level(ms_level) if ms_level else self
        
        for spec_idx, spec in enumerate(spectra):
            mz, intensity = spec.peaks
            
            for peak_idx in range(len(mz)):
                row = {
                    'spectrum_index': spec_idx,
                    'peak_index': peak_idx,
                    'mz': mz[peak_idx],
                    'intensity': intensity[peak_idx],
                    'retention_time': spec.retention_time,
                    'ms_level': spec.ms_level,
                }
                
                if spec.ms_level >= 2:
                    row['precursor_mz'] = spec.precursor_mz
                    row['precursor_charge'] = spec.precursor_charge
                
                data.append(row)
        
        return pd.DataFrame(data)
    
    # ==================== Filtering (Method Chaining) ====================
    def filter_by_ms_level(self, ms_level: int) -> 'Py_MSExperiment':
        """
        Filter spectra by MS level.
        
        Args:
            ms_level: MS level to keep
            
        Returns:
            New MSExperiment with filtered spectra
            
        Example:
            >>> ms1_exp = exp.filter_by_ms_level(1)
        """
        new_exp = oms.MSExperiment()
        for spec in self.spectra_by_level(ms_level):
            new_exp.addSpectrum(spec.native)
        return Py_MSExperiment(new_exp)
    
    def filter_by_spectrum_index(self, min_index, max_index, step=1) -> 'Py_MSExperiment':
        """
        Filter spectra by their index in the experiment.
        
        Args:
            min_index: Minimum spectrum index (inclusive)
            max_index: Maximum spectrum index (inclusive)
            step: Step size for slicing (default is 1)
            
        Returns:
            New Py_MSExperiment with filtered spectra
            
        Example:
            >>> subset = exp.filter_by_spectrum_index(0, 5, step=2)
        """

        if min_index < 0 or min_index >= len(self):
            raise IndexError(f"Start index {min_index} out of range [0, {len(self)})")
        if max_index > len(self):
            raise ValueError(f"Stop index {max_index} out of range [0, {len(self)})")

        new_exp = oms.MSExperiment()
        indices = range(min_index, max_index, step)

        for i in indices:
            new_exp.addSpectrum(self._experiment.getSpectrum(i))

        return Py_MSExperiment(new_exp)

    
    def get_by_spectrum_index(self, index: int) -> Py_MSSpectrum:
        """
        Get a single spectrum by its index.
        
        Args:
            index: Spectrum index
            
        Returns:
            Spectrum object
            
        Example:
            >>> spec = exp.get_by_spectrum_index(10)
            >>> print(f"RT: {spec.retention_time:.2f}")
        """
        # Handle integer index - return single Spectrum
        if index < 0:
            index = len(self) + index  # Convert negative index

        if index < 0 or index >= len(self):
            raise IndexError(f"Spectrum index {index} out of range [0, {len(self)})")

        return Py_MSSpectrum(self._experiment.getSpectrum(index))

    def filter_by_mz_range(self, mz_min: float, mz_max: float) -> 'Py_MSExperiment':
        """
        Filter spectra to keep only peaks within m/z range.
        
        This method filters the peaks in each spectrum to keep only those
        within the specified m/z range. The filtering is done on each
        spectrum individually.
        
        Args:
            mz_min: Minimum m/z value to keep
            mz_max: Maximum m/z value to keep
            
        Returns:
            New Py_MSExperiment with filtered peaks
            
        Example:
            >>> # Filter all spectra to keep m/z 400-500
            >>> filtered = exp.filter_by_mz_range(400.0, 500.0)
            >>> 
            >>> # Chain with slicing: get spectra 1-5, then filter peaks
            >>> subset = exp[1:5].filter_by_mz_range(400.0, 500.0)
            >>> 
            >>> # Or use tuple indexing
            >>> subset = exp[1:5][(400.0, 500.0)]
        """
        new_exp = oms.MSExperiment()
        
        for i in range(len(self)):
            spec = self[i]
            # Get the native spectrum and apply m/z filtering
            native_spec = spec.native
            
            # Get current peaks
            mz_array, intensity_array = native_spec.get_peaks()
            
            # Filter by m/z range
            mask = (mz_array >= mz_min) & (mz_array <= mz_max)
            filtered_mz = mz_array[mask]
            filtered_intensity = intensity_array[mask]
            
            # Create new spectrum with filtered peaks
            new_spec = oms.MSSpectrum()
            new_spec.setRT(native_spec.getRT())
            new_spec.setMSLevel(native_spec.getMSLevel())
            new_spec.set_peaks((filtered_mz, filtered_intensity))
            
            # Copy meta data if available
            if hasattr(native_spec, 'getPrecursors') and len(native_spec.getPrecursors()) > 0:
                new_spec.setPrecursors(native_spec.getPrecursors())
            
            new_exp.addSpectrum(new_spec)
        
        return Py_MSExperiment(new_exp)
    
    def filter_by_rt_range(self, rt_min: float, rt_max: float) -> 'Py_MSExperiment':
        """
        Filter spectra to keep only those within retention time range.
        
        This method filters the experiment to keep only spectra that have
        retention times within the specified range.
        
        Args:
            rt_min: Minimum retention time (seconds)
            rt_max: Maximum retention time (seconds)
            
        Returns:
            New MSExperiment with filtered spectra
            
        Example:
            >>> # Filter spectra by retention time 2010-2010.5 seconds
            >>> filtered = exp.filter_by_rt_range(2010.0, 2010.5)
            >>> 
            >>> # Chain with spectrum slicing: get spectra 1-5, then filter by RT
            >>> subset = exp[1:5].filter_by_rt_range(2010.0, 2010.5)
            >>> 
            >>> # Or use tuple indexing with explicit RT marker
            >>> subset = exp[(2010.0, 2010.5, 'rt')]
            >>> 
            >>> # Or use property slicing
            >>> subset = exp[1:5].rt[2010.0:2010.5]
        """
        new_exp = oms.MSExperiment()
        for spec in self.spectra_in_rt_range(rt_min, rt_max):
            new_exp.addSpectrum(spec.native)
        return Py_MSExperiment(new_exp)

    def filter_top_n_peaks(self, n: int) -> 'Py_MSExperiment':
        """
        Keep only top N peaks in each spectrum.
        
        Args:
            n: Number of peaks to keep per spectrum
            
        Returns:
            New MSExperiment with filtered spectra
            
        Example:
            >>> top_peaks = exp.filter_top_n_peaks(100)
        """
        new_exp = oms.MSExperiment()
        for spec in self:
            filtered_spec = spec.top_n_peaks(n)
            new_exp.addSpectrum(filtered_spec.native)
        return Py_MSExperiment(new_exp)

    def pick_peaks(
        self,
        method: str = "HiRes",
        params: Optional[Dict[str, Any]] = None,
        ms_levels: Optional[Union[int, Sequence[int]]] = 1,
        inplace: bool = False,
    ) -> 'Py_MSExperiment':
        """
        Convenience interface for pyOpenMS peak pickers.

        Args:
            method: Name of the peak picking algorithm ("HiRes", "CWT", ...).
            params: Optional dictionary of parameters passed to the picker.
            ms_levels: Single MS level, list of MS levels, or None for all levels.
            inplace: If True, modify the current experiment and return self.

        Returns:
            A Py_MSExperiment with picked spectra (self when inplace=True).

        Example:
            >>> picked = exp.pick_peaks(method="HiRes", ms_levels=1)
        """

        picker_cls = PEAK_PICKER_REGISTRY.get(method.lower())
        if picker_cls is None:
            raise ValueError(
                f"Unknown peak picking method '{method}'. "
                f"Available options: {sorted(PEAK_PICKER_REGISTRY.keys())}"
            )

        picker = picker_cls()

        if params:
            picker_params = picker.getParameters()
            for key, value in params.items():
                picker_params.setValue(key, value)
            picker.setParameters(picker_params)

        if ms_levels is None:
            target_levels = None
        elif isinstance(ms_levels, int):
            target_levels = {int(ms_levels)}
        else:
            target_levels = {int(level) for level in ms_levels}

        working_exp = oms.MSExperiment(self._experiment)
        working_exp.clear(True)

        source_exp = self._experiment

        for idx in range(source_exp.getNrSpectra()):
            source_spec = source_exp.getSpectrum(idx)
            target_spec = oms.MSSpectrum(source_spec)

            if target_levels is None or source_spec.getMSLevel() in target_levels:
                picker.pick(source_spec, target_spec)

            working_exp.addSpectrum(target_spec)

        if inplace:
            self._experiment = working_exp
            return self

        return Py_MSExperiment(working_exp)

    def smooth_gaussian(
        self,
        ms_levels: Optional[Union[int, Sequence[int]]] = None,
        inplace: bool = False,
        params: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> 'Py_MSExperiment':
        """Apply a Gaussian smoothing filter to spectra.

        Args:
            ms_levels: Optional single MS level or iterable of MS levels to smooth.
                When ``None`` (default) all spectra are smoothed.
            inplace: When ``True`` the current experiment is modified and returned.
            params: Optional dictionary of GaussFilter parameters.
            **kwargs: Additional GaussFilter parameters passed as keyword arguments.

        Returns:
            A Py_MSExperiment with smoothed spectra (``self`` when ``inplace=True``).

        Example:
            >>> smoothed = exp.smooth_gaussian(gaussian_width=0.1)
            >>> exp.smooth_gaussian(ms_levels=1, inplace=True)
        """

        smoother = oms.GaussFilter()
        configured = self._configure_filter(smoother, params, kwargs)
        return self._apply_spectrum_transform(
            configured.filter,
            ms_levels=ms_levels,
            inplace=inplace,
        )

    def smooth_savitzky_golay(
        self,
        ms_levels: Optional[Union[int, Sequence[int]]] = None,
        inplace: bool = False,
        params: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> 'Py_MSExperiment':
        """Apply a Savitzky-Golay smoothing filter to spectra.

        Args:
            ms_levels: Optional single MS level or iterable of MS levels to smooth.
            inplace: When ``True`` the current experiment is modified and returned.
            params: Optional dictionary of filter parameters.
            **kwargs: Additional filter parameters (e.g. ``frame_length``).

        Returns:
            A Py_MSExperiment with smoothed spectra (``self`` when ``inplace=True``).

        Example:
            >>> exp.smooth_savitzky_golay(frame_length=7)
            >>> exp.smooth_savitzky_golay(ms_levels=[2], inplace=True)
        """

        smoother = oms.SavitzkyGolayFilter()
        configured = self._configure_filter(smoother, params, kwargs)
        return self._apply_spectrum_transform(
            configured.filter,
            ms_levels=ms_levels,
            inplace=inplace,
        )


    def __getitem__(self, key) -> Union['Py_MSExperiment', Py_MSSpectrum]:
        """Return spectra using Python's indexing semantics."""

        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            new_experiment = oms.MSExperiment()

            for idx in range(start, stop, step):
                new_experiment.addSpectrum(oms.MSSpectrum(self._experiment.getSpectrum(idx)))

            return Py_MSExperiment(new_experiment)

        if isinstance(key, int):
            return self.get_by_spectrum_index(key)

        raise TypeError(f"Invalid index type: {type(key)}")

    def append(self, spectrum: Union[Py_MSSpectrum, oms.MSSpectrum]) -> 'Py_MSExperiment':
        """Append a spectrum to the experiment."""

        native = self._coerce_spectrum(spectrum)
        self._experiment.addSpectrum(native)
        return self

    def extend(self, spectra: Iterable[Union[Py_MSSpectrum, oms.MSSpectrum]]) -> 'Py_MSExperiment':
        """Append multiple spectra to the experiment."""

        for spectrum in spectra:
            self.append(spectrum)
        return self

    def remove(self, index: int) -> 'Py_MSExperiment':
        """Remove the spectrum at *index* and return ``self`` for chaining."""

        normalized = self._normalize_index(index)
        self._delete_indices([normalized])
        return self

    def __delitem__(self, key) -> None:
        """Delete spectra using integer indices or slices."""

        if isinstance(key, int):
            self.remove(key)
            return

        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            indices = list(range(start, stop, step))
            self._delete_indices(indices)
            return

        raise TypeError(f"Invalid deletion index type: {type(key)}")
    
    # ==================== Analysis Methods ====================
    def summary(self) -> dict:
        """
        Get summary statistics of the experiment.
        
        Returns:
            Dictionary with summary information
            
        Example:
            >>> stats = exp.summary()
            >>> print(f"Total spectra: {stats['nr_spectra']}")
        """
        ms1_count = sum(1 for _ in self.ms1_spectra())
        ms2_count = sum(1 for _ in self.ms2_spectra())
        rt_min, rt_max = self.rt_range
        
        total_peaks = sum(len(spec) for spec in self)
        avg_peaks = total_peaks / len(self) if len(self) > 0 else 0
        
        return {
            'nr_spectra': len(self),
            'n_ms1_spectra': ms1_count,
            'n_ms2_spectra': ms2_count,
            'ms_levels': sorted(self.ms_levels),
            'rt_range': (rt_min, rt_max),
            'total_peaks': total_peaks,
            'avg_peaks_per_spectrum': avg_peaks,
        }
    
    def print_summary(self):
        """
        Print a formatted summary of the experiment.
        
        Example:
            >>> exp.print_summary()
            MSExperiment Summary
            ===================
            Total Spectra: 1000
            MS1 Spectra: 200
            MS2 Spectra: 800
            ...
        """
        stats = self.summary()
        print("MSExperiment Summary")
        print("=" * 50)
        print(f"Total Spectra: {stats['nr_spectra']}")
        print(f"MS1 Spectra: {stats['n_ms1_spectra']}")
        print(f"MS2 Spectra: {stats['n_ms2_spectra']}")
        print(f"MS Levels: {stats['ms_levels']}")
        print(f"RT Range: {stats['rt_range'][0]:.2f} - {stats['rt_range'][1]:.2f} seconds")
        print(f"Total Peaks: {stats['total_peaks']}")
        print(f"Avg Peaks/Spectrum: {stats['avg_peaks_per_spectrum']:.1f}")
    
    @property
    def native(self) -> oms.MSExperiment:
        """
        Get the underlying pyOpenMS MSExperiment object.
        
        Use this when you need to access pyOpenMS-specific methods
        not wrapped by this class.
        """
        return self._experiment

    # ==================== Private Helpers ====================

    def _coerce_spectrum(self, spectrum: Union[Py_MSSpectrum, oms.MSSpectrum]) -> oms.MSSpectrum:
        if isinstance(spectrum, Py_MSSpectrum):
            return oms.MSSpectrum(spectrum.native)
        if isinstance(spectrum, oms.MSSpectrum):
            return oms.MSSpectrum(spectrum)
        raise TypeError(
            "append expects pyopenms.MSSpectrum or Py_MSSpectrum instances"
        )

    def _normalize_index(self, index: int) -> int:
        length = len(self)
        if length == 0:
            raise IndexError("MSExperiment is empty")
        if index < 0:
            index += length
        if index < 0 or index >= length:
            raise IndexError(f"Spectrum index {index} out of range [0, {length})")
        return index

    def _delete_indices(self, indices: Iterable[int]) -> None:
        drop = sorted(set(indices))
        if not drop:
            return

        source_exp = self._experiment
        length = source_exp.getNrSpectra()
        drop_set = set(drop)

        new_exp = oms.MSExperiment(source_exp)
        new_exp.clear(False)

        for idx in range(length):
            if idx in drop_set:
                continue
            new_exp.addSpectrum(oms.MSSpectrum(source_exp.getSpectrum(idx)))

        self._experiment = new_exp

    def _configure_filter(
        self,
        filter_obj: Any,
        params: Optional[Dict[str, Any]],
        extra_params: Dict[str, Any],
    ) -> Any:
        """Return the filter object after applying the provided parameters."""

        combined: Dict[str, Any] = {}
        if params:
            combined.update(params)
        if extra_params:
            combined.update(extra_params)

        if combined:
            filter_params = filter_obj.getParameters()
            for key, value in combined.items():
                filter_params.setValue(str(key), value)
            filter_obj.setParameters(filter_params)

        return filter_obj

    def _apply_spectrum_transform(
        self,
        transform: Callable[[oms.MSSpectrum], None],
        ms_levels: Optional[Union[int, Sequence[int]]],
        inplace: bool,
    ) -> 'Py_MSExperiment':
        target_levels = self._normalize_ms_levels(ms_levels)

        source_exp = self._experiment
        working_exp = oms.MSExperiment(source_exp)
        working_exp.clear(True)

        for idx in range(source_exp.getNrSpectra()):
            spectrum = oms.MSSpectrum(source_exp.getSpectrum(idx))
            if target_levels is None or spectrum.getMSLevel() in target_levels:
                transform(spectrum)
            working_exp.addSpectrum(spectrum)

        if inplace:
            self._experiment = working_exp
            return self

        return Py_MSExperiment(working_exp)

    def _normalize_ms_levels(
        self, ms_levels: Optional[Union[int, Sequence[int]]]
    ) -> Optional[Set[int]]:
        if ms_levels is None:
            return None
        if isinstance(ms_levels, int):
            return {int(ms_levels)}
        return {int(level) for level in ms_levels}
