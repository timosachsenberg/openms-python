# openms-python

**CAUTION: this package is under heavy development and is largely LLM generated. Breaking changes are expected and documentation might be out of date.** 

**A Pythonic wrapper around pyOpenMS for mass spectrometry data analysis**

`openms-python` provides an intuitive, Python-friendly interface to OpenMS, making mass spectrometry data analysis feel natural for Python developers and data scientists.

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-BSD--3-green)](LICENSE)

## Why openms-python?

[pyOpenMS](https://pyopenms.readthedocs.io/) is a Python binding for the powerful OpenMS C++ library. However, being a direct C++ binding, it doesn't always feel "Pythonic". This package wraps pyOpenMS to provide:

✅ **Pythonic properties** instead of verbose getters/setters
✅ **Intuitive iteration** with smart filtering
✅ **Identification helpers** for protein/peptide results with idXML IO
✅ **pandas DataFrame integration** for data analysis
✅ **Method chaining** for processing pipelines  
✅ **Type hints** for better IDE support  
✅ **Clean, documented API** with examples  

### Before (pyOpenMS)
```python
import pyopenms as oms

exp = oms.MSExperiment()
oms.MzMLFile().load("data.mzML", exp)

n_spectra = exp.getNrSpectra()
for i in range(n_spectra):
    spec = exp.getSpectrum(i)
    if spec.getMSLevel() == 1:
        rt = spec.getRT()
        peaks = spec.get_peaks()
	mz = peaks[0]
        intensity = peaks[1]
	print(f"RT: {spec.retention_time:.2f}s, Peaks: {spec.mz}, intens: {intensity}")
```

### After (openms-python)
```python
from openms_python import MSExperiment

exp = Py_MSExperiment.from_file("data.mzML")

print(f"Loaded {len(exp)} spectra")
for spec in exp.ms1_spectra():
    print(f"RT: {spec.retention_time:.2f}s, mz: {spec.mz}, intens: {spec.intensity}")


## or convert to pandas dataframe    
df = spec.to_dataframe()  # Get peaks as DataFrame
```

### Reading mzML Files

```python
from openms_python import Py_MSExperiment

# Load experiment
exp = Py_MSExperiment.from_file('data.mzML')

# Get basic info
print(f"Total spectra: {len(exp)}")
print(f"RT range: {exp.rt_range}")
print(f"MS levels: {exp.ms_levels}")

# Print summary
exp.print_summary()
```

# Identification workflows

Protein- and peptide-identification results often originate from search
engines that export **idXML** files. The `Identifications` helper reads these
files into convenient containers for downstream processing.

```python
from openms_python import Identifications

# Load both protein and peptide identifications from idXML
ids = Identifications.from_idxml("search_results.idXML")

print(ids.summary())
# {'proteins': 12, 'peptides': 42, 'protein_hits': 12, 'peptide_hits': 42}

# Filter peptides by their top-hit score while keeping the protein context
high_conf = ids.filter_peptides_by_score(0.05)

# Look up peptides that match a particular protein accession
matches = high_conf.peptides_for_protein("P01234")
for pep in matches:
    hit = pep.getHits()[0]
    print(hit.getSequence(), hit.getScore())

# Persist the curated results back to disk
high_conf.to_idxml("curated_results.idXML")
```

### Working with Spectra

```python
# Access individual spectra
spec = exp[0]

# Access multiple spectra with slicing
first_10 = exp[0:10]           # First 10 spectra
last_5 = exp[-5:]              # Last 5 spectra  
every_other = exp[::2]         # Every other spectrum
ms1_only = exp[1:4]            # Spectra 2-4 (0 indexing)

print(f"First spectrum: {spec}")
print(f"First 10 spectra: {len(first_10)} spectra")
print(f"Last 5 spectra: {len(last_5)} spectra")

# Pythonic properties
print(f"Retention time: {spec.retention_time} seconds")
print(f"MS level: {spec.ms_level}")
print(f"Number of peaks: {len(spec)}")
print(f"Total ion current: {spec.total_ion_current}")

# Boolean helpers
if spec.is_ms1:
    print("This is an MS1 spectrum")

# Get peaks as NumPy arrays
mz, intensity = spec.peaks

# Or as a DataFrame
peaks_df = spec.to_dataframe()
print(peaks_df.head())
```

### Peak Picking in One Line

#### Before (pyOpenMS)
```python
import pyopenms as oms

picker = oms.PeakPickerHiRes()
params = picker.getDefaults()
params.setValue("signal_to_noise", 3.0)
picker.setParameters(params)

centroided = oms.MSExperiment(exp)
picker.pickExperiment(exp, centroided, True)
```

#### After (openms-python)
```python
from openms_python import Py_MSExperiment

centroided = exp.pick_peaks(method="HiRes", params={"signal_to_noise": 3.0})
# or modify in-place
exp.pick_peaks(inplace=True)
```

### Smoothing Spectra

```python
# Apply a GaussFilter with a custom width
smoothed = exp.smooth_gaussian(gaussian_width=0.1)

# Smooth only MS2 spectra in-place using Savitzky-Golay
exp.smooth_savitzky_golay(ms_levels=2, inplace=True, frame_length=7)
```

### Smart Iteration

```python
# Iterate over MS1 spectra only
for spec in exp.ms1_spectra():
    print(f"MS1 at RT={spec.retention_time:.2f}s")

# Iterate over MS2 spectra only
for spec in exp.ms2_spectra():
    print(f"MS2: precursor m/z = {spec.precursor_mz:.4f}")

# Filter by retention time
for spec in exp.rt_filter[100:200]:
    print(f"Spectrum at RT={spec.retention_time:.2f}s")
```

### DataFrame Integration

```python
# Convert entire experiment to DataFrame
df = exp.to_dataframe(include_peaks=True)
print(df.head())

# Spectrum-level DataFrame
df_spectra = exp.to_dataframe(include_peaks=False)

# MS2 peaks only
df_ms2 = exp.to_dataframe(include_peaks=True, ms_level=2)
```

### Method Chaining

```python
# Filter and process in a pipeline
filtered_exp = (exp
    .filter_by_ms_level(1)
    .filter_by_rt(100, 500)
    .filter_by_mz(400, 500)
    .filter_top_n_peaks(100))

# OR
filtered_exp = (exp
    .filter_by_ms_level(1)
    .rt_filter[100:500]
    .mz_filter[400:500]
    .filter_top_n_peaks(100))

print(f"After filtering: {len(filtered_exp)} spectra")
```

### Data Manipulation

```python
# Filter peaks by m/z
filtered_spec = spec.filter_by_mz(100, 500)
# OR
filtered_spec = spec.mz_filter[100:500]

# Get top N peaks
top_10 = spec.top_n_peaks(10)
```

### Writing Files

```python
# Save to mzML
exp.to_mzml('output.mzML')

# Or use convenience function
from openms_python import write_mzml
write_mzml(exp, 'output.mzML')
```

### Wrapper-Aware Load/Store

All high-level containers expose `load`/`store` helpers that infer the correct
pyOpenMS reader from the file extension (including `.gz`).

#### Before (pyOpenMS)
```python
import pyopenms as oms

exp = oms.MSExperiment()
oms.MzMLFile().load("input.mzML", exp)
oms.MzMLFile().store("output.mzML", exp)

fmap = oms.FeatureMap()
oms.FeatureXMLFile().load("input.featureXML", fmap)
oms.FeatureXMLFile().store("output.featureXML", fmap)
```

#### After (openms-python)
```python
from openms_python import Py_MSExperiment, Py_FeatureMap, Py_ConsensusMap

exp = Py_MSExperiment().load("input.mzML")
exp.store("output.mzML")

feature_map = Py_FeatureMap().load("input.featureXML")
feature_map.store("output.featureXML")

cons_map = Py_ConsensusMap().load("input.consensusXML")
cons_map.store("output.consensusXML")
```

### Context Managers

```python
from openms_python.io import MzMLReader, MzMLWriter

# Reading with context manager
with MzMLReader('data.mzML') as exp:
    print(f"Loaded {len(exp)} spectra")
    for spec in exp.ms1_spectra():
        print(spec)

# Writing with context manager
with MzMLWriter('output.mzML') as writer:
    writer.write(exp)
```

### Streaming Large mzML Files

#### Before (pyOpenMS)
```python
import pyopenms as oms

class SpectrumCounter(oms.MSExperimentConsumer):
    def __init__(self):
        super().__init__()
        self.ms2 = 0

    def consumeSpectrum(self, spec):
        if spec.getMSLevel() == 2:
            self.ms2 += 1

consumer = SpectrumCounter()
oms.MzMLFile().transform("big.mzML", consumer)
print(f"Processed {consumer.ms2} MS2 spectra")
```

#### After (openms-python)
```python
from openms_python import stream_mzml

with stream_mzml("big.mzML") as spectra:
    ms2 = sum(1 for spec in spectra if spec.ms_level == 2)
print(f"Processed {ms2} MS2 spectra")
```

## High-level Workflows

### Feature Detection as a Method

```python
from openms_python import Py_MSExperiment

exp = Py_MSExperiment.from_file("data.mzML")
features = exp.detect_features(min_intensity=1_000.0)
print(f"Detected {len(features)} coarse features")
```

### Streaming FASTA Digestion and Theoretical Spectra

```python
from openms_python import FastaWorkflow

workflow = FastaWorkflow.from_file("proteins.fasta")
for entry in workflow.theoretical_spectra(charges=(2, 3)):
    spectrum = entry["spectrum"]
    print(entry["protein_id"], entry["peptide_sequence"], spectrum.getMSLevel())
```

### High-level FDR Filtering

```python
from openms_python import Identifications

ids = Identifications.from_idxml("search_results.idXML")
high_conf = ids.filter_by_fdr(max_q_value=0.01)
print(f"{len(high_conf.peptide_identifications)} peptides at 1% FDR")
```

### ID Mapping, Alignment and Linking

```python
from openms_python import Identifications, Py_MSExperiment, Py_FeatureMap

exp = Py_MSExperiment.from_file("data.mzML")
ids = Identifications.from_idxml("search_results.idXML")
high_conf = ids.filter_by_fdr()

features = exp.detect_features()
annotated = features.annotate_with_identifications(high_conf)
aligned = Py_FeatureMap.align([annotated])
consensus = Py_FeatureMap.link(aligned)
```

### Quantitation Table Export

```python
quant_table = consensus.to_quant_table()
print(quant_table.head())
```

### Pythonic Mutation Helpers

All wrappers behave like mutable Python sequences.

#### MSExperiment
```python
# Append and extend
exp.append(new_spec).extend(other_exp)

# Remove by index or slice
exp.remove(-1)
del exp[::2]
```

#### FeatureMap / ConsensusMap
```python
feature_map.append(feature)
feature_map.extend(iter_of_features)
feature_map.remove(0)     # delete by index

cons_map.append(cons_feature)
del cons_map[-3:]
```

Behind the scenes the wrappers copy the retained entries back into the
underlying pyOpenMS container, preserving meta data while exposing the
expected Python semantics. By contrast, pyOpenMS requires manually creating a
new container, copying every element except the ones you wish to remove, and
reassigning the result.

## Advanced Examples

### Creating Spectra from Scratch

```python
import pandas as pd
from openms_python import Py_MSSpectrum

# From DataFrame
df = pd.DataFrame({
    'mz': [100.0, 200.0, 300.0],
    'intensity': [50.0, 100.0, 75.0]
})

spec = Spectrum.from_dataframe(
    df, 
    retention_time=120.5,
    ms_level=1,
    native_id='spectrum=1'
)
```

### Creating Experiments from DataFrames

```python
# Create experiment from grouped data
df = pd.DataFrame({
    'spectrum_id': [0, 0, 1, 1, 2, 2],
    'mz': [100, 200, 150, 250, 120, 220],
    'intensity': [50, 100, 60, 110, 55, 105],
    'retention_time': [10.0, 10.0, 20.0, 20.0, 30.0, 30.0],
    'ms_level': [1, 1, 1, 1, 1, 1]
})

exp = Py_MSExperiment.from_dataframe(df)
```

### Analysis Workflow

```python
from openms_python import Py_MSExperiment
import pandas as pd
import matplotlib.pyplot as plt

# Load data
exp = Py_MSExperiment.from_file('data.mzML')

# Get MS2 spectra as DataFrame
df_ms2 = exp.to_dataframe(include_peaks=True, ms_level=2)

# Analyze precursor distribution
precursor_stats = df_ms2.groupby('precursor_mz').agg({
    'intensity': 'sum',
    'spectrum_index': 'count'
}).rename(columns={'spectrum_index': 'n_spectra'})

print(precursor_stats.head())

# Plot TIC over time
df_spectra = exp.to_dataframe(include_peaks=False)
plt.figure(figsize=(10, 4))
plt.plot(df_spectra['retention_time'], df_spectra['total_ion_current'])
plt.xlabel('Retention Time (s)')
plt.ylabel('Total Ion Current')
plt.title('TIC over time')
plt.show()
```

## API Reference

### Py_MSExperiment

**Properties:**
- `n_spectra`: Number of spectra
- `rt_range`: Tuple of (min_rt, max_rt)
- `ms_levels`: Set of MS levels present

**Methods:**
- `from_file(filepath)`: Load from mzML file (class method)
- `from_dataframe(df, group_by)`: Create from DataFrame (class method)
- `to_file(filepath)`: Save to mzML file
- `to_dataframe(include_peaks, ms_level)`: Convert to DataFrame
- `ms1_spectra()`: Iterator over MS1 spectra
- `ms2_spectra()`: Iterator over MS2 spectra
- `spectra_by_level(level)`: Iterator over specific MS level
- `spectra_in_rt_range(min_rt, max_rt)`: Iterator over RT range
- `filter_by_ms_level(level)`: Filter by MS level
- `filter_by_rt(min_rt, max_rt)`: Filter by RT range
- `filter_top_n_peaks(n)`: Keep top N peaks per spectrum
- `summary()`: Get summary statistics
- `print_summary()`: Print formatted summary

### Py_MSSpectrum

**Properties:**
- `retention_time`: RT in seconds
- `ms_level`: MS level (1, 2, etc.)
- `is_ms1`, `is_ms2`: Boolean helpers
- `precursor_mz`: Precursor m/z (MS2+)
- `precursor_charge`: Precursor charge (MS2+)
- `native_id`: Native spectrum ID
- `total_ion_current`: Sum of intensities
- `base_peak_mz`: m/z of most intense peak
- `base_peak_intensity`: Intensity of base peak
- `peaks`: Tuple of (mz_array, intensity_array)

**Methods:**
- `from_dataframe(df, **metadata)`: Create from DataFrame (class method)
- `to_dataframe()`: Convert to DataFrame
- `filter_by_mz(min_mz, max_mz)`: Filter peaks by m/z
- `filter_by_intensity(min_intensity)`: Filter peaks by intensity
- `top_n_peaks(n)`: Keep top N peaks
- `normalize_intensity(max_value)`: Normalize intensities

- `normalize_intensity(max_value)`: Normalize intensities

### Identifications, ProteinIdentifications & PeptideIdentifications

**Identifications** combines protein and peptide search results and keeps the
two collections synchronized.

**Constructors / IO:**
- `Identifications.from_idxml(path)`: Load protein & peptide IDs from an idXML file
- `Identifications.store(path)` / `to_idxml(path)`: Save to idXML

**Convenience helpers:**
- `Identifications.summary()`: Quick counts of IDs and hits
- `Identifications.find_protein(...)`, `find_peptide(...)`: Search by identifier
- `Identifications.find_protein_by_accession(...)`: Look up proteins by accession
- `Identifications.find_peptide_by_sequence(...)`: Case-insensitive peptide search
- `Identifications.filter_peptides_by_score(threshold)`: Return a copy with filtered peptides
- `Identifications.peptides_for_protein(accession)`: Retrieve peptides linked to a protein

Underlying containers (`ProteinIdentifications`, `PeptideIdentifications`) behave
like Python sequences: they support slicing, appending, iterating, and provide
convenience methods such as `find_by_identifier`, `find_by_accession`, and
`filter_by_score` for quickly triaging search hits.

## Development

### Setup Development Environment

```bash
git clone https://github.com/openms/openms-python.git
cd openms-python
pip install -e ".[dev]"
```

## Comparison with pyOpenMS

| Feature | pyOpenMS | openms-python |
|---------|----------|---------------|
| Get spectrum count | `exp.getNrSpectra()` | `len(exp)` |
| Get retention time | `spec.getRT()` | `spec.retention_time` |
| Check MS1 | `spec.getMSLevel() == 1` | `spec.is_ms1` |
| Load file | `MzMLFile().load(path, exp)` | `exp = MSExperiment.from_file(path)` |
| Iterate MS1 | Manual loop + level check | `for spec in exp.ms1_spectra():` |
| Peak data | `peaks = spec.get_peaks(); mz = peaks[0]` | `mz, intensity = spec.peaks` |
| DataFrame | Not available | `df = exp.to_dataframe()` |

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the BSD-3-Clause License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built on top of the excellent [pyOpenMS](https://pyopenms.readthedocs.io/) library
- Part of the [OpenMS](https://www.openms.de/) ecosystem

## Citation

If you use openms-python in your research, please cite:

```
Röst HL, Sachsenberg T, Aiche S, et al. OpenMS: a flexible open-source software platform 
for mass spectrometry data analysis. Nat Methods. 2016;13(9):741-748.
```

## Support

- **Documentation**: [https://openms-python.readthedocs.io](https://openms-python.readthedocs.io)
- **Issues**: [GitHub Issues](https://github.com/openms/openms-python/issues)
- **Discussions**: [GitHub Discussions](https://github.com/openms/openms-python/discussions)
