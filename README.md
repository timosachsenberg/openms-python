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

### Built-in example data

New to OpenMS or don't have data handy? `openms_python` ships with a tiny
``small.mzML`` example that is perfect for quick experiments.

```python
from openms_python import Py_MSExperiment, get_example

example_path = get_example("small.mzML")
exp = Py_MSExperiment.from_file(example_path)
print(f"Loaded {len(exp)} spectra from the example file")

# Or load the raw bytes directly
example_bytes = get_example("small.mzML", load=True)
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

## Consensus alignment and linking

Multiple `FeatureMap` runs can be aligned and converted into a single
`ConsensusMap` directly from Python. The `Py_ConsensusMap.align_and_link`
helper performs three steps:

1. Copies the incoming feature maps to avoid mutating your data
2. Aligns the feature maps with your choice of OpenMS alignment algorithm
3. Links the aligned runs using `FeatureGroupingAlgorithmQT`

```python
from openms_python import Py_FeatureMap, Py_ConsensusMap

feature_maps = [
    Py_FeatureMap.from_dataframe(run_a_df),
    Py_FeatureMap.from_dataframe(run_b_df),
]

consensus = Py_ConsensusMap.align_and_link(
    feature_maps,
    alignment_method="pose_clustering",  # or "identification" / "identity"
    alignment_params={"max_rt_shift": 15.0},
)

print(f"Consensus contains {len(consensus)} features")
```

The helper returns a fresh `Py_ConsensusMap` instance that can be exported,
converted to a pandas DataFrame, or iterated for downstream analysis.

## Protein inference and rollups

Recent wrappers expose multiple entry points for inferring proteins directly
from the Python API—either by starting from identification files, feature maps,
or full consensus maps.

```python
from openms_python import Identifications, Py_FeatureMap, Py_ConsensusMap

# 1) Run inference straight from an idXML file
ids = Identifications.from_idxml("search_results.idXML")
protein_summary = ids.infer_proteins(algorithm="bayesian")
print(protein_summary.summary())

# 2) Trigger inference on a feature map (assigned + unassigned peptides)
fmap = Py_FeatureMap().load("sample.featureXML")
proteins = fmap.infer_proteins(include_unassigned=True)
proteins.to_idxml("sample_proteins.idXML")

# 3) Operate directly on a consensus map
consensus = Py_ConsensusMap().load("merged.consensusXML")
consensus.infer_proteins(algorithm="basic")

# Optionally compute quantitative protein ratios in place
consensus.infer_protein_quantities(reference_map=1)
consensus.store("merged_with_proteins.consensusXML")
```

All helpers share the same ergonomic parameter handling, accept native
`pyopenms` parameters (`oms.Param`) or plain dictionaries, and return
`Identifications` or the map instance itself for easy method chaining.

## Identification performance showcase

Looking for a larger end-to-end example? `tests/test_idperformance.py` ships with
the repository as a miniature-yet-realistic identification workflow that ties
many wrapper conveniences together. The test builds a tiny FASTA database,
simulates MS2 spectra, runs a Hyperscore-style search (complete with target/decoy
competition and q-value estimation), and even records how long the round-trip
through mzML takes. It demonstrates how concise—and how performant—a simple
search engine can be when built with the high-level helpers in
`openms_python`. Reuse the test as inspiration for bespoke pipelines or as a
regression harness when experimenting with search-related utilities.

### Iterate over containers and metadata

All sequence-like wrappers (feature maps, consensus maps, identification containers,
and experiments) support Python's iteration protocol. Metadata-aware wrappers such
as :class:`Py_Feature` and :class:`Py_MSSpectrum` expose their meta values like a
regular mapping, so you can loop over keys or call ``len()``.

```python
from openms_python import Py_Feature, Py_FeatureMap, Identifications

# Assume ``feature_df`` is a pandas DataFrame with feature columns
feature_df = ...
fmap = Py_FeatureMap.from_dataframe(feature_df)
for feature in fmap:
    print("Feature UID", feature.getUniqueId())

ids = Identifications.from_idxml("results.idXML")
for peptide in ids:  # equivalent to iterating over ids.peptide_identifications
    print(peptide.getIdentifier())

feature = Py_Feature()
feature["label"] = "sample_a"
for key in feature:
    print(key, feature[key])
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

## Workflow helpers

`openms_python` now exposes opinionated utilities that combine the primitive
wrappers into streaming pipelines and end-to-end quantitation helpers:

```python
from openms_python import (
    Py_MSExperiment,
    Identifications,
    ProteinStream,
    map_identifications_to_features,
    align_feature_maps,
    link_features,
    export_quant_table,
)

# 1) Detect features in an experiment
experiment = Py_MSExperiment.from_file("run.mzML")
feature_map = experiment.detect_features()

# 2) Map identifications and filter them by FDR
identifications = Identifications.from_idxml("search_results.idXML")
filtered = identifications.filter_by_fdr(threshold=0.01)
annotated = map_identifications_to_features(feature_map, filtered)

# 3) Align multiple maps and link them into a consensus representation
aligned = align_feature_maps([annotated, second_run])
consensus = link_features(aligned)

# 4) Export a tidy quantitation table (per-sample intensities)
quant_df = export_quant_table(consensus)

# Bonus: streaming FASTA digestion & theoretical spectra generation
for record in ProteinStream.from_fasta("proteins.fasta").digest().theoretical_spectra():
    print(record.protein.identifier, record.peptide.toString(), len(record.spectrum))
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

```python
# DataFrame round-trip
df = feature_map.to_dataframe()
df["mz"] += 0.01  # manipulate with pandas
feature_map = Py_FeatureMap.from_dataframe(df)

cons_df = cons_map.get_df()
cons_df["quality"] = 0.95
cons_map = Py_ConsensusMap.from_df(cons_df)

peaks_df = exp.get_df()
peaks_df["intensity"] *= 1.1
exp = Py_MSExperiment.from_df(peaks_df)
```

Behind the scenes the wrappers copy the retained entries back into the
underlying pyOpenMS container, preserving meta data while exposing the
expected Python semantics. By contrast, pyOpenMS requires manually creating a
new container, copying every element except the ones you wish to remove, and
reassigning the result.

### Dictionary-Style Meta Data Access

Any pyOpenMS type that derives from `MetaInfoInterface` (features, spectra,
consensus features, etc.) now behaves like a standard Python mapping for its
meta annotations:

```python
import pyopenms as oms
from openms_python import Py_Feature, Py_MSSpectrum

feature = Py_Feature()
feature["label"] = "Sample1"
feature.update(condition="control", replicate="R1")
assert feature["label"] == "Sample1"
assert feature.get("missing", "n/a") == "n/a"

spectrum = Py_MSSpectrum(oms.MSSpectrum())
spectrum["IonInjectTime"] = 12.3
assert "IonInjectTime" in spectrum
spectrum.pop("IonInjectTime")
```

Internally this syntax delegates to the familiar `setMetaValue`,
`getMetaValue`, and `removeMetaValue` calls on the wrapped pyOpenMS object, so
no information is lost compared to the C++ interface.

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
