# gDNA barcode fastq → counts pipeline

A fork of [`SydShafferLab/BarcodeAnalysis_v2_gDNA`](https://github.com/SydShafferLab/BarcodeAnalysis_v2_gDNA)
(gDNA-only barcode extraction + starcode consensus clustering), with several
correctness bugs fixed, support added for a tag field adjacent to the barcode
(needed for multi-tag lineage-barcode constructs), and a single config-driven
driver script (`run_barcode_pipeline.py`) that replaces copy-pasting a new
one-off script per sequencing run.

Originally developed for the eGRK1059 project (RajLab); see `CHANGES.md` for
exactly what differs from upstream and why.

## Quick start

```bash
conda create -n barcode_analysis_gDNA python=3.9 -c bioconda -c conda-forge starcode
conda activate barcode_analysis_gDNA
pip install -r BarcodeAnalysis_v2_gDNA/ExtractBarcodes/requirements.txt

cp -r run_configs/TEMPLATE run_configs/<your_run_name>
# edit run_configs/<your_run_name>/config.json -- see run_configs/TEMPLATE/README.md
python run_barcode_pipeline.py run_configs/<your_run_name>/config.json
```

Output: `<out_dir>/starcode_outputs/sc_output_counts_<sample>.txt`, tab-separated
`barcode\ttag\tcount` (one row per consensus-barcode/tag pair per sample).

## Repo layout

- `run_barcode_pipeline.py` — the pipeline driver. Takes one JSON config and
  runs: (optional) fastq staging → per-read extraction (anchor match, tag
  call, quality filters) → per-group starcode clustering → per-sample count
  tables.
- `run_configs/TEMPLATE/` — a documented config template (`config.json`) and
  field-by-field `README.md`. Copy this per new run/dataset.
- `run_configs/EXAMPLE_single_tag/` — a worked example for the simpler case
  where every sample carries the same single tag (vs. TEMPLATE's example,
  which shows multiple tags distinguishing different sample families).
- `analysis/` — post-processing scripts (operate on `sc_output_counts_*.txt`,
  not raw fastqs): diagnosing whether an inflated barcode count is real
  cross-contamination from another tagged library vs. sequencing/PCR error,
  inferring the abundance cutoff below which even a "real"-tagged barcode is
  statistically indistinguishable from that contamination, and a pairwise
  frequency scatter for sanity-checking it across samples. See its own
  `README.md` for the full workflow.
- `BarcodeAnalysis_v2_gDNA/` — the vendored, patched extraction pipeline.
  `ExtractBarcodes/Scripts/Functions/gDNA_starcode_prep.py` is the core
  per-read extraction logic `run_barcode_pipeline.py` calls directly;
  `ExtractBarcodes/Scripts/Step3.py` is the original upstream driver (adapted
  for 10x CellRanger custom feature-barcoding output), kept here for
  reference/for anyone who does need that CellRanger output format, patched
  with the same bug fixes.
  - `starcode/` is vendored from [gui11aume/starcode](https://github.com/gui11aume/starcode)
    and ships a prebuilt Linux binary — doesn't run on macOS. Get starcode via
    conda/bioconda as in Quick start above instead of relying on this copy.

## Requirements

- Python packages: `BarcodeAnalysis_v2_gDNA/ExtractBarcodes/requirements.txt`
  (biopython, rapidfuzz, jstyleson, numpy, pandas, matplotlib, ...)
- [`starcode`](https://github.com/gui11aume/starcode) on `PATH` (install via
  bioconda, not the vendored copy above).

## See also

- `CHANGES.md` — detailed diff against upstream: bugs fixed and tag/quality
  features added.
