# Config fields for run_barcode_pipeline.py

Copy this folder to `run_configs/<your_run_name>/`, edit `config.json` (plain JSON,
no comments), then:

```
conda activate barcode_analysis_gDNA   # for starcode on PATH
python run_barcode_pipeline.py run_configs/<your_run_name>/config.json
```

## Paths

- `fastq_dir`: where per-sample, decompressed R1 fastq files live (or will be
  staged to). Layout the pipeline expects: `<fastq_dir>/<sample>/<sample>_L001_R1_001.fastq`
  (the `_L001_R1_001` Illumina-style suffix matters -- the extraction code parses
  the sample name by splitting on `_L0`).
- `sc_in_dir`, `filt_dir`: scratch directories for intermediate per-read files.
  `filt_dir` is a required parameter of the underlying extraction function but
  unused downstream in this script -- just give it an empty directory.
- `out_dir`: final output location. Per-sample count tables land in
  `<out_dir>/starcode_outputs/sc_output_counts_<sample>.txt`.

## `stage_fastqs` (optional)

`{sample_name: path_to_raw_R1_fastq_or_fastq.gz}`. If present, the pipeline copies
(decompressing if `.gz`) each entry into the `fastq_dir` layout above before
extraction. **Omit this key entirely** if you've already staged fastqs yourself --
useful when your sequencer's demux output needs its own lane/index-picking logic
first (see e.g. `setup_20260821_new_samples_fastq.py` in this repo for an example
of that kind of one-off staging step; write a similar small script for a new run's
idiosyncratic raw layout, then point `fastq_dir` at its output and skip `stage_fastqs`).

## Extraction parameters

These are per-read settings, applied before any clustering:

- `strtseq`: the constant anchor sequence immediately before the barcode (reads
  are reverse-complemented first, so in search orientation the tag sits between
  the anchor and the barcode). Must match your construct's backbone sequence.
- `bclen`: barcode length in nt to slice out after the anchor (+ tag window, if
  present).
- `startseq_match`: % match required to call the anchor present; reads failing
  this are dropped.
- `tag_seqs`: `{tag_name: forward-orientation sequence immediately after the
  barcode}`. Leave `{}` if your construct has no tag field. A read's tag window
  is fuzzy-matched against all entries; reads with none confidently found keep a
  0nt offset (so untagged reads like ladder/spike-ins stay in register) and are
  labeled tag `"None"`.
- `tag_match`: % match required to call a tag.
- `include_tag_in_clustering`: `false` (default) strips the tag window before
  the barcode goes to starcode, so clustering never sees it (tag is still called
  and reported per-read separately). `true` carries the tag window through as
  part of what gets clustered -- only use this if you specifically want
  tag+barcode treated as one unit.
- `min_qual`: minimum mean Phred quality over the tag window and (separately)
  the barcode window. Low barcode-window quality drops the read entirely; low
  tag-window quality only downgrades the tag call to `"None"`. Set to `null` to
  disable.
- `exclude_tag_seqs`: `{tag_name: forward-orientation sequence}` for tags known
  to NOT belong to this experiment (e.g. cross-contamination from another
  library sharing sequencing capacity). Reads whose best tag match is one of
  these are dropped entirely, not just relabeled `"None"`. Leave `{}` if not
  applicable.

## Clustering parameters

- `sc_mm`: starcode mismatch tolerance (max 8). `4` has been the validated
  choice for this barcode library; re-derive this per construct (e.g. via a
  Hamming-distance-between-true-barcodes check) rather than assuming it carries
  over.
- `groups`: list of lists of sample names, e.g. `[["s1","s2"], ["s3"]]`. Samples
  within a group are pooled into one starcode run before being split back out
  per-sample -- use this when samples should share barcode-consensus calling
  (e.g. same lineage pool / same tag family). A group of one sample is fine if
  you don't want any pooling across samples. Every sample referenced here must
  have a corresponding entry (directly, or via `stage_fastqs`) in `fastq_dir`.

## Output format

`<out_dir>/starcode_outputs/sc_output_counts_<sample>.txt`, tab-separated:

```
barcode	tag	count
```

One row per (consensus barcode, tag) pair observed in that sample.
