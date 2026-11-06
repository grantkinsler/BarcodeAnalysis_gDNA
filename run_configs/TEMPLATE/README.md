# Config fields for run_barcode_pipeline.py

Copy this folder to `run_configs/<your_run_name>/`, edit `config.json` (plain JSON,
no comments), then:

```
conda activate barcode_analysis_gDNA   # for starcode on PATH
python run_barcode_pipeline.py run_configs/<your_run_name>/config.json
```

## Paths

- `fastq_dir`: where per-sample, decompressed barcode-read fastq files live (or
  will be staged to). Layout the pipeline expects:
  `<fastq_dir>/<sample>/<sample>_L001_R<N>_001.fastq`, `<N>` being 1 or 2
  depending on which read carries the barcode for that sample (see
  `reverse_complement` below) -- the `_L001_R<N>_001` Illumina-style suffix
  matters, the extraction code parses the sample name by splitting on `_L0`.
  Each sample folder should contain exactly one such file.
- `sc_in_dir`, `filt_dir`: scratch directories for intermediate per-read files.
  `filt_dir` is a required parameter of the underlying extraction function but
  unused downstream in this script -- just give it an empty directory.
- `out_dir`: final output location. Per-sample count tables land in
  `<out_dir>/starcode_outputs/sc_output_counts_<sample>.txt`.

## `stage_fastqs` (optional)

If present, the pipeline copies (decompressing if `.gz`) each entry into the
`fastq_dir` layout above before extraction. Two forms per sample -- a whole
run is either all-R1 or all-R2 (never mixed), so you'll use one form for
every sample in a given config:

- `"sample_name": "/path/to/raw_R1.fastq.gz"` -- the common case (barcode
  read is R1).
- `"sample_name": {"src": "/path/to/raw_R2.fastq.gz", "read": "R2"}` -- for a
  run where read orientation was swapped, so the barcode-bearing read is R2
  instead (see `reverse_complement` below). `"read"` must be `"R1"` or `"R2"`.

**Omit this key entirely** if you've already staged fastqs yourself -- useful
when your sequencer's demux output needs its own lane/index-picking logic
first (see e.g. `setup_20260821_new_samples_fastq.py` in this repo for an
example of that kind of one-off staging step; write a similar small script
for a new run's idiosyncratic raw layout, then point `fastq_dir` at its
output -- using the `_R1_`/`_R2_` naming above for whichever read actually
carries the barcode -- and skip `stage_fastqs`).

## Extraction parameters

These are per-read settings, applied before any clustering:

- `strtseq`: the constant anchor sequence immediately before the barcode, in
  the read's orientation *after* the `reverse_complement` transform below is
  applied (so in the normal R1/reverse-complement=true case, the tag sits
  between the anchor and the barcode in that post-transform orientation).
  Must match your construct's backbone sequence.
- `bclen`: barcode length in nt to slice out after the anchor (+ tag window, if
  present).
- `reverse_complement` (optional, omit for the default): controls whether the
  raw read is reverse-complemented before searching for `strtseq`. Default
  (key omitted, or `null` for a given sample) **auto-detects per sample from
  its fastq filename**: an R1 file gets reverse-complemented, an R2 file
  doesn't -- this matches both the normal R1 setup for this construct and the
  confirmed behavior for a read-orientation-swapped R2 sample.

  `strtseq` and `tag_seqs` describe one fixed "canonical" orientation --
  the one you get *after* this reverse-complement step (if any) is applied --
  and that canonical form doesn't change based on which physical read it came
  from or whether getting there needed a flip. Concretely: `strtseq` is
  matched directly against the (possibly reverse-complemented) read, and a
  tag is matched against the reverse complement of its `tag_seqs` value,
  regardless of `reverse_complement`'s value for that sample. So for a
  swapped R2 sample needing no flip, its raw read should already read
  directly as `strtseq + revcomp(tag) + barcode` -- the exact same string an
  R1 read produces only after being reverse-complemented. This has been
  verified against synthetic reads built that way; **spot-check your first
  real swapped-orientation run** (a suspiciously low read-pass rate, or tag
  calls coming back mostly `"None"`, both mean this assumption doesn't hold
  for your construct -- come back and we'll adjust the tag-matching logic to
  match what your real reads actually look like).

  Leave `reverse_complement` out entirely for the normal case -- the
  filename-based default already does the right thing uniformly across a
  whole run (all-R1 or all-R2; samples are never mixed within one run). Only
  set it explicitly (a single bool, applied to every sample) if a specific
  construct/swap needs something other than that default.
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
