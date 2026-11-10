# Config fields for run_barcode_pipeline.py

Copy this folder to `run_configs/<your_run_name>/`, edit `config.json` (plain JSON,
no comments), then:

```
conda activate barcode_analysis_gDNA   # for starcode on PATH
python run_barcode_pipeline.py run_configs/<your_run_name>/config.json
```

## Paths

- `fastq_dir`: where per-sample barcode-read fastq files live (or will be
  staged to). Plain `.fastq` or gzipped `.fastq.gz` both work (gzipped files
  are decompressed on the fly, so symlinks to the original `.fastq.gz` files
  are enough -- no need to decompress to disk). Layout the pipeline expects:
  `<fastq_dir>/<sample>/<sample>_L001_R<N>_001.fastq[.gz]`, `<N>` set for the whole
  run by `barcode_read` below -- the `_L001_R<N>_001` Illumina-style suffix
  matters, the extraction code parses the sample name by splitting on `_L0`.
  Each sample folder should contain exactly one such file.
- `sc_in_dir`, `filt_dir`: scratch directories for intermediate per-read files.
  `filt_dir` is a required parameter of the underlying extraction function but
  unused downstream in this script -- just give it an empty directory.
- `out_dir`: final output location. Per-sample count tables land in
  `<out_dir>/starcode_outputs/sc_output_counts_<sample>.txt`.

## `barcode_read` (optional, omit for the default)

`"R1"` (default) or `"R2"`: which physical read carries the barcode, for
**every sample in this run** -- a run is always uniformly R1 or R2, never
mixed within one config. Set this to `"R2"` for a run where read orientation
was deliberately swapped in library prep. It only controls the `_R1_`/`_R2_`
tag in the filenames `stage_fastqs` writes (see below) and which file gets
read -- **it does not change the reverse-complement step**: both R1 and R2
reads get reverse-complemented before the `strtseq` search by default (see
`reverse_complement` under Extraction parameters).

If you've already staged fastqs yourself (skipping `stage_fastqs`), just name
your files with the matching `_R1_`/`_R2_` tag and set `barcode_read` to match.

## `stage_fastqs` (optional)

If present, the pipeline copies (decompressing if `.gz`) each entry into the
`fastq_dir` layout above before extraction, using the read set by
`barcode_read`:

```json
"stage_fastqs": {
  "sample1": "/path/to/raw_data/sample1_R1.fastq.gz",
  "sample2": "/path/to/raw_data/sample2_R1.fastq.gz"
}
```

(An entry can instead be `{"src": "...", "read": "R2"}` to override
`barcode_read` for one sample specifically, if that's ever genuinely needed --
not the normal case, since a run's samples are never mixed.)

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
  applied (so the tag sits between the anchor and the barcode in that
  post-transform orientation). Must match your construct's backbone sequence.
- `bclen`: barcode length in nt to slice out after the anchor (+ tag window, if
  present).
- `reverse_complement` (optional, omit for the default): controls whether the
  raw read is reverse-complemented before searching for `strtseq`. **Defaults
  to `true` for every sample, R1 or R2 alike** -- the raw read, whichever
  physical read it is, needs flipping before `strtseq -> tag -> barcode`
  appears in that order. `barcode_read` does *not* change this default (it
  only picks which file gets read); only set `reverse_complement` to `false`
  if a specific construct/read genuinely reads directly into
  `strtseq -> tag -> barcode` without needing a flip -- either a single bool
  applied to every sample, or a per-sample object, e.g.
  `{"some_sample": false}`, with `true` as the default for any sample not
  listed.

  `strtseq` and `tag_seqs` describe one fixed "canonical" orientation -- the
  one you get *after* this reverse-complement step is applied -- so `strtseq`
  is matched directly against the (reverse-complemented) read, and a tag is
  matched against the reverse complement of its `tag_seqs` value.
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
- `assume_tag_present` (optional, `false` by default): `false` means the
  tag_len-nt skip (whether the barcode window starts after the tag) only
  happens when a read's tag window confidently fuzzy-matches a declared tag
  -- correct when some reads genuinely have no tag at that position (e.g. a
  tag-less ladder/spike-in mixed in with tagged samples). Set to `true` when
  every real read in this sample is known to carry a tag at that fixed
  position, with no genuinely tag-less population to protect (e.g. a
  single-tag library). **This matters a lot**: with the default `false`, a
  read whose tag is genuinely present but has enough sequencing error to miss
  `tag_match` falls back to a 0nt offset instead of the correct tag_len-nt
  one -- this doesn't drop the read, it silently reports
  `[the noisy tag][first (bclen-tag_len) nt of the real barcode]` as a fake
  "barcode" instead of the true one. Because the noise pattern differs read
  to read, this fabricates a different spurious barcode per affected read
  rather than merging them, which can inflate the apparent distinct-barcode
  count by an order of magnitude or more for a construct where the tag is
  always truly present. A telltale sign you need `true`: a cluster of
  "barcodes" in your output that are shorter than `bclen` nt of real
  diversity (i.e. share a long common prefix) at tag `"None"`.
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
