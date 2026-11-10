# What's different from upstream (`SydShafferLab/BarcodeAnalysis_v2_gDNA`)

Verified by diffing `BarcodeAnalysis_v2_gDNA/` here directly against a fresh
clone of the upstream repo. Only two files differ in actual logic:
`ExtractBarcodes/Scripts/Functions/gDNA_starcode_prep.py` and
`ExtractBarcodes/Scripts/Step3.py`.

## Correctness bugs fixed (affect any run, even without using tags)

1. **Last-group read-count truncation** (`Step3.py`) — the final read-index
   boundary for each sample group was computed from `len(cat_fastq)`, a
   variable left over from whichever group the earlier "modifying fastq
   files" loop happened to process *last* — not necessarily the group
   currently being split out. Every group except that one got silently
   truncated or overextended, losing most of the last sample's reads in all
   but one group. Fixed to read each group's own combined starcode-input file
   length directly instead.
2. **Off-by-one dropped reads for multi-sample groups** (`Step3.py`) — a
   stray `+1` in the fastq-splitting index arithmetic dropped the first read
   of every sample after the first, for any group of ≥2 samples. (Upstream's
   own shipped example config only ever groups one sample per starcode run,
   so this never surfaced there.)
3. **WSN (homopolymer/N) filter never ran in gDNA mode** (`Step3.py`) — gated
   on `barcodeSource == 'cDNA'`, a typo; should have been `'gDNA'`. In
   gDNA-only mode (the only mode this fork supports) the filter silently
   never executed as shipped upstream.
4. **starcode PATH handling** (`Step3.py`) — upstream appends the vendored
   `starcode/` source directory to `PATH` and assumes a prebuilt Linux binary
   lives there. That binary doesn't run on macOS. Fixed to require `starcode`
   already resolvable on `PATH` (e.g. via a conda env with it installed from
   bioconda) and fail loudly with a clear message if it isn't.

## Features added (for constructs with a tag field between barcode and anchor)

The naive upstream pipeline's only notion of read structure is
**anchor sequence → barcode**. Our barcode constructs have a short tag field
immediately after the barcode (in construct orientation), so
`gDNA_starcode_prep` gained:

- **Tag calling** (`tag_seqs`, `tagMatch`) — fuzzy-matches a tag window right
  after the anchor (in the reverse-complemented search orientation this sits
  between the anchor and the barcode) and reports it per read. Upstream's
  output is `barcode\tcount`; this fork's is `barcode\ttag\tcount`.
- **Conditional barcode offset** — the barcode-start offset only skips past
  the tag window when a tag *confidently matches*; reads with no tag at all
  (e.g. ladder/spike-in controls with no tag field) keep a 0nt offset instead
  of a fixed skip that would misalign and corrupt their barcode. This is the
  key fix that lets tag-less control reads survive extraction at all.
- **Tag exclusion** (`exclude_tag_seqs`) — reads whose tag window
  best-matches a *known-foreign* tag (e.g. cross-contamination from another
  multiplexed library sharing sequencing capacity) are dropped entirely
  rather than merely mislabeled.
- **Per-base quality gating** (`minQual`) — upstream never looks at Phred
  quality, only the anchor's fuzzy-match percentage. Added a mean-quality
  threshold on the tag window (downgrades the tag call to `"None"`, doesn't
  drop the read — we're just no longer confident *which* tag) and separately
  on the barcode window (drops the read entirely, since an unreliable
  barcode call shouldn't enter starcode clustering).
- **Per-read anchor position** — upstream forces every read in a file to a
  single global modal stagger position (`statistics.mode` over all reads).
  This fork finds each read's own exact anchor position first, falling back
  to the file mode only when a read has no exact match — needed because a
  secondary read population (e.g. a ladder/spike-in construct) can have a
  different, internally-consistent stagger of its own that the global mode
  would otherwise misalign.
- **`include_tag_in_clustering` option** — normally the tag window is
  stripped before the barcode goes to starcode (so clustering never sees it,
  though the tag is still called and reported separately). Setting this to
  `True` instead carries the tag window through as part of what gets
  clustered, for cases where you want tag+barcode treated as one unit.
- **R2 / read-orientation-swap support** (`barcode_read`, `reverse_complement`)
  — upstream only ever reads `*_R1*.fastq`. This fork also globs
  `*_R2*.fastq`, driven by `run_barcode_pipeline.py`'s top-level
  `barcode_read` config field (`"R1"`/`"R2"`) -- which `_R1_`/`_R2_`-tagged
  file gets staged/read for every sample in a run (always uniformly one or
  the other, never mixed within one config). Reads are reverse-complemented
  before the `strtseq` search by default regardless of `barcode_read` --
  upstream's original, only-ever-R1 behavior carries over unchanged to R2 --
  controlled by `reverse_complement` (bool, default `True`, or
  `{sample: bool}` for an override) if a specific construct/read genuinely
  needs something else. `strtseq`/`tag_seqs` describe one fixed canonical
  orientation reached after that flip, so tag matching (always against the
  reverse complement of a `tag_seqs` value) doesn't change with
  `reverse_complement` — confirmed against synthetic R1 and R2 reads built
  that way (see git history for the test).
- **`assume_tag_present` option** — by default, the tag_len-nt skip (whether
  the barcode window starts after the tag) only happens when a read's tag
  window confidently fuzzy-matches a declared tag; a read that fails that
  match keeps a 0nt offset, correct for a construct where some reads
  genuinely have no tag there (e.g. a tag-less ladder/spike-in). For a
  single-tag library where the tag is always truly physically present, that
  default is actively harmful: a read whose real tag has enough sequencing
  error to miss the match threshold doesn't just get mislabeled, it silently
  reports `[the noisy tag][first (bclen-tag_len) nt of the real barcode]` as
  a fake barcode instead of the true one — a different fabricated barcode per
  affected read, which inflated an observed distinct-barcode count by ~20x
  before being diagnosed and fixed. `assume_tag_present=True` skips the tag
  window unconditionally regardless of match confidence, for exactly this
  single-tag-library case (the tag LABEL reported is unaffected, still driven
  by the fuzzy match, just no longer gates the offset).

- **Gzipped fastq input** — upstream only globs for and opens plain
  `*.fastq`. `gDNA_starcode_prep` now also finds `*.fastq.gz` (R1/R2/I1/I2)
  and decompresses on the fly for both read passes, so per-sample fastqs can
  be symlinks to the original compressed files instead of ~20 GB of
  decompressed copies. Verified byte-identical `sc_input_*`/`tag_calls_*`
  output for the same reads given as `.fastq` vs `.fastq.gz`. `filt_dir`
  outputs are always written as plain text (`.gz` stripped from the name).

## New on top of both

`run_barcode_pipeline.py` is a new, generalized driver (not part of upstream
at all). It skips the parts of `Step3.py` that exist only to feed 10x
CellRanger custom feature-barcoding (`FeatureReference.csv` generation,
modified-fastq rewriting, CellRanger input directory layout) — none of that
is relevant if you just want flat per-sample count tables — and instead calls
extraction (`gDNA_starcode_prep`) and starcode clustering directly, driven by
a single JSON config per run rather than a hand-edited copy of the driver
script. `Step3.py` is retained for anyone who does need the CellRanger output
format.
