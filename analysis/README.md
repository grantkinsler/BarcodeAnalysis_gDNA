# Post-processing analysis: tag-contamination diagnosis and real-signal cutoff

These operate on already-extracted `sc_output_counts_<sample>.txt` files (the
pipeline's output), not raw fastqs. They don't belong in the extraction
pipeline itself -- see the design rationale below.

## Why this exists

If a `sc_output_counts_*.txt` file has way more distinct barcodes than you
expected, two unrelated things can cause it, and they need different fixes:

1. **A pipeline bug**, not a data problem -- see `assume_tag_present` in the
   main README/CHANGES.md. Symptom: barcodes at tag `"None"` that share a
   long common prefix (shorter real diversity than `bclen` nt). Fixed in the
   extraction step itself, not here.
2. **Real cross-contamination** from another tagged library sharing the same
   sequencing run (documented precedent: eGRK1059's Tag3-6). This is a
   genuine property of the data, not a bug -- the tool here is for
   confirming it's real (vs. a 3rd possibility, sequencing/PCR error) and
   then deciding how much of even the nominally-correct-tag population is
   still unreliable because of it.

## Step 1: is there a secondary tag signal at all?

Tally the raw 8(or however long your tag is)-nt window immediately after the
anchor (`strtseq`) across a few hundred thousand raw reads, independent of
the pipeline/tag_seqs config -- i.e. don't trust the pipeline's own tag
calls for this, go back to the fastq directly:

```python
import gzip
from Bio.Seq import Seq
from collections import Counter

STRTSEQ = "..."       # your strtseq
TAG_LEN = 8            # your tag length

window_counts = Counter()
n_checked = 0
with gzip.open("sample_R1_or_R2.fastq.gz", "rt") as f:
    for i, line in enumerate(f):
        if i % 4 != 1:
            continue
        rc = str(Seq(line.strip()).reverse_complement())  # match your reverse_complement setting
        idx = rc.find(STRTSEQ)
        if idx == -1:
            continue
        window = rc[idx+len(STRTSEQ):idx+len(STRTSEQ)+TAG_LEN]
        if len(window) == TAG_LEN:
            window_counts[window] += 1
            n_checked += 1
        if n_checked >= 300000:
            break

for seq, cnt in window_counts.most_common(20):
    print(seq, cnt, f"{100*cnt/n_checked:.2f}%")
```

A single declared tag should dominate (one sharp peak). If a *second* sharp,
high-frequency peak shows up at a sequence you didn't declare, go to Step 2
before assuming anything about it.

## Step 2: contamination, or sequencing/PCR error? Four checks + one control

Don't just assume a secondary peak is contamination -- verify it isn't
explainable by error, which would need a different fix (e.g. a looser/
tighter `tag_match`, not `exclude_tag_seqs`):

1. **Edit distance** from your real tag(s) to the candidate. Distance 1-2
   (reachable by a couple of substitution errors) is consistent with error.
   Distance 4+ essentially rules out single-read sequencing error producing
   it at any meaningful frequency (a specific k-substitution pattern on one
   read has probability ~(per-base error rate)^k).
2. **Base quality** in that window for reads matching the candidate, vs. for
   reads matching your real tag. Should be statistically indistinguishable
   if it's a genuine, correctly-sequenced different sequence; degraded if
   it's noise.
3. **Read-level mixing within individual barcodes**: for barcodes whose
   dominant tag is your real one, what fraction have *any* reads under the
   candidate tag? If it's real contamination (separate molecules), this
   should be low (contamination lives in its own barcodes, not smeared
   across your real population). If it's widespread at low level across
   most of your real barcodes, that's more consistent with per-read error.
4. **Cross-sample barcode identity**, if you have >1 sample from the same
   sequencing run: do the *specific barcodes* under the candidate tag
   overlap between samples? Near-zero overlap despite the *same* candidate
   tag identity recurring in both (classic index-hopping signature: same
   foreign tag, different random draw of that foreign library's barcodes
   each time) is strong evidence of real contamination from an external
   source. Your real tag's barcodes should show *meaningful* overlap across
   replicate/related samples, for contrast.
5. **Random-sequence negative control**: generate a few random 8-mers at the
   *same* Hamming distance from your real tag as the candidate, and tally
   them the same way (Step 1). If they come back at ~0% while the candidate
   is at several percent, that rules out "any sufficiently-distant sequence
   would show up this way" as an explanation -- it's specifically that
   sequence, not the distance, driving the signal.

If multiple candidates check out, add them all to `exclude_tag_seqs` in your
pipeline config and re-run extraction (this drops the contaminating reads
entirely before clustering, which is the structurally correct fix -- don't
try to filter them out after the fact in these barcode-count files, since
by then they've already potentially merged into starcode clusters with real
reads at `sc_mm` distance).

## Step 3: even within the real tag, where's the real-signal cutoff?

Contamination tags are also useful as an independent reference for "what does
a spurious/background barcode's abundance look like" -- even among barcodes
confidently called your real tag, the lowest-abundance ones may be
statistically indistinguishable from that background rate.

Re-run extraction with **every** candidate tag (real + contamination)
declared in `tag_seqs` and `exclude_tag_seqs` left empty, so every read gets
labeled by whichever it best matches instead of being dropped. Then:

```bash
python analysis/plot_barcode_abundance_by_tag.py \
    --samples sc_output_counts_<sample>.txt:<label> \
    --out out.png \
    --metric frequency \
    --tags <YourRealTag> \
    --infer-cutoff --noise-reference-tags <Contam1>,<Contam2>,...
```

Per abundance bin, this compares your real tag's distinct-barcode count to
the *mean* count of the contamination tags in that same bin, and reports the
abundance above which contamination drops to a sustained zero -- barcodes at
or below that cutoff should be treated as unreliable even though they're
nominally your real tag. Accepts multiple `--samples` for side-by-side small
multiples with shared bins/axes.

## `plot_pairwise_frequency.py`

Pairwise frequency scatter between two samples' *corrected* (contamination-
excluded) counts, colored by the Step 3 real/contamination classification
(needs both samples' corrected AND all-tags counts files). Useful as a sanity
check: the "likely real in both" population should show a real correlation
structure; "likely contamination in both" should sit at the axis floor with
no cross-sample relationship.

```bash
python analysis/plot_pairwise_frequency.py \
    --corrected sampleA_corrected.txt:LabelA sampleB_corrected.txt:LabelB \
    --alltags sampleA_alltags.txt sampleB_alltags.txt \
    --primary-tag <YourRealTag> --noise-reference-tags <Contam1>,<Contam2>,... \
    --out out.png
```

## Why this lives here, not in the pipeline

The extraction pipeline's job is mechanical: given reads and known tag
sequences, produce accurate counts (`exclude_tag_seqs`/`assume_tag_present`
belong there because they change what counts as "the barcode" or whether a
read is clustered at all). Whether to *trust* an already-correctly-counted
low-abundance barcode as real biology is an analysis-stage judgment -- same
category as a limit-of-detection threshold -- that needs the full
all-tags breakdown (which a routine clean run shouldn't have to produce) and
a user-chosen reference-tag set and threshold that's worth being able to
re-examine and adjust without re-running extraction.
