"""
Generalized gDNA barcode fastq -> counts pipeline.

This replaces the one-off setup_<run>_fastq.py + extract_and_cluster_<run>_new_samples.py
script pairs (see e.g. extract_and_cluster_20260821_new_samples.py) with a single script
driven by a JSON config, so a new sequencing run/dataset means writing a new config, not a
new copy of the pipeline logic.

Pipeline (matches the validated/corrected eGRK1059 20260723+ approach -- see
load_corrected_clustering.py and the eGRK1059_lod_and_hamming_clustering memory):

  1. (optional) stage raw per-sample fastq.gz files into the uniform layout
     gDNA_starcode_prep expects: <fastq_dir>/<sample>/<sample>_L001_R<N>_001.fastq
     (decompressed -- gDNA_starcode_prep globs for *_R1*.fastq / *_R2*.fastq, not
     .fastq.gz). <N> is set for the whole run by the config's "barcode_read"
     ("R1" by default, or "R2" for a run where read orientation was swapped --
     a run is always uniformly one or the other, never mixed). Skip staging
     (omit "stage_fastqs" from the config) if you've already staged fastqs
     yourself, e.g. because your sequencer's demux layout needs run-specific
     lane/index logic beyond a flat sample->file mapping -- just name the files
     you stage with the matching _R1_/_R2_ tag.

  2. per-read extraction (gDNA_starcode_prep): reverse-complement the barcode-bearing
     read -- True by default for every sample, R1 or R2 alike (confirmed correct for
     a read-orientation-swapped library too; override via "reverse_complement" if a
     specific construct/read genuinely doesn't need it -- see
     run_configs/TEMPLATE/README.md), anchor on strtseq, call an optional tag window
     immediately after the anchor (conditional offset -- only skipped when a tag
     confidently matches, so tag-less reads like ladder/spike-ins stay in register),
     slice the barcode, apply startseq/quality filters, and optionally drop reads
     matching exclude_tag_seqs entirely.

  3. starcode clustering, run once per group in "groups" (pool samples that should share
     barcode-consensus calling -- e.g. same physical lineage pool / same tag family; a
     group of one is fine if you don't want any pooling).

  4. per-sample count tables: starcode's --seq-id output is split back out by sample and
     cross-tabulated against each read's tag call, written to
     <out_dir>/starcode_outputs/sc_output_counts_<sample>.txt as
     barcode\ttag\tcount

Requires `starcode` on PATH (e.g. `conda activate barcode_analysis_gDNA`).

Usage:
    python run_barcode_pipeline.py path/to/config.json

See run_configs/TEMPLATE/config.json for a documented config template.
"""
import argparse
import gzip
import json
import os
import shutil
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "BarcodeAnalysis_v2_gDNA", "ExtractBarcodes", "Scripts", "Functions"))
from gDNA_starcode_prep import gDNA_starcode_prep  # noqa: E402


def stage_fastqs(stage_map, fastq_dir, default_read="R1"):
    """Copy/decompress each sample's raw barcode-read fastq(.gz) into the uniform
    per-sample layout gDNA_starcode_prep expects.

    stage_map: {sample_name: path_to_source_fastq} for the common case, where
    every sample in this run uses the same read (set via default_read, itself
    driven by the config's top-level "barcode_read" -- a whole run is always
    uniformly R1 or R2, never mixed). A sample entry can instead be
    {"src": path_to_source_fastq, "read": "R2"} to override default_read for
    that one sample, if ever needed."""
    print("=== staging fastqs ===")
    for sample, entry in stage_map.items():
        if isinstance(entry, dict):
            src, read = entry["src"], entry.get("read", default_read)
        else:
            src, read = entry, default_read
        assert read in ("R1", "R2"), f"{sample}: read must be 'R1' or 'R2', got {read!r}"
        dest_dir = os.path.join(fastq_dir, sample)
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, f"{sample}_L001_{read}_001.fastq")
        if src.endswith(".gz"):
            with gzip.open(src, "rb") as f_in, open(dest, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copyfile(src, dest)
        size_mb = os.path.getsize(dest) / 1e6
        print(f"  {sample:<30} <- {src}  ({read}, {size_mb:,.1f} MB)")
    print("Done staging.\n")


def run_starcode_and_split(groups, sc_in_stage, sc_out, sc_mm):
    """For each group: concatenate its samples' per-read sc_input/tag_calls files,
    run starcode once on the pooled group, then split the --seq-id output back out
    per sample, cross-tabulated against each read's tag call."""
    starcode_cmd = shutil.which("starcode")
    if starcode_cmd is None:
        raise RuntimeError("starcode not found on PATH -- activate the conda env that has it "
                            "installed (e.g. `conda activate barcode_analysis_gDNA`) first.")

    print("=== pooling by group + running starcode ===")
    import subprocess

    for gi, group in enumerate(groups, start=1):
        group_sc_in = os.path.join(sc_out, f"sc_in_group{gi}_comb.txt")
        boundaries = {}
        idx = 1
        tag_lines_by_sample = {}
        with open(group_sc_in, "w") as out_sc:
            for sample in group:
                with open(os.path.join(sc_in_stage, f"sc_input_{sample}.txt")) as f:
                    lines = [l for l in f.read().split("\n") if l]
                with open(os.path.join(sc_in_stage, f"tag_calls_{sample}.txt")) as f:
                    tag_lines = [l for l in f.read().split("\n") if l]
                assert len(lines) == len(tag_lines), (
                    f"{sample}: {len(lines)} sc_input lines vs {len(tag_lines)} tag_calls lines")
                tag_lines_by_sample[sample] = tag_lines
                boundaries[sample] = (idx, idx + len(lines) - 1)
                out_sc.write("\n".join(lines) + "\n")
                idx += len(lines)

        out_path = os.path.join(sc_out, f"sc_output_group{gi}_comb.txt")
        print(f"  group{gi} ({len(group)} samples, {idx - 1:,} total reads): running starcode -d {sc_mm} ...")
        subprocess.run([starcode_cmd, "-d", str(sc_mm), "-t", "20", "-i", group_sc_in,
                         "-o", out_path, "--seq-id"], check=True)

        all_tag_lines = []
        for sample in group:
            all_tag_lines.extend(tag_lines_by_sample[sample])

        for sample in group:
            lo, hi = boundaries[sample]
            with open(os.path.join(sc_out, f"sc_output_counts_{sample}.txt"), "w") as outfile:
                outfile.write("barcode\ttag\tcount\n")
                with open(out_path) as f:
                    for line in f:
                        bc, cnt, ids = line.rstrip("\n").split("\t")
                        ids = [int(x) for x in ids.split(",")]
                        in_range = [r for r in ids if lo <= r <= hi]
                        if not in_range:
                            continue
                        tag_counts = defaultdict(int)
                        for r in in_range:
                            tag_counts[all_tag_lines[r - 1]] += 1
                        for tag in sorted(tag_counts):
                            outfile.write(f"{bc}\t{tag}\t{tag_counts[tag]}\n")
        print(f"  wrote per-sample count tables for: {group}")


def main(config_path):
    with open(config_path) as f:
        cfg = json.load(f)

    fastq_dir = cfg["fastq_dir"]
    sc_in_stage = cfg["sc_in_dir"]
    filt_dir = cfg["filt_dir"]  # required positional param of gDNA_starcode_prep, unused downstream
    out_dir = cfg["out_dir"]
    sc_out = os.path.join(out_dir, "starcode_outputs")
    for d in (fastq_dir, sc_in_stage, filt_dir, sc_out):
        os.makedirs(d, exist_ok=True)

    # barcode_read: which physical read carries the barcode for every sample in this
    # run (a run is always uniformly R1 or R2, never mixed -- see
    # run_configs/TEMPLATE/README.md). Only controls the destination filename used
    # when staging fastqs and which *_R1*/*_R2* file gets globbed -- both R1 and R2
    # reads get reverse-complemented before the strtseq search by default (see
    # "reverse_complement" below), so this has no effect on that transform.
    barcode_read = cfg.get("barcode_read", "R1")
    assert barcode_read in ("R1", "R2"), f"barcode_read must be 'R1' or 'R2', got {barcode_read!r}"

    if cfg.get("stage_fastqs"):
        stage_fastqs(cfg["stage_fastqs"], fastq_dir, default_read=barcode_read)

    groups = cfg["groups"]
    samples = [s for group in groups for s in group]
    assert len(samples) == len(set(samples)), "a sample name appears in more than one group"

    print("=== extracting reads ===")
    gDNA_starcode_prep(
        fastq_dir, sc_in_stage, groups,
        cfg["bclen"], filt_dir, cfg["strtseq"], cfg["startseq_match"], cfg["sc_mm"],
        cfg.get("tag_seqs", {}), cfg.get("tag_match", 75),
        cfg.get("include_tag_in_clustering", False),
        cfg.get("min_qual"), cfg.get("exclude_tag_seqs", {}),
        cfg.get("reverse_complement", True),
        cfg.get("assume_tag_present", False),
    )
    print("Done extracting.\n")

    run_starcode_and_split(groups, sc_in_stage, sc_out, cfg["sc_mm"])

    print(f"\nDone. Per-sample count tables in: {sc_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="path to a config JSON (see run_configs/TEMPLATE/config.json)")
    args = parser.parse_args()
    main(args.config)
