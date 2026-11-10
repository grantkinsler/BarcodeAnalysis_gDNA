"""
Pairwise barcode-frequency scatter between two samples, following the house
convention from eGRK1059's plot_pairwise_correlation.py: log-log, a barcode
absent from one sample is placed at a shared floor (half the global min
nonzero frequency) rather than dropped, a y=x reference line, and a Pearson r
on log-transformed values for barcodes present in both.

Frequencies come from each sample's Tag3-6-excluded ("corrected") counts file
(the real experimental signal, not diluted by cross-contamination reads).
Points are colored by the real-vs-contamination classification from the
all-tags diagnostic run (see plot_barcode_abundance_by_tag.py --infer-cutoff):
a barcode's raw Tag1 read count in the ALL-TAGS run is compared to that
sample's inferred absolute-count cutoff (independent of which frequency
denominator a downstream plot happens to use).

Usage:
    python plot_pairwise_frequency.py \
        --corrected p1_corrected.txt:1P_p1_Init_Rep2 p2_corrected.txt:1P_p2_Init \
        --alltags p1_alltags.txt p2_alltags.txt \
        --primary-tag Tag1 --noise-reference-tags Tag3,Tag4,Tag5,Tag6 \
        --out out.png
"""
import argparse
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr


def load_corrected_freqs(path):
    barcode_total = defaultdict(int)
    with open(path) as f:
        next(f)
        for line in f:
            bc, tag, cnt = line.rstrip("\n").split("\t")
            barcode_total[bc] += int(cnt)
    total = sum(barcode_total.values())
    return {bc: cnt / total for bc, cnt in barcode_total.items()}


def load_alltags_dominant_counts(path, primary_tag):
    """Returns {barcode: raw read count of its dominant tag call}, and the
    count-space cutoff (see infer_cutoff in plot_barcode_abundance_by_tag.py)."""
    reads_by_bc_tag = defaultdict(lambda: defaultdict(int))
    with open(path) as f:
        next(f)
        for line in f:
            bc, tag, cnt = line.rstrip("\n").split("\t")
            reads_by_bc_tag[bc][tag] += int(cnt)
    return reads_by_bc_tag


def infer_count_cutoff(reads_by_bc_tag, primary_tag, reference_tags):
    barcodes = list(reads_by_bc_tag)
    totals = np.array([sum(reads_by_bc_tag[bc].values()) for bc in barcodes])
    dominant = np.array([max(reads_by_bc_tag[bc], key=reads_by_bc_tag[bc].get) for bc in barcodes])
    bins = np.logspace(0, np.ceil(np.log10(totals.max())), 50)
    last_contaminated_edge = None
    for lo, hi in zip(bins[:-1], bins[1:]):
        in_bin = (totals >= lo) & (totals < hi)
        ref_counts = [int(((dominant == t) & in_bin).sum()) for t in reference_tags]
        if np.mean(ref_counts) > 0:
            last_contaminated_edge = hi
    return last_contaminated_edge


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--corrected", nargs=2, required=True, metavar="FILE:LABEL",
                         help="two 'path:Label' corrected (contamination-excluded) counts files")
    parser.add_argument("--alltags", nargs=2, required=True, metavar="FILE",
                         help="the two matching all-tags (nothing excluded) counts files, same order as --corrected")
    parser.add_argument("--primary-tag", default="Tag1")
    parser.add_argument("--noise-reference-tags", default="Tag3,Tag4,Tag5,Tag6")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    reference_tags = [t.strip() for t in args.noise_reference_tags.split(",")]

    specs = [s.partition(":") for s in args.corrected]
    (path_a, _, label_a), (path_b, _, label_b) = specs
    freq_a = load_corrected_freqs(path_a)
    freq_b = load_corrected_freqs(path_b)

    rbt_a = load_alltags_dominant_counts(args.alltags[0], args.primary_tag)
    rbt_b = load_alltags_dominant_counts(args.alltags[1], args.primary_tag)
    cutoff_a = infer_count_cutoff(rbt_a, args.primary_tag, reference_tags)
    cutoff_b = infer_count_cutoff(rbt_b, args.primary_tag, reference_tags)

    def is_real(bc, reads_by_bc_tag, cutoff):
        if bc not in reads_by_bc_tag:
            return False
        total = sum(reads_by_bc_tag[bc].values())
        return cutoff is None or total > cutoff

    all_bc = sorted(set(freq_a) | set(freq_b))
    all_vals = [v for v in list(freq_a.values()) + list(freq_b.values())]
    floor = min(all_vals) / 2
    lims = [floor * 0.5, 1.0]

    xs = np.array([freq_a.get(bc, floor) for bc in all_bc])
    ys = np.array([freq_b.get(bc, floor) for bc in all_bc])
    real_a = np.array([is_real(bc, rbt_a, cutoff_a) for bc in all_bc])
    real_b = np.array([is_real(bc, rbt_b, cutoff_b) for bc in all_bc])

    both_real = real_a & real_b
    one_real = (real_a | real_b) & ~both_real
    neither_real = ~real_a & ~real_b

    fig, ax = plt.subplots(figsize=(6.5, 7.2))
    ax.scatter(xs[neither_real], ys[neither_real], s=8, color="#9E9E9E", alpha=0.35, linewidths=0,
               zorder=2, label=f"likely contamination in both (n={neither_real.sum():,})")
    ax.scatter(xs[one_real], ys[one_real], s=10, color="#FF7F0E", alpha=0.6, linewidths=0,
               zorder=3, label=f"likely real in only one (n={one_real.sum():,})")
    ax.scatter(xs[both_real], ys[both_real], s=14, color="#1F77B4", alpha=0.8, linewidths=0,
               zorder=4, label=f"likely real in both (n={both_real.sum():,})")
    ax.plot(lims, lims, color="0.6", linestyle="--", linewidth=1, zorder=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel(f"{label_a} barcode frequency (log scale)")
    ax.set_ylabel(f"{label_b} barcode frequency (log scale)")
    ax.spines[["top", "right"]].set_visible(False)

    both_present = np.array([bc in freq_a and bc in freq_b for bc in all_bc])
    if both_present.sum() > 1:
        r, _ = pearsonr(np.log10(xs[both_present]), np.log10(ys[both_present]))
        ax.text(0.05, 0.95, f"r={r:.2f} (all barcodes present in both, log-log)",
                transform=ax.transAxes, fontsize=9, va="top")
    if both_real.sum() > 1:
        r_real, _ = pearsonr(np.log10(xs[both_real]), np.log10(ys[both_real]))
        ax.text(0.05, 0.90, f"r={r_real:.2f} (likely-real-in-both only)",
                transform=ax.transAxes, fontsize=9, va="top", color="#1F77B4")

    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    ax.set_title(f"{label_a} vs {label_b}: pairwise barcode frequency\n"
                 f"(Tag3-6 excluded; colored by real-signal cutoff\nfrom all-tags contamination analysis)",
                 fontsize=10)

    fig.tight_layout()
    fig.savefig(args.out, dpi=170)
    plt.close(fig)
    print(f"Saved {args.out}")
    print(f"cutoff ({label_a}): {cutoff_a} reads; cutoff ({label_b}): {cutoff_b} reads")
    print(f"n_total_barcodes={len(all_bc)}, both_real={both_real.sum()}, one_real={one_real.sum()}, neither_real={neither_real.sum()}")


if __name__ == "__main__":
    main()
