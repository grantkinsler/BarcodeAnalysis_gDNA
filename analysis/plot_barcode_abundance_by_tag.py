"""
Plot the distribution of barcode abundance for one or more sc_output_counts_<sample>.txt
files, split by each barcode's dominant tag call. One or more samples are rendered as
side-by-side small multiples sharing the same x/y axes and bins, for direct comparison.

Each distinct barcode is classified by whichever tag accounts for the most of its
reads (ties broken by dict order); its TOTAL read count (summed across all tag
calls for that barcode) determines which histogram bin it falls into.

Usage:
    python plot_barcode_abundance_by_tag.py --samples <counts_file>[:Label] [...] \
        --out <out_png> [--metric count|frequency] [--tags Tag1,Tag2,...] [--title-suffix "..."]

--samples: one or more counts files. Each may be "path" (sample name inferred
    from the sc_output_counts_<name>.txt filename) or "path:Label" to set the
    label explicitly. One file -> a single plot; multiple -> small multiples
    sharing bins/axes so panels are directly comparable.
--tags: which tag names to break out as their own category (stacking order,
    bottom to top); any barcode whose dominant call isn't in this list falls
    into a single catch-all "Other" category. Defaults to "Tag1".

--infer-cutoff: if set, treats --noise-reference-tags (e.g. known cross-
    contamination tags) as an independent estimate of the background rate at
    which ANY tag -- including the primary one (first entry in --tags) --
    gets spuriously called at a given abundance, purely from library/sequencing
    noise. Per abundance bin, compares the primary tag's distinct-barcode count
    to the mean per-tag count across the reference tags; the inferred cutoff is
    the upper edge of the highest-abundance bin where that reference rate is
    still nonzero. Below the cutoff, the primary tag is no more numerous than
    an independent contaminant tag at the same abundance, so those calls can't
    be trusted as real signal even though they're nominally "Tag1" (etc.) --
    above it, the reference tags disappear entirely and the primary tag stands
    alone. Drawn as a vertical dashed line per panel, annotated with the value
    and the count of primary-tag barcodes it flags as likely-still-contamination.
"""
import argparse
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

# fixed categorical palette, assigned by position in --tags (not cycled) -- Tag1
# (or whichever is listed first) gets the "hero" blue, "Other" is always neutral gray
PALETTE = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2"]
OTHER_COLOR = "#9E9E9E"


def load_and_classify(counts_file, tag_names):
    reads_by_bc_tag = defaultdict(lambda: defaultdict(int))
    with open(counts_file) as f:
        next(f)  # header
        for line in f:
            bc, tag, cnt = line.rstrip("\n").split("\t")
            reads_by_bc_tag[bc][tag] += int(cnt)
    barcodes = list(reads_by_bc_tag)
    totals = np.array([sum(reads_by_bc_tag[bc].values()) for bc in barcodes])
    dominant_tag = np.array([max(reads_by_bc_tag[bc], key=reads_by_bc_tag[bc].get) for bc in barcodes])
    return totals, dominant_tag


def infer_cutoff(values, dominant_tag, bins, primary_tag, reference_tags):
    """Upper edge of the highest-abundance bin where the mean per-reference-tag
    distinct-barcode count is still > 0 -- i.e. the point above which the
    reference (contamination) tags vanish entirely. Returns (cutoff, n_flagged),
    n_flagged = count of primary_tag barcodes at or below that cutoff."""
    last_contaminated_edge = None
    for lo, hi in zip(bins[:-1], bins[1:]):
        in_bin = (values >= lo) & (values < hi)
        ref_counts = [int(((dominant_tag == t) & in_bin).sum()) for t in reference_tags]
        if np.mean(ref_counts) > 0:
            last_contaminated_edge = hi
    if last_contaminated_edge is None:
        return None, 0
    n_flagged = int(((dominant_tag == primary_tag) & (values <= last_contaminated_edge)).sum())
    return last_contaminated_edge, n_flagged


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", required=True, help="counts file(s), each 'path' or 'path:Label'")
    parser.add_argument("--out", required=True, dest="out_png")
    parser.add_argument("--metric", choices=["count", "frequency"], default="count")
    parser.add_argument("--tags", default="Tag1", help="comma-separated tag names to break out individually")
    parser.add_argument("--title-suffix", default="")
    parser.add_argument("--infer-cutoff", action="store_true",
                         help="annotate the inferred real-signal cutoff (see --noise-reference-tags)")
    parser.add_argument("--noise-reference-tags", default="",
                         help="comma-separated tags used as the contamination-rate reference for --infer-cutoff")
    args = parser.parse_args()

    tag_names = [t.strip() for t in args.tags.split(",")]

    samples = []
    for spec in args.samples:
        path, _, label = spec.partition(":")
        if not label:
            label = path.split("sc_output_counts_")[-1].removesuffix(".txt")
        totals, dominant_tag = load_and_classify(path, tag_names)
        samples.append((label, totals, dominant_tag))

    # one shared metric + bin range across all samples, so panels are directly comparable
    all_totals = np.concatenate([s[1] for s in samples])
    all_total_reads_per_sample = [t.sum() for _, t, _ in samples]

    values_per_sample = []
    for _, totals, _ in samples:
        if args.metric == "frequency":
            values_per_sample.append(totals / totals.sum())
        else:
            values_per_sample.append(totals.astype(float))
    all_values = np.concatenate(values_per_sample)
    xlabel = "Barcode frequency within sample (log scale)" if args.metric == "frequency" else "Reads per barcode (log scale)"
    bins = np.logspace(np.floor(np.log10(all_values.min())), np.ceil(np.log10(all_values.max())), 50)

    n = len(samples)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 5.5), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)

    for ax, (label, totals, dominant_tag), values, total_reads in zip(
        axes, samples, values_per_sample, all_total_reads_per_sample
    ):
        n_barcodes = len(totals)
        groups, colors, labels, counts_per_group = [], [], [], []
        for i, tag in enumerate(tag_names):
            mask = dominant_tag == tag
            groups.append(values[mask])
            colors.append(PALETTE[i % len(PALETTE)])
            labels.append(tag)
            counts_per_group.append(int(mask.sum()))
        other_mask = ~np.isin(dominant_tag, tag_names)
        groups.append(values[other_mask])
        colors.append(OTHER_COLOR)
        labels.append("Other")
        counts_per_group.append(int(other_mask.sum()))

        ax.hist(groups, bins=bins, stacked=True, color=colors, label=labels, edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.spines[["top", "right"]].set_visible(False)

        breakdown = ", ".join(f"{c:,} {l}" for l, c in zip(labels, counts_per_group))
        title = f"{label}\n{n_barcodes:,} distinct barcodes, {total_reads:,} total reads\n({breakdown})"

        print(f"{label}: n_barcodes={n_barcodes}, total_reads={total_reads}")
        for l, c in zip(labels, counts_per_group):
            print(f"  {l}: {c}")

        if args.infer_cutoff and args.noise_reference_tags:
            reference_tags = [t.strip() for t in args.noise_reference_tags.split(",")]
            primary_tag = tag_names[0]
            cutoff, n_flagged = infer_cutoff(values, dominant_tag, bins, primary_tag, reference_tags)
            if cutoff is not None:
                n_real = counts_per_group[0] - n_flagged
                ax.axvline(cutoff, color="black", linewidth=1.3, linestyle="--", zorder=5)
                cutoff_str = f"{cutoff:.1e}" if args.metric == "frequency" else f"{cutoff:.0f}"
                title += (f"\ncutoff={cutoff_str}: {n_real:,} {primary_tag} likely real (above), "
                          f"{n_flagged:,} of {counts_per_group[0]:,} likely still contamination (at/below)")
                print(f"  inferred cutoff={cutoff_str} ({primary_tag} vs mean[{','.join(reference_tags)}]); "
                      f"{n_real}/{counts_per_group[0]} {primary_tag} barcodes likely real (above cutoff), "
                      f"{n_flagged}/{counts_per_group[0]} at/below it")

                ax.annotate(f"{n_flagged:,} likely\ncontamination", xy=(cutoff, 0.96), xycoords=("data", "axes fraction"),
                            ha="right", va="top", fontsize=8, color="#555555",
                            xytext=(-6, 0), textcoords="offset points")
                ax.annotate(f"{n_real:,} likely\nreal", xy=(cutoff, 0.96), xycoords=("data", "axes fraction"),
                            ha="left", va="top", fontsize=8, color="#1F77B4", fontweight="bold",
                            xytext=(6, 0), textcoords="offset points")

        ax.set_title(title, fontsize=9, wrap=True)

    axes[0].set_ylabel("Number of distinct barcodes (log scale)")
    handles, labels_ = axes[-1].get_legend_handles_labels()
    suptitle = args.title_suffix if args.title_suffix else None
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    if suptitle:
        fig.suptitle(suptitle, y=0.99, fontsize=10)
    fig.legend(handles, labels_, frameon=False, fontsize=9, loc="upper center",
               ncol=len(labels_), bbox_to_anchor=(0.5, 0.94 if suptitle else 0.99))
    fig.savefig(args.out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {args.out_png}")


if __name__ == "__main__":
    main()
