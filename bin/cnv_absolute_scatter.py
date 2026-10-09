#!/usr/bin/env python3
"""Absolute-copy-number scatter + BAF panel for one region (a whole
chromosome, or a chrom:start-end window), reading a CNVkit .cnr/.cns pair
directly - no dependency on cnvkit's own (log2-only) `scatter` command or its
internal cnvlib package. Renders one PNG per invocation, driven by one row
of generate_report.py's scatter manifest (CREATE_CNVKIT_SCATTER module).

Per-bin and per-segment log2 is converted to absolute copy number using the
same formula as CNVkit's own `cnvlib/call.py:_log2_ratio_to_absolute` (for
autosomes, ref_copies == expect_copies == ploidy):
    CN = (ploidy * 2**log2 - ploidy * (1 - purity)) / purity
A segment's own `cn`/`cn1`/`cn2` columns (already computed by the upstream
`cnvkit.py call --purity --ploidy` step) are used as-is for segment lines -
NOT re-derived from the segment's own log2, which would double-apply the
purity/ploidy rescaling already baked into that log2 by `call`.
"""
import argparse
import csv
import gzip
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

POINT_COLOR = "#606060"
NEUTRAL_COLOR = "#A0A0A0"      # not a reportable event - background/baseline
EVENT_COLOR = "#d62728"        # this plot's own reported/target event - red
OTHER_EVENT_COLOR = "darkorange"  # a different reportable event incidentally in view
GENE_COLOR = "gold"
CENTROMERE_COLOR = "#8B4513"

LABEL_STAGGER_LEVELS = 8  # see _stagger_y


def absolute_cn(log2, ploidy, purity):
    return (ploidy * 2.0 ** log2 - ploidy * (1 - purity)) / purity


def _opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_region(region, cnr_path):
    """'chr7' (whole chromosome) or 'chr3:100000-200000'. Whole-chromosome
    bounds come from the .cnr's own bins (its last bin's end on that
    chromosome), so no extra reference file is needed."""
    if ":" in region:
        chrom, coords = region.split(":", 1)
        start, end = (int(x) for x in coords.split("-", 1))
        return chrom, start, end
    chrom = region
    end = 0
    with _opener(cnr_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["chromosome"] == chrom:
                end = max(end, int(row["end"]))
    if end == 0:
        sys.exit(f"ERROR: chromosome {chrom!r} not found in {cnr_path}")
    return chrom, 0, end


def load_cnr(path, chrom, start, end, ploidy, purity):
    pos, cn = [], []
    with _opener(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["chromosome"] != chrom:
                continue
            s, e = int(row["start"]), int(row["end"])
            if e < start or s > end:
                continue
            pos.append(0.5 * (s + e) / 1e6)
            cn.append(absolute_cn(float(row["log2"]), ploidy, purity))
    return np.array(pos), np.array(cn)


def load_cns(path, chrom, start, end):
    segs = []
    with _opener(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["chromosome"] != chrom:
                continue
            s, e = int(row["start"]), int(row["end"])
            if e < start or s > end:
                continue
            cn1_raw, cn2_raw = row.get("cn1"), row.get("cn2")
            segs.append({
                "start": max(s, start) / 1e6, "end": min(e, end) / 1e6,
                "cn": int(round(float(row["cn"]))),
                "cn1": int(round(float(cn1_raw))) if cn1_raw not in (None, "") else None,
                "cn2": int(round(float(cn2_raw))) if cn2_raw not in (None, "") else None,
                "gene": row.get("gene", ""),
            })
    return segs


def load_vcf_vaf(path, chrom, start, end, sample_id, normal_id,
                  zygosity_freq, min_depth):
    """Tumour VAF at SNVs whose normal-sample AF looks heterozygous - same
    het-site selection --zygosity-freq drives in cnvkit's own scatter
    --vcf, read directly from the FORMAT/AF field."""
    pos, vaf = [], []
    with _opener(path) as f:
        header = None
        tumour_col = normal_col = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                tumour_col = header.index(sample_id) if sample_id in header else None
                normal_col = header.index(normal_id) if normal_id in header else None
                continue
            if tumour_col is None or normal_col is None:
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] != chrom:
                continue
            p = int(cols[1])
            if not (start <= p <= end):
                continue
            if len(cols[3]) != 1 or len(cols[4]) != 1:
                continue  # SNVs only
            fmt = cols[8].split(":")
            try:
                dp_i, af_i = fmt.index("DP"), fmt.index("AF")
            except ValueError:
                continue
            try:
                t_fields = cols[tumour_col].split(":")
                n_fields = cols[normal_col].split(":")
                t_dp, t_af = int(t_fields[dp_i]), float(t_fields[af_i])
                n_af = float(n_fields[af_i])
            except (ValueError, IndexError):
                continue
            if t_dp < min_depth:
                continue
            if not (zygosity_freq <= n_af <= 1 - zygosity_freq):
                continue
            pos.append(p / 1e6)
            vaf.append(t_af)
    return np.array(pos), np.array(vaf)


def load_centromere(bed_path, chrom):
    """p-arm end == q-arm start, from the 4-column chrom_arm.bed
    (chrom, start, end, 'p'/'q'). None if the chromosome/bed is unavailable
    (e.g. chrX/chrY, or no bed given) - acrocentric-p chromosomes still have
    a real centromere position here even though ACROCENTRIC_P_ARMS excludes
    their p-arm from CNV calling elsewhere in the report."""
    if not bed_path:
        return None
    with open(bed_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 4 and fields[0] == chrom and fields[3].strip().lower() == "p":
                return int(fields[2])
    return None


def _stagger_y(index, y_top, y_bottom, step_frac, cycle=LABEL_STAGGER_LEVELS):
    """Vertical position for the index-th label in a group sharing one plot:
    1st at the top, each subsequent one a bit lower, wrapping back to the
    top every `cycle` labels (default 8) so a busy plot doesn't run out of
    headroom - avoids a dense cluster of same-height labels overlapping."""
    level = index % cycle
    return y_top - level * step_frac * (y_top - y_bottom)


def find_gene_spans(cnr_path, chrom, start, end, genes):
    if not genes:
        return []
    wanted = set(genes)
    spans = {}
    with _opener(cnr_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["chromosome"] != chrom:
                continue
            s, e = int(row["start"]), int(row["end"])
            if e < start or s > end:
                continue
            for gene in row["gene"].split(","):
                gene = gene.strip()
                if gene not in wanted:
                    continue
                lo, hi = spans.get(gene, (s, e))
                spans[gene] = (min(lo, s), max(hi, e))
    return [(lo / 1e6, hi / 1e6, name) for name, (lo, hi) in sorted(
        spans.items(), key=lambda kv: kv[1][0])]


def highlight_genes(ax, gene_spans, y_frac=0.78, step_frac=0.07):
    if not gene_spans:
        return
    ylo, yhi = ax.get_ylim()
    y_top = ylo + y_frac * (yhi - ylo)
    n = len(gene_spans)
    rot = "horizontal" if n <= 3 else (30 if n <= 6 else (45 if n <= 10 else (60 if n <= 20 else "vertical")))
    size = "small" if n <= 6 else "x-small"
    for i, (gstart, gend, name) in enumerate(gene_spans):
        y_posn = _stagger_y(i, y_top, ylo, step_frac)
        ax.axvspan(gstart, gend, alpha=0.35, color=GENE_COLOR, zorder=-1)
        ax.text(0.5 * (gstart + gend), y_posn, name, ha="center", rotation=rot, size=size)


def _is_event(cn, cn1, cn2, ploidy):
    """A segment worth calling out on a plot: genuine CN-LOH (always,
    regardless of ploidy), or a total CN that differs from ploidy and isn't
    itself a plain neutral CN=2 - same rule as generate_report.py's
    _is_reportable_segment with no arm baseline (see its module docstring)."""
    if cn == 2 and cn1 == 2 and cn2 == 0:
        return True
    if cn == ploidy:
        return False
    return cn != 2


def segment_color(seg, ploidy, target_start=None, target_end=None):
    """green: this is the specific segment the plot was generated for (focal
    plots), or any reportable event on a whole-chromosome/arm plot (no single
    target there). orange: a DIFFERENT reportable event that's merely
    incidentally captured in a focal plot's view window. grey: not a
    reportable event at all - background/baseline."""
    is_event = _is_event(seg["cn"], seg["cn1"], seg["cn2"], ploidy)
    if target_start is not None:
        is_target = abs(seg["start"] - target_start) < 1e-6 and abs(seg["end"] - target_end) < 1e-6
        if is_target:
            return EVENT_COLOR
        return OTHER_EVENT_COLOR if is_event else NEUTRAL_COLOR
    return EVENT_COLOR if is_event else NEUTRAL_COLOR


def segment_vaf_lines(vaf_pos, vaf_val, seg_start, seg_end):
    mask = (vaf_pos >= seg_start) & (vaf_pos <= seg_end)
    vals = vaf_val[mask]
    lines = []
    for sel in (vals > 0.5, vals <= 0.5):
        if sel.sum() > 1:
            lines.append(np.median(vals[sel]))
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cnr", required=True)
    ap.add_argument("--cns", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--normal-id", required=True)
    ap.add_argument("--region", required=True, help="'chr7' (whole chromosome) or 'chr3:100000-200000'")
    ap.add_argument("--purity", type=float, required=True)
    ap.add_argument("--ploidy", type=int, required=True)
    ap.add_argument("--genes", default="", help="Comma-separated gene symbols to highlight (blank: none)")
    ap.add_argument("--zygosity-freq", type=float, default=0.35)
    ap.add_argument("--min-variant-depth", type=int, default=30)
    ap.add_argument("--y-min", type=float, default=-1.0)
    ap.add_argument("--y-max", type=float, default=10.0,
                     help="Fixed y-axis ceiling - never expanded. A segment whose CN exceeds it is "
                          "drawn clipped at the ceiling with its true CN labelled above the plot")
    ap.add_argument("--chrom-arm-bed", default=None,
                     help="chrom_arm.bed (params.chrom_arm_bed) - draws a centromere marker. Optional.")
    ap.add_argument("--target-start", type=int, default=None,
                     help="Exact (unpadded) start of the segment this plot was generated for - "
                          "omitted for whole-chromosome/arm plots, which have no single target")
    ap.add_argument("--target-end", type=int, default=None)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    chrom, start, end = parse_region(args.region, args.cnr)
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    target_start_mb = args.target_start / 1e6 if args.target_start is not None else None
    target_end_mb = args.target_end / 1e6 if args.target_end is not None else None

    bin_pos, bin_cn = load_cnr(args.cnr, chrom, start, end, args.ploidy, args.purity)
    segs = load_cns(args.cns, chrom, start, end)
    vaf_pos, vaf_val = load_vcf_vaf(
        args.vcf, chrom, start, end, args.sample_id, args.normal_id,
        args.zygosity_freq, args.min_variant_depth)
    gene_spans = find_gene_spans(args.cnr, chrom, start, end, genes)
    centromere = load_centromere(args.chrom_arm_bed, chrom)
    centromere_mb = centromere / 1e6 if centromere is not None and start <= centromere <= end else None

    print(f"{len(bin_pos)} bins, {len(segs)} segments, {len(vaf_pos)} informative SNVs, "
          f"{len(gene_spans)} genes in {chrom}:{start:,}-{end:,}", file=sys.stderr)

    y_min, y_max = args.y_min, args.y_max

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 2], "hspace": 0.08})

    # ---- top: absolute copy number ----
    ax1.scatter(bin_pos, bin_cn, s=10, alpha=0.35, color=POINT_COLOR, edgecolors="none")
    ax1.axhline(args.ploidy, color="k", lw=0.8, ls=":", alpha=0.6, zorder=-2)
    if centromere_mb is not None:
        ax1.axvline(centromere_mb, color=CENTROMERE_COLOR, lw=1.2, ls="--", alpha=0.7, zorder=-2)
    off_scale_i = 0
    for seg in segs:
        c = segment_color(seg, args.ploidy, target_start_mb, target_end_mb)
        if seg["cn"] > y_max:
            # Off-scale: draw clipped at the ceiling, marked with a triangle
            # at the edge and its true CN labelled just inside the plot.
            # Labels stagger downward (wrapping every 8) so several off-scale
            # segments in one plot don't collide with each other.
            mid = 0.5 * (seg["start"] + seg["end"])
            label_y = _stagger_y(off_scale_i, y_max - 0.04 * (y_max - y_min), y_min, 0.05)
            off_scale_i += 1
            ax1.plot([seg["start"], seg["end"]], [y_max] * 2, color=c, lw=8,
                      solid_capstyle="round", zorder=3)
            ax1.plot(mid, y_max, marker="^", color=c, markersize=9, zorder=4, clip_on=False)
            ax1.text(mid, label_y, f"CN={seg['cn']}",
                      ha="center", va="top", fontsize=9, fontweight="bold", color=c, zorder=4)
        else:
            ax1.plot([seg["start"], seg["end"]], [seg["cn"]] * 2, color=c, lw=8,
                      solid_capstyle="round", zorder=3)
    ax1.set_ylim(y_min, y_max)
    tick_step = max(1, round((y_max - y_min) / 10))
    ax1.yaxis.set_major_locator(MultipleLocator(tick_step))
    ax1.set_ylabel("Absolute copy number")
    highlight_genes(ax1, gene_spans)
    title = f"{args.sample_id}  {chrom}:{start:,}-{end:,}  (purity={args.purity}, ploidy={args.ploidy})"
    ax1.set_title(title, fontsize=11)

    # ---- bottom: BAF ----
    ax2.scatter(vaf_pos, vaf_val, s=14, alpha=0.35, color=POINT_COLOR, edgecolors="none")
    ax2.axhline(0.5, color="k", lw=0.8, alpha=0.5, zorder=-2)
    if centromere_mb is not None:
        ax2.axvline(centromere_mb, color=CENTROMERE_COLOR, lw=1.2, ls="--", alpha=0.7, zorder=-2)
    for seg in segs:
        c = segment_color(seg, args.ploidy, target_start_mb, target_end_mb)
        for v in segment_vaf_lines(vaf_pos, vaf_val, seg["start"], seg["end"]):
            ax2.plot([seg["start"], seg["end"]], [v, v], color=c, lw=5,
                      solid_capstyle="round", zorder=3)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("VAF")
    ax2.set_xlabel("Position (Mb)")

    fig.tight_layout()
    fig.savefig(args.output, dpi=130)
    print(f"wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
