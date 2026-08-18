#!/usr/bin/env python3
"""Per-patient report: SNV findings, CNV findings, MSI/TMB summary, tool versions."""

import argparse
import csv
import datetime
import os
import sys
from collections import defaultdict


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--patient-id',    required=True)
    p.add_argument('--case-id',       required=True)
    p.add_argument('--germline-vcf',  required=True, help='VCF_PRIORITISE germline output (already curated: ClinVar P/LP + eviQ)')
    p.add_argument('--somatic-vcf',   required=True, help='VCF_PRIORITISE somatic output (already curated: OncoKB/COSMIC-tier1/hotspot/CNS5-whitelist)')
    p.add_argument('--cns',             required=True, help='CNVkit _call.cns file')
    p.add_argument('--cnr-genemetrics', required=True, help='CNVkit _cnr_genemetrics.txt file (CNR only)')
    p.add_argument('--msi',           required=True, help='MSIsensor output file')
    p.add_argument('--tmb',           required=True, help='pyTMB TSV output')
    p.add_argument('--markers',       required=True, help='Marker gene CSV (CNV section only - SNV section reads pre-curated VCF_PRIORITISE output instead)')
    p.add_argument('--chrom-arm',     required=True, help='chrom_arm.bed file')
    p.add_argument('--mane-db-name',  default='',    help='SnpEff MANE DB directory name')
    p.add_argument('--clinvar-name',  default='',    help='ClinVar VCF filename')
    p.add_argument('--vep-species',   default='homo_sapiens')
    p.add_argument('--vep-assembly',  default='GRCh38')
    p.add_argument('--vep-version',   default='')
    p.add_argument('--chrom-frac-thr', type=float, default=0.9,
                    help='Fraction of chromosome length that must be AMP/DEL for a chromosome-level call')
    p.add_argument('--arm-frac-thr',   type=float, default=0.9,
                    help='Fraction of arm length that must be AMP/DEL for a chromosome-arm-level call')
    p.add_argument('--log2-thr',       type=float, default=0.4,
                    help='Absolute log2 ratio threshold for a segment to count as AMP/DEL (chromosome, arm, and segment level)')
    p.add_argument('--out',           default='report.md')
    return p.parse_args()


# ── Markers ───────────────────────────────────────────────────────────────────

def load_markers(marker_csv):
    markers = set()
    with open(marker_csv) as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) >= 2 and row[1].strip():
                markers.add(row[1].strip())
    if not markers:
        sys.exit(f"ERROR: No markers loaded from {marker_csv}")
    return markers


# ── VCF parsing ───────────────────────────────────────────────────────────────

def _parse_info(info_str):
    d = {}
    for part in info_str.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            d[k] = v
        else:
            d[part] = True
    return d


def _extract_ann(ann_str):
    """Return (gene, transcript, hgvs_c) from the first ANN record."""
    first = ann_str.split(',')[0]
    parts = first.split('|')
    gene = parts[3]  if len(parts) > 3  else ''
    tx   = parts[6]  if len(parts) > 6  else ''
    cdna = parts[9]  if len(parts) > 9  else ''
    return gene, tx, cdna


def _iter_vcf_records(vcf_path):
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue

            chrom, pos, _, ref, alt = cols[:5]
            info = _parse_info(cols[7])

            try:
                af = float(info.get('AF', 0))
                dp = int(  info.get('DP', 0))
            except (ValueError, TypeError):
                af, dp = 0.0, 0

            gene, tx, cdna = _extract_ann(info.get('ANN', '')) if info.get('ANN') else ('', '', '')

            yield chrom, pos, ref, alt, info, af, dp, gene, tx, cdna


def parse_germline_vcf(vcf_path):
    """VCF_PRIORITISE's germline output is already curated (ClinVar P/LP,
    high-confidence review status, eviQ germline gene) - report every row."""
    out = []
    for chrom, pos, ref, alt, info, af, dp, gene, tx, cdna in _iter_vcf_records(vcf_path):
        out.append({
            'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
            'af': af, 'dp': dp,
            'gene': gene, 'tx': tx, 'cdna': cdna, 'clnsig': info.get('CLNSIG', ''),
        })
    return out


def parse_somatic_vcf(vcf_path):
    """VCF_PRIORITISE's somatic output is already curated (OncoKB-actionable +
    COSMIC CGC tier 1, COSMIC CMC significance tier, cancerhotspots, or a WHO
    CNS5 whitelist match) - report every row, not just marker-gene rows."""
    out = []
    for chrom, pos, ref, alt, info, af, dp, gene, tx, cdna in _iter_vcf_records(vcf_path):
        out.append({
            'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
            'status': info.get('STATUS', ''), 'af': af, 'dp': dp,
            'gene': gene, 'tx': tx, 'cdna': cdna,
        })
    return out


# ── CNV analysis ──────────────────────────────────────────────────────────────

VALID_CHROMS = {f"chr{i}" for i in range(1, 23)}


def load_gene_cn_from_cns(cns_path):
    """Return {gene: cn} by reading the cn column from _call.cns segments."""
    d = {}
    with open(cns_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            cn_val = row.get('cn', '')
            if cn_val == '':
                continue
            for g in [x.strip() for x in row.get('gene', '').split(',') if x.strip()]:
                if g not in d:
                    try:
                        d[g] = int(round(float(cn_val)))
                    except (ValueError, TypeError):
                        pass
    return d


def load_chrom_arm_info(arms_bed):
    """Return (chr_lengths, arms_by_chr) from a 4-column BED of chromosome arms."""
    chr_lengths = {}
    arms_by_chr = defaultdict(list)
    with open(arms_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            ch, s, e, label = fields[:4]
            if ch not in VALID_CHROMS:
                continue
            try:
                s, e = int(s), int(e)
            except ValueError:
                continue
            arm = label[0].lower()
            chr_lengths[ch] = max(chr_lengths.get(ch, 0), e)
            arms_by_chr[ch].append({'arm': arm, 'start': s, 'end': e})
    return chr_lengths, dict(arms_by_chr)


def gene_cnv_from_cnr_genemetrics(cnr_genemetrics_path, markers):
    """Return marker gene AMP/DEL rows from the CNR-only genemetrics file.
    The file is already filtered at |log2| >= 0.2 by cnvkit genemetrics."""
    rows = []
    with open(cnr_genemetrics_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            gene = row.get('gene', '').strip()
            if gene not in markers:
                continue
            ch = row.get('chromosome', '')
            if ch not in VALID_CHROMS:
                continue
            try:
                l2 = float(row['log2'])
                s  = int(row['start'])
                e  = int(row['end'])
            except (ValueError, KeyError):
                continue
            event = 'AMP' if l2 > 0 else 'DEL'
            rows.append({'chromosome': ch, 'start': s, 'end': e,
                         'gene': gene, 'event': event, 'log2': l2})
    rows.sort(key=lambda r: (_chr_order(r['chromosome']), r['start']))
    return rows


def _gene_arm(gene_row, arms_by_chr):
    ch  = gene_row['chromosome']
    mid = (gene_row['start'] + gene_row['end']) // 2
    for ai in arms_by_chr.get(ch, []):
        if ai['start'] <= mid < ai['end']:
            return ai['arm']
    return None


def _chr_order(chrom):
    try:
        return int(chrom.replace("chr", ""))
    except ValueError:
        return 10 ** 9


def chr_cnv_from_cns(cns_path, chr_lengths, frac_thr=0.9, log2_thr=0.4):
    per = defaultdict(lambda: {'amp': 0, 'del': 0})
    with open(cns_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ch = row.get('chromosome', '')
            if ch not in VALID_CHROMS:
                continue
            try:
                s, e, l2 = int(row['start']), int(row['end']), float(row['log2'])
            except (ValueError, KeyError):
                continue
            n = max(0, e - s)
            if l2 >= log2_thr:
                per[ch]['amp'] += n
            elif l2 <= -log2_thr:
                per[ch]['del'] += n

    calls = []
    for ch in sorted(per, key=_chr_order):
        total = chr_lengths.get(ch, 0)
        af = per[ch]['amp'] / total if total else 0
        df = per[ch]['del'] / total if total else 0
        call = ('AMP' if af >= frac_thr and af >= df
                else 'DEL' if df >= frac_thr and df > af
                else 'NEUTRAL')
        if call != 'NEUTRAL':
            calls.append({'chromosome': ch, 'call': call,
                          'amp_fraction': af, 'del_fraction': df})
    return calls


def arm_cnv_from_cns(cns_path, arms_bed, frac_thr=0.9, log2_thr=0.4):
    arms_by_chr = defaultdict(list)
    with open(arms_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            ch, s, e, label = fields[:4]
            if ch not in VALID_CHROMS:
                continue
            try:
                arms_by_chr[ch].append({'arm': label[0].lower(), 'start': int(s), 'end': int(e)})
            except ValueError:
                continue

    per = {}
    with open(cns_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ch = row.get('chromosome', '')
            if ch not in VALID_CHROMS or ch not in arms_by_chr:
                continue
            try:
                ss, se, l2 = int(row['start']), int(row['end']), float(row['log2'])
            except (ValueError, KeyError):
                continue
            for ai in arms_by_chr[ch]:
                ovl = max(0, min(se, ai['end']) - max(ss, ai['start']))
                if ovl == 0:
                    continue
                key = (ch, ai['arm'])
                if key not in per:
                    per[key] = {'total': 0, 'amp': 0, 'del': 0,
                                'arm_length': ai['end'] - ai['start']}
                per[key]['total'] += ovl
                if l2 >= log2_thr:
                    per[key]['amp'] += ovl
                elif l2 <= -log2_thr:
                    per[key]['del'] += ovl

    calls = []
    for (ch, arm), st in sorted(per.items(), key=lambda x: (_chr_order(x[0][0]), x[0][1])):
        al = st['arm_length']
        af = st['amp'] / al if al else 0
        df = st['del'] / al if al else 0
        call = ('AMP' if af >= frac_thr and af >= df
                else 'DEL' if df >= frac_thr and df > af
                else 'NEUTRAL')
        if call != 'NEUTRAL':
            calls.append({'chromosome': ch, 'arm': arm, 'call': call,
                          'amp_fraction': af, 'del_fraction': df,
                          'arm_length': al})
    return calls


def seg_cnv_from_cns(cns_path, log2_thr=0.4):
    rows = []
    with open(cns_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ch = row.get('chromosome', '')
            if ch not in VALID_CHROMS:
                continue
            try:
                s, e, l2 = int(row['start']), int(row['end']), float(row['log2'])
            except (ValueError, KeyError):
                continue
            if e <= s:
                continue
            if l2 >= log2_thr:
                grp = 'AMP'
            elif l2 <= -log2_thr:
                grp = 'DEL'
            else:
                continue
            rows.append({'chromosome': ch, 'start': s, 'end': e,
                         'gene': row.get('gene', ''), 'log2': l2, 'group': grp})

    rows.sort(key=lambda r: (
        _chr_order(r['chromosome']),
        0 if r['group'] == 'AMP' else 1,
        -r['log2'] if r['group'] == 'AMP' else r['log2'],
        r['start'],
    ))
    return rows


def gene_cnv_from_seg(seg_rows, markers):
    out = []
    for r in seg_rows:
        for g in [x.strip() for x in r['gene'].split(',') if x.strip()]:
            if g in markers:
                out.append({'chromosome': r['chromosome'], 'start': r['start'],
                            'end': r['end'], 'gene': g, 'event': r['group'], 'log2': r['log2']})
    out.sort(key=lambda r: (_chr_order(r['chromosome']), r['start']))
    return out


def filter_seg_by_markers(seg_rows, markers):
    """Keep segments that contain at least one marker gene. Show all genes; bold the markers."""
    out = []
    for r in seg_rows:
        all_genes = [g.strip() for g in r['gene'].split(',') if g.strip()]
        if not any(g in markers for g in all_genes):
            continue
        displayed = ', '.join(f"**{g}**" if g in markers else g for g in all_genes)
        out.append({**r, 'gene_display': displayed})
    return out


# ── MSI / TMB ─────────────────────────────────────────────────────────────────

def parse_msi(msi_path):
    try:
        with open(msi_path) as f:
            lines = [l.strip() for l in f if l.strip()]
        if len(lines) < 2:
            return {}
        header = lines[0].split('\t')
        values = lines[1].split('\t')
        d = dict(zip(header, values))
        pct = float(d.get('%', d.get('Percentage', 0)))
        return {
            'total_sites':    d.get('Total_Number_of_Sites', '?'),
            'unstable_sites': d.get('Number_of_Unstable_Sites', '?'),
            'pct':    pct,
            'status': 'MSI-H' if pct >= 3.5 else 'MSS',
        }
    except Exception as e:
        return {'error': str(e)}


def parse_tmb(tmb_path):
    """Parse pyTMB key=value report format."""
    try:
        d = {}
        with open(tmb_path) as f:
            for line in f:
                line = line.strip()
                if '=' in line and not line.startswith('#'):
                    k, v = line.split('=', 1)
                    d[k.strip()] = v.strip()
        return {
            'tmb':             d.get('TMB', '?'),
            'total_variants':  d.get('Total number of variants', '?'),
            'coding_variants': d.get('Variants after filters', '?'),
            'sample':          d.get('Sample', '?'),
        }
    except Exception as e:
        return {'error': str(e)}


# ── Markdown helpers ──────────────────────────────────────────────────────────

def _format_clinvar_date(filename):
    """Extract YYYYMMDD from clinvar filename and return DD/MM/YYYY."""
    import re
    m = re.search(r'(\d{8})', filename)
    if m:
        d = m.group(1)
        return f"{d[6:8]}/{d[4:6]}/{d[:4]}"
    return filename


def _wrap(s, width=25):
    if not s or len(s) <= width:
        return s
    return '<br>'.join(s[i:i+width] for i in range(0, len(s), width))


# ── Report writer ─────────────────────────────────────────────────────────────

def write_report(args, somatic, germline, chr_calls, arm_calls, gene_rows, seg_rows, msi, tmb, gene_cn, arms_by_chr):
    today = datetime.date.today().isoformat()

    with open(args.out, 'w') as f:

        # ── Header ──
        f.write(f"# Patient Report – {args.patient_id}\n\n")
        f.write(f"**Case ID:** {args.case_id}  \n")
        f.write(f"**Date:** {today}\n\n")
        f.write("---\n\n")

        # ── 1. SNV Findings ──
        f.write("## 1. SNV Findings\n\n")

        f.write("### Somatic mutations (actionable / hotspot / WHO CNS5 whitelist)\n\n")
        if somatic:
            f.write("| Gene | Location | Transcript | Variant (HGVS) | AF |\n")
            f.write("| --- | --- | --- | --- | --- |\n")
            for v in somatic:
                f.write(f"| {v['gene']} | {v['chrom']}:{v['pos']} | {v['tx']} "
                        f"| {_wrap(v['cdna'])} | {v['af']:.4f} |\n")
        else:
            f.write("_No actionable/hotspot/whitelisted somatic variants detected._\n")
        f.write("\n")

        f.write("### Germline pathogenic variants\n\n")
        if germline:
            f.write("| Gene | Location | Transcript | Variant (HGVS) | ClinVar |\n")
            f.write("| --- | --- | --- | --- | --- |\n")
            for v in germline:
                f.write(f"| {v['gene']} | {v['chrom']}:{v['pos']} | {v['tx']} "
                        f"| {_wrap(v['cdna'])} | {v['clnsig']} |\n")
        else:
            f.write("_No germline pathogenic/likely-pathogenic variants detected._\n")
        f.write("\n---\n\n")

        # ── 2. CNV Findings ──
        f.write("## 2. CNV Findings\n\n")

        chr_with_call = {c['chromosome'] for c in chr_calls}
        arm_calls_display = [a for a in arm_calls if a['chromosome'] not in chr_with_call]
        arm_with_call = {(a['chromosome'], a['arm']) for a in arm_calls_display}
        gene_rows_display = [
            g for g in gene_rows
            if g['chromosome'] not in chr_with_call
            and (g['chromosome'], _gene_arm(g, arms_by_chr)) not in arm_with_call
        ]

        f.write("### Chromosome-level events\n\n")
        if chr_calls:
            f.write("| Chromosome | Event | AMP (%) | DEL (%) |\n")
            f.write("| ---: | --- | ---: | ---: |\n")
            for c in chr_calls:
                f.write(f"| {c['chromosome']} | {c['call']} "
                        f"| {c['amp_fraction']*100:.1f} | {c['del_fraction']*100:.1f} |\n")
        else:
            f.write("_No whole-chromosome AMP/DEL events detected._\n")
        f.write("\n")

        f.write("### Chromosome-arm events\n\n")
        if arm_calls_display:
            f.write("| Chromosome | Arm | Event | AMP (%) | DEL (%) | Arm length (bp) |\n")
            f.write("| ---: | :---: | --- | ---: | ---: | ---: |\n")
            for a in arm_calls_display:
                f.write(f"| {a['chromosome']} | {a['arm']} | {a['call']} "
                        f"| {a['amp_fraction']*100:.1f} | {a['del_fraction']*100:.1f} | {a['arm_length']} |\n")
        else:
            f.write("_No chromosome-arm AMP/DEL events detected (or all arms covered by chromosome-level calls above)._\n")
        f.write("\n")

        f.write("### Gene-level events (marker genes*)\n\n")
        if gene_rows_display:
            f.write("| Chromosome | Start | End | Gene | Event | CN |\n")
            f.write("| ---: | ---: | ---: | --- | --- | ---: |\n")
            for g in gene_rows_display:
                cn_str = str(gene_cn.get(g['gene'], '?'))
                f.write(f"| {g['chromosome']} | {g['start']} | {g['end']} "
                        f"| {g['gene']} | {g['event']} | {cn_str} |\n")
            f.write(f"\n\\*Marker genes are genes found in WHO CNS5, OncoKB (Actionable and Cancer Genes), "
                    f"and all sources listed in `marker.csv`.\n")
        else:
            f.write("_No marker genes with AMP/DEL events detected._\n")
        f.write("\n---\n\n")

        # ── 3. MSI and TMB Summary ──
        f.write("## 3. MSI and TMB Summary\n\n")
        f.write("| Metric | Value |\n")
        f.write("| --- | --- |\n")
        if msi and 'error' not in msi:
            f.write(f"| MSI total sites | {msi.get('total_sites', '?')} |\n")
            f.write(f"| MSI unstable sites | {msi.get('unstable_sites', '?')} |\n")
            f.write(f"| MSI score (%) | {msi.get('pct', 0):.2f} |\n")
            f.write(f"| MSI status | **{msi.get('status', '?')}** |\n")
        else:
            f.write("| MSI | _N/A_ |\n")
        if tmb and 'error' not in tmb:
            f.write(f"| TMB (mut/Mb) | {tmb.get('tmb', tmb.get('TMB', '?'))} |\n")
            f.write(f"| Total variants | {tmb.get('total_variants', '?')} |\n")
            f.write(f"| Coding variants | {tmb.get('coding_variants', '?')} |\n")
        else:
            f.write("| TMB | _N/A_ |\n")
        f.write("\n---\n\n")

        # ── Appendix A: Segments with marker genes ──
        f.write("## Appendix A – Segments with marker genes\n\n")
        if seg_rows:
            for s in seg_rows:
                loc = f"{s['chromosome']}:{s['start']:,}–{s['end']:,}"
                f.write(f"**{loc}** &nbsp; {s['group']} &nbsp; log2 = {s['log2']:.4f}\n\n")
                f.write(f"{s['gene_display']}\n\n")
        else:
            f.write("_No AMP/DEL segments overlapping marker genes._\n")
        f.write("\n---\n\n")

        # ── Supplementary: Tool and annotation versions ──
        f.write("## Supplementary – Tool and Annotation Versions\n\n")
        f.write("| Tool | Version / Build |\n")
        f.write("| --- | --- |\n")
        if args.mane_db_name:
            f.write(f"| SnpEff MANE | `{args.mane_db_name}` |\n")
        if args.clinvar_name:
            f.write(f"| ClinVar access date | {_format_clinvar_date(args.clinvar_name)} |\n")
        if args.vep_version:
            f.write(f"| VEP | `v{args.vep_version} – {args.vep_species} {args.vep_assembly}` |\n")
        f.write("\n")

        # ── Manual comment box ──
        f.write("---\n\n")
        f.write("## Laboratory Comment\n\n")
        f.write("> _(Please enter your comments here)_\n\n")

        f.write("### Filtering parameters\n\n")
        f.write("_SNV germline/somatic curation (ClinVar P/LP + eviQ; OncoKB/COSMIC-tier1/hotspot/CNS5-whitelist) "
                "is applied upstream by VCF_PRIORITISE - see that module for the exact bcftools criteria._\n\n")
        f.write("| Parameter | Value |\n")
        f.write("| --- | --- |\n")
        f.write(f"| Chromosome AMP/DEL fraction threshold | {args.chrom_frac_thr} |\n")
        f.write(f"| Chromosome-arm AMP/DEL fraction threshold | {args.arm_frac_thr} |\n")
        f.write(f"| CNV log2 ratio threshold | {args.log2_thr} |\n")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    markers = load_markers(args.markers)

    somatic  = parse_somatic_vcf(args.somatic_vcf)
    germline = parse_germline_vcf(args.germline_vcf)

    chr_lengths, arms_by_chr = load_chrom_arm_info(args.chrom_arm)
    chr_calls    = chr_cnv_from_cns(args.cns, chr_lengths, frac_thr=args.chrom_frac_thr, log2_thr=args.log2_thr)
    arm_calls    = arm_cnv_from_cns(args.cns, args.chrom_arm, frac_thr=args.arm_frac_thr, log2_thr=args.log2_thr)
    seg_rows_all = seg_cnv_from_cns(args.cns, log2_thr=args.log2_thr)
    gene_rows    = gene_cnv_from_cnr_genemetrics(args.cnr_genemetrics, markers)
    seg_rows     = filter_seg_by_markers(seg_rows_all, markers)
    gene_cn      = load_gene_cn_from_cns(args.cns)

    msi = parse_msi(args.msi)
    tmb = parse_tmb(args.tmb)

    # Use basename of db paths for display in the version table
    args.mane_db_name = os.path.basename(args.mane_db_name) if args.mane_db_name else ''
    args.clinvar_name = os.path.basename(args.clinvar_name) if args.clinvar_name else ''

    write_report(args, somatic, germline, chr_calls, arm_calls, gene_rows, seg_rows, msi, tmb, gene_cn, arms_by_chr)
    print(f"Report written to {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
