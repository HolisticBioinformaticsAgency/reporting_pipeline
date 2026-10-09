#!/usr/bin/env python3
"""Per-patient report: SNV findings, CNV findings, MSI/TMB summary, tool versions."""

import argparse
import csv
import datetime
import json
import os
import sys
from collections import defaultdict

# Tumour-purity estimates below this are withheld from the report (shown as
# "not reported"); findings should then be read with extra caution.
PURITY_MIN_REPORTABLE = 0.2
DEFAULT_PLOIDY = 2  # fallback when FACETS ploidy is missing, "NA", or unparseable


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--patient-id',    required=True)
    p.add_argument('--case-id',       required=True)
    p.add_argument('--germline-vcf',  required=True, help='VCF_PRIORITISE germline output (already curated: ClinVar P/LP + eviQ)')
    p.add_argument('--somatic-vcf',   required=True, help='VCF_PRIORITISE somatic output (already curated: whitelist / marker gene / COSMIC-tier1-oncogene, CMC-tier / cancerhotspots / HIGH-MODERATE consequence)')
    p.add_argument('--hotspot-counts', default=None, help='WHITELIST module TSV - tumour/normal REF/ALT read counts at whitelisted hotspot positions, pulled straight from the raw alignment (BAM). Optional: omitted when no whitelist is configured (params.whitelist unset).')
    p.add_argument('--facets-purity', default=None, help='FACETS tumour-purity text file - a single line: a fraction in (0, 1], or the literal "NA". Optional: missing / "NA" renders as "not available".')
    p.add_argument('--facets-ploidy', default=None, help='FACETS ploidy text file - a single line: a positive number, or the literal "NA". Rounded and used as the CNV section\'s baseline "neutral" CN. Optional: missing / "NA" / unparseable falls back to the historical assumption of 2.')
    p.add_argument('--whitelist', default=None, help='Same whitelist CSV VCF_PRIORITISE used (params.whitelist), or the assets/NO_WHITELIST placeholder - recovers the "Whitelist match" reason for the somatic table\'s Selection reason column. Optional.')
    p.add_argument('--patient-metadata', default=None, help='Per-patient metadata TSV (params.patient_metadata), keyed by case_id in MOLBI###### form. Fills the report front matter; an absent file or unmatched patient leaves those fields blank.')
    p.add_argument('--cns',             required=True, help='CNVkit _call.cns file')
    p.add_argument('--msi',           required=True, help='MSIsensor output file')
    p.add_argument('--tmb',           required=True, help='pyTMB TSV output')
    p.add_argument('--tumour-hla',    required=True, help='Polysolver tumour sample _hla.txt file')
    p.add_argument('--germline-hla',  required=True, help='Polysolver germline/normal sample _hla.txt file')
    p.add_argument('--markers',       required=True, help='Marker gene CSV (CNV section only - SNV section reads pre-curated VCF_PRIORITISE output instead)')
    p.add_argument('--chrom-arm',     required=True, help='chrom_arm.bed file')
    p.add_argument('--supplementary-json', default=None,
                    help='Shared cohort-wide manifest of tool/database versions (two top-level '
                         'keys, "Tools" and "Annotation Databases", each name -> version string) '
                         'for the Supplementary section. Optional: missing file or absent key '
                         'renders that table as not provided.')
    p.add_argument('--pct-alt-thr',   type=float, default=90.0,
                    help='Chromosome/arm %% length AMP, DEL, or CN-LOH required to call that chromosome/arm')
    p.add_argument('--amp-seg-cn',    type=int, default=3,
                    help='CN at or above which a segment counts as amplified')
    p.add_argument('--del-seg-cn',    type=int, default=1,
                    help='CN at or below which a segment counts as deleted')
    p.add_argument('--out',              default='report.md')
    p.add_argument('--scatter-manifest', default='scatter_manifest.tsv',
                    help='Output TSV listing chromosomes/segments that need a CNVkit scatter plot')
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


def load_whitelist_keys(whitelist_csv):
    """{(CHROM, POS, REF, ALT)} exact-match keys, mirroring VCF_PRIORITISE's
    own WL_EXPR construction. None when no whitelist was configured (the
    assets/NO_WHITELIST placeholder, or path not given)."""
    if not whitelist_csv or os.path.basename(whitelist_csv) == 'NO_WHITELIST':
        return None
    keys = set()
    with open(whitelist_csv) as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) >= 4 and row[0].startswith('chr'):
                keys.add((row[0], row[1], row[2], row[3]))
    return keys


def load_supplementary(path):
    """{"Tools": {name: version}, "Annotation Databases": {name: version}} -
    a manually-maintained, cohort-wide manifest (see supplementary.json): the
    tool/database versions behind the report's inputs can't be recovered
    reliably from filenames or container tags alone (many annotation DBs
    carry no version in their filename at all, and a container tag is
    sometimes a build hash, not a semantic version), so this is asserted by
    whoever ran the pipeline rather than inferred. Missing file, invalid
    JSON, or an absent top-level key all degrade to an empty dict for that
    table, not an error - the report then states that section wasn't
    provided."""
    empty = {"Tools": {}, "Annotation Databases": {}}
    if not path or not os.path.isfile(path):
        return empty
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return empty
    return {
        "Tools": data.get("Tools") or {},
        "Annotation Databases": data.get("Annotation Databases") or {},
    }


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
    """Return (gene, transcript, hgvs_c, hgvs_p, consequence) from the first
    ANN record. `consequence` is SnpEff's own Annotation field (e.g.
    missense_variant, stop_gained, upstream_gene_variant) - as-is, including
    any '&'-joined multiple consequences on one variant."""
    first = ann_str.split(',')[0]
    parts = first.split('|')
    consequence = parts[1]  if len(parts) > 1  else ''
    gene        = parts[3]  if len(parts) > 3  else ''
    tx          = parts[6]  if len(parts) > 6  else ''
    cdna        = parts[9]  if len(parts) > 9  else ''
    prot        = parts[10] if len(parts) > 10 else ''
    return gene, tx, cdna, prot, consequence


def _ann_all_genes(ann_str):
    """Every Gene_Name across ALL ANN blocks (not just the first) - same
    multi-gene scan VCF_PRIORITISE's gate B uses, needed here because a
    variant's marker-gene match can come from a non-first block (see
    vcf_prioritise.nf's own header comment for the EPM2AIP1/MLH1 example)."""
    genes = []
    for block in ann_str.split(','):
        parts = block.split('|')
        if len(parts) > 3 and parts[3] and parts[3] not in genes:
            genes.append(parts[3])
    return genes


# Same putative HIGH/MODERATE-impact vocabulary as VCF_PRIORITISE's
# CSQ_EVIDENCE regex (modules/vcf_prioritise.nf) - kept in sync by hand.
_HIGH_MODERATE_CONSEQUENCES = {
    'bidirectional_gene_fusion', 'chromosome_number_variation', 'conservative_inframe_deletion',
    'conservative_inframe_insertion', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion',
    'duplication', 'exon_loss_variant', 'exon_loss', 'feature_ablation', 'frameshift_variant',
    'gene_fusion', 'inversion', 'missense_variant', 'protein_protein_contact',
    'rare_amino_acid_variant', 'rearranged_at_DNA_level', 'sequence_feature&exon_loss_variant',
    'splice_acceptor_variant', 'splice_branch_variant_U12', 'splice_donor_variant', 'start_lost',
    'stop_gained', 'stop_lost', 'structural_interaction_variant', 'transcript_ablation',
    '3_prime_UTR_truncation&exon_loss', '5_prime_UTR_truncation&exon_loss_variant',
}


def _is_high_moderate_ann(ann_str):
    for block in ann_str.split(','):
        parts = block.split('|')
        if len(parts) > 1 and parts[1] in _HIGH_MODERATE_CONSEQUENCES:
            return True
    return False


def selection_reason_somatic(chrom, pos, ref, alt, info, ann_str, markers, whitelist_keys):
    """Human-readable reason(s) this somatic variant passed VCF_PRIORITISE
    (mirrors its gate A / gate B exactly - see modules/vcf_prioritise.nf).
    Lists every qualifying piece of evidence, not just the minimum sufficient
    one, since a variant can clear the filter multiple ways at once."""
    reasons = []
    if whitelist_keys is not None and (chrom, pos, ref, alt) in whitelist_keys:
        reasons.append('Whitelist match')

    all_genes = _ann_all_genes(ann_str) if ann_str else []
    marker_hits = [g for g in all_genes if g in markers]
    if marker_hits:
        reasons.append(f"Marker gene ({', '.join(marker_hits)})")

    onc_tsg = info.get('COSMIC_ONC_TSG', '')
    tier = info.get('COSMIC_CGC_TIER', '')
    if tier == '1' and 'oncogene' in onc_tsg.split(','):
        reasons.append('COSMIC CGC tier 1 oncogene')

    cmc_tier = info.get('COSMIC_MUTATION_SIGNIFICANCE_TIER', '.')
    if cmc_tier not in ('', '.'):
        reasons.append(f'COSMIC CMC tier {cmc_tier}')

    hs_versions = [v for v in ('1', '2', '3') if info.get(f'cancerhotspots_v{v}') == 'T']
    if hs_versions:
        reasons.append(f"cancerhotspots (v{'/'.join(hs_versions)})")

    if ann_str and _is_high_moderate_ann(ann_str):
        reasons.append('HIGH/MODERATE-impact consequence')

    return '; '.join(reasons) if reasons else 'Unknown'


_CLNREVSTAT_LABELS = {
    'reviewed_by_expert_panel': 'reviewed by expert panel',
    'criteria_provided,_multiple_submitters,_no_conflicts': 'multiple submitters, no conflicts',
}


def selection_reason_germline(info):
    """VCF_PRIORITISE's germline gate is one fixed AND (ClinVar P/LP with a
    high-confidence review status, on the eviQ germline gene list) - every
    row satisfies the same structure, so this just surfaces which specific
    ClinVar significance/review status applied for that row."""
    clnsig = info.get('CLNSIG', '').replace('_', ' ')
    revstat_raw = info.get('CLNREVSTAT', '')
    revstat = _CLNREVSTAT_LABELS.get(revstat_raw, revstat_raw.replace('_', ' '))
    parts = [f'ClinVar {clnsig}' + (f' ({revstat})' if revstat else '')]
    if info.get('eviQ_germline') == 'T':
        parts.append('eviQ germline gene')
    return '; '.join(parts)


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

            gene, tx, cdna, prot, consequence = (
                _extract_ann(info.get('ANN', '')) if info.get('ANN') else ('', '', '', '', ''))

            yield chrom, pos, ref, alt, info, af, dp, gene, tx, cdna, prot, consequence


def parse_germline_vcf(vcf_path):
    """VCF_PRIORITISE's germline output is already curated (ClinVar P/LP,
    high-confidence review status, eviQ germline gene) - report every row."""
    out = []
    for chrom, pos, ref, alt, info, af, dp, gene, tx, cdna, prot, consequence in _iter_vcf_records(vcf_path):
        out.append({
            'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
            'af': af, 'dp': dp,
            'gene': gene, 'tx': tx, 'cdna': cdna, 'prot': prot, 'consequence': consequence,
            'clnsig': info.get('CLNSIG', ''),
            'reason': selection_reason_germline(info),
        })
    return out


def parse_somatic_vcf(vcf_path, markers, whitelist_keys):
    """VCF_PRIORITISE's somatic output is already curated (gate A: whitelist
    match, OR gate B: marker gene / COSMIC CGC tier-1 oncogene AND a COSMIC
    CMC tier, cancerhotspots match, or HIGH/MODERATE consequence - see
    modules/vcf_prioritise.nf) - report every row, not just marker-gene rows.
    `markers`/`whitelist_keys` are only used to recover *why* each row
    already passed, for the report's "Selection reason" column."""
    out = []
    for chrom, pos, ref, alt, info, af, dp, gene, tx, cdna, prot, consequence in _iter_vcf_records(vcf_path):
        out.append({
            'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
            'status': info.get('STATUS', ''), 'af': af, 'dp': dp,
            'gene': gene, 'tx': tx, 'cdna': cdna, 'prot': prot, 'consequence': consequence,
            'reason': selection_reason_somatic(chrom, pos, ref, alt, info, info.get('ANN', ''),
                                                markers, whitelist_keys),
        })
    return out


# ── CNV analysis ──────────────────────────────────────────────────────────────
#
# Ploidy (from FACETS, rounded - see parse_facets_ploidy) is the report's
# baseline "neutral" CN, replacing the old fixed assumption of 2 everywhere
# below except in the CN-LOH pattern itself (CN-LOH is a fixed biological/
# CNVkit-calling definition - cn == 2, cn1 == 2, cn2 == 0 - independent of
# tumour ploidy, and is never suppressed by a ploidy match; see point 3).
#
# Segment classification (per segment, from `cn`/`cn1`/`cn2`):
#   CN-LOH : cn == 2 and cn1 == 2 and cn2 == 0   (copy-neutral loss of heterozygosity)
#   GAIN   : cn >= --amp-seg-cn
#   LOSS   : cn <= --del-seg-cn
#   (checked in that order - CN-LOH first, since it's itself a cn == 2 segment
#   that would otherwise fall through both GAIN/LOSS absolute-threshold checks)
#
# Chromosome/arm calling is ARM-FIRST:
#   1. For each arm (skipping acrocentric p-arms - see ACROCENTRIC_P_ARMS),
#      pool the segment-length covered by each CATEGORY (GAIN / LOSS / CN-LOH
#      - not each exact CN value) and check whether one category covers
#      >= --pct-alt-thr% of the arm's total segment length. An arm that's
#      50% CN=1 and 50% CN=0 IS now a LOSS arm (both pool into the LOSS
#      category), unlike the old exact-CN-value grouping. See
#      _arm_majority_category.
#   2. A GAIN/LOSS-qualifying arm is then dropped entirely (treated as
#      non-qualifying, its segments falling through to the multigene/focal
#      section like any other unreported region) if its own rounded weighted
#      CN equals ploidy - it's not a real event, just the tumour's own
#      background. CN-LOH-qualifying arms are never dropped this way: a
#      whole-arm CN-LOH pattern is itself the finding, regardless of ploidy.
#   3. Surviving arms compute their "interesting segments" (see below)
#      immediately, against that arm's own weighted CN.
#   4. The sibling arm is then checked: if it ALSO survived and qualified for
#      the SAME call (GAIN/LOSS/CN-LOH - not requiring the same CN any more),
#      the two arm-level reports are replaced by one chromosome-level report,
#      whose displayed CN is a fresh weight-pooled mean across BOTH arms'
#      segments (see compute_chrom_weighted_cn) - not either arm's own value.
#      If that pooled CN itself equals ploidy (GAIN/LOSS only - can happen
#      even when neither arm's own CN did), promotion is skipped and the two
#      arms are reported individually instead, since each already
#      individually differs from ploidy.
#   Weighted CN (CNVkit `weight`-column weighted mean of `cn`) is rounded to
#   the nearest integer with Python's round() (banker's rounding) - this is
#   what's displayed and used as CNc/CNa below, a separate quantity from the
#   category-majority value that decided qualification.
#
# Within a qualifying arm, a segment is "interesting" (and, if it contains
# >=1 marker gene, surfaced in the multigene/focal events section) when:
#   CNs != CNa   AND   (is_cnloh(CNs, cn1, cn2)  OR  CNs != ploidy)   AND
#   (CNs != 2   OR   is_cnloh(CNs, cn1, cn2))
# CNa is always the segment's own ARM's weighted CN (even when that arm was
# promoted to a chromosome-level call) - a scatter plot of the segment will
# be zoomed to its arm, so the local arm baseline is the relevant comparison,
# not the coarser whole-chromosome pooled mean. A true CN-LOH segment
# overrides the ploidy-equality check (so CN-LOH is surfaced even in a
# diploid-ploidy tumour) but never the CNa check (so a segment matching its
# own arm's already-reported CN-LOH pattern isn't redundantly re-listed).
# See _is_reportable_segment.

VALID_CHROMS = {f"chr{i}" for i in range(1, 23)}

# Acrocentric p-arms: short-arm satellite/rDNA repeat content, not
# meaningfully callable as a CNV event - excluded from arm-level reporting
# (and therefore can never trigger chromosome-level promotion either, since
# that requires both arms to independently qualify).
ACROCENTRIC_P_ARMS = {('chr13', 'p'), ('chr14', 'p'), ('chr15', 'p'), ('chr21', 'p'), ('chr22', 'p')}


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


def _to_num(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def read_cns_segments(cns_path):
    """Parse a CNVkit _call.cns file into segment dicts: chromosome, start,
    end, cn (int), cn1/cn2 (int or None - allele-specific CN, absent when
    `call` ran without a --vcf), weight (float, CNVkit's own per-segment
    confidence), and genes (list of gene symbols overlapping the segment)."""
    rows = []
    with open(cns_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ch = row.get('chromosome', '')
            if ch not in VALID_CHROMS:
                continue
            try:
                s = int(row['start'])
                e = int(row['end'])
                cn = int(round(float(row['cn'])))
                w = float(row['weight'])
            except (ValueError, KeyError, TypeError):
                continue
            if e <= s:
                continue
            cn1_raw = _to_num(row.get('cn1'))
            cn2_raw = _to_num(row.get('cn2'))
            cn1 = int(round(cn1_raw)) if cn1_raw is not None else None
            cn2 = int(round(cn2_raw)) if cn2_raw is not None else None
            genes = [g.strip() for g in row.get('gene', '').split(',') if g.strip()]
            rows.append({'chromosome': ch, 'start': s, 'end': e, 'cn': cn,
                         'cn1': cn1, 'cn2': cn2, 'weight': w, 'genes': genes})
    return rows


def _is_cnloh(cn, cn1, cn2):
    return cn == 2 and cn1 == 2 and cn2 == 0


def _is_reportable_segment(cn, cn1, cn2, ploidy, baseline=None):
    """Whether a segment is a genuine, worth-surfacing event for the
    multigene/focal section (see module docstring for the full rule):
      - never reportable if it matches `baseline` (the segment's own arm's
        weighted CN) - this IS the arm's own already-reported pattern, not a
        deviation from it. `baseline=None` when the segment isn't inside any
        qualifying arm (nothing to compare against here).
      - always reportable if it's genuinely CN-LOH (cn==2, cn1==2, cn2==0) -
        overrides the ploidy check, but not the baseline check above.
      - otherwise not reportable if it sits at ploidy (the tumour's own
        background - not an event), or if it's a plain neutral CN=2 segment
        (cn<=1 or cn>=3 is required otherwise, same as the old absolute
        GAIN/LOSS thresholds - cn != 2 is equivalent to that for integer CN).
    """
    if baseline is not None and cn == baseline:
        return False
    if _is_cnloh(cn, cn1, cn2):
        return True
    if cn == ploidy:
        return False
    return cn != 2


def _event_label(cn, cn1=None, cn2=None):
    """GAIN / LOSS / CN-LOH for a segment already selected for reporting -
    always defined (unlike _is_reportable_segment, which decides whether a
    segment is worth reporting at all)."""
    if _is_cnloh(cn, cn1, cn2):
        return 'CN-LOH'
    return 'GAIN' if cn > 2 else 'LOSS'


def _weighted_cn(rows):
    total_w = sum(r['weight'] for r in rows)
    if total_w <= 0:
        return None
    return sum(r['cn'] * r['weight'] for r in rows) / total_w


def compute_arm_weighted_cn(segments, arms_by_chr):
    """{(chromosome, arm): weighted_cn} across segments overlapping that arm.
    This is the displayed "Copy Number" for a standalone arm-level report,
    and is always what "interesting segment" comparisons use as CNa (even
    for a promoted chromosome's arms) - a separate quantity from the
    category-majority qualification check below."""
    per = defaultdict(list)
    for r in segments:
        for ai in arms_by_chr.get(r['chromosome'], []):
            ovl = max(0, min(r['end'], ai['end']) - max(r['start'], ai['start']))
            if ovl > 0:
                per[(r['chromosome'], ai['arm'])].append(r)
    return {k: _weighted_cn(rows) for k, rows in per.items()}


def compute_chrom_weighted_cn(segments, arms_by_chr, chrom):
    """Weighted CN pooled across every segment overlapping BOTH arms of
    `chrom` - the displayed "Copy Number" for a promoted chromosome-level
    report (distinct from either arm's own weighted CN, which is what
    "interesting segment" comparisons keep using instead - see
    compute_arm_weighted_cn)."""
    rows = []
    for ai in arms_by_chr.get(chrom, []):
        rows.extend(r for _, r in _arm_segments(segments, ai, chrom))
    return _weighted_cn(rows)


def _arm_segments(segments, ai, chrom):
    """(length, segment) pairs for segments overlapping one arm, clipped to
    the arm's own bounds (so a segment spanning the centromere only
    contributes its true overlap)."""
    out = []
    for r in segments:
        if r['chromosome'] != chrom:
            continue
        ovl = max(0, min(r['end'], ai['end']) - max(r['start'], ai['start']))
        if ovl > 0:
            out.append((ovl, r))
    return out


def _arm_majority_category(arm_segs, amp_cn, del_cn, pct_thr):
    """An arm is GAIN/LOSS/CN-LOH if that CATEGORY's segments - pooled
    together regardless of their individual CN, not grouped by exact CN
    value - cover >= pct_thr% of the arm's total segment length. E.g. an arm
    that's 50% CN=1 and 50% CN=0 IS a LOSS arm (both pool into LOSS).
    Denominator is every segment overlapping the arm. Returns None if
    nothing qualifies."""
    total = sum(length for length, _ in arm_segs)
    if total <= 0:
        return None

    loh_len = sum(length for length, r in arm_segs
                  if r['cn'] == 2 and r['cn1'] == 2 and r['cn2'] == 0)
    if 100 * loh_len / total >= pct_thr:
        return 'CN-LOH'

    gain_len = sum(length for length, r in arm_segs if r['cn'] >= amp_cn)
    if 100 * gain_len / total >= pct_thr:
        return 'GAIN'

    loss_len = sum(length for length, r in arm_segs if r['cn'] <= del_cn)
    if 100 * loss_len / total >= pct_thr:
        return 'LOSS'

    return None


def _chr_order(chrom):
    try:
        return int(chrom.replace("chr", ""))
    except ValueError:
        return 10 ** 9


def _segment_arm(row, arms_by_chr):
    """Which arm a segment's midpoint falls in."""
    ch = row['chromosome']
    mid = (row['start'] + row['end']) // 2
    for ai in arms_by_chr.get(ch, []):
        if ai['start'] <= mid < ai['end']:
            return ai['arm']
    return None


def cnv_calls(segments, arms_by_chr, amp_cn, del_cn, pct_thr, ploidy):
    """Arm-first chromosome/arm calling (see module docstring). Each arm is
    qualified independently via _arm_majority_category (a category's pooled
    segments must cover >= pct_thr% of the arm's length), then a GAIN/LOSS
    qualification is dropped if its own rounded weighted CN equals ploidy -
    it's just the tumour's background, not a real event (CN-LOH is exempt).
    Only once BOTH arms of a chromosome independently survive qualification
    for the SAME call (GAIN/LOSS/CN-LOH - CN need not match any more) are
    they promoted to a chromosome-level report, whose CN is a fresh
    weight-pooled mean across both arms' segments; if that pooled CN itself
    equals ploidy (GAIN/LOSS only), promotion is skipped and the two arms
    are kept as individual reports instead.

    Returns (chr_calls, arm_calls, qualifying_arms):
      chr_calls: [{'chromosome', 'call', 'weighted_cn'}] - promoted whole chromosomes
      arm_calls: [{'chromosome', 'arm', 'call', 'weighted_cn'}] - individually reported arms
      qualifying_arms: {(chromosome, arm): {'call', 'weighted_cn'}} for EVERY
        arm that survived qualification, whether or not it ended up promoted
        - this is the basis for "interesting segment" comparisons (always
        arm-scoped, using the arm's own weighted CN - see module docstring).
    `weighted_cn` is always CNVkit's `weight`-column weighted mean, rounded.
    """
    arm_wcn = compute_arm_weighted_cn(segments, arms_by_chr)

    chr_calls, arm_calls = [], []
    qualifying_arms = {}

    for ch in sorted(arms_by_chr, key=_chr_order):
        arm_list = arms_by_chr[ch]
        arm_result = {}
        for ai in arm_list:
            arm = ai['arm']
            key = (ch, arm)
            if key in ACROCENTRIC_P_ARMS:
                continue
            wcn = arm_wcn.get(key)
            if wcn is None:
                continue
            cat = _arm_majority_category(_arm_segments(segments, ai, ch), amp_cn, del_cn, pct_thr)
            if not cat:
                continue
            rounded = round(wcn)
            if cat != 'CN-LOH' and rounded == ploidy:
                continue  # just the tumour's background at this ploidy, not a real event
            arm_result[arm] = {'call': cat, 'weighted_cn': rounded}

        qualifying_arms.update({(ch, arm): res for arm, res in arm_result.items()})

        promoted = False
        if len(arm_list) == 2:
            a1, a2 = (ai['arm'] for ai in arm_list)
            r1, r2 = arm_result.get(a1), arm_result.get(a2)
            if r1 is not None and r2 is not None and r1['call'] == r2['call']:
                chrom_wcn = compute_chrom_weighted_cn(segments, arms_by_chr, ch)
                if chrom_wcn is not None:
                    chrom_rounded = round(chrom_wcn)
                    if r1['call'] == 'CN-LOH' or chrom_rounded != ploidy:
                        chr_calls.append({'chromosome': ch, 'call': r1['call'], 'weighted_cn': chrom_rounded})
                        promoted = True

        if not promoted:
            for arm, res in arm_result.items():
                arm_calls.append({'chromosome': ch, 'arm': arm, **res})

    chr_calls.sort(key=lambda c: _chr_order(c['chromosome']))
    arm_calls.sort(key=lambda a: (_chr_order(a['chromosome']), a['arm']))
    return chr_calls, arm_calls, qualifying_arms


def build_multigene_focal_rows(segments, markers, chr_calls, arm_calls, qualifying_arms,
                                arms_by_chr, ploidy):
    """Multigene/focal events: one row per qualifying segment (all its marker
    genes listed together), combining two sources (see _is_reportable_segment
    and the module docstring for the shared reportability rule):
      (a) marker-gene-containing segments reportable against ploidy alone,
          outside any reported chromosome/arm - this is also where a
          non-qualifying arm's segments end up, including one that qualified
          by category but was then dropped for sitting at ploidy.
      (b) marker-gene-containing segments reportable against their own
          qualifying arm's weighted CN (as well as against ploidy) - checked
          for EVERY arm in qualifying_arms, whether or not it ended up
          promoted to chromosome level (a promoted chromosome's interesting
          segments are just the union of what its two arms already found,
          each still compared to its own arm's CN - see module docstring).
    Deduplicated by (chromosome, start, end); sorted by genomic coordinate,
    then event, then CN.
    """
    reported_chroms = {c['chromosome'] for c in chr_calls}
    reported_arms = {(a['chromosome'], a['arm']) for a in arm_calls}

    combined = {}

    def _add(r):
        marker_genes = [g for g in r['genes'] if g in markers]
        if not marker_genes:
            return
        key = (r['chromosome'], r['start'], r['end'])
        if key in combined:
            return
        combined[key] = {
            'chromosome': r['chromosome'], 'start': r['start'], 'end': r['end'],
            'cn': r['cn'], 'event': _event_label(r['cn'], r['cn1'], r['cn2']),
            'marker_genes': sorted(set(marker_genes)),
        }

    # (a) segments outside any reported chromosome/arm, reportable against ploidy
    for r in segments:
        ch = r['chromosome']
        arm = _segment_arm(r, arms_by_chr)
        if ch in reported_chroms or (ch, arm) in reported_arms:
            continue
        if _is_reportable_segment(r['cn'], r['cn1'], r['cn2'], ploidy):
            _add(r)

    # (b) segments within every qualifying arm, reportable against that arm's
    # own weighted CN (regardless of whether the arm ended up promoted)
    for (ch, arm), res in qualifying_arms.items():
        ai = next((a for a in arms_by_chr.get(ch, []) if a['arm'] == arm), None)
        if ai is None:
            continue
        for length, r in _arm_segments(segments, ai, ch):
            if _is_reportable_segment(r['cn'], r['cn1'], r['cn2'], ploidy, baseline=res['weighted_cn']):
                _add(r)

    rows = list(combined.values())
    rows.sort(key=lambda x: (_chr_order(x['chromosome']), x['start'], x['event'], x['cn']))
    return rows


def appendix_marker_genes_in_reported(segments, chr_calls, arm_calls, arms_by_chr, markers):
    """Appendix A: one row per reported chromosome/arm feature, listing every
    marker gene physically located in its segments (regardless of that
    segment's own CN) alongside the feature's own AMP/DEL/CN-LOH call."""
    region_event = {c['chromosome']: c['call'] for c in chr_calls}
    region_event.update({f"{a['chromosome']}{a['arm']}": a['call'] for a in arm_calls})

    reported_chroms = {c['chromosome'] for c in chr_calls}
    reported_arms = {(a['chromosome'], a['arm']) for a in arm_calls}

    region_genes = defaultdict(set)
    for r in segments:
        ch = r['chromosome']
        marker_genes = [g for g in r['genes'] if g in markers]
        if not marker_genes:
            continue

        if ch in reported_chroms:
            region_genes[ch].update(marker_genes)
        for ai in arms_by_chr.get(ch, []):
            key = (ch, ai['arm'])
            if key not in reported_arms:
                continue
            ovl = max(0, min(r['end'], ai['end']) - max(r['start'], ai['start']))
            if ovl > 0:
                region_genes[f"{ch}{ai['arm']}"].update(marker_genes)

    rows = [{'region': region, 'event': region_event[region], 'genes': sorted(genes)}
            for region, genes in region_genes.items()]
    rows.sort(key=lambda x: (_chr_order(x['region']), x['region']))
    return rows


def write_scatter_manifest(path, chrom_keys, segment_rows):
    """TSV consumed by the CREATE_CNVKIT_SCATTER module: one row per
    chromosome-wide scatter plot (chromosome/arm-level events) plus one row
    per multigene/focal event segment (region ± 5 Mb, labelled with that
    segment's marker gene(s)). The appendices carry no scatter plots of
    their own - all plots are consolidated in Appendix B.

    `target_start`/`target_end` are the segment's own (unpadded) bounds -
    blank for chrom-type rows, where there's no single "target" segment.
    cnv_absolute_scatter.py uses this to colour the one segment a focal plot
    was actually generated for green, versus orange for any other segment
    that also happens to be reportable but is merely incidentally in view."""
    seen = set()
    with open(path, 'w') as f:
        f.write("key\ttype\tregion\tgenes\ttarget_start\ttarget_end\n")
        for ch in sorted(set(chrom_keys), key=_chr_order):
            if ch in seen:
                continue
            seen.add(ch)
            f.write(f"{ch}\tchrom\t{ch}\t\t\t\n")
        for r in segment_rows:
            key = f"{r['chromosome']}_{r['start']}_{r['end']}"
            if key in seen:
                continue
            seen.add(key)
            region = f"{r['chromosome']}:{max(0, r['start'] - 5_000_000)}-{r['end'] + 5_000_000}"
            genes = ','.join(sorted(set(r['marker_genes'])))
            f.write(f"{key}\tsegment\t{region}\t{genes}\t{r['start']}\t{r['end']}\n")


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
            'status': 'Microsatellite Instability' if pct >= 3.5 else 'Microsatellite Stable',
        }
    except Exception as e:
        return {'error': str(e)}


def _format_hla_allele(raw):
    """polysolver alleles look like 'hla_a_32_01_01' or 'hla_a_02_01_01_01' -
    gene letter + numeric fields; display convention is 2-field resolution:
    'hla_a_32_01_01' -> 'A*32:01'."""
    parts = raw.strip().split('_')
    if len(parts) < 4 or parts[0] != 'hla':
        return raw.strip()
    gene = parts[1].upper()
    f1, f2 = parts[2], parts[3]
    return f"{gene}*{f1}:{f2}"


def parse_hla(hla_path):
    """Polysolver _hla.txt output embeds a long alignment/calling log; the
    result is just the 'winners1'/'winners2' lines (tab-separated: label,
    HLA-A, HLA-B, HLA-C allele). winners1 is written twice in practice -
    take the first occurrence of each label."""
    try:
        winners = {}
        with open(hla_path) as f:
            for line in f:
                if line.startswith('winners1') or line.startswith('winners2'):
                    cols = line.rstrip('\n').split('\t')
                    if len(cols) < 4:
                        continue
                    label = cols[0]
                    if label not in winners:
                        winners[label] = (cols[1], cols[2], cols[3])
        if 'winners1' not in winners or 'winners2' not in winners:
            return {'error': 'winners1/winners2 lines not found'}
        a1, b1, c1 = winners['winners1']
        a2, b2, c2 = winners['winners2']
        return {
            'A': (_format_hla_allele(a1), _format_hla_allele(a2)),
            'B': (_format_hla_allele(b1), _format_hla_allele(b2)),
            'C': (_format_hla_allele(c1), _format_hla_allele(c2)),
        }
    except Exception as e:
        return {'error': str(e)}


def parse_hotspot_counts(tsv_path):
    """Rows from the WHITELIST module's TSV - one row per WHO CNS5 hotspot
    (exact CHROM/POS/REF/ALT match), with tumour/normal read depths pulled
    straight from the raw VCF, or 0/0 if that exact position+allele wasn't
    called. See modules/whitelist.nf for how the match is made."""
    rows = []
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.append(row)
    return rows


def parse_patient_metadata(path, patient_id, case_id):
    """Look up this patient's row in the metadata TSV (--patient-metadata).
    The TSV is keyed by 'case_id' in MOLBI###### form; the pipeline's
    patient_id is 'M######' (e.g. 'M24007' -> 'MOLBI24007'). Falls back to an
    exact patient_id / case_id / tumour-label match. Returns the row as a
    dict, or {} when the file is missing or no row matches - front-matter
    fields then render as an em dash."""
    if not path or not os.path.exists(path):
        return {}
    digits = ''.join(c for c in patient_id if c.isdigit())
    keys = {f"MOLBI{digits}", patient_id, case_id, case_id.split('_')[0]}
    try:
        with open(path, newline='') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                if (row.get('case_id') or '').strip() in keys:
                    return {k: (v or '').strip() for k, v in row.items()}
    except OSError:
        return {}
    return {}


def parse_facets_purity(path):
    """FACETS writes one line: a tumour-purity fraction in (0, 1], or the
    literal 'NA' when the fit failed or there was no tumour sample. Returns a
    float, or None when unavailable (no path, missing file, 'NA', empty,
    unparseable, or out of the (0, 1] range)."""
    if not path:
        return None
    try:
        with open(path) as f:
            raw = f.read().strip()
        val = float(raw)
    except (OSError, ValueError):
        return None
    return val if 0.0 < val <= 1.0 else None


def parse_facets_ploidy(path):
    """FACETS writes one line: a positive ploidy estimate, or the literal
    'NA' when unavailable. Returns the rounded int ploidy - unlike purity,
    ploidy is used as a hard baseline throughout the CNV calling below, so
    unlike parse_facets_purity this always returns a usable value: it falls
    back to DEFAULT_PLOIDY (2, the historical fixed assumption) when the
    file is absent, 'NA', empty, unparseable, or not positive."""
    if not path:
        return DEFAULT_PLOIDY
    try:
        with open(path) as f:
            raw = f.read().strip()
        val = float(raw)
    except (OSError, ValueError):
        return DEFAULT_PLOIDY
    return round(val) if val > 0.0 else DEFAULT_PLOIDY


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


def _wrap(s, width=25):
    if not s or len(s) <= width:
        return s
    return '<br>'.join(s[i:i+width] for i in range(0, len(s), width))


def _comma(v):
    """Format a genomic position with thousands separators (208248388 ->
    208,248,388). Values that aren't integer-like (missing-data placeholders
    such as '' or '?') pass through unchanged."""
    try:
        return f"{int(v):,}"
    except (TypeError, ValueError):
        return v


# ── Report writer ─────────────────────────────────────────────────────────────

def write_report(args, patient_meta, purity, ploidy, somatic, germline, hotspot_counts, chr_calls, arm_calls,
                  multigene_focal_rows, appendix_regions_rows, chrom_plot_keys,
                  msi, tmb, tumour_hla, germline_hla, supplementary):
    today = datetime.date.today().isoformat()
    m = patient_meta or {}
    def mv(key):
        return (m.get(key) or '').strip() or '—'

    with open(args.out, 'w') as f:

        # ── Front matter ──
        f.write("# CGenWES · WES-based genomic profiling report\n\n")
        f.write("Clinical Genomics and Bioinformatics (CGen) Group  \n")
        f.write("Department of Cancer Medicine, School of Translational Medicine, Monash University\n\n")
        f.write(f"**Date of report:** {today}\n\n")
        f.write("_Research-grade report. Produced by the Clinical Genomics and Bioinformatics "
                "(CGen) Group, Monash University. It is not issued by an accredited pathology "
                "service, is not validated for clinical use, and does not constitute medical "
                "advice. The result is one line of evidence for multidisciplinary consideration "
                "and must be interpreted by a clinician in the context of the patient's complete "
                "clinical picture._\n\n")

        f.write("### Patient\n\n")
        f.write("|  |  |\n| --- | --- |\n")
        f.write(f"| Case ID | {args.patient_id} |\n")
        f.write(f"| Year of birth | {mv('year_of_birth')} |\n")
        f.write(f"| Sex | {mv('sex')} |\n\n")

        f.write("### Tumour\n\n")
        f.write("|  |  |\n| --- | --- |\n")
        f.write(f"| Tumour event | {mv('tumour_event')} |\n")
        f.write(f"| Tumour type | {mv('tumour_type')} |\n")
        f.write(f"| WHO grade | {mv('who_grade')} |\n\n")

        f.write("### Assay & specimen\n\n")
        f.write("|  |  |\n| --- | --- |\n")
        f.write(f"| Assay | {mv('assay')} |\n")
        f.write(f"| Platform | {mv('sequencing_platform')} |\n")
        f.write(f"| Tissue collection date | {mv('tissue_collected')} |\n")
        f.write(f"| Normal collection date | {mv('normal_collected')} |\n\n")

        # ── Filtering parameters (moved to page 1) ──
        f.write("### Filtering parameters\n\n")
        f.write("| Parameter | Value |\n")
        f.write("| --- | --- |\n")
        f.write(f"| Chromosome/arm GAIN/LOSS/CN-LOH threshold (% length) | ≥ {args.pct_alt_thr:.0f}% |\n")
        f.write(f"| Gained segment threshold (CN ≥) | {args.amp_seg_cn} |\n")
        f.write(f"| Lost segment threshold (CN ≤) | {args.del_seg_cn} |\n")
        f.write(f"| CN-LOH segment definition | CN = 2, allele-specific CN 2 / 0 |\n")
        f.write("\n")

        # ── Tumour Purity (moved below Filtering parameters, above Limitations) ──
        f.write("### Tumour Purity\n\n")
        f.write("| Metric | Value |\n")
        f.write("| --- | --- |\n")
        if purity is None:
            f.write("| Tumour purity | _Not available_ |\n")
        elif purity < PURITY_MIN_REPORTABLE:
            f.write(f"| Tumour purity | _Not reported (estimate < {PURITY_MIN_REPORTABLE:.1f})_ |\n")
        else:
            f.write(f"| Tumour purity | {purity:.2f} |\n")
        f.write("\n")
        f.write("_Tumour purity is the fraction of tumour cells in the biopsy. It affects how "
                "credible the molecular findings below are - higher purity means the findings "
                f"are more reliable. A purity estimate below {PURITY_MIN_REPORTABLE:.1f} is not "
                "reported, and the findings should then be interpreted with consideration. "
                "Where a tumour-purity estimate was available it was used to adjust the CNV "
                "findings._\n")
        f.write("\n")

        # ── Limitations (moved to page 1) ──
        f.write("### Limitations\n\n")
        f.write("- This genomic test does not identify gene fusions nor somatic variants "
                "with an allele frequency lower than 0.10.\n")
        f.write("- CNV calling in this report is based on total copy number only. Segments "
                "with the same total copy number but a different allele-specific split are "
                "reported as the same event - e.g. CN=4 composed of 3 copies of one allele "
                "and 1 of the other (allelic imbalance) is reported identically to CN=4 "
                "composed of 4 copies of one allele and 0 of the other (loss of the other "
                "allele), even though the underlying biology differs.\n")
        f.write("- Sensitivity of the test is limited by tumour purity, which was not "
                "assessed by a histopathologist.\n")
        f.write("- This report is the outcome of a research test conducted under strict "
                "scientific protocols. It is not a diagnostic report.\n")
        f.write("\n")

        # ── Manual comment box (moved to page 1, last) ──
        f.write("## Laboratory Comment\n\n")
        f.write("> _(Please enter your comments here)_\n\n")

        # ── 1. SNV Findings ──
        f.write("## 1. SNV Findings\n\n")

        f.write("### Somatic mutations (actionable / hotspot)\n\n")
        if somatic:
            f.write("| Gene | Location | Transcript | HGVSC | HGVSP | AF | Depth | Consequence | Selection reason |\n")
            f.write("| --- | --- | --- | --- | --- | --- | ---: | --- | --- |\n")
            for v in somatic:
                consequence = '' if v['consequence'] == 'whitelist_hotspot' else v['consequence']
                f.write(f"| {v['gene']} | {v['chrom']}:{_comma(v['pos'])} | {v['tx']} "
                        f"| {_wrap(v['cdna'])} | {v['prot']} | {v['af']:.2f} | {v['dp']} "
                        f"| {consequence} | {v['reason']} |\n")
        else:
            f.write("_No actionable/hotspot somatic variants detected._\n")
        f.write("\n")
        f.write("_Variants with Selection reason \"Whitelist match\" are included from the WHO "
                "CNS5 hotspot whitelist by exact genomic position, regardless of this pipeline's "
                "depth/allele-frequency filtering._\n")
        f.write("\n")

        f.write("### Germline pathogenic variants\n\n")
        if germline:
            f.write("| Gene | Location | Transcript | HGVSC | HGVSP | CLNSIG | Depth | Consequence | Selection reason |\n")
            f.write("| --- | --- | --- | --- | --- | --- | ---: | --- | --- |\n")
            for v in germline:
                f.write(f"| {v['gene']} | {v['chrom']}:{_comma(v['pos'])} | {v['tx']} "
                        f"| {_wrap(v['cdna'])} | {v['prot']} | {v['clnsig']} | {v['dp']} "
                        f"| {v['consequence']} | {v['reason']} |\n")
        else:
            f.write("_No germline pathogenic/likely-pathogenic variants detected._\n")
        f.write("\n_Only variants in genes on the eviQ germline testing panel are included._\n")
        f.write("\n")

        f.write("### WHO CNS5 Hotspot Read Counts\n\n")
        if hotspot_counts:
            f.write("| Gene | Chrom | Position | Hotspot | Tumour REF depth | Tumour ALT depth "
                    "| Germline REF depth | Germline ALT depth |\n")
            f.write("| --- | --- | ---: | --- | ---: | ---: | ---: | ---: |\n")
            for h in hotspot_counts:
                f.write(f"| {h.get('GENE','')} | {h.get('CHROM','')} | {_comma(h.get('POS',''))} | {_wrap(h.get('HOTSPOT',''))} "
                        f"| {h.get('TUMOUR_REF_DEPTH','?')} | {h.get('TUMOUR_ALT_DEPTH','?')} "
                        f"| {h.get('GERMLINE_REF_DEPTH','?')} | {h.get('GERMLINE_ALT_DEPTH','?')} |\n")
        elif hotspot_counts is None:
            f.write("_No whitelist was included in this run, so no hotspot read counts were "
                    "computed. If a whitelist is provided (via `params.whitelist`), this section "
                    "would show a table of tumour and germline REF/ALT read depths at each "
                    "whitelisted position, pulled directly from the raw alignment (BAM)._\n")
        else:
            f.write("_No hotspot read-count data available._\n")
        f.write("\n")

        # ── 2. CNV Findings ──
        f.write("## 2. CNV Findings\n\n")

        f.write("| Metric | Value |\n")
        f.write("| --- | --- |\n")
        f.write(f"| Tumour ploidy (rounded) | {ploidy} |\n\n")
        f.write("_Tumour ploidy is the estimated whole-genome average copy number, rounded to "
                "the nearest integer. It is the baseline \"neutral\" copy number below - a "
                "chromosome/arm or segment whose copy number equals ploidy reflects the "
                "tumour's own background, not a distinct event, and is not reported as one "
                "(copy-neutral LOH is the exception: it is always reported regardless of "
                "ploidy, since it is a finding about allelic imbalance, not total copy "
                "number)._\n\n")

        f.write("### Chromosome-level events\n\n")
        if chr_calls:
            f.write("| Chromosome | Event | Copy Number |\n")
            f.write("| ---: | --- | ---: |\n")
            for c in chr_calls:
                f.write(f"| {c['chromosome']} | {c['call']} | {c['weighted_cn']} |\n")
        else:
            f.write("_No whole-chromosome GAIN/LOSS/CN-LOH events detected._\n")
        f.write("\n")

        f.write("### Chromosome-arm events\n\n")
        if arm_calls:
            f.write("| Chromosome | Arm | Event | Copy Number |\n")
            f.write("| ---: | :---: | --- | ---: |\n")
            for a in arm_calls:
                f.write(f"| {a['chromosome']} | {a['arm']} | {a['call']} | {a['weighted_cn']} |\n")
        else:
            f.write("_No chromosome-arm GAIN/LOSS/CN-LOH events detected (or promoted to "
                    "chromosome-level above)._\n")
        acrocentric_labels = sorted((f"{ch}{arm}" for ch, arm in ACROCENTRIC_P_ARMS),
                                     key=lambda s: _chr_order(s))
        f.write(f"\n_Acrocentric p-arms ({', '.join(acrocentric_labels)}) are never reported. "
                f"Consequently, these chromosomes can at most have q-arm level events._\n")
        f.write("\n")

        f.write("### Multigene and focal events (cancer genes¹ only)\n\n")
        if multigene_focal_rows:
            f.write("| Chromosome | Start | End | Event | Copy Number | Cancer genes |\n")
            f.write("| ---: | ---: | ---: | --- | ---: | --- |\n")
            for r in multigene_focal_rows:
                f.write(f"| {r['chromosome']} | {_comma(r['start'])} | {_comma(r['end'])} | {r['event']} | {r['cn']} "
                        f"| {', '.join(r['marker_genes'])} |\n")
        else:
            f.write("_No multigene/focal cancer-gene events detected._\n")
        f.write("\n_¹ Cancer genes include oncogenes from OncoKB and COSMIC._\n")
        f.write("\n")

        # ── 3. MSI  ──
        f.write("## 3. Microsatellite Instability Assessment \n\n")
        f.write("| Metric | Value |\n")
        f.write("| --- | --- |\n")
        if msi and 'error' not in msi:
            f.write(f"| MSI total sites | {msi.get('total_sites', '?')} |\n")
            f.write(f"| MSI unstable sites | {msi.get('unstable_sites', '?')} |\n")
            f.write(f"| MSI score (%) | {msi.get('pct', 0):.2f} |\n")
            f.write(f"| MSI status | **{msi.get('status', '?')}** |\n")
        else:
            f.write("| MSI | _N/A_ |\n")
        f.write("\n")
        f.write("_Tumours are classified as MSI if the MSI score is ≥ 3.5 and MSS if < 3.5._\n")
        f.write("\n")

        # ── 4. HLA Typing ──
        f.write("## 4. HLA Typing\n\n")
        if tumour_hla and 'error' not in tumour_hla:
            f.write("| Locus | Tumour Allele 1 | Tumour Allele 2 |\n")
            f.write("| --- | --- | --- |\n")
            for locus in ('A', 'B', 'C'):
                t1, t2 = tumour_hla.get(locus, ('?', '?'))
                f.write(f"| HLA-{locus} | {t1} | {t2} |\n")
        else:
            f.write("_No HLA typing data available._\n")
        f.write("\n")

        # ── 5. Tumour Molecular Burden  ──
        f.write("## 5. Tumour Molecular Burden Summary\n\n")
        f.write("| TMB (mut/Mb) | Total variants | Variants after filters |\n")
        f.write("| ---: | ---: | ---: |\n")
        if tmb and 'error' not in tmb:
            f.write(f"| {tmb.get('tmb', tmb.get('TMB', '?'))} | {tmb.get('total_variants', '?')} | {tmb.get('coding_variants', '?')} |\n")
        else:
            f.write("| _N/A_ | _N/A_ | _N/A_ \n")
        f.write("\n")

        # ── Appendix A: cancer genes within reported chromosome/arm events ──
        f.write("## Appendix A – Cancer genes¹ within reported chromosome/arm events\n\n")
        if appendix_regions_rows:
            f.write("| Feature | Event | Cancer genes |\n")
            f.write("| --- | --- | --- |\n")
            for r in appendix_regions_rows:
                f.write(f"| {r['region']} | {r['event']} | {', '.join(r['genes'])} |\n")
        else:
            f.write("_No reported chromosome/arm events, or no cancer genes within them._\n")
        f.write("\n_¹ Cancer genes include oncogenes from OncoKB and COSMIC._\n")
        f.write("\n")

        # ── Appendix B: scatter plots (consolidated - none appear inline above) ──
        f.write("## Appendix B – Scatter Plots\n\n")

        f.write("### Chromosome/arm-level events\n\n")
        region_genes_lookup = {r['region']: r['genes'] for r in appendix_regions_rows}
        chrom_events = defaultdict(list)
        for c in chr_calls:
            chrom_events[c['chromosome']].append(f"{c['chromosome']} {c['call']}")
        for a in arm_calls:
            chrom_events[a['chromosome']].append(f"{a['chromosome']}{a['arm']} {a['call']}")
        if chrom_events:
            for ch in sorted(chrom_events, key=_chr_order):
                # Genes for every reported feature (chromosome and/or its arm(s))
                # on this line, deduplicated and combined across all of them.
                feature_keys = [ch] + [f"{ch}{a['arm']}" for a in arm_calls if a['chromosome'] == ch]
                genes = sorted({g for k in feature_keys for g in region_genes_lookup.get(k, [])})
                gene_str = f"  \nCancer genes: {', '.join(genes)}" if genes else ""
                f.write(f"**{ch}** — {', '.join(chrom_events[ch])}{gene_str}\n\n")
                f.write(f"<!--SCATTER:{ch}-->\n\n")
        else:
            f.write("_No chromosome/arm-level events to plot._\n\n")

        f.write("### Multigene/focal events\n\n")
        if multigene_focal_rows:
            for r in multigene_focal_rows:
                loc = f"{r['chromosome']}:{r['start']:,}–{r['end']:,}"
                key = f"{r['chromosome']}_{r['start']}_{r['end']}"
                f.write(f"**{loc}** &nbsp; {r['event']} &nbsp; CN = {r['cn']} &nbsp; "
                        f"Cancer genes: {', '.join(r['marker_genes'])}\n\n")
                f.write(f"<!--SCATTER:{key}-->\n\n")
        else:
            f.write("_No multigene/focal events to plot._\n\n")

        # ── Supplementary: Tool and annotation versions - manually asserted,
        #    shared across the cohort (see supplementary.json / load_supplementary) ──
        f.write("## Supplementary – Tool and Annotation Versions\n\n")
        f.write(f"_As asserted in `supplementary.json` on {today}._\n\n")

        f.write("### Tools\n\n")
        tools = [(n, v) for n, v in (supplementary.get("Tools") or {}).items() if v]
        if tools:
            f.write("| Tool | Version |\n")
            f.write("| --- | --- |\n")
            for name, version in tools:
                f.write(f"| {name} | `{version}` |\n")
        else:
            f.write("_Tools not provided by user._\n")
        f.write("\n")

        f.write("### Annotation Databases\n\n")
        dbs = [(n, v) for n, v in (supplementary.get("Annotation Databases") or {}).items() if v]
        if dbs:
            f.write("| Database | Version / Build |\n")
            f.write("| --- | --- |\n")
            for name, version in dbs:
                f.write(f"| {name} | `{version}` |\n")
        else:
            f.write("_Annotation databases not provided by user._\n")
        f.write("\n")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    markers = load_markers(args.markers)
    whitelist_keys = load_whitelist_keys(args.whitelist)

    somatic        = parse_somatic_vcf(args.somatic_vcf, markers, whitelist_keys)
    germline       = parse_germline_vcf(args.germline_vcf)
    hotspot_counts = parse_hotspot_counts(args.hotspot_counts) if args.hotspot_counts else None

    _, arms_by_chr = load_chrom_arm_info(args.chrom_arm)
    segments = read_cns_segments(args.cns)

    ploidy = parse_facets_ploidy(args.facets_ploidy)

    chr_calls, arm_calls, qualifying_arms = cnv_calls(
        segments, arms_by_chr, args.amp_seg_cn, args.del_seg_cn, args.pct_alt_thr, ploidy)

    multigene_focal_rows = build_multigene_focal_rows(
        segments, markers, chr_calls, arm_calls, qualifying_arms, arms_by_chr, ploidy)

    appendix_regions_rows = appendix_marker_genes_in_reported(
        segments, chr_calls, arm_calls, arms_by_chr, markers)

    chrom_plot_keys = {c['chromosome'] for c in chr_calls} | {a['chromosome'] for a in arm_calls}
    write_scatter_manifest(args.scatter_manifest, chrom_plot_keys, multigene_focal_rows)

    patient_meta = parse_patient_metadata(args.patient_metadata, args.patient_id, args.case_id)
    purity = parse_facets_purity(args.facets_purity)
    msi = parse_msi(args.msi)
    tmb = parse_tmb(args.tmb)
    tumour_hla   = parse_hla(args.tumour_hla)
    germline_hla = parse_hla(args.germline_hla)

    supplementary = load_supplementary(args.supplementary_json)

    write_report(args, patient_meta, purity, ploidy, somatic, germline, hotspot_counts, chr_calls, arm_calls,
                 multigene_focal_rows, appendix_regions_rows, chrom_plot_keys,
                 msi, tmb, tumour_hla, germline_hla, supplementary)
    print(f"Report written to {args.out}", file=sys.stderr)
    print(f"Scatter manifest written to {args.scatter_manifest}", file=sys.stderr)


if __name__ == '__main__':
    main()
