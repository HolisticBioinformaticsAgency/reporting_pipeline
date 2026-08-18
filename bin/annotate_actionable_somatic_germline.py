#!/usr/bin/env python3
"""
Annotate a VCF directly from raw sources - OncoKB, COSMIC Cancer Mutation
Census (CMC), cancerhotspots.org v1/v2/v3 (kept separate, not merged), and
the eviQ germline testing gene list - with no intermediate reference table.

Every annotation is its own separate INFO tag, semicolon-delimited (the
normal VCF INFO separator - no packed '|'-delimited sub-values), named
[Database]_[Column_Name] using each source's own real column name (spaces
replaced with underscores). A tag is only written when its match succeeds;
absence of the tag is what "false" means for OncoKB/COSMIC, since the
presence of real data already implies a positive match. The two exceptions
that stay explicit T/F are cancerhotspots (no informative value to attach,
per version) and eviQ_germline (a boolean testing-protocol flag):

  OncoKB_Gene_Type=<value>
      Only present if the gene (ANN Gene_Name) is in OncoKB's gene list.

  COSMIC_ONC_TSG=<value>;COSMIC_CGC_TIER=<value>
      Only present if the gene has a non-blank ONC_TSG/CGC_TIER in CMC (the
      ~768 Cancer-Gene-Census-equivalent genes, not CMC's full ~20,000-gene
      universe).

  COSMIC_Mutation_CDS=<value>;COSMIC_Mutation_AA=<value>;
  COSMIC_Mutation_genome_position_GRCh38=<value>;
  COSMIC_MUTATION_SIGNIFICANCE_TIER=<value>
      Only present if the gene already qualified above AND this variant's
      CHROM:POS (chr-prefix stripped) and protein change both exactly match
      the same cataloged CMC entry. Unlike earlier versions of this script,
      MUTATION_SIGNIFICANCE_TIER is NOT restricted to 1/2/3 anymore - 'Other'
      tier mutations are included and annotated too, so this now reflects
      every CMC-cataloged mutation with a matching position + protein
      change, not just the significance-tiered ones. Note this means the
      mutation lookup now holds CMC's full ~5.8M rows in memory rather than
      the ~60k tier 1/2/3 subset used previously - slower to build, more
      memory.

  cancerhotspots_v1=T, cancerhotspots_v2=T, cancerhotspots_v3=T
      Checked independently per version (not merged, since v1/v2/v3 have
      genuinely different schemas). Only present when this variant's
      gene+codon position matches that version's own hotspot table -
      omitted entirely (not an explicit F) otherwise.

  eviQ_germline=T
      Only present when the gene is on the eviQ germline testing gene list
      (gene-symbol match only - see load_eviq_germline_genes() for the
      qualifier-suffix caveat) - omitted entirely (not an explicit F)
      otherwise.

Usage:
    python3 annotate_actionable_somatic_germline.py \
        --vcf M24087T_M24087N_vardict_MANE_clinvar_gnomad_dbNSFP_spliceAI.vcf \
        --oncokb-gene-list ../data/6.cancerGeneList.tsv \
        --cosmic-cmc ../data/7.CancerMutationCensus_extracted/CancerMutationCensus_AllData_v104_GRCh37.tsv.gz \
        --hotspots-v1 ../data/9.hotspots.xls \
        --hotspots-v2 ../data/9.hotspots_v2.xls \
        --hotspots-v3 ../data/9.hotspots_v3.xlsx \
        --eviq-germline ../data/10.germline_eviq.csv \
        --out M24087T_M24087N_vardict_MANE_clinvar_gnomad_dbNSFP_spliceAI_oncorecur.vcf
"""
import argparse
import csv
import datetime
import gzip
import re
import sys

import pandas as pd
import pysam

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

HGVSP_3LETTER_RE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
HGVSP_POS_RE = re.compile(r"^p\.[A-Z][a-z]{2}(\d+)")
COSMIC_AA_1LETTER_RE = re.compile(r"^p\.([A-Z])(\d+)([A-Z*])$")
GENOME_POS_RE = re.compile(r"^(\S+):(\d+)-(\d+)$")
CODON_POS_RE = re.compile(r"(\d+)")


# ---------------------------------------------------------------- loaders

def load_oncokb(path):
    """Gene -> Gene Type."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["Hugo Symbol"].strip()
            gt = row["Gene Type"].strip()
            if gene and gt:
                out[gene] = gt
    return out


def load_cosmic_cmc(path):
    """Single pass over CMC building two lookups:
      - gene_role_tier: GENE_NAME -> (ONC_TSG, CGC_TIER), gene-level.
      - mutation_lookup: (Gene, chrom, start, ref, aa_pos, alt) -> (Mutation
        CDS, Mutation AA, GRCh38 position string, MUTATION_SIGNIFICANCE_TIER)
        for EVERY tier including 'Other' (not filtered to 1/2/3)."""
    gene_role_tier = {}
    mutation_lookup = {}
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["GENE_NAME"].strip()
            if not gene:
                continue

            if gene not in gene_role_tier:
                # ONC_TSG uses ", " for compound categories (e.g. "oncogene, fusion"),
                # but VCF INFO values can't contain literal whitespace - strip it,
                # keeping the comma (e.g. "oncogene,fusion").
                role = row["ONC_TSG"].strip().replace(" ", "")
                tier = row["CGC_TIER"].strip()
                if role or tier:
                    gene_role_tier[gene] = (role, tier)

            mut_aa = row["Mutation AA"].strip()
            genome_pos = row["Mutation genome position GRCh38"].strip()
            mut_tier = row["MUTATION_SIGNIFICANCE_TIER"].strip()
            if mut_aa and genome_pos:
                aa_m = COSMIC_AA_1LETTER_RE.match(mut_aa)
                pos_m = GENOME_POS_RE.match(genome_pos)
                if aa_m and pos_m:
                    ref, aa_pos, alt = aa_m.group(1), aa_m.group(2), aa_m.group(3)
                    chrom, start = pos_m.group(1), pos_m.group(2)
                    mutation_lookup[(gene, chrom, start, ref, aa_pos, alt)] = (
                        row["Mutation CDS"].strip(),
                        mut_aa,
                        genome_pos,
                        mut_tier,
                    )
    return gene_role_tier, mutation_lookup


def load_hotspots_v1(path):
    """Hugo Symbol -> set of positions, parsed from Codon (e.g. 'V600').
    Real data starts at row 3 (title/license preamble above it); 3 trailing
    footnote rows with no real gene/codon are naturally excluded since they
    have no parseable position."""
    df = pd.read_excel(path, engine="xlrd", header=2)
    out = {}
    for gene, codon in zip(df["Hugo Symbol"], df["Codon"]):
        m = CODON_POS_RE.search(str(codon))
        gene = str(gene).strip()
        if gene and gene != "nan" and m:
            out.setdefault(gene, set()).add(m.group(1))
    return out


def load_hotspots_v2(path):
    """Hugo_Symbol -> set of positions, from the explicit
    Amino_Acid_Position column (kept as-is, not regex re-parsed, since a
    handful of entries are non-numeric splice-site notations like
    'X187_splice')."""
    df = pd.read_excel(path, engine="xlrd")
    out = {}
    for gene, pos in zip(df["Hugo_Symbol"], df["Amino_Acid_Position"]):
        gene = str(gene).strip()
        pos = str(pos).strip()
        if gene and gene != "nan" and pos and pos != "nan":
            out.setdefault(gene, set()).add(pos)
    return out


def load_hotspots_v3(path):
    """Hugo_Symbol -> set of positions, parsed from Codon (e.g. 'E40')."""
    df = pd.read_excel(path)
    out = {}
    for gene, codon in zip(df["Hugo_Symbol"], df["Codon"]):
        m = CODON_POS_RE.search(str(codon))
        gene = str(gene).strip()
        if gene and gene != "nan" and m:
            out.setdefault(gene, set()).add(m.group(1))
    return out


def load_eviq_germline_genes(path):
    """Gene column is usually a plain symbol, but 10/63 entries append a
    qualifier describing which variant type the protocol applies to (e.g.
    'CHEK2 (truncating variants)', 'MUTYH (biallelic variants)', 'HOXB13
    G84E (c.251G>A)'). The qualifier is ignored here - matching is on gene
    symbol alone (first whitespace-separated token). Known limitation: this
    reports eviQ_germline=T for any variant in these genes, not just ones
    meeting the qualifier's actual condition (e.g. a single heterozygous
    MUTYH variant, even though its protocol is specifically for biallelic
    carriers; any HOXB13 variant, even though its protocol is specifically
    for the G84E hotspot)."""
    genes = set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            gene_field = row["Gene"].strip()
            if gene_field:
                genes.add(gene_field.split()[0])
    return genes


# ------------------------------------------------------------- VCF parsing

def parse_protein_change_1letter(hgvsp):
    """'p.Val600Glu' -> ('V', '600', 'E'); None if not a simple substitution."""
    if not hgvsp:
        return None
    m = HGVSP_3LETTER_RE.match(hgvsp)
    if not m:
        return None
    ref3, pos, alt3 = m.group(1), m.group(2), m.group(3)
    ref1, alt1 = AA3TO1.get(ref3), AA3TO1.get(alt3)
    if ref1 is None or alt1 is None:
        return None
    return ref1, pos, alt1


def parse_protein_position(hgvsp):
    if not hgvsp:
        return None
    m = HGVSP_POS_RE.match(hgvsp)
    return m.group(1) if m else None


def iter_ann_blocks(record):
    """Yield (gene, hgvsp) for each ANN block on this record."""
    ann = record.info.get("ANN")
    if not ann:
        return
    if isinstance(ann, str):
        ann = (ann,)
    for block in ann:
        fields = block.split("|")
        gene = fields[3] if len(fields) > 3 else ""
        hgvsp = fields[10] if len(fields) > 10 else ""
        yield gene, hgvsp


# ------------------------------------------------------------- annotation

def annotate_record(record, oncokb, gene_role_tier, mutation_lookup,
                     hotspots_v1, hotspots_v2, hotspots_v3, eviq_germline_genes):
    result = {}  # INFO tag -> value, only for tags that actually apply

    record_chrom = record.chrom[3:] if record.chrom.startswith("chr") else record.chrom
    record_pos = str(record.pos)

    oncokb_done = False
    cosmic_gene_done = False
    hs_v1 = hs_v2 = hs_v3 = False
    eviq = False

    seen_genes = set()
    for gene, hgvsp in iter_ann_blocks(record):
        if not gene or gene in seen_genes:
            continue
        seen_genes.add(gene)

        if not oncokb_done:
            gt = oncokb.get(gene)
            if gt:
                result["OncoKB_Gene_Type"] = gt
                oncokb_done = True

        if not cosmic_gene_done:
            role_tier = gene_role_tier.get(gene)
            if role_tier:
                role, tier = role_tier
                # Number=. field: pass a tuple so pysam writes it as a proper
                # comma-separated multi-value list (e.g. "oncogene,fusion" ->
                # 2 elements), not a single opaque string.
                result["COSMIC_ONC_TSG"] = tuple(role.split(","))
                result["COSMIC_CGC_TIER"] = tier
                cosmic_gene_done = True

                protein_change = parse_protein_change_1letter(hgvsp)
                if protein_change:
                    ref1, pos, alt1 = protein_change
                    mut = mutation_lookup.get((gene, record_chrom, record_pos, ref1, pos, alt1))
                    if mut:
                        cds, aa, genome_pos, mut_tier = mut
                        result["COSMIC_Mutation_CDS"] = cds
                        result["COSMIC_Mutation_AA"] = aa
                        result["COSMIC_Mutation_genome_position_GRCh38"] = genome_pos
                        result["COSMIC_MUTATION_SIGNIFICANCE_TIER"] = mut_tier

        if not hs_v1:
            positions = hotspots_v1.get(gene)
            if positions:
                pos = parse_protein_position(hgvsp)
                if pos and pos in positions:
                    hs_v1 = True
        if not hs_v2:
            positions = hotspots_v2.get(gene)
            if positions:
                pos = parse_protein_position(hgvsp)
                if pos and pos in positions:
                    hs_v2 = True
        if not hs_v3:
            positions = hotspots_v3.get(gene)
            if positions:
                pos = parse_protein_position(hgvsp)
                if pos and pos in positions:
                    hs_v3 = True

        if not eviq and gene in eviq_germline_genes:
            eviq = True

    if hs_v1:
        result["cancerhotspots_v1"] = "T"
    if hs_v2:
        result["cancerhotspots_v2"] = "T"
    if hs_v3:
        result["cancerhotspots_v3"] = "T"
    if eviq:
        result["eviQ_germline"] = "T"
    return result


INFO_HEADER_LINES = [
    ("OncoKB_Gene_Type", 1, "String", "OncoKB Gene Type (ONCOGENE/TSG/ONCOGENE_AND_TSG), present only if the gene is in OncoKB's gene list (gene-symbol match)"),
    ("COSMIC_ONC_TSG", ".", "String", "COSMIC CMC ONC_TSG, present only if the gene has a non-blank role/tier in CMC (i.e. is a Cancer Gene Census-equivalent gene). Compound categories (e.g. oncogene,fusion) are a comma-separated list of 1-3 elements, not a single string - hence Number=."),
    ("COSMIC_CGC_TIER", 1, "String", "COSMIC CMC CGC_TIER, present only if the gene has a non-blank role/tier in CMC"),
    ("COSMIC_Mutation_CDS", 1, "String", "COSMIC CMC Mutation CDS, present only if this variant's CHROM:POS and protein change both match a cataloged CMC mutation for the gene"),
    ("COSMIC_Mutation_AA", 1, "String", "COSMIC CMC Mutation AA, present only if this variant's CHROM:POS and protein change both match a cataloged CMC mutation for the gene"),
    ("COSMIC_Mutation_genome_position_GRCh38", 1, "String", "COSMIC CMC Mutation genome position GRCh38, present only if this variant's CHROM:POS and protein change both match a cataloged CMC mutation for the gene"),
    ("COSMIC_MUTATION_SIGNIFICANCE_TIER", 1, "String", "COSMIC CMC MUTATION_SIGNIFICANCE_TIER (1/2/3/Other - all tiers included), present only if this variant's CHROM:POS and protein change both match a cataloged CMC mutation for the gene"),
    ("cancerhotspots_v1", 1, "String", "T, present only if this variant's gene+codon position matches cancerhotspots.org v1"),
    ("cancerhotspots_v2", 1, "String", "T, present only if this variant's gene+codon position matches cancerhotspots.org v2"),
    ("cancerhotspots_v3", 1, "String", "T, present only if this variant's gene+codon position matches cancerhotspots.org v3"),
    ("eviQ_germline", 1, "String", "T, present only if the gene is on the eviQ germline testing gene list, gene-symbol match only"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--oncokb-gene-list", required=True)
    ap.add_argument("--cosmic-cmc", required=True)
    ap.add_argument("--hotspots-v1", required=True)
    ap.add_argument("--hotspots-v2", required=True)
    ap.add_argument("--hotspots-v3", required=True)
    ap.add_argument("--eviq-germline", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print("Loading OncoKB...")
    oncokb = load_oncokb(args.oncokb_gene_list)
    print("Loading COSMIC CMC (all tiers - this is the slow one)...")
    gene_role_tier, mutation_lookup = load_cosmic_cmc(args.cosmic_cmc)
    print("Loading cancerhotspots v1/v2/v3...")
    hotspots_v1 = load_hotspots_v1(args.hotspots_v1)
    hotspots_v2 = load_hotspots_v2(args.hotspots_v2)
    hotspots_v3 = load_hotspots_v3(args.hotspots_v3)
    print("Loading eviQ germline gene list...")
    eviq_germline_genes = load_eviq_germline_genes(args.eviq_germline)

    vin = pysam.VariantFile(args.vcf)
    for tag, number, type_, desc in INFO_HEADER_LINES:
        vin.header.info.add(tag, number, type_, desc)
    command_line = "annotate_actionable_somatic_germline.py " + " ".join(sys.argv[1:])
    vin.header.add_line(
        f'##annotate_actionable_somatic_germline_command="{command_line}"; '
        f'Date={datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")}'
    )
    vout = pysam.VariantFile(args.out, "w", header=vin.header)

    n_total = 0
    counts = {tag: 0 for tag, *_ in INFO_HEADER_LINES}
    for record in vin:
        n_total += 1
        result = annotate_record(record, oncokb, gene_role_tier, mutation_lookup,
                                  hotspots_v1, hotspots_v2, hotspots_v3, eviq_germline_genes)
        for tag, value in result.items():
            record.info[tag] = value
            counts[tag] += 1
        vout.write(record)

    vin.close()
    vout.close()
    print(f"Total variants: {n_total}")
    for tag, *_ in INFO_HEADER_LINES:
        print(f"{tag}: {counts[tag]}")


if __name__ == "__main__":
    main()
