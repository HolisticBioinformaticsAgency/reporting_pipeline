#!/usr/bin/env python3
"""Build a bcftools -i expression from cns5_whitelist.csv, OR-ing together
one clause per usable row. Printed to stdout for a Nextflow process to
capture into a shell variable.

Only rows where confidence != "needs_verification" and match_type is one of
exact_position / codon_hotspot / gene_any_variant are used - cnv_arm_level,
wildtype_status, and pathway_level rows aren't representable as a VCF variant
filter (see cns5_whitelist.csv's own notes column per row).

Gene matching uses unescaped "|GENE|" bounding against SnpEff's ANN field
(Allele|Effect|Impact|Gene_Name|...) - empirically confirmed against real
annotated VCF data that bcftools' regex engine treats a *backslash-escaped*
pipe ("\\|GENE\\|") as matching every record (always-true), while the bare,
unescaped "|GENE|" matches correctly and rejects substring false positives
(e.g. "FAP126" wrongly matching "CFAP126"). Do not "fix" this by escaping.
"""
import csv
import sys

USABLE_MATCH_TYPES = {"exact_position", "codon_hotspot", "gene_any_variant"}


def escape_hgvsp_dot(pattern):
    """Escape only the literal 'p.' prefix dot; leave any trailing wildcard
    (e.g. the '.*' in 'p.Arg132.*') as a real regex wildcard."""
    if pattern.startswith("p."):
        return "p\\." + pattern[2:]
    return pattern


def build_clauses(csv_path):
    clauses = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("confidence", "").strip() == "needs_verification":
                continue
            match_type = row.get("match_type", "").strip()
            if match_type not in USABLE_MATCH_TYPES:
                continue

            gene      = row.get("gene_symbol", "").strip()
            chrom     = row.get("chrom", "").strip()
            pos_start = row.get("pos_start", "").strip()
            pos_end   = row.get("pos_end", "").strip()
            hgvsp     = row.get("hgvsp_pattern", "").strip()

            if match_type == "exact_position" and chrom and pos_start and pos_end:
                clauses.append(f'(CHROM="{chrom}" && POS>={pos_start} && POS<={pos_end})')
            elif match_type == "codon_hotspot" and gene and hgvsp:
                clauses.append(
                    f'(INFO/ANN~"|{gene}|" && INFO/ANN~"{escape_hgvsp_dot(hgvsp)}")'
                )
            elif match_type == "gene_any_variant" and gene:
                clauses.append(f'INFO/ANN~"|{gene}|"')

    return clauses


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} cns5_whitelist.csv")
    clauses = build_clauses(sys.argv[1])
    if not clauses:
        sys.exit("ERROR: no usable whitelist rows found")
    print(" || ".join(clauses))


if __name__ == "__main__":
    main()
