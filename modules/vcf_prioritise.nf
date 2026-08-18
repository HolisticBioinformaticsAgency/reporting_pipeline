// Prioritise reportable variants from the fully-annotated VCF, via bcftools:
//   - germline: STATUS=Germline, ClinVar P/LP with a high-confidence review
//     status, and on the eviQ germline testing gene list
//   - somatic: STATUS!=Germline, and either OncoKB-actionable + COSMIC CGC
//     tier 1, a COSMIC CMC significance-tiered mutation, a cancerhotspots
//     match (any version), OR a WHO CNS5 whitelist match (cns5_whitelist.csv)
//
// CLNREVSTAT and gene-bounded ANN matching both had real bcftools gotchas
// caught empirically before shipping this (see build_whitelist_filter.py's
// own docstring for the ANN gene-bounding one):
//   - CLNREVSTAT is Number=. in the VCF header, so bcftools comma-splits it
//     into an array. A pattern containing commas itself (e.g. the literal
//     ClinVar string "criteria_provided,_multiple_submitters,_no_conflicts")
//     never matches any single split token - always a silent false negative.
//     Fixed by matching on the unique post-split token ("no_conflicts")
//     instead of the full comma-containing string.
//
// Building the whitelist filter expression needs Python (csv module, quote-
// aware parsing of cns5_whitelist.csv's notes column) but applying it needs
// bcftools - kept as two processes rather than one, since mixing toolchains
// in a single containerized process doesn't work (the bcftools container has
// no python3, and /apps isn't bind-mounted into it per nextflow.config).

process BUILD_WHITELIST_FILTER {
    tag "whitelist"

    input:
    path whitelist_csv

    output:
    path "whitelist_filter.txt", emit: filter

    script:
    """
    set -euo pipefail
    /apps/anaconda/2024.02-1/bin/python3 build_whitelist_filter.py "${whitelist_csv}" > whitelist_filter.txt
    """
}


process PRIORITISE_FILTER {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'
    publishDir { "${params.outdir}/${patient_id}/vcfs/prioritised" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf), path(whitelist_filter_txt)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_prioritised_germline.vcf"), emit: germline
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_prioritised_somatic.vcf"),  emit: somatic

    script:
    """
    set -euo pipefail

    WHITELIST_FILTER=\$(cat "${whitelist_filter_txt}")

    GERMLINE_FILTER='
        (INFO/STATUS="Germline") &&
        (INFO/CLNSIG="Pathogenic" || INFO/CLNSIG="Likely_pathogenic" || INFO/CLNSIG="Pathogenic/Likely_pathogenic" || INFO/CLNSIG="Pathogenic/Likely_pathogenic|other") &&
        (INFO/CLNREVSTAT~"no_conflicts" || INFO/CLNREVSTAT~"reviewed_by_expert_panel") &&
        INFO/eviQ_germline="T"
    '

    SOMATIC_FILTER="
        INFO/STATUS!=\\"Germline\\" &&
        ((INFO/OncoKB_Gene_Type!=\\".\\" && INFO/OncoKB_Gene_Type!=\\"INSUFFICIENT_EVIDENCE\\" && INFO/COSMIC_CGC_TIER=\\"1\\")
            || INFO/COSMIC_MUTATION_SIGNIFICANCE_TIER!=\\".\\"
            || INFO/cancerhotspots_v1=\\"T\\" || INFO/cancerhotspots_v2=\\"T\\" || INFO/cancerhotspots_v3=\\"T\\"
            || (\$WHITELIST_FILTER))
    "

    bcftools view -i "\$GERMLINE_FILTER" -o "${case_id}_prioritised_germline.vcf" "${vcf}"
    bcftools view -i "\$SOMATIC_FILTER"  -o "${case_id}_prioritised_somatic.vcf"  "${vcf}"

    for f in "${case_id}_prioritised_germline.vcf" "${case_id}_prioritised_somatic.vcf"; do
        if ! grep -q '^#CHROM' "\$f"; then
            echo "[PRIORITISE_FILTER] ERROR: \$f missing #CHROM header" >&2
            head -50 "\$f" >&2
            exit 1
        fi
    done
    """
}


workflow VCF_PRIORITISE {
    take:
    ch_vcf        // (patient_id, case_id, mode, vcf)
    ch_whitelist  // value channel: path to cns5_whitelist.csv

    main:
    BUILD_WHITELIST_FILTER( ch_whitelist )
    PRIORITISE_FILTER( ch_vcf.combine( BUILD_WHITELIST_FILTER.out.filter ) )

    emit:
    germline = PRIORITISE_FILTER.out.germline
    somatic  = PRIORITISE_FILTER.out.somatic
}
