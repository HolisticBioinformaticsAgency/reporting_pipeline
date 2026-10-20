process GENERATE_REPORT {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0'
    publishDir { "${params.outdir}/${patient_id}/report" }, mode: 'copy', pattern: '*_scatter_manifest.tsv'

    input:
    tuple val(patient_id), val(case_id), path(germline_vcf), path(somatic_vcf), path(hotspot_counts), path(facets_purity), path(facets_ploidy), path(cns), path(msi), path(tmb), path(tumour_hla), path(germline_hla)
    path marker_csv
    path chrom_arm_bed
    path patient_metadata
    path whitelist_csv
    path supplementary_json

    output:
    tuple val(patient_id), val(case_id), path("${patient_id}_report.draft.md"),      emit: report
    tuple val(patient_id), val(case_id), path("${patient_id}_scatter_manifest.tsv"), emit: manifest

    script:
    // hotspot_counts is a placeholder (assets/NO_WHITELIST) when no whitelist
    // was configured - drop the flag entirely so the report renders the
    // "no whitelist included" note instead of a table.
    def hotspot_arg = hotspot_counts.name == 'NO_WHITELIST' ? '' : "--hotspot-counts \"${hotspot_counts}\""
    """
    generate_report.py \\
        --patient-id    "${patient_id}" \\
        --case-id       "${case_id}" \\
        --germline-vcf  "${germline_vcf}" \\
        --somatic-vcf   "${somatic_vcf}" \\
        ${hotspot_arg} \\
        --facets-purity   "${facets_purity}" \\
        --facets-ploidy   "${facets_ploidy}" \\
        --cns             "${cns}" \\
        --msi             "${msi}" \\
        --tmb           "${tmb}" \\
        --tumour-hla    "${tumour_hla}" \\
        --germline-hla  "${germline_hla}" \\
        --markers       "${marker_csv}" \\
        --chrom-arm     "${chrom_arm_bed}" \\
        --patient-metadata "${patient_metadata}" \\
        --whitelist     "${whitelist_csv}" \\
        --supplementary-json "${supplementary_json}" \\
        --pct-alt-thr    ${params.cn_pct_alt_threshold} \\
        --amp-seg-cn     ${params.amp_segment_cn} \\
        --del-seg-cn     ${params.del_segment_cn} \\
        --out             "${patient_id}_report.draft.md" \\
        --scatter-manifest "${patient_id}_scatter_manifest.tsv"
    """
}
