process GENERATE_REPORT {
    tag "${patient_id}:${case_id}"
    container 'python:3.11-slim'
    publishDir { "${params.outdir}/${patient_id}/report" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), path(germline_vcf), path(somatic_vcf), path(cns), path(cnr_genemetrics), path(msi), path(tmb)
    path marker_csv
    path chrom_arm_bed

    output:
    tuple val(patient_id), val(case_id), path("${patient_id}_report.md"), emit: report

    script:
    def mane_db_name = params.mane_db_dir  ? new File(params.mane_db_dir  as String).name : ''
    def clinvar_name = params.clinvar_vcf  ? new File(params.clinvar_vcf  as String).name : ''
    """
    generate_report.py \\
        --patient-id    "${patient_id}" \\
        --case-id       "${case_id}" \\
        --germline-vcf  "${germline_vcf}" \\
        --somatic-vcf   "${somatic_vcf}" \\
        --cns             "${cns}" \\
        --cnr-genemetrics "${cnr_genemetrics}" \\
        --msi             "${msi}" \\
        --tmb           "${tmb}" \\
        --markers       "${marker_csv}" \\
        --chrom-arm     "${chrom_arm_bed}" \\
        --mane-db-name  "${mane_db_name}" \\
        --clinvar-name  "${clinvar_name}" \\
        --vep-species   "${params.vep_species}" \\
        --vep-assembly  "${params.vep_assembly}" \\
        --vep-version   "${params.vep_version}" \\
        --chrom-frac-thr ${params.chromosomal_percentage_decision} \\
        --arm-frac-thr   ${params.arm_percentage_decision} \\
        --log2-thr       ${params.cnv_log2_threshold} \\
        --out           "${patient_id}_report.md"
    """
}
