process ADD_SCATTER_TO_REPORT {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0'
    publishDir { "${params.outdir}/${patient_id}/report" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), path(draft_report), path(manifest), path(scatter_dir)

    output:
    tuple val(patient_id), val(case_id), path("${patient_id}_report.md"), emit: report

    script:
    """
    add_scatter_to_report.py \\
        --report      "${draft_report}" \\
        --manifest    "${manifest}" \\
        --scatter-dir "${scatter_dir}" \\
        --out         "${patient_id}_report.md"
    """
}
