/*
 * CGen's standard somatic gate (recipe A0-A9 + germline arm) on the full paired VarDict norm VCF, before annotation,
 * so the CGen branch annotates the same call set as CGen's local runs. Script lives in the CGenReporter install.
 */
process CGEN_GATE {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'
    publishDir { "${params.outdir}/${patient_id}/vcfs/cgen_gate" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id.split('_')[0]}_reportable.vcf.gz"), emit: vcf
    path "${case_id.split('_')[0]}_reportable.vcf.gz.tbi"

    script:
    def (t, n) = case_id.split('_')
    """
    set -euo pipefail
    bash ${params.cgen_reporter_dir}/scripts/somatic_gate.sh "${vcf}" "${t}_reportable.vcf.gz" "${t}" "${n}"
    """
}
