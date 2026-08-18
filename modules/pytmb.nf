// TMB needs the full somatic mutation burden, not the curated
// actionable/hotspot set that VCF_PRIORITISE now produces for reporting -
// so this consumes the fully-annotated VCF directly (same input as
// VCF_PRIORITISE), restricted only to non-Germline STATUS calls, and lets
// pyTMB's own --minDepth/--minAltDepth/--vaf/--filterSyn/--filterNonCoding/
// --filterLowQual flags do the quality filtering.

process SOMATIC_FILTER_FOR_TMB {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'
    publishDir { "${params.outdir}/${patient_id}/vcfs/tmb_somatic" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_somatic_for_tmb.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    if [ "${mode}" = "paired" ]; then
        bcftools view -i 'INFO/STATUS!="Germline"' -o "${case_id}_somatic_for_tmb.vcf" "${vcf}"
    else
        bcftools view -o "${case_id}_somatic_for_tmb.vcf" "${vcf}"
    fi

    if ! grep -q '^#CHROM' "${case_id}_somatic_for_tmb.vcf"; then
        echo "[SOMATIC_FILTER_FOR_TMB] ERROR: output VCF missing #CHROM header" >&2
        head -50 "${case_id}_somatic_for_tmb.vcf" >&2
        exit 1
    fi
    """
}


process RUN_PYTMB {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/tmb:1.5.0--pyhdfd78af_1'
    publishDir { "${params.outdir}/${patient_id}/tmb" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(tumor_label), path(vcf)
    path db_config
    path var_config

    output:
    tuple val(patient_id), val(case_id), path("${case_id}_pytmb.tsv"), emit: tmb

    script:
    """
    set -euo pipefail

    cp ${db_config} db.yml
    cp ${var_config} var.yml

    PLAIN_VCF="${case_id}_input.vcf"

    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c "${vcf}" > "\${PLAIN_VCF}"
    else
        cp "${vcf}" "\${PLAIN_VCF}"
    fi

    if ! grep -q '^##' "\${PLAIN_VCF}"; then
        echo "ERROR: VCF has no header lines" >&2
        exit 1
    fi

    VARIANT_COUNT=\$(grep -v "^#" "\${PLAIN_VCF}" | wc -l)
    echo "INFO: \${VARIANT_COUNT} variants" >&2

    if [ "\${VARIANT_COUNT}" -eq 0 ]; then
        echo "WARNING: No variants — writing empty TMB output" >&2
        printf "sample\\ttmb\\ttotal_variants\\tcoding_variants\\n" > "${case_id}_pytmb.tsv"
        printf "%s\\t0\\t0\\t0\\n" "${tumor_label}"                 >> "${case_id}_pytmb.tsv"
    else
        pyTMB.py \\
            -i "\${PLAIN_VCF}" \\
            --bed ${params.pytmb_capture_bed} \\
            --sample "${tumor_label}" \\
            --dbConfig db.yml \\
            --varConfig var.yml \\
            --vaf ${params.pytmb_vaf} \\
            --minDepth ${params.pytmb_min_depth} \\
            --minAltDepth ${params.pytmb_min_alt} \\
            ${params.pytmb_filter_syn ? "--filterSyn" : ""} \\
            ${params.pytmb_filter_nc  ? "--filterNonCoding" : ""} \\
            ${params.pytmb_filter_lq  ? "--filterLowQual" : ""} \\
            > "${case_id}_pytmb.tsv"
    fi
    """
}


workflow PYTMB {
    take:
    ch_vcf          // (patient_id, case_id, mode, vcf) - fully-annotated VCF
    ch_tumor_label  // (patient_id, case_id, tumor_label)
    ch_db_config    // value channel: path to pyTMB db config
    ch_var_config   // value channel: path to pyTMB var config

    main:
    SOMATIC_FILTER_FOR_TMB( ch_vcf )

    ch_tmb_in = SOMATIC_FILTER_FOR_TMB.out.vcf
        .map { pid, cid, mode, vcf -> tuple(pid, cid, vcf) }
        .join( ch_tumor_label, by: [0, 1] )
        .map { pid, cid, vcf, tlabel -> tuple(pid, cid, tlabel, vcf) }

    RUN_PYTMB( ch_tmb_in, ch_db_config, ch_var_config )

    emit:
    tmb = RUN_PYTMB.out.tmb
}
