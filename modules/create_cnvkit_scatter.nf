process CREATE_CNVKIT_SCATTER {
    tag "${patient_id}:${case_id}"
    // Not a CNVkit call: bin/cnv_absolute_scatter.py is our own standalone
    // script, reading the .cnr/.cns directly. This container is reused only
    // as a ready-made numpy+matplotlib Python runtime (both are CNVkit's own
    // dependencies) - no cnvkit.py command, internal cnvlib import, or file
    // override is involved anywhere below.
    container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
    publishDir { "${params.outdir}/${patient_id}/report" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), path(cnr), path(cns), path(manifest), path(vcf), path(vcf_tbi), path(facets_purity), path(facets_ploidy)
    path chrom_arm_bed

    output:
    tuple val(patient_id), val(case_id), path("scatter"), emit: scatter_dir

    script:
    // case_id is "<tumour_id>_<normal_id>" (e.g. M24007T_M24007N) - same
    // convention main.nf itself uses to split tumor_label/normal_label.
    def mainId = case_id.split('_')[0]
    def ctrlId = case_id.split('_')[1]
    """
    set -euo pipefail
    export MPLCONFIGDIR='/tmp'
    export MPLBACKEND='Agg'
    export XDG_CACHE_HOME='/tmp'
    export HOME='/tmp'

    mkdir -p scatter

    # Purity/ploidy: same fallback as the .cns this module was actually
    # handed (main.nf prefers the purity+ploidy call.cns; falls back to the
    # plain one, which was called with neither correction, when FACETS
    # purity is unavailable) - so mirror that here: purity=1.0/ploidy=2
    # (a no-op in the absolute-CN formula) unless a usable purity exists.
    PURITY_RAW=\$(tr -d '[:space:]' < "${facets_purity}" 2>/dev/null || echo NA)
    if echo "\${PURITY_RAW}" | grep -qE '^[0-9]*\\.?[0-9]+\$' \\
       && awk -v p="\${PURITY_RAW}" 'BEGIN{exit !(p > 0 && p <= 1)}'; then
        PURITY="\${PURITY_RAW}"
        PLOIDY_RAW=\$(tr -d '[:space:]' < "${facets_ploidy}" 2>/dev/null || echo NA)
        if echo "\${PLOIDY_RAW}" | grep -qE '^[0-9]*\\.?[0-9]+\$' \\
           && awk -v p="\${PLOIDY_RAW}" 'BEGIN{exit !(p > 0)}'; then
            PLOIDY=\$(awk -v p="\${PLOIDY_RAW}" 'BEGIN{printf "%.0f", p}')
        else
            PLOIDY=2
        fi
    else
        PURITY=1.0
        PLOIDY=2
    fi
    echo "[CREATE_CNVKIT_SCATTER] ${case_id}: purity=\${PURITY} ploidy=\${PLOIDY}"

    tail -n +2 "${manifest}" | while IFS=\$'\\t' read -r key type region genes target_start target_end; do
        [ -z "\${key}" ] && continue
        TARGET_ARGS=""
        if [ -n "\${target_start}" ] && [ -n "\${target_end}" ]; then
            TARGET_ARGS="--target-start \${target_start} --target-end \${target_end}"
        fi
        cnv_absolute_scatter.py \\
            --cnr "${cnr}" --cns "${cns}" --vcf "${vcf}" \\
            --sample-id ${mainId} --normal-id ${ctrlId} \\
            --zygosity-freq 0.35 --min-variant-depth 30 \\
            --region "\${region}" --genes "\${genes}" \\
            --purity "\${PURITY}" --ploidy "\${PLOIDY}" \\
            --y-max ${params.scatter_y_max} \\
            --chrom-arm-bed "${chrom_arm_bed}" \\
            \${TARGET_ARGS} \\
            -o "scatter/\${key}.png" || true
    done
    """
}
