// High-confidence variant selection, applied before annotation.
//   1. Split off calls that overlap repeat regions (union of microsatellite,
//      RepeatMasker, segmental duplication, simple tandem repeat tracks) and
//      rescue the top `repetitive_regions_rescue` fraction of those by tumour
//      mapping quality (FMT/MQ, sample 0) — a uniquely, cleanly-mapped read
//      set in a repeat region is still trustworthy evidence.
//   2. Apply per-STATUS depth/AF cutoffs (paired mode only — these compare
//      tumour vs normal DP/AF, which isn't meaningful for a single-sample call).

process VCF_FILTER {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'
    publishDir { "${params.outdir}/${patient_id}/vcfs/highconf" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf)
    path repeat_bed

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_highconf.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    NONREPEAT_VCF="nonrepeat.vcf.gz"
    REPEAT_VCF="repeat.vcf.gz"
    RESCUE_VCF="repeat_rescued.vcf.gz"
    PREFILTER_VCF="prefilter.vcf.gz"

    # ── Step 1: split repeat-region calls from everything else ─────────────
    bcftools view -T '^${repeat_bed}' -Oz -o "\$NONREPEAT_VCF" "${vcf}"
    bcftools view -T '${repeat_bed}'  -Oz -o "\$REPEAT_VCF"    "${vcf}"
    bcftools index -t -f "\$NONREPEAT_VCF"

    # ── Step 2: rescue the top ${params.repetitive_regions_rescue} fraction
    #    of repeat-region calls by tumour mapping quality (FMT/MQ[0]) ───────
    N_REPEAT=\$(bcftools view -H "\$REPEAT_VCF" | wc -l)
    if [ "\$N_REPEAT" -gt 0 ]; then
        MQ_THRESH=\$(bcftools query -f '[%MQ\\t]\\n' "\$REPEAT_VCF" \\
            | awk '{print \$1}' \\
            | sort -n \\
            | awk -v frac=${params.repetitive_regions_rescue} '{a[NR]=\$1} END{
                  idx = int(NR * (1 - frac));
                  if (idx < 1) idx = 1;
                  print a[idx]
              }')
        bcftools view -i "FMT/MQ[0]>=\$MQ_THRESH" -Oz -o "\$RESCUE_VCF" "\$REPEAT_VCF"
    else
        MQ_THRESH="NA"
        cp "\$REPEAT_VCF" "\$RESCUE_VCF"
    fi
    bcftools index -t -f "\$RESCUE_VCF"

    bcftools concat -a "\$NONREPEAT_VCF" "\$RESCUE_VCF" -Oz -o "\$PREFILTER_VCF"
    bcftools index -t -f "\$PREFILTER_VCF"

    echo "[VCF_FILTER] repeat-region calls: \$N_REPEAT | MQ rescue threshold (top ${params.repetitive_regions_rescue} fraction): \$MQ_THRESH" >&2

    # ── Step 3: per-STATUS depth/AF cutoffs (paired mode only) ─────────────
    if [ "${mode}" = "paired" ]; then
        CUTOFF_FILTER='
           (INFO/STATUS="LikelySomatic" && FMT/DP[0]>=${params.min_likelysomatic_tumour_depth} && FMT/DP[1]>=${params.min_likelysomatic_normal_depth} && FMT/AF[0:0]>=${params.min_universal_tumour_af})
        || (INFO/STATUS="StrongSomatic" && FMT/DP[0]>=${params.min_strongsomatic_tumour_depth} && FMT/DP[1]>=${params.min_strongsomatic_normal_depth} && FMT/AF[0:0]>=${params.min_universal_tumour_af})
        || (INFO/STATUS="Germline"      && FMT/DP[0]>=${params.min_germline_tumour_depth} && FMT/DP[1]>=${params.min_germline_normal_depth} && FMT/AF[1:0]>=${params.min_universal_normal_af})
        || (INFO/STATUS="StrongLOH"     && FMT/DP[0]>=${params.min_strongloh_tumour_depth} && FMT/DP[1]>=${params.min_strongloh_normal_depth} && FMT/AF[1:0]>=${params.min_universal_normal_af})
        || (INFO/STATUS="LikelyLOH" && FMT/AF[0:0]<0.05  && FMT/DP[0]>=${params.min_likelyloh005_tumour_depth} && FMT/DP[1]>=${params.min_likelyloh005_normal_depth} && FMT/AF[1:0]>=${params.min_universal_normal_af})
        || (INFO/STATUS="LikelyLOH" && FMT/AF[0:0]>0.95  && FMT/DP[0]>=${params.min_likelyloh095_tumour_depth} && FMT/DP[1]>=${params.min_likelyloh095_normal_depth} && FMT/AF[1:0]>=${params.min_universal_normal_af})
        '
        bcftools view -i "\$CUTOFF_FILTER" -o "${case_id}_highconf.vcf" "\$PREFILTER_VCF"
    else
        bcftools view -o "${case_id}_highconf.vcf" "\$PREFILTER_VCF"
    fi

    if ! grep -q '^#CHROM' "${case_id}_highconf.vcf"; then
        echo "[VCF_FILTER] ERROR: output VCF missing #CHROM header" >&2
        head -50 "${case_id}_highconf.vcf" >&2
        exit 1
    fi
    """
}
