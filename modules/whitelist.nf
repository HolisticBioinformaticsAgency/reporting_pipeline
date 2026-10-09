// For each hotspot in whitelist.csv, pull tumour and germline REF/ALT read
// depths straight from the BAMs via `bcftools mpileup -r <all hotspots at
// once>` piped into `bcftools query -f '%CHROM\t%POS\t%ALT\t[%AD]\n'`, rather
// than from a called VCF. `bcftools mpileup`'s AD field is (REF depth,
// depth-per-observed-ALT-in-column-order), computed directly from the
// pileup - so unlike a VCF lookup, it always reports a REF depth even when
// no variant was called there. That's what disambiguates a genuinely
// wildtype hotspot (REF depth high, ALT depth 0) from one with no coverage
// at all (REF depth 0 too) - no separate raw-depth columns needed.
//
// whitelist.csv's REF/ALT are curated to genomic (+strand) convention by
// hand. GENE/STRAND/TRANSCRIPT/CDNA/PROTEIN_OR_NAME/CONTEXT are read for
// context but only GENE and a CDNA+PROTEIN "hotspot" label are carried into
// the output - see /fs04/scratch2/vh83/liem/areas/whitelist/whitelist.csv
// for the source data.
//
// If a hotspot's specific ALT allele was never observed in the pileup at
// all (no reads support it), its ALT depth is 0 - not dropped.

process WHITELIST {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'
    publishDir { "${params.outdir}/${patient_id}/whitelist" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), path(tumour_bam), path(tumour_bai), path(germline_bam), path(germline_bai)
    path whitelist_csv
    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(patient_id), val(case_id), path("${case_id}_hotspot_counts.tsv"), emit: counts

    script:
    """
    set -euo pipefail

    # One comma-joined region list ("chr:pos-pos,chr:pos-pos,...") covering
    # every hotspot, so each BAM only needs a single mpileup call.
    REGIONS=\$(tail -n +2 "${whitelist_csv}" | awk -F',' '\$1!="" {printf "%s%s:%s-%s", (NR>1?",":""), \$1, \$2, \$2}')

    # bcftools query extracts CHROM/POS/ALT/AD as clean tab fields directly
    # from the mpileup stream - no manual FORMAT-column parsing needed.
    bcftools mpileup -r "\$REGIONS" --fasta-ref "${fasta}" "${tumour_bam}"   2>/dev/null \\
        | bcftools query -f '%CHROM\\t%POS\\t%ALT\\t[%AD]\\n' - > tumour_pileup.tsv
    bcftools mpileup -r "\$REGIONS" --fasta-ref "${fasta}" "${germline_bam}" 2>/dev/null \\
        | bcftools query -f '%CHROM\\t%POS\\t%ALT\\t[%AD]\\n' - > germline_pileup.tsv

    printf 'GENE\\tCHROM\\tPOS\\tHOTSPOT\\tTUMOUR_REF_DEPTH\\tTUMOUR_ALT_DEPTH\\tGERMLINE_REF_DEPTH\\tGERMLINE_ALT_DEPTH\\n' > "${case_id}_hotspot_counts.tsv"

    tail -n +2 "${whitelist_csv}" | while IFS=, read -r CHROM POS REF ALT GENE STRAND TRANSCRIPT CDNA PROTEIN CONTEXT; do
        [[ -z "\$CHROM" ]] && continue

        T_LINE=\$(awk -F'\\t' -v c="\$CHROM" -v p="\$POS" '\$1==c && \$2==p {print; exit}' tumour_pileup.tsv)
        N_LINE=\$(awk -F'\\t' -v c="\$CHROM" -v p="\$POS" '\$1==c && \$2==p {print; exit}' germline_pileup.tsv)

        # AD is REF depth followed by one depth per comma-listed ALT allele,
        # in the same order as the ALT (\$3) column.
        if [ -n "\$T_LINE" ]; then
            read -r T_REF T_ALT <<< "\$(printf '%s' "\$T_LINE" | awk -F'\\t' -v alt="\$ALT" '
                {
                    n = split(\$3, alts, ",");
                    split(\$4, ad, ",");
                    alt_d = 0;
                    for (i = 1; i <= n; i++) if (alts[i] == alt) alt_d = ad[i+1] + 0;
                    print ad[1] + 0, alt_d;
                }')"
        else
            T_REF=0; T_ALT=0
        fi

        if [ -n "\$N_LINE" ]; then
            read -r N_REF N_ALT <<< "\$(printf '%s' "\$N_LINE" | awk -F'\\t' -v alt="\$ALT" '
                {
                    n = split(\$3, alts, ",");
                    split(\$4, ad, ",");
                    alt_d = 0;
                    for (i = 1; i <= n; i++) if (alts[i] == alt) alt_d = ad[i+1] + 0;
                    print ad[1] + 0, alt_d;
                }')"
        else
            N_REF=0; N_ALT=0
        fi

        HOTSPOT="\${CDNA} \${PROTEIN}"

        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
            "\$GENE" "\$CHROM" "\$POS" "\$HOTSPOT" "\$T_REF" "\$T_ALT" "\$N_REF" "\$N_ALT" >> "${case_id}_hotspot_counts.tsv"
    done
    """
}
