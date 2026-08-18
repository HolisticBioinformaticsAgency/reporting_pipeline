// Sequential annotation: SnpEff MANE → ClinVar (SnpSift) → VEP → gnomAD →
// dbNSFP → SpliceAI → OncoKB/COSMIC/cancerhotspots/eviQ

process SNPEFF {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_1'
    publishDir { "${params.outdir}/${patient_id}/vcfs/snpeff" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf), path(mane_db_dir)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_vardict_MANE.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    BASE_CFG="/snpeff_data/snpEff.config"
    CFG="snpeff.mane.config"

    if [[ -s "\$BASE_CFG" ]]; then
        cp "\$BASE_CFG" "\$CFG"
    else
        echo "data.dir = \$PWD" > "\$CFG"
    fi

    {
        echo ""
        echo "data.dir = \$PWD"
        echo "GRCh38.mane.1.2.refseq.genome : Homo sapiens (GRCh38 MANE v1.2 RefSeq)"
    } >> "\$CFG"

    [[ -d "GRCh38.mane.1.2.refseq" ]] || ln -s "${mane_db_dir}" "GRCh38.mane.1.2.refseq"

    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c "${vcf}" | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} \$5 ~ /^</ {next} {print}' > "${case_id}_snpeff_input.vcf"
    else
        awk 'BEGIN{OFS="\\t"} /^#/ {print; next} \$5 ~ /^</ {next} {print}' "${vcf}" > "${case_id}_snpeff_input.vcf"
    fi

    snpEff ann \\
        -c "\$CFG" \\
        -dataDir "\$PWD" \\
        -noStats \\
        "GRCh38.mane.1.2.refseq" \\
        "${case_id}_snpeff_input.vcf" > "${case_id}_vardict_MANE.vcf"
    """
}


process CLINVAR {
    tag "${patient_id}:${case_id}"
    container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'
    publishDir { "${params.outdir}/${patient_id}/vcfs/clinvar" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(mane_vcf), path(clinvar_vcf)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_vardict_MANE_*.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    if [[ "${mane_vcf}" == *.gz ]]; then
        gunzip -c "${mane_vcf}" > "clinvar_input.tmp.vcf"
    else
        cp "${mane_vcf}" "clinvar_input.tmp.vcf"
    fi

    clin_stem=\$(basename "${clinvar_vcf}")
    clin_stem="\${clin_stem%.gz}"
    clin_stem="\${clin_stem%.vcf}"
    OUT_VCF="${case_id}_vardict_MANE_\${clin_stem}.vcf"

    if [[ "${clinvar_vcf}" == *.gz ]]; then
        gunzip -c "${clinvar_vcf}" > "clinvar_db.tmp.vcf"
    else
        cp "${clinvar_vcf}" "clinvar_db.tmp.vcf"
    fi
    SnpSift annotate -sorted "clinvar_db.tmp.vcf" "clinvar_input.tmp.vcf" > "\${OUT_VCF}"

    if ! grep -q '^#CHROM' "\${OUT_VCF}"; then
        echo "[CLINVAR] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "\${OUT_VCF}" >&2
        exit 1
    fi
    """
}


process VEP {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/ensembl-vep:115--pl5321h2a3209d_0'
    publishDir { "${params.outdir}/${patient_id}/vcfs/vep" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf)

    output:
    tuple val(patient_id), val(case_id), val(mode),
          path("${case_id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}.vcf"),
          emit: vcf
    path("${case_id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}_stats.txt"),
          emit: stats

    script:
    def out_vcf   = "${case_id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}.vcf"
    def out_stats = "${case_id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}_stats.txt"
    """
    set -euo pipefail

    vep \\
        --offline --cache \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --dir_cache ${params.vep_cache_dir} \\
        --fasta ${params.vep_fasta} \\
        --vcf --force_overwrite --format vcf \\
        --buffer_size 50 \\
        --symbol --numbers --mane --hgvs --total_length \\
        --canonical --biotype --protein --ccds --domains --uniprot \\
        --af --af_gnomad --max_af --variant_class \\
        --pick \\
        --input_file ${vcf} \\
        --output_file ${out_vcf} \\
        --stats_text --stats_file ${out_stats}
    """
}


process GNOMAD {
    tag "${patient_id}:${case_id}"
    container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'
    publishDir { "${params.outdir}/${patient_id}/vcfs/gnomad" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf), path(gnomad_vcf)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_gnomad.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    SnpSift annotate -name ${params.gnomad_info_prefix} "${gnomad_vcf}" "${vcf}" > "${case_id}_gnomad.vcf"

    if ! grep -q '^#CHROM' "${case_id}_gnomad.vcf"; then
        echo "[GNOMAD] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "${case_id}_gnomad.vcf" >&2
        exit 1
    fi
    """
}


process DBNSFP {
    tag "${patient_id}:${case_id}"
    container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'
    publishDir { "${params.outdir}/${patient_id}/vcfs/dbnsfp" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf)
    tuple path(dbnsfp_db), path(dbnsfp_tbi)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_dbNSFP.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    # dbNSFP is ~50GB - SnpSift's container-baked default heap (-Xmx1g) OOMs
    # partway through without surfacing an error, silently truncating the
    # variant set. 12g heap (needs >=16G task memory - see nextflow.config
    # withName:'DBNSFP') is what actually got through the full file previously.
    export _JAVA_OPTIONS="-Xmx12g"

    SnpSift dbnsfp -db "${dbnsfp_db}" -f ${params.dbnsfp_fields} \\
        "${vcf}" > "${case_id}_dbNSFP.vcf.new"

    # SnpSift dbnsfp writes its whole header block a second time at the very
    # end of the output (trailing, spurious duplicate #CHROM line, no data
    # after it) - strip it, keeping the real header once plus every
    # non-header data line in original order.
    awk '
      seen==0 { print; if (/^#CHROM/) seen=1; next }
      !/^#/ { print }
    ' "${case_id}_dbNSFP.vcf.new" > "${case_id}_dbNSFP.vcf"
    rm -f "${case_id}_dbNSFP.vcf.new"

    IN_COUNT=\$(grep -vc "^#" "${vcf}")
    OUT_COUNT=\$(grep -vc "^#" "${case_id}_dbNSFP.vcf")
    if [ "\$IN_COUNT" != "\$OUT_COUNT" ]; then
        echo "[DBNSFP] ERROR: variant count changed (\$IN_COUNT -> \$OUT_COUNT) - dbNSFP annotation likely truncated" >&2
        exit 1
    fi
    """
}


process SPLICEAI {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'
    publishDir { "${params.outdir}/${patient_id}/vcfs/spliceai" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf)
    tuple path(spliceai_vcf), path(spliceai_tbi)

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_spliceAI.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    bgzip -f -c "${vcf}" > "${case_id}_input.vcf.gz"
    tabix -f -p vcf "${case_id}_input.vcf.gz"

    bcftools annotate -a "${spliceai_vcf}" -c INFO/SpliceAI \\
        "${case_id}_input.vcf.gz" \\
        -o "${case_id}_spliceAI.vcf"

    if ! grep -q '^#CHROM' "${case_id}_spliceAI.vcf"; then
        echo "[SPLICEAI] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "${case_id}_spliceAI.vcf" >&2
        exit 1
    fi
    """
}


process ONCOKB_COSMIC_HOTSPOTS_EVIQ {
    tag "${patient_id}:${case_id}"
    publishDir { "${params.outdir}/${patient_id}/vcfs/oncorecur" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf)
    path oncokb_gene_list
    path cosmic_cmc
    path hotspots_v1
    path hotspots_v2
    path hotspots_v3
    path eviq_germline

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_oncorecur.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    /apps/anaconda/2024.02-1/bin/python3 annotate_actionable_somatic_germline.py \\
        --vcf "${vcf}" \\
        --oncokb-gene-list "${oncokb_gene_list}" \\
        --cosmic-cmc "${cosmic_cmc}" \\
        --hotspots-v1 "${hotspots_v1}" \\
        --hotspots-v2 "${hotspots_v2}" \\
        --hotspots-v3 "${hotspots_v3}" \\
        --eviq-germline "${eviq_germline}" \\
        --out "${case_id}_oncorecur.vcf"

    if ! grep -q '^#CHROM' "${case_id}_oncorecur.vcf"; then
        echo "[ONCOKB_COSMIC_HOTSPOTS_EVIQ] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "${case_id}_oncorecur.vcf" >&2
        exit 1
    fi
    """
}


workflow VCF_ANNOTATE {
    take:
    ch_vcf          // (patient_id, case_id, mode, vcf)
    ch_mane_db      // value channel: path to MANE DB dir
    ch_clinvar      // value channel: path to ClinVar VCF
    ch_gnomad       // value channel: path to gnomAD VCF
    ch_dbnsfp       // value channel: tuple(dbNSFP db, dbNSFP .tbi)
    ch_spliceai     // value channel: tuple(SpliceAI VCF, SpliceAI .tbi)
    ch_oncokb       // value channel: path to OncoKB gene list
    ch_cosmic       // value channel: path to COSMIC CMC
    ch_hotspots_v1  // value channel: path to cancerhotspots v1
    ch_hotspots_v2  // value channel: path to cancerhotspots v2
    ch_hotspots_v3  // value channel: path to cancerhotspots v3
    ch_eviq         // value channel: path to eviQ germline gene list

    main:
    SNPEFF(   ch_vcf.combine(ch_mane_db) )
    CLINVAR(  SNPEFF.out.vcf.combine(ch_clinvar) )
    VEP(      CLINVAR.out.vcf )
    GNOMAD(   VEP.out.vcf.combine(ch_gnomad) )
    DBNSFP(   GNOMAD.out.vcf, ch_dbnsfp )
    SPLICEAI( DBNSFP.out.vcf, ch_spliceai )
    ONCOKB_COSMIC_HOTSPOTS_EVIQ(
        SPLICEAI.out.vcf,
        ch_oncokb, ch_cosmic, ch_hotspots_v1, ch_hotspots_v2, ch_hotspots_v3, ch_eviq
    )

    emit:
    vcf   = ONCOKB_COSMIC_HOTSPOTS_EVIQ.out.vcf
    stats = VEP.out.stats
}
