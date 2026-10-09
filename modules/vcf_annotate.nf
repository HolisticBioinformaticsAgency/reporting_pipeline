// Sequential annotation: SnpEff MANE → ClinVar (SnpSift) → VEP → gnomAD →
// dbNSFP → SpliceAI → OncoKB/COSMIC/cancerhotspots/eviQ
//
// Each step's output filename is built by appending its own suffix onto
// whatever filename came in, so by the end of the chain the file carries its
// full annotation history, e.g.:
//   M24087T_M24087N_vardict_norm_hiconf_MANE_clinvar_20260615_VEP_gnomad_dbNSFP_spliceAI_oncorecur.vcf
// (the "_hiconf" link comes from VCF_FILTER, upstream of this workflow).
//
// NOTE: `output:` filename expressions can only reference process INPUT
// variables (and params/globals) - NOT `def`s created inside `script:`.
// Those two blocks are evaluated separately; a script:-scoped def is
// invisible to output:. So the filename expression is written out fully in
// both places (output: for Nextflow's own bookkeeping, script: for the
// shell body) rather than computed once and shared.

process SNPEFF {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_1'
    publishDir { "${params.outdir}/${patient_id}/vcfs/snpeff" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf), path(mane_db_dir)

    output:
    tuple val(patient_id), val(case_id), val(mode),
          path("${vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_MANE.vcf"),
          emit: vcf

    script:
    def base    = vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def out_vcf = "${base}_MANE.vcf"
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
        "${case_id}_snpeff_input.vcf" > "${out_vcf}"
    """
}


process CLINVAR {
    tag "${patient_id}:${case_id}"
    container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'
    publishDir { "${params.outdir}/${patient_id}/vcfs/clinvar" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(mane_vcf), path(clinvar_vcf)

    output:
    tuple val(patient_id), val(case_id), val(mode),
          path("${mane_vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_${clinvar_vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}.vcf"),
          emit: vcf

    script:
    def base      = mane_vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def clin_stem = clinvar_vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def out_vcf   = "${base}_${clin_stem}.vcf"
    """
    set -euo pipefail

    if [[ "${mane_vcf}" == *.gz ]]; then
        gunzip -c "${mane_vcf}" > "clinvar_input.tmp.vcf"
    else
        cp "${mane_vcf}" "clinvar_input.tmp.vcf"
    fi

    if [[ "${clinvar_vcf}" == *.gz ]]; then
        gunzip -c "${clinvar_vcf}" > "clinvar_db.tmp.vcf"
    else
        cp "${clinvar_vcf}" "clinvar_db.tmp.vcf"
    fi
    SnpSift annotate -sorted "clinvar_db.tmp.vcf" "clinvar_input.tmp.vcf" > "${out_vcf}"

    if ! grep -q '^#CHROM' "${out_vcf}"; then
        echo "[CLINVAR] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "${out_vcf}" >&2
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
          path("${vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_VEP.vcf"),
          emit: vcf
    path("${vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_VEP_stats.txt"),
          emit: stats

    script:
    def base      = vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def out_vcf   = "${base}_VEP.vcf"
    def out_stats = "${base}_VEP_stats.txt"
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
    tuple val(patient_id), val(case_id), val(mode),
          path("${vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_gnomad.vcf"),
          emit: vcf

    script:
    def base    = vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def out_vcf = "${base}_gnomad.vcf"
    """
    set -euo pipefail

    # SnpSift annotate doesn't auto-decompress a .gz database - passing it
    # directly makes it read the raw bgzip bytes as VCF text (visible as
    # "Impropper VCF entry" garbage containing the BGZF block's own "BC"
    # header marker). Same fix CLINVAR already needed for its own DB.
    gunzip -c "${gnomad_vcf}" > "gnomad_db.tmp.vcf"

    # gnomAD is ~14.6GB uncompressed. -sorted does NOT avoid building a full
    # in-memory index - confirmed empirically (DbVcfSorted.open() still calls
    # VcfIndex.index()/loadIntervals() internally) - kept anyway since both
    # inputs are genuinely sorted and it's still the right mode; the actual
    # fix for OOM is heap size. Reproduced a real OutOfMemoryError at the
    # container's baked-in default -Xmx1g (same class of bug as DBNSFP).
    # See nextflow.config withName:'GNOMAD' for why 36g/48G is a reasoned
    # estimate, not an empirically-confirmed minimum.
    export _JAVA_OPTIONS="-Xmx36g"

    SnpSift annotate -name ${params.gnomad_info_prefix} "gnomad_db.tmp.vcf" "${vcf}" > "${out_vcf}"
    rm -f "gnomad_db.tmp.vcf"

    if ! grep -q '^#CHROM' "${out_vcf}"; then
        echo "[GNOMAD] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "${out_vcf}" >&2
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
    tuple val(patient_id), val(case_id), val(mode),
          path("${vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_dbNSFP.vcf"),
          emit: vcf

    script:
    def base    = vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def out_vcf = "${base}_dbNSFP.vcf"
    """
    set -euo pipefail

    # dbNSFP is ~50GB - SnpSift's container-baked default heap (-Xmx1g) OOMs
    # partway through without surfacing an error, silently truncating the
    # variant set. 12g heap (needs >=16G task memory - see nextflow.config
    # withName:'DBNSFP') is what actually got through the full file previously.
    export _JAVA_OPTIONS="-Xmx12g"

    SnpSift dbnsfp -db "${dbnsfp_db}" -f ${params.dbnsfp_fields} \\
        "${vcf}" > "${out_vcf}.new"

    # SnpSift dbnsfp writes its whole header block a second time at the very
    # end of the output (trailing, spurious duplicate #CHROM line, no data
    # after it) - strip it, keeping the real header once plus every
    # non-header data line in original order.
    awk '
      seen==0 { print; if (/^#CHROM/) seen=1; next }
      !/^#/ { print }
    ' "${out_vcf}.new" > "${out_vcf}"
    rm -f "${out_vcf}.new"

    IN_COUNT=\$(grep -vc "^#" "${vcf}")
    OUT_COUNT=\$(grep -vc "^#" "${out_vcf}")
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
    tuple val(patient_id), val(case_id), val(mode),
          path("${vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_spliceAI.vcf"),
          emit: vcf

    script:
    def base    = vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def out_vcf = "${base}_spliceAI.vcf"
    """
    set -euo pipefail

    bgzip -f -c "${vcf}" > "${case_id}_input.vcf.gz"
    tabix -f -p vcf "${case_id}_input.vcf.gz"

    bcftools annotate -a "${spliceai_vcf}" -c INFO/SpliceAI \\
        "${case_id}_input.vcf.gz" \\
        -o "${out_vcf}"

    if ! grep -q '^#CHROM' "${out_vcf}"; then
        echo "[SPLICEAI] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "${out_vcf}" >&2
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
    tuple val(patient_id), val(case_id), val(mode),
          path("${vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')}_oncorecur.vcf"),
          emit: vcf

    script:
    def base    = vcf.name.replaceAll(/\.vcf\.gz$|\.vcf$/, '')
    def out_vcf = "${base}_oncorecur.vcf"
    """
    set -euo pipefail

    /apps/anaconda/2024.02-1/bin/python3 "${projectDir}/bin/annotate_actionable_somatic_germline.py" \\
        --vcf "${vcf}" \\
        --oncokb-gene-list "${oncokb_gene_list}" \\
        --cosmic-cmc "${cosmic_cmc}" \\
        --hotspots-v1 "${hotspots_v1}" \\
        --hotspots-v2 "${hotspots_v2}" \\
        --hotspots-v3 "${hotspots_v3}" \\
        --eviq-germline "${eviq_germline}" \\
        --out "${out_vcf}"

    if ! grep -q '^#CHROM' "${out_vcf}"; then
        echo "[ONCOKB_COSMIC_HOTSPOTS_EVIQ] ERROR: annotated VCF missing #CHROM header" >&2
        head -50 "${out_vcf}" >&2
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
