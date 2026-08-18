// -------------------- Includes --------------------
include { VCF_ANNOTATE }    from './modules/vcf_annotate.nf'
include { VCF_PRIORITISE }  from './modules/vcf_prioritise.nf'
include { PYTMB }           from './modules/pytmb.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

// -------------------- Workflow --------------------
workflow {

    def excluded = ['multiqc', 'reference', 'pipeline_info'] as Set

    // ── Collect preprocessing outputs per patient ──────────────────────────
    patient_data_ch = channel
        .fromList( file(params.preprocessing_outdir).listFiles()?.toList() ?: [] )
        .filter  { it.isDirectory() && !excluded.contains(it.name) }
        .map { patient_dir ->
            def patient_id = patient_dir.name

            def norm_vcf = files("${patient_dir}/vcfs/norm/*_vardict_norm.vcf.gz")[0]
            if (!norm_vcf) error "No norm VCF found for patient: ${patient_id}"

            // case_id from VCF filename; mode from T_N pattern; tumor = first sample
            def case_id     = norm_vcf.name.replace('_vardict_norm.vcf.gz', '')

            def cns_file             = file("${patient_dir}/cnvkit/${case_id}_call.cns")
            def cnr_genemetrics_file = file("${patient_dir}/cnvkit/${case_id}_cnr_genemetrics.txt")
            def msi_file             = files("${patient_dir}/msisensor/*.msi")[0]

            if (!cns_file.exists())             error "No call CNS file found for patient: ${patient_id}"
            if (!cnr_genemetrics_file.exists()) error "No CNR genemetrics file found for patient: ${patient_id}"
            if (!msi_file)                      error "No MSI file found for patient: ${patient_id}"
            def mode        = (case_id ==~ /.*T_.*N/) ? 'paired' : 'single'
            def tumor_label = case_id.split('_')[0]

            tuple(patient_id, case_id, mode, tumor_label, norm_vcf, cns_file, cnr_genemetrics_file, msi_file)
        }

    // ── Split: VCF pipeline channel vs metadata ────────────────────────────
    ch_vcf_in = patient_data_ch.map { pid, cid, mode, tlabel, vcf, cns, cnr_gmet, msi ->
        tuple(pid, cid, mode, vcf)
    }

    ch_meta = patient_data_ch.map { pid, cid, mode, tlabel, vcf, cns, cnr_gmet, msi ->
        tuple(pid, cid, tlabel, cns, cnr_gmet, msi)
    }

    // ── Resource channels ──────────────────────────────────────────────────
    ch_mane_db     = channel.value( file(params.mane_db_dir) )
    ch_clinvar     = channel.value( file(params.clinvar_vcf)  )
    ch_gnomad      = channel.value( file(params.gnomad_vcf) )
    ch_dbnsfp      = channel.value( tuple( file(params.dbnsfp_db),   file("${params.dbnsfp_db}.tbi") ) )
    ch_spliceai    = channel.value( tuple( file(params.spliceai_vcf), file("${params.spliceai_vcf}.tbi") ) )
    ch_oncokb      = channel.value( file(params.oncokb_gene_list) )
    ch_cosmic      = channel.value( file(params.cosmic_cmc) )
    ch_hotspots_v1 = channel.value( file(params.cancerhotspots_v1) )
    ch_hotspots_v2 = channel.value( file(params.cancerhotspots_v2) )
    ch_hotspots_v3 = channel.value( file(params.cancerhotspots_v3) )
    ch_eviq        = channel.value( file(params.eviq_germline) )
    ch_whitelist   = channel.value( file(params.whitelist) )

    // ── Annotation: SnpEff MANE → ClinVar → VEP → gnomAD → dbNSFP → SpliceAI → OncoKB/COSMIC/cancerhotspots/eviQ ──
    VCF_ANNOTATE(
        ch_vcf_in, ch_mane_db, ch_clinvar,
        ch_gnomad, ch_dbnsfp, ch_spliceai,
        ch_oncokb, ch_cosmic, ch_hotspots_v1, ch_hotspots_v2, ch_hotspots_v3, ch_eviq
    )

    // ── Prioritise reportable germline/somatic variants ────────────────────
    VCF_PRIORITISE( VCF_ANNOTATE.out.vcf, ch_whitelist )

    // ── PyTMB (full somatic burden, independent of the curated VCF_PRIORITISE output) ──
    PYTMB(
        VCF_ANNOTATE.out.vcf,
        ch_meta.map { pid, cid, tlabel, cns, cnr_gmet, msi -> tuple(pid, cid, tlabel) },
        channel.value( file(params.pytmb_db_config)  ),
        channel.value( file(params.pytmb_var_config) )
    )

    // ── Per-patient report — reads VCF_PRIORITISE's curated germline/somatic
    //    outputs (whitelist-matched somatic variants included), not the raw
    //    annotated VCF ────────────────────────────────────────────────────
    ch_report_in = VCF_PRIORITISE.out.germline
        .map   { pid, cid, mode, gvcf -> tuple(pid, cid, gvcf) }
        .join( VCF_PRIORITISE.out.somatic.map { pid, cid, mode, svcf -> tuple(pid, cid, svcf) }, by: [0, 1] )
        .join( ch_meta,         by: [0, 1] )
        .join( PYTMB.out.tmb,   by: [0, 1] )
        .map   { pid, cid, gvcf, svcf, tlabel, cns, cnr_gmet, msi, tmb ->
            tuple(pid, cid, gvcf, svcf, cns, cnr_gmet, msi, tmb)
        }

    GENERATE_REPORT(
        ch_report_in,
        channel.value( file(params.hotspot_markers) ),
        channel.value( file(params.chrom_arm_bed)   )
    )
}
