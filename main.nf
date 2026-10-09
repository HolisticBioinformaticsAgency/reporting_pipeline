// -------------------- Includes --------------------
include { VCF_FILTER }              from './modules/vcf_filter.nf'
include { VCF_ANNOTATE }            from './modules/vcf_annotate.nf'
include { VCF_ANNOTATE as VCF_ANNOTATE_CGEN } from './modules/vcf_annotate.nf'
include { CGEN_GATE }               from './modules/cgen_gate.nf'
include { VCF_PRIORITISE }          from './modules/vcf_prioritise.nf'
include { WHITELIST }               from './modules/whitelist.nf'
include { PYTMB }                   from './modules/pytmb.nf'
include { GENERATE_REPORT }         from './modules/generate_report.nf'
include { CREATE_CNVKIT_SCATTER }   from './modules/create_cnvkit_scatter.nf'
include { ADD_SCATTER_TO_REPORT }   from './modules/add_scatter_to_report.nf'

// A paired tumour/normal norm VCF directly under `dir` (T_..N pattern only -
// excludes normal-only VCFs like a bare "M24017N_vardict_norm.vcf.gz").
def findPairedNormVcf(dir) {
    def hits = files("${dir}/vcfs/norm/*T_*N_vardict_norm.vcf.gz")
    hits ? hits[0] : null
}

// -------------------- Workflow --------------------
workflow {

    def excluded = ['multiqc', 'reference', 'pipeline_info'] as Set

    // ── Collect preprocessing outputs per patient ──────────────────────────
    // Most patient dirs are one tumour/normal case directly. Two wrinkles:
    //   - multi-tumour patients (e.g. M24012, M24020) have no VCF directly
    //     under the patient dir - each tumour instead sits one level deeper,
    //     in its own "<tumour>T_<normal>N/" case subdirectory (vcfs/cnvkit/
    //     facets_purity/msisensor/polysolver all nested there; bams/fastqc/
    //     hsmetrics stay shared at the patient-dir level) - these become
    //     separate patients (M24012A, M24012B, ...) here.
    //   - normal-only patients (e.g. M24017: only an "N" BAM, a germline-only
    //     norm VCF, no tumour) have no paired T_..N VCF anywhere and are
    //     skipped entirely - the reporting pipeline needs a tumour/normal pair.
    patient_data_ch = channel
        .fromList( file(params.input_dir).listFiles()?.toList() ?: [] )
        .filter  { it.isDirectory() && !excluded.contains(it.name) }
        .flatMap { patient_dir ->
            def direct_vcf = findPairedNormVcf(patient_dir)
            def cases = []
            if (direct_vcf) {
                cases << [patient_dir: patient_dir, case_dir: patient_dir, norm_vcf: direct_vcf]
            } else {
                (patient_dir.listFiles()?.toList() ?: []).each { sub ->
                    if (sub.isDirectory()) {
                        def sub_vcf = findPairedNormVcf(sub)
                        if (sub_vcf) cases << [patient_dir: patient_dir, case_dir: sub, norm_vcf: sub_vcf]
                    }
                }
            }
            if (!cases) {
                log.warn "Skipping ${patient_dir.name}: no paired tumour/normal norm VCF found (normal-only sample?)"
            }
            cases
        }
        .map { c ->
            def patient_dir = c.patient_dir
            def case_dir    = c.case_dir
            def norm_vcf    = c.norm_vcf

            def norm_vcf_tbi = files("${case_dir}/vcfs/norm/${norm_vcf.name}.tbi")[0]
            if (!norm_vcf_tbi) error "No norm VCF index (.tbi) found for: ${norm_vcf.name}"

            // case_id from VCF filename; tumor/normal labels from the T_N pattern;
            // patient_id is the tumour label with its trailing 'T' stripped, so a
            // multi-tumour case dir "M24012AT_M24012N" yields patient_id "M24012A".
            def case_id      = norm_vcf.name.replace('_vardict_norm.vcf.gz', '')
            def tumor_label  = case_id.split('_')[0]
            def normal_label = case_id.split('_')[1]
            def patient_id   = tumor_label.replaceAll(/T$/, '')
            def mode         = 'paired'

            def cnr_file              = files("${case_dir}/cnvkit/${case_id}*.cnr")[0]
            // Prefer the FACETS-purity+ploidy-adjusted call.cns when the preprocessing
            // run produced one (*purity.ploidy.call.cns); fall back to the plain one.
            def cns_purity            = files("${case_dir}/cnvkit/${case_id}*purity.ploidy.call.cns")
            def cns_plain             = files("${case_dir}/cnvkit/${case_id}*.call.cns")
            def cns_file              = cns_purity ? cns_purity[0] : (cns_plain ? cns_plain[0] : null)
            def msi_file              = files("${case_dir}/msisensor/*.msi")[0]

            // Optional: FACETS tumour-purity estimate (one line, a fraction or "NA").
            // Fall back to the assets/NO_PURITY placeholder ("NA") when absent.
            def fp_glob               = files("${case_dir}/facets_purity/*.facets_purity.txt")
            def facets_purity_file    = fp_glob ? fp_glob[0] : file("${projectDir}/assets/NO_PURITY")

            // Optional: FACETS ploidy estimate (one line, a number or "NA") -
            // the CNV section's new baseline "neutral" CN (rounded), replacing
            // the old fixed assumption of 2. Fall back to assets/NO_PLOIDY.
            def fpl_glob              = files("${case_dir}/facets_purity/*.facets_ploidy.txt")
            def facets_ploidy_file    = fpl_glob ? fpl_glob[0] : file("${projectDir}/assets/NO_PLOIDY")

            if (!cnr_file) error "No CNR file found for patient: ${patient_id}"
            if (!cns_file) error "No call CNS file found for patient: ${patient_id}"
            if (!msi_file) error "No MSI file found for patient: ${patient_id}"

            def tumour_hla_file   = files("${case_dir}/polysolver/${tumor_label}_hla.txt")[0]
            def germline_hla_file = files("${case_dir}/polysolver/${normal_label}_hla.txt")[0]
            if (!tumour_hla_file)   error "No tumour HLA file found for patient: ${patient_id}"
            if (!germline_hla_file) error "No germline HLA file found for patient: ${patient_id}"

            // bams/ is always at the patient-dir level (shared across a
            // multi-tumour patient's cases), never inside the case subdir.
            def tumour_bam   = files("${patient_dir}/bams/${tumor_label}_sorted_mkdup.bam")[0]
            def tumour_bai   = files("${patient_dir}/bams/${tumor_label}_sorted_mkdup.bam.bai")[0]
            def germline_bam = files("${patient_dir}/bams/${normal_label}_sorted_mkdup.bam")[0]
            def germline_bai = files("${patient_dir}/bams/${normal_label}_sorted_mkdup.bam.bai")[0]
            if (!tumour_bam || !tumour_bai)     error "No tumour BAM/BAI found for patient: ${patient_id}"
            if (!germline_bam || !germline_bai) error "No germline BAM/BAI found for patient: ${patient_id}"

            tuple(patient_id, case_id, mode, tumor_label, norm_vcf, norm_vcf_tbi, cnr_file, cns_file, msi_file,
                  tumour_hla_file, germline_hla_file, tumour_bam, tumour_bai, germline_bam, germline_bai,
                  facets_purity_file, facets_ploidy_file)
        }

    // ── Split: VCF pipeline channel vs metadata ────────────────────────────
    ch_vcf_in = patient_data_ch.map { pid, cid, mode, tlabel, vcf, vcf_tbi, cnr, cns, msi, thla, ghla, tbam, tbai, gbam, gbai, fpur, fploidy ->
        tuple(pid, cid, mode, vcf)
    }

    ch_meta = patient_data_ch.map { pid, cid, mode, tlabel, vcf, vcf_tbi, cnr, cns, msi, thla, ghla, tbam, tbai, gbam, gbai, fpur, fploidy ->
        tuple(pid, cid, tlabel, cnr, cns, msi, thla, ghla, fpur, fploidy)
    }

    // ── WHITELIST reads tumour/germline REF+ALT depth straight off the BAMs
    //    (bcftools mpileup) - no VCF needed ──────────────────────────────────
    ch_whitelist_in = patient_data_ch.map { pid, cid, mode, tlabel, vcf, vcf_tbi, cnr, cns, msi, thla, ghla, tbam, tbai, gbam, gbai, fpur, fploidy ->
        tuple(pid, cid, tbam, tbai, gbam, gbai)
    }

    // ── norm VCF (+ index), for CREATE_CNVKIT_SCATTER's BAF panel ──────────
    ch_norm_vcf = patient_data_ch.map { pid, cid, mode, tlabel, vcf, vcf_tbi, cnr, cns, msi, thla, ghla, tbam, tbai, gbam, gbai, fpur, fploidy ->
        tuple(pid, cid, vcf, vcf_tbi)
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
    ch_repeat_bed  = channel.value( file(params.repetitive_regions_bed) )
    ch_ref_fasta   = channel.value( tuple( file(params.vep_fasta), file("${params.vep_fasta}.fai") ) )
    ch_markers     = channel.value( file(params.hotspot_markers) )
    ch_patient_meta = channel.value( file(params.patient_metadata) )

    // params.whitelist is optional - stage the real file when set, else the
    // assets/NO_WHITELIST placeholder (VCF_PRIORITISE and WHITELIST both
    // detect that name and skip their whitelist-dependent behaviour).
    ch_whitelist_file = params.whitelist
        ? channel.value( file(params.whitelist) )
        : channel.value( file("${projectDir}/assets/NO_WHITELIST") )

    // params.supplementary_json is optional - a manually-maintained, cohort-
    // wide manifest of tool/database versions (see supplementary.json and
    // generate_report.py's load_supplementary()). Missing/invalid JSON
    // degrades gracefully to "not provided" in the report, so the
    // assets/NO_SUPPLEMENTARY placeholder needs no special detection logic.
    ch_supplementary_json = params.supplementary_json
        ? channel.value( file(params.supplementary_json) )
        : channel.value( file("${projectDir}/assets/NO_SUPPLEMENTARY") )

    // ── High-confidence variant selection, applied before annotation ───────
    VCF_FILTER( ch_vcf_in, ch_repeat_bed )

    // ── Annotation: SnpEff MANE → ClinVar → VEP → gnomAD → dbNSFP → SpliceAI → OncoKB/COSMIC/cancerhotspots/eviQ ──
    VCF_ANNOTATE(
        VCF_FILTER.out.vcf, ch_mane_db, ch_clinvar,
        ch_gnomad, ch_dbnsfp, ch_spliceai,
        ch_oncokb, ch_cosmic, ch_hotspots_v1, ch_hotspots_v2, ch_hotspots_v3, ch_eviq
    )

    // ── CGen branch (2026-10-09): CGen's somatic gate on the full norm VCF, then the same annotation chain on the
    //    gated calls, for CGenReporter and CGenGx. Paired cases only; the highconf path above is unchanged ──────
    CGEN_GATE( ch_vcf_in.filter { pid, cid, mode, raw -> mode == 'paired' } )
    VCF_ANNOTATE_CGEN(
        CGEN_GATE.out.vcf, ch_mane_db, ch_clinvar,
        ch_gnomad, ch_dbnsfp, ch_spliceai,
        ch_oncokb, ch_cosmic, ch_hotspots_v1, ch_hotspots_v2, ch_hotspots_v3, ch_eviq
    )

    // ── Prioritise reportable germline/somatic variants. Somatic list is
    //    marker-gene-restricted and whitelist-rescued (see vcf_prioritise.nf).
    //    The raw (pre-VCF_FILTER) norm VCF rides along too, so a whitelisted
    //    variant dropped by high-confidence filtering can still be rescued
    //    straight from it there (gate A's raw-VCF fallback) ──────────────────
    ch_prioritise_in = VCF_ANNOTATE.out.vcf
        .join( ch_vcf_in.map { pid, cid, mode, raw -> tuple(pid, cid, raw) }, by: [0, 1] )
        .map  { pid, cid, mode, avcf, raw -> tuple(pid, cid, mode, avcf, raw) }

    VCF_PRIORITISE( ch_prioritise_in, ch_whitelist_file, ch_markers )

    // ── WHO CNS5 hotspot REF/ALT read counts, off the BAMs directly
    //    (independent of VCF_PRIORITISE — feeds the report's supplementary
    //    section, not the main SNV findings).
    //
    //    Optional: only runs when params.whitelist is set. When it isn't,
    //    every patient gets a placeholder (assets/NO_WHITELIST) instead, and
    //    the report renders a "no whitelist included" note in place of the
    //    hotspot read-count table. ───────────────────────────────────────────
    if( params.whitelist ) {
        WHITELIST( ch_whitelist_in, ch_whitelist_file, ch_ref_fasta )
        ch_hotspot_counts = WHITELIST.out.counts
    } else {
        def no_whitelist = file("${projectDir}/assets/NO_WHITELIST")
        ch_hotspot_counts = ch_whitelist_in.map { pid, cid, tbam, tbai, gbam, gbai ->
            tuple(pid, cid, no_whitelist)
        }
    }

    // ── PyTMB (full somatic burden, independent of the curated VCF_PRIORITISE output) ──
    PYTMB(
        VCF_ANNOTATE.out.vcf,
        ch_meta.map { pid, cid, tlabel, cnr, cns, msi, thla, ghla, fpur, fploidy -> tuple(pid, cid, tlabel) },
        channel.value( file(params.pytmb_db_config)  ),
        channel.value( file(params.pytmb_var_config) )
    )

    // ── Per-patient report — reads VCF_PRIORITISE's curated germline/somatic
    //    outputs plus WHITELIST's hotspot read-count table (supplementary
    //    section), not the raw annotated VCF ────────────────────────────────
    ch_report_in = VCF_PRIORITISE.out.germline
        .map   { pid, cid, mode, gvcf -> tuple(pid, cid, gvcf) }
        .join( VCF_PRIORITISE.out.somatic.map { pid, cid, mode, svcf -> tuple(pid, cid, svcf) }, by: [0, 1] )
        .join( ch_hotspot_counts,      by: [0, 1] )
        .join( ch_meta,                by: [0, 1] )
        .join( PYTMB.out.tmb,          by: [0, 1] )
        .map   { pid, cid, gvcf, svcf, counts, tlabel, cnr, cns, msi, thla, ghla, fpur, fploidy, tmb ->
            tuple(pid, cid, gvcf, svcf, counts, fpur, fploidy, cns, msi, tmb, thla, ghla)
        }

    GENERATE_REPORT(
        ch_report_in,
        ch_markers,
        channel.value( file(params.chrom_arm_bed)   ),
        ch_patient_meta,
        ch_whitelist_file,
        ch_supplementary_json
    )

    // ── CNVkit scatter plots — one per chromosome with a chr/arm-level call,
    //    one per appendix segment (region ± 5 Mb, labelled with its marker
    //    gene(s)) — driven by GENERATE_REPORT's scatter manifest ───────────
    ch_scatter_in = ch_meta.map { pid, cid, tlabel, cnr, cns, msi, thla, ghla, fpur, fploidy -> tuple(pid, cid, cnr, cns, fpur, fploidy) }
        .join( GENERATE_REPORT.out.manifest, by: [0, 1] )
        .join( ch_norm_vcf,                  by: [0, 1] )
        .map  { pid, cid, cnr, cns, fpur, fploidy, manifest, vcf, vcf_tbi ->
            tuple(pid, cid, cnr, cns, manifest, vcf, vcf_tbi, fpur, fploidy)
        }

    CREATE_CNVKIT_SCATTER( ch_scatter_in, channel.value( file(params.chrom_arm_bed) ) )

    // ── Embed the scatter plots into the draft report's <!--SCATTER:key-->
    //    placeholders to produce the final published report ────────────────
    ch_final_report_in = GENERATE_REPORT.out.report
        .join( GENERATE_REPORT.out.manifest,        by: [0, 1] )
        .join( CREATE_CNVKIT_SCATTER.out.scatter_dir, by: [0, 1] )

    ADD_SCATTER_TO_REPORT( ch_final_report_in )
}


// Run summary for the copy-number record and any report built on these outputs (CGen, 2026-10-09):
// code version, Nextflow version, dates, outcome, command line and every parameter.
workflow.onComplete {
    def rev = workflow.commitId ?: "git -C ${projectDir} rev-parse HEAD".execute().text.trim()
    def dirty = workflow.commitId ? "" : "git -C ${projectDir} status --porcelain".execute().text.trim()
    def info = [revision: rev, uncommitted_changes: dirty ? true : false, nextflow: nextflow.version.toString(),
                started: workflow.start.toString(), completed: workflow.complete.toString(), success: workflow.success,
                command: workflow.commandLine, params: params]
    file("${params.outdir}/pipeline_info").mkdirs()
    file("${params.outdir}/pipeline_info/run_info_${workflow.start.format('yyyyMMdd_HHmmss')}.json").text =
        groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
}
