// CGen final reports (2026-10-14, both HTMLs added 2026-10-18): the last step in the
// pipeline. Runs CGen's own cohort_run.py chain (gate + TMB, then CGenReporter
// ingest/annotate/classify/match/report, then CGenGx) and publishes both HTML reports -
// CGenGx's (the canonical, clinician-facing molecular record) and CGenReporter's (the
// full variant-by-variant classification detail; see CGenGx/README.md: "a molecular
// record, not a variant list; the full variant table stays in the CGenReporter report").
// No PDF: neither tool writes one from this cohort_run.py call. No container:
// /fs04/vh83/local_software/cgen_tools/conda_env is a real, already-built conda env on
// the shared filesystem, called by absolute path directly, same recipe already proven by
// hand via sbatch.
process CGEN_REPORT {
    tag "${patient_id}:${case_id}"
    beforeScript 'module load singularity; module load R/4.4.0-mkl; module load samtools/1.19.3'
    publishDir { "${params.outdir}/${patient_id}/report" }, mode: 'copy'

    input:
    // vcf/mode here are taken purely to force ordering onto VCF_ANNOTATE_CGEN's own
    // terminal output (2026-10-15 fix) - CGEN_REPORT reads everything it needs from
    // params.outdir/params.input_dir directly, not from this tuple's files. Was
    // ADD_SCATTER_TO_REPORT.out.report, which only guarantees the OLD highconf branch
    // (VCF_FILTER->...->ADD_SCATTER_TO_REPORT) finished - a separate, independently-
    // scheduled branch from the CGen one (CGEN_GATE->VCF_ANNOTATE_CGEN) this process
    // actually needs published. On the real 41-patient cohort, ADD_SCATTER_TO_REPORT
    // finished before VCF_ANNOTATE_CGEN:SPLICEAI had even run for 2 patients
    // (M24026/M24046) - cohort_run.py's glob for the reportable spliceai VCF found
    // zero files. The M24087 single-sample test had worked only by scheduling luck.
    tuple val(patient_id), val(case_id), val(mode), path(vcf)
    path metadata_tsv
    path run_config_json

    output:
    tuple val(patient_id), val(case_id),
          path("${patient_id}_CGenGx_report.html"),
          path("${patient_id}_CGenReporter_report.html"), emit: html

    script:
    """
    export PATH="/fs04/vh83/local_software/cgen_tools/bcftools_bin:\$PATH"
    export APPTAINER_BINDPATH=/fs04
    export CGEN_REF_MANIFEST=/fs04/vh83/local_software/cgen_tools/CGenReporter/config/references_manifest_m3.yaml
    export ONCOKB_API_TOKEN="${params.oncokb_api_token}"

    # params.input_dir/outdir are resolved relative to the launch dir, not this
    # task's own scratch cwd - file() here resolves them to absolute paths, same
    # as how the rest of main.nf already reads params.input_dir.
    #
    # --out must be absolute, not "." - cohort_run.py stores case.json's
    # inputs.reporter_classified as str(out/case_id/.../classified_*.json); when
    # out is relative that string starts with the case_id again, and cgengx's own
    # (case_dir / that_string).resolve() then doubles the case_id segment
    # (case_dir is already absolute by the time cgengx resolves it) - confirmed
    # 2026-10-14 (FileNotFoundError: .../M24087/M24087/CGenReporter/classified_*.json).
    /fs04/vh83/local_software/cgen_tools/conda_env/bin/python3 -m cgen_reporter.cli cohort run \\
        --metadata "${metadata_tsv}" \\
        --case "${patient_id}" \\
        --preprocessing "${file(params.input_dir)}" \\
        --reporting "${file(params.outdir)}" \\
        --out "\$PWD" \\
        --config "${run_config_json}" \\
        --force

    cp "${patient_id}"/CGenGx/"${patient_id}"_cgengx_report_*.html "${patient_id}_CGenGx_report.html"
    cp "${patient_id}"/CGenReporter/"${patient_id}"_Comprehensive_Variant_Classification_Report_*.html "${patient_id}_CGenReporter_report.html"
    """
}
