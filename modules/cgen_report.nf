// CGen final reports (2026-10-14, both HTMLs added 2026-10-18): the last step in the
// pipeline. Runs CGen's own cohort_run.py chain (gate + TMB, then CGenCurator
// ingest/annotate/classify/match/report, then CGenGx) and publishes both HTML reports -
// CGenGx's (the canonical, clinician-facing molecular record) and CGenCurator's (the
// full variant-by-variant classification detail; see CGenGx/README.md: "a molecular
// record, not a variant list; the full variant table stays in the CGenCurator report").
// No PDF: neither tool writes one from this cohort_run.py call. No container:
// params.cgen_python is a real, already-built conda env on the shared filesystem,
// called by absolute path directly, same recipe already proven by hand via sbatch.
//
// 2026-10-20 (CGen): the case's full run now goes to <outdir>/<patient>/cgen_run/ and is
// kept on M3 (cohort collect and the diagnostic-report comparison read it); before, it was
// written to the task's scratch folder and discarded. Tools come from the cgen-tools
// repository (CGenCurator, renamed from CGenReporter 2026-10-17). Layout:
// cgen-tools/CGenCurator/docs/INSTALL_M3.md.
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
          path("${patient_id}_CGenCurator_report.html"),
          path("${patient_id}_cgengx_results.json"), emit: html

    script:
    def tools = params.cgen_tools_dir
    // params.outdir is resolved relative to the launch dir, not this task's own scratch
    // cwd - file() resolves it to an absolute path, as the rest of main.nf already does.
    // --out must be absolute: cohort_run.py stores absolute input paths in case.json
    // (confirmed 2026-10-14, FileNotFoundError: .../M24087/M24087/CGenReporter/classified_*.json).
    def results = file(params.outdir)
    """
    set -euo pipefail
    export PATH="/fs04/vh83/local_software/cgen_tools/bcftools_bin:\$PATH"
    export APPTAINER_BINDPATH=/fs04
    export CGEN_REF_MANIFEST="${tools}/CGenCurator/config/references_manifest_m3.yaml"
    export ONCOKB_API_TOKEN="${params.oncokb_api_token}"

    # stop at once if the install is not where INSTALL_M3.md says, rather than run on another copy
    "${params.cgen_python}" -c "import sys, cgen_curator; f = cgen_curator.__file__; sys.exit(0 if f.startswith('${tools}/CGenCurator/') else 'cgen_curator is installed from ' + f + ', not ${tools}/CGenCurator (see CGenCurator/docs/INSTALL_M3.md)')"
    test -f "${tools}/CGenGx/cgengx/__init__.py" || { echo "CGenGx not found at ${tools}/CGenGx (see CGenCurator/docs/INSTALL_M3.md)" >&2; exit 1; }

    "${params.cgen_python}" -m cgen_curator.cli cohort run \\
        --metadata "${metadata_tsv}" \\
        --case "${patient_id}" \\
        --preprocessing "${file(params.input_dir)}" \\
        --reporting "${results}" \\
        --out "${results}" --case-subfolder cgen_run \\
        --config "${run_config_json}" \\
        --force

    # the kept folder holds every run's dated files, so take the newest of each
    run="${results}/${patient_id}/cgen_run"
    cp "\$(ls -t "\$run"/CGenGx/${patient_id}_cgengx_report_*.html | head -1)" "${patient_id}_CGenGx_report.html"
    cp "\$(ls -t "\$run"/CGenCurator/${patient_id}_Comprehensive_Variant_Classification_Report_*.html | head -1)" "${patient_id}_CGenCurator_report.html"
    cp "\$(ls -t "\$run"/CGenGx/${patient_id}_cgengx_results_*.json | head -1)" "${patient_id}_cgengx_results.json"
    """
}
