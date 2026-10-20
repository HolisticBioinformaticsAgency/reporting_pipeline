// CGen (2026-10-20): each case classifies against its own copy of the OncoKB cache (<outdir>/<patient>/cgen_run/oncokb_cache/),
// so five cases can run at once without sharing a SQLite file. This step runs once, after every CGEN_REPORT, and folds the
// copies back into the shared cache, row by row (the more recently fetched answer stays).
process CGEN_CACHE_MERGE {
    tag "oncokb cache"

    input:
    val patient_ids

    script:
    def results = file(params.outdir)
    def copies = patient_ids.collect { pid -> "${results}/${pid}/cgen_run/oncokb_cache/oncokb_byvariant.sqlite" }.join(' ')
    """
    set -euo pipefail
    export CGEN_REF_MANIFEST="${params.cgen_tools_dir}/CGenCurator/config/references_manifest_m3.yaml"
    "${params.cgen_python}" -m cgen_curator.cli cache merge ${copies}
    """
}
