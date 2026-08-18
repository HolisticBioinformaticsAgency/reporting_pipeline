process CLINVAR_ANNOTATE {
  
  tag "${subject}:${id}"
  container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'
  
  input:
  // (pub_base, subject, id, mode, mane_vcf, clinvar_vcf)
  tuple val(pub_base), val(subject), val(id), val(mode), path(mane_vcf), path(clinvar_vcf)

  output:
  // emit ONLY the final annotated VCF (the one we want bgzipped)
  tuple val(pub_base), val(subject), val(id), val(mode), path("${id}_*.vcf"), emit: snpeffvcf

  script:
  """
  set -euo pipefail

  # helper names MUST NOT start with \${id}_ or they will match the output glob
  MANE_PLAIN="clinvar_input.tmp.vcf"
  CLINVAR_PLAIN="clinvar_db.tmp.vcf"

  # 1) normalize MANE input to plain VCF
  zcat -f "${mane_vcf}" > "\${MANE_PLAIN}"

  # 2) figure out final name
  clin_stem=\$(basename "${clinvar_vcf}")
  clin_stem="\${clin_stem%.gz}"
  clin_stem="\${clin_stem%.vcf}"

  OUT_VCF="${id}_vardict_MANE_\${clin_stem}.vcf"

  # 3) fast path if tabix exists
  if [[ -s "${clinvar_vcf}.tbi" ]]; then
    echo "[CLINVAR_ANNOTATE] Using tabix mode with ${clinvar_vcf}" >&2
    SnpSift annotate -tabix "${clinvar_vcf}" "\${MANE_PLAIN}" > "\${OUT_VCF}"
  else
    echo "[CLINVAR_ANNOTATE] No .tbi for ${clinvar_vcf} — using sorted fallback" >&2
    gunzip -c "${clinvar_vcf}" > "\${CLINVAR_PLAIN}"
    SnpSift annotate -sorted "\${CLINVAR_PLAIN}" "\${MANE_PLAIN}" > "\${OUT_VCF}"
  fi

  # 4) sanity
  if ! grep -q '^#CHROM' "\${OUT_VCF}"; then
    echo "[CLINVAR_ANNOTATE] ERROR: annotated VCF missing #CHROM" >&2
    head -50 "\${OUT_VCF}" >&2
    exit 1
  fi
  """
}
