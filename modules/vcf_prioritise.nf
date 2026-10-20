// Prioritise reportable variants from the fully-annotated VCF, via bcftools:
//   - germline: STATUS=Germline, ClinVar P/LP with a high-confidence review
//     status, and on the eviQ germline testing gene list
//   - somatic: STATUS!=Germline AND either
//       * (gate A) an exact CHROM/POS/REF/ALT match to a params.whitelist
//         row - unconditional: a whitelisted variant is reported no matter
//         what else is or isn't true about it (this is how e.g. a TERT
//         promoter call, consequence `upstream_gene_variant`, still gets
//         through). Checked first against the annotated VCF; a whitelist
//         row not found there (because VCF_FILTER's earlier high-confidence
//         filtering dropped it before it ever reached annotation) is
//         rescued a second time straight from the RAW VCF, with a synthetic
//         ANN block built from the whitelist's own GENE/TRANSCRIPT/CDNA/
//         PROTEIN_OR_NAME columns standing in for the SnpEff annotation it
//         never received - see the raw-VCF rescue step near the end of the
//         script below, or
//       * (gate B) any of its ANN blocks names a marker gene (col 2 of
//         params.hotspot_markers) OR it's a COSMIC oncogene or TSG
//         (COSMIC_ONC_TSG contains "oncogene" or "TSG") at COSMIC_CGC_TIER=1,
//         AND it has a COSMIC CMC significance tier, a cancerhotspots match
//         (any version), or a putative HIGH/MODERATE-impact ANN consequence
//         (the ANN regex in CSQ_EVIDENCE below).
//     Gate B's second half (CMC/hotspot/impact), plus gate A itself, is
//     applied first via `bcftools view`; gate B's first half (marker gene /
//     COSMIC tier-1 oncogene) is applied last, as a post-filter awk pass that
//     scans every ANN block's Gene_Name (4th subfield) - not just the first -
//     plus the COSMIC_CGC_TIER/COSMIC_ONC_TSG INFO tags. Scanning all blocks,
//     rather than only the first, matters because SnpEff sorts ANN blocks by
//     consequence severity (impact, then effect type, then canonical
//     transcript, then genomic coordinate - never by gene identity): a
//     variant that's merely upstream/intronic of a marker gene (MODIFIER
//     impact) but exonic-and-synonymous in some unrelated overlapping gene
//     (LOW impact) sorts the irrelevant gene first, so checking only ANN
//     block 1 would miss the marker gene entirely (e.g. chr3:36991971 is
//     synonymous in EPM2AIP1 but upstream of MLH1 - MLH1 is block 2, not 1).
//     Note generate_report.py's displayed Gene_Name still comes from block 1
//     only, so the gene that qualified a variant here is not always the one
//     shown for it there. A whitelist row is carried into this same awk pass
//     so it survives even when neither half of gate B's first clause is true.
//     OncoKB-actionable status is no longer part of this filter.
//
//     A marker-gene hit only counts via its own HIGH/MODERATE consequence,
//     and the CMC-tier/cancerhotspots rescue only counts via block 1's own
//     gene - see the in-line comments on gene_is_marker_with_consequence /
//     block1_gene_is_marker below for why (FFAR1/CD22 at chr19:35351756,
//     NRAS/CSDE1 at chr1:114713908 - both in 261008.md).
//
//     Separately, a compound ANN consequence that spans an entire intron -
//     both splice_acceptor_variant AND splice_donor_variant on the same
//     block at once (e.g. TYRO3 c.1107+2_1108del: frameshift_variant&
//     splice_acceptor_variant&splice_donor_variant&splice_region_variant&
//     intron_variant) - no longer qualifies as a HIGH/MODERATE consequence,
//     even though frameshift_variant/splice_acceptor_variant/
//     splice_donor_variant are each individually HIGH. A genuine point
//     mutation at a splice site only ever disrupts one side of its intron
//     (splice_acceptor_variant&intron_variant OR splice_donor_variant&
//     intron_variant alone); seeing both together means a single event
//     destroyed the whole intron, which on WES is almost always a large
//     VarDict soft-clip-realigned deletion (INFO_REALIGN_3_PRIME) in
//     repetitive intronic sequence that exome capture was never designed to
//     confidently call - not a real coding/splice-disrupting mutation. This
//     does NOT affect real canonical splice-site point mutations (e.g. ATM
//     c.73-1G>A, RB1 c.607+1G>A), which are exon-adjacent, reliably
//     WES-captured, and only ever carry one splice-site term. See
//     spans_whole_intron() below and 261008.md.
//
// The whitelist (params.whitelist) is OPTIONAL. When unset the workflow
// stages the assets/NO_WHITELIST placeholder instead; the process detects
// that name and skips the whitelist rescue. Separately, WHITELIST in
// modules/whitelist.nf still consumes the same file for the report's
// supplementary REF/ALT hotspot-coverage table - a different code path.
//
// CLNREVSTAT had a real bcftools gotcha caught empirically before shipping
// this: it's Number=. in the VCF header, so bcftools comma-splits it into an
// array. A pattern containing commas itself (e.g. the literal ClinVar string
// "criteria_provided,_multiple_submitters,_no_conflicts") never matches any
// single split token - always a silent false negative. Fixed by matching on
// the unique post-split token ("no_conflicts") instead of the full
// comma-containing string.
//
// Second bcftools gotcha, relevant to the two regexes below: `~` uses POSIX
// *basic* regular expressions (BRE), not extended. So alternation is `\|`,
// grouping is `\( \)`, "optional" is `\?`, and a bare `|` / `(` / `)` is a
// literal character. `~` on a Number=. tag (CLNSIG, ANN) is true if *any*
// comma-split element matches, and `^`/`$` anchor each element - but `!~` is
// "some element does NOT match", not the negation of `~`, so neither regex is
// ever negated here.
//
// The ANN consequence check (one arm of gate B's second clause, not an
// unconditional requirement - a CMC tier or hotspot match alone is enough
// without it): ANN's 2nd subfield (pipe-delimited: Allele|Annotation|
// Impact|...) is SnpEff's effect, itself a `&`-joined list.
// `^[^|]*|[^|]*\(<terms>\)[^|]*|` pins the match to that 2nd field of any
// transcript annotation. <terms> is SnpEff's putative HIGH+MODERATE impact
// vocabulary (frameshift/stop/start/splice-site/inframe-indel/missense/
// fusion/exon-loss/... ). `exon_loss` on its own covers the
// "3_prime_UTR_truncation & exon_loss" compound; `exon_loss_variant` covers
// the "sequence_feature & exon_loss_variant" and "5_prime_UTR_truncation &
// exon_loss_variant" compounds.
//
// COSMIC_ONC_TSG (used in gate B's first clause) is Number=. like
// CLNSIG/ANN above, so the same BRE gotcha applies:
// `~"\(^\|,\)\(oncogene\|TSG\)\(,\|$\)"` anchors on the comma-split token,
// since the tag can be a compound like "oncogene,TSG" or "oncogene,fusion"
// or "TSG,fusion" - either role alone is enough to qualify.

process PRIORITISE_FILTER {
    tag "${patient_id}:${case_id}"
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'
    publishDir { "${params.outdir}/${patient_id}/vcfs/prioritised" }, mode: 'copy'

    input:
    tuple val(patient_id), val(case_id), val(mode), path(vcf), path(raw_vcf)
    path whitelist_csv
    path marker_csv

    output:
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_prioritised_germline.vcf"), emit: germline
    tuple val(patient_id), val(case_id), val(mode), path("${case_id}_prioritised_somatic.vcf"),  emit: somatic

    script:
    """
    set -euo pipefail

    # ---- marker gene symbols (col 2 of marker_csv) for the final post-filter ----
    awk -F, 'NR>1 && \$2!="" {print \$2}' "${marker_csv}" | sort -u > markers.txt

    # ---- gate A: whitelist exact CHROM/POS/REF/ALT rescue. WL_EXPR (a
    #      bcftools OR of one clause per row) is added to the somatic
    #      selection; wl_keys.tsv carries the same rows into the gate-B
    #      post-filter so a whitelisted variant survives even in a non-marker,
    #      non-CGC-tier-1-oncogene gene. Both empty when no whitelist was
    #      configured (assets/NO_WHITELIST placeholder). ----
    if [ "${whitelist_csv}" = "NO_WHITELIST" ]; then
        : > wl_keys.tsv
        WL_EXPR=""
    else
        awk -F, 'NR>1 && \$1 ~ /^chr/ {print \$1"\\t"\$2"\\t"\$3"\\t"\$4}' "${whitelist_csv}" > wl_keys.tsv
        WL_EXPR=\$(awk -F, 'NR>1 && \$1 ~ /^chr/ {printf "%s(CHROM=\\"%s\\" && POS=%s && REF=\\"%s\\" && ALT=\\"%s\\")", sep, \$1, \$2, \$3, \$4; sep=" || "}' "${whitelist_csv}")
    fi

    GERMLINE_FILTER='
        (INFO/STATUS="Germline") &&
        (INFO/CLNSIG~"^\\(Likely_pathogenic\\|Pathogenic\\(/Likely_pathogenic\\(|other\\)\\?\\)\\?\\)\$") &&
        (INFO/CLNREVSTAT~"no_conflicts" || INFO/CLNREVSTAT~"reviewed_by_expert_panel") &&
        INFO/eviQ_germline="T"
    '

    # Gate B, second clause: a COSMIC CMC significance tier, a cancerhotspots
    # match (any version), or a putative HIGH/MODERATE-impact ANN consequence
    # - any one is sufficient, this is an OR not an AND. Gate B's first clause
    # (marker gene / COSMIC tier-1 oncogene) is applied afterwards, in the awk
    # pass below.
    CSQ_EVIDENCE='INFO/COSMIC_MUTATION_SIGNIFICANCE_TIER!="." || INFO/cancerhotspots_v1="T" || INFO/cancerhotspots_v2="T" || INFO/cancerhotspots_v3="T" || INFO/ANN~"^[^|]*|[^|]*\\(bidirectional_gene_fusion\\|chromosome_number_variation\\|conservative_inframe_deletion\\|conservative_inframe_insertion\\|disruptive_inframe_deletion\\|disruptive_inframe_insertion\\|duplication\\|exon_loss_variant\\|exon_loss\\|feature_ablation\\|frameshift_variant\\|gene_fusion\\|inversion\\|missense_variant\\|protein_protein_contact\\|rare_amino_acid_variant\\|rearranged_at_DNA_level\\|sequence_feature&exon_loss_variant\\|splice_acceptor_variant\\|splice_branch_variant_U12\\|splice_donor_variant\\|start_lost\\|stop_gained\\|stop_lost\\|structural_interaction_variant\\|transcript_ablation\\|3_prime_UTR_truncation&exon_loss\\|5_prime_UTR_truncation&exon_loss_variant\\)[^|]*|"'

    if [ -n "\$WL_EXPR" ]; then
        SOMATIC_FILTER='INFO/STATUS!="Germline" && (( '"\$WL_EXPR"' ) || ( '"\$CSQ_EVIDENCE"' ))'
    else
        SOMATIC_FILTER='INFO/STATUS!="Germline" && ( '"\$CSQ_EVIDENCE"' )'
    fi

    bcftools view -i "\$GERMLINE_FILTER" -o "${case_id}_prioritised_germline.vcf" "${vcf}"

    bcftools view -i "\$SOMATIC_FILTER" "${vcf}" > somatic.prefilt.vcf

    # ---- final somatic list (gate B, first clause): keep a record only if
    #      ANY of its ANN blocks' Gene_Name is a marker gene, OR it is a
    #      COSMIC oncogene or TSG at CGC tier 1, OR it is an exact whitelist
    #      CHROM/POS/REF/ALT match (gate A) ----
    #
    #      LOH exclusion (gate B only - never overrides a whitelist match,
    #      gate A stays unconditional): STrongLOH and LikelyLOH both include
    #      a near-zero-tumour-AF class (StrongLOH only ever shows tAF=0 in
    #      this cohort; LikelyLOH splits cleanly into a <0.05 AND a >0.95
    #      cluster, each matching VCF_FILTERs own thresholds) - that near-0
    #      AF reflects near-total loss of the ALT allele, not a real mutation
    #      to report, so StrongLOH and low-AF LikelyLOH are excluded from
    #      gate B. LikelyLOH with high AF (ALT allele retained) is unaffected.
    # Same HIGH/MODERATE vocabulary as CSQ_EVIDENCE above, reused here to tie
    # a marker-gene match to that SAME gene's own consequence - see the "why
    # per-block, not per-variant" note just below.
    HM_TERMS='(bidirectional_gene_fusion|chromosome_number_variation|conservative_inframe_deletion|conservative_inframe_insertion|disruptive_inframe_deletion|disruptive_inframe_insertion|duplication|exon_loss_variant|exon_loss|feature_ablation|frameshift_variant|gene_fusion|inversion|missense_variant|protein_protein_contact|rare_amino_acid_variant|rearranged_at_DNA_level|sequence_feature&exon_loss_variant|splice_acceptor_variant|splice_branch_variant_U12|splice_donor_variant|start_lost|stop_gained|stop_lost|structural_interaction_variant|transcript_ablation|3_prime_UTR_truncation&exon_loss|5_prime_UTR_truncation&exon_loss_variant)'

    awk -v mf=markers.txt -v wf=wl_keys.tsv -v hm="\$HM_TERMS" '
        function spans_whole_intron(cons) {
            # A single real splice-site point mutation (the overwhelming
            # majority of splice calls on WES - exome capture kits reliably
            # cover the canonical +-1/+-2 bases right at the exon boundary)
            # only ever disrupts ONE side of its intron, so SnpEff annotates
            # it splice_acceptor_variant&intron_variant OR
            # splice_donor_variant&intron_variant - never both at once.
            # Seeing BOTH splice_acceptor_variant and splice_donor_variant on
            # the same ANN block means the call spans the ENTIRE intron (both
            # splice sites destroyed by one event), which for WES is almost
            # always a large VarDict soft-clip-realigned deletion in
            # repetitive intronic sequence - depth there is incidental, not
            # part of the capture design, so breakpoint calls spanning
            # hundreds of bases of deep intron are far less trustworthy than
            # an exon-adjacent point mutation (see 261008.md).
            return (cons ~ /splice_acceptor_variant/) && (cons ~ /splice_donor_variant/)
        }
        BEGIN {
            while ((getline l < mf) > 0) M[l] = 1
            while ((getline l < wf) > 0) { split(l, p, "\\t"); W[p[1] SUBSEP p[2] SUBSEP p[3] SUBSEP p[4]] = 1 }
        }
        /^#/ { print; next }
        {
            # A marker-gene hit only "counts" on its own if that SAME gene is
            # the one with the HIGH/MODERATE consequence - a marker gene that
            # is merely co-located with an unrelated gene carrying the real
            # evidence (e.g. FFAR1 missense annotated right next to CD22,
            # CD22 itself only downstream_gene_variant/MODIFIER) no longer
            # qualifies this way. COSMIC CMC tier / cancerhotspots are keyed
            # by gene+protein-change (COSMIC_Mutation_AA etc.), i.e. tied to
            # whichever gene actually has a protein-coding consequence here -
            # in practice that is always the block-1 gene (SnpEff sorts the
            # protein-altering/most-severe annotation first). So the rescue
            # is restricted to block 1 gene too, otherwise an unrelated
            # marker gene merely co-located nearby (e.g. CSDE1 next to an
            # NRAS hotspot at chr1:114713908) could ride along on evidence
            # that is not actually its own.
            #
            # A compound consequence spanning a whole intron (both splice
            # sites at once - see spans_whole_intron() above) no longer
            # qualifies here even if frameshift_variant/splice_acceptor_
            # variant/splice_donor_variant are individually HIGH, since that
            # pattern only arises from a large, low-confidence intronic
            # deletion call on WES, not a genuine point mutation.
            gene_is_marker_with_consequence = 0
            block1_gene_is_marker = 0
            if (match(\$8, /(^|;)ANN=[^;]*/)) {
                ann = substr(\$8, RSTART, RLENGTH); sub(/^;?ANN=/, "", ann)
                n = split(ann, aa, ",")
                for (i = 1; i <= n; i++) {
                    split(aa[i], f, "|")
                    if (f[4] in M) {
                        if ((f[2] ~ hm) && !spans_whole_intron(f[2])) gene_is_marker_with_consequence = 1
                        if (i == 1) block1_gene_is_marker = 1
                    }
                }
            }
            tier = ""
            if (match(\$8, /(^|;)COSMIC_CGC_TIER=[^;]*/)) {
                t = substr(\$8, RSTART, RLENGTH); sub(/^;?COSMIC_CGC_TIER=/, "", t); tier = t
            }
            onc_tsg = ""
            if (match(\$8, /(^|;)COSMIC_ONC_TSG=[^;]*/)) {
                o = substr(\$8, RSTART, RLENGTH); sub(/^;?COSMIC_ONC_TSG=/, "", o); onc_tsg = o
            }
            is_cgc1_rescue = (tier == "1") && (onc_tsg ~ /(^|,)(oncogene|TSG)(,|\$)/)

            cmc_tier = ""
            if (match(\$8, /(^|;)COSMIC_MUTATION_SIGNIFICANCE_TIER=[^;]*/)) {
                c = substr(\$8, RSTART, RLENGTH); sub(/^;?COSMIC_MUTATION_SIGNIFICANCE_TIER=/, "", c); cmc_tier = c
            }
            has_cmc_or_hotspot = (cmc_tier != "" && cmc_tier != ".") \\
                || (\$8 ~ /(^|;)cancerhotspots_v1=T(;|\$)/) \\
                || (\$8 ~ /(^|;)cancerhotspots_v2=T(;|\$)/) \\
                || (\$8 ~ /(^|;)cancerhotspots_v3=T(;|\$)/)

            marker_qualifies = gene_is_marker_with_consequence || (block1_gene_is_marker && has_cmc_or_hotspot)

            status = ""
            if (match(\$8, /(^|;)STATUS=[^;]*/)) {
                st = substr(\$8, RSTART, RLENGTH); sub(/^;?STATUS=/, "", st); status = st
            }
            taf = -1
            nf9 = split(\$9, fmt, ":"); af_col = 0
            for (fi = 1; fi <= nf9; fi++) if (fmt[fi] == "AF") af_col = fi
            if (af_col > 0) { split(\$10, tum, ":"); taf = tum[af_col] + 0 }
            is_low_loh = (status == "StrongLOH") || (status == "LikelyLOH" && taf < 0.05)

            key = \$1 SUBSEP \$2 SUBSEP \$4 SUBSEP \$5
            is_whitelisted = (key in W)
            if (is_whitelisted || ((marker_qualifies || is_cgc1_rescue) && !is_low_loh)) print
        }
    ' somatic.prefilt.vcf > "${case_id}_prioritised_somatic.vcf"

    # ---- raw-VCF whitelist rescue: a whitelist row whose record never
    #      reached this point (dropped upstream by VCF_FILTER's
    #      high-confidence filtering, so it was never annotated) is pulled
    #      straight from the RAW (pre-filter) VCF instead. Since there's no
    #      SnpEff ANN to read (by definition - this variant was never
    #      annotated), a synthetic ANN block is built from the whitelist's
    #      own curated GENE/TRANSCRIPT/CDNA/PROTEIN_OR_NAME columns in its
    #      place. generate_report.py's _extract_ann() reads this exactly
    #      like a real ANN block (field order: Allele|Annotation|Impact|
    #      Gene_Name|Gene_ID|Feature_Type|Feature_ID|...|HGVS.c|HGVS.p), and
    #      selection_reason_somatic() already labels any whitelist-key match
    #      "Whitelist match" regardless of whether ANN is real or synthetic
    #      - no changes needed in generate_report.py. Skipped when no
    #      whitelist was configured. ----
    if [ "${whitelist_csv}" != "NO_WHITELIST" ] && [ -s wl_keys.tsv ]; then
        bcftools view -H "${raw_vcf}" | awk -F'\\t' -v OFS='\\t' -v wcsv="${whitelist_csv}" -v somvcf="${case_id}_prioritised_somatic.vcf" '
            BEGIN {
                while ((getline line < wcsv) > 0) {
                    n = split(line, f, ",")
                    if (f[1] !~ /^chr/) continue
                    key = f[1] SUBSEP f[2] SUBSEP f[3] SUBSEP f[4]
                    GENE[key] = f[5]; TX[key] = f[7]; CDNA[key] = f[8]; PROT[key] = f[9]
                }
                while ((getline line < somvcf) > 0) {
                    if (line ~ /^#/) continue
                    split(line, c, "\\t")
                    seen[c[1] SUBSEP c[2] SUBSEP c[4] SUBSEP c[5]] = 1
                }
            }
            {
                key = \$1 SUBSEP \$2 SUBSEP \$4 SUBSEP \$5
                if ((key in GENE) && !(key in seen)) {
                    ann = ".|whitelist_hotspot|.|" GENE[key] "|.|.|" TX[key] "|.|.|" CDNA[key] "|" PROT[key]
                    \$8 = \$8 ";ANN=" ann
                    print
                }
            }
        ' >> "${case_id}_prioritised_somatic.vcf"
    fi

    for f in "${case_id}_prioritised_germline.vcf" "${case_id}_prioritised_somatic.vcf"; do
        if ! grep -q '^#CHROM' "\$f"; then
            echo "[PRIORITISE_FILTER] ERROR: \$f missing #CHROM header" >&2
            head -50 "\$f" >&2
            exit 1
        fi
    done
    """
}


workflow VCF_PRIORITISE {
    take:
    ch_vcf        // (patient_id, case_id, mode, annotated_vcf, raw_vcf)
    ch_whitelist  // value: params.whitelist file, or the assets/NO_WHITELIST placeholder
    ch_markers    // value: params.hotspot_markers file

    main:
    PRIORITISE_FILTER( ch_vcf, ch_whitelist, ch_markers )

    emit:
    germline = PRIORITISE_FILTER.out.germline
    somatic  = PRIORITISE_FILTER.out.somatic
}
