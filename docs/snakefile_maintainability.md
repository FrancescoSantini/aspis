# Snakefile Maintainability Notes

The root `Snakefile` remains the only supported user entry point. At roughly
5,900 lines, it is large enough that future edits should be guided by explicit
module boundaries rather than opportunistic refactors.

## Safe To Do Before Full Real-Data Validation

- Keep this document as the split map while real-data validation is still
  running.
- Add or update comments around major rule groups only when touching those
  groups for another reason.
- Keep smoke and contract tests close to the behavior they protect.
- Avoid moving rules across files until the real G100 run has validated output
  paths and report contracts.

## Candidate Module Boundaries

These boundaries are based on the current rule layout and should be treated as
a future mechanical split plan, not as an instruction to refactor immediately.

| Candidate module | Current line range | Contents |
| --- | ---: | --- |
| `workflow/rules/common.smk` | 1-1989 | config loading, helpers, final targets, run dashboard, environment and execution reports, provenance |
| `workflow/rules/materialization.smk` | 2047-2238 | library materialization, manifest generation, branch readiness, branch design |
| `workflow/rules/common_qc.smk` | 2239-2354 | raw FASTQ inspection, per-file FastQC, branch FastQC/MultiQC |
| `workflow/rules/smallrna.smk` | 2355-3749 | smallRNA reference prep, preprocessing, depletion, alignment, quantification, differential, targets, reports |
| `workflow/rules/rnaseq_preprocess_alignment.smk` | 1995-2046 and 3751-4129 | RNA-seq index building, preprocessing, alignment, alignment QC, strandedness inference |
| `workflow/rules/rnaseq_quantification.smk` | 4138-4499 | RNA-seq quantification planning, featureCounts, StringTie, gffcompare, transcript matrix, sample QC |
| `workflow/rules/rnaseq_differential_reports.smk` | 4513-5530 | gene/transcript DESeq2, isoform switch, DTU planning, enrichment, summaries, reports |
| `workflow/rules/smoke_targets.smk` | 5554-5932 | DESeq2/report smoke helper targets |

## Fast And Easy Next Steps

1. Add `workflow/rules/README.md` with this map when the first real G100 run is
   complete.
2. Move only the smoke helper targets first; they are the least coupled to the
   production branch DAG.
3. Run local smoke dry-runs after every module move.
4. Move common QC next, then materialization, then assay-specific modules.

## Not Fast Or Safe Yet

- Moving the config/helper block before real validation. Many downstream input
  functions and final target calculations depend on it.
- Splitting smallRNA or RNA-seq report rules while report readability is still
  being assessed.
- Renaming outputs or changing manifest contracts as part of a module split.
