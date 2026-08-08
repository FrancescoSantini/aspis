# ASPIS configuration reference

ASPIS loads `config/aspis.yaml` as base defaults and extends it with the YAML
passed through `--configfile`. Real projects should therefore copy
`config/templates/aspis_project.template.yaml` into `config/local/` and change
only the values relevant to the project.

## Intake table

Each row is one sequencing library.

| Column | Purpose |
| --- | --- |
| `library_id` | Unique, stable library identifier. |
| `biospecimen_id` | Biological specimen identifier; required to match assays. |
| `project` | Project identifier used to group analyses and reports. |
| `input_1` | First local FASTQ path or public SRA/ENA accession. |
| `input_2` | Second FASTQ path for paired-end libraries; leave empty for single-end. |
| `assay_hint` | `rnaseq` or `smallrna`. |
| `condition` | Group used for contrasts; must include `control_label`. |
| `treatment`, `dose`, `dose_unit`, `time_h` | Study metadata for contrasts and reporting. |
| `replicate` | Biological or technical replicate label. |
| `batch` | Independent technical/batch label, if present. |

The metadata columns can be extended. Columns listed in `contrast_by`,
covariates, or model formulas must be present and correctly populated.

## Run namespace: `paths`

All paths in this block must use one unique namespace. The template uses
`MY_PROJECT` consistently:

- `raw_dir`, `scratch_dir`: materialized inputs and temporary work;
- `metadata_dir`, `manifest`, `analysis_plan`, `environment_report`: audit
  metadata;
- `branch_dir`: assay-specific outputs;
- `run_dashboard`, `run_dashboard_done`: principal final targets;
- `project_report_dir`: integrated project reports.

Never reuse a namespace for a materially different configuration.

## Reference resources: `resources`

Provide the resources appropriate to enabled branches.

| Section | Required when |
| --- | --- |
| `genome.fasta`, `genome.gtf` | RNA-seq alignment/quantification, residual smallRNA genome annotation, or DTU. |
| `genome.star_genome_dir` or `genome.hisat2_index_prefix` | The selected RNA-seq aligner uses a prebuilt index. |
| `mirbase.mature_fasta`, `mirbase.species_prefix` | smallRNA miRNA alignment. |
| `smallrna_contaminants.fasta` | smallRNA contaminant depletion. |
| `residual_genome.*` | `smallrna.residual_run: true`. |
| `rnaseq_feature_sets.*` | RNA-seq feature-set enrichment. |
| `smallrna_targets.*` | miRNA target enrichment or cross-assay integration. |

Record the provider, release, stable local path, and checksum whenever the
resource is public or redistributed. ASPIS does not bundle biological
databases or licence-restricted resources.

## Experimental design: `design`

`condition_col` and `control_label` define the basic differential comparison.
Use `model_formula` for an explicit design such as `~ batch + condition` or a
paired design such as `~ subject + condition`. Only include variables with
independent variation; do not add a batch or positional factor that is fully
confounded with treatment. `covariates`, `batch_factors`, `blocking_factors`,
and `interaction_terms` feed diagnostics and planning metadata.

Assay-level `design_formula` values override `design.model_formula` for the
corresponding differential analysis. Leave them empty to use the default
condition-based model.

## RNA-seq sections

- `rnaseq_alignment`: select `star` or `hisat2`, provide the matching index,
  genome/annotation paths, strandness settings, and compute threads.
- `rnaseq_quantification`: configures featureCounts and StringTie based
  quantification; specify `read_length` and annotation/genome paths.
- `rnaseq_differential`: gene and transcript DESeq2 analysis. `contrast_by`
  defines metadata strata for separate contrasts; keep it empty for a single
  condition comparison.
- `rnaseq_dtu`: optional DTU/splicing methods. Configure only methods whose
  dependencies are installed. Native rMATS requires an explicit positive
  `rmats_read_length` and the correct library type.

## smallRNA section

Set `smallrna.run: true` for a smallRNA branch. The default branch includes
adapter trimming, contaminant depletion, miRNA alignment/quantification,
differential analysis, and reports. Enable `residual_run` only when residual
genome reference/index and annotation paths are ready. Residual results are a
diagnostic classification of miRNA-unassigned reads, not a replacement for
miRNA quantification.

## Matched miRNA-mRNA integration

Enable `mirna_mrna_integration.run` only when both assays use matched
`biospecimen_id` values and valid miRNA target resources are configured. The
integration layer is evidence aggregation; it does not establish causality.

## Compute and provenance

Use `execution` for portable defaults only. Set site-specific account,
partition, and resource requests through the local YAML, profile, or command
line. Keep `provenance.run: true` for normal production analyses.

Before execution, always use a dry run. After execution, review the run
dashboard, environment report, execution log, and report-inventory validation
before distributing outputs.
