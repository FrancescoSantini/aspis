# ASPIS

ASPIS is a Snakemake 9 workflow for reproducible RNA-seq, smallRNA/miRNA-seq,
and integrated multi-assay sequencing projects. It materializes sequencing
libraries into a canonical project manifest, routes libraries through
assay-specific branches, and produces both human-readable reports and
machine-readable result tables.

The workflow is designed for project-scale analyses where the same biological
samples may have RNA-seq and smallRNA-seq evidence. ASPIS keeps those layers
separate while also building matched cross-assay summaries for review.

## Main Capabilities

- Materialize local FASTQ files or public SRA/ENA run accessions into a
  reproducible run namespace.
- Plan assay-aware branches from one intake table, with support for `rnaseq`
  and `smallrna` libraries.
- Run FASTQ inspection, preprocessing, alignment, quantification, QC summaries,
  and design diagnostics.
- Analyze RNA-seq at gene and transcript level with DESeq2.
- Run optional RNA-seq feature-set enrichment when explicit GMT or table
  resources are configured.
- Run isoform-switch reports and differential transcript usage/splicing
  methods, including DRIMSeq, DEXSeq transcript-feature usage, DEXSeq exon-bin
  usage, SUPPA2, and rMATS when the required tools and inputs are available.
- Analyze smallRNA/miRNA differential expression and optional miRNA target
  enrichment.
- Integrate matched RNA-seq and smallRNA-seq contrasts through miRNA-target
  evidence tables.
- Render project dashboards, branch reports, contrast summaries, plot atlases,
  technical PDFs, status manifests, and audit tables.

## Repository Layout

```text
Snakefile                         Workflow rules
config/aspis.yaml                 Base configuration
config/*project.example.yaml      Real-project configuration templates
config/*project.example.tsv       Intake table templates
envs/aspis-snakemake.yaml         Conda/Mamba environment
workflow/scripts/                 Python and R helper scripts
profiles/slurm/                   Snakemake SLURM profile
docs/                             Detailed usage and resource documentation
tests/                            Smoke tests and contract tests
```

## Installation

ASPIS is distributed as a Snakemake workflow. The recommended environment is
defined in `envs/aspis-snakemake.yaml`.

```bash
git clone https://github.com/FrancescoSantini/aspis.git
cd aspis

mamba env create -f envs/aspis-snakemake.yaml
conda activate aspis-smk9

snakemake --version
```

If the environment already exists:

```bash
mamba env update -n aspis-smk9 -f envs/aspis-snakemake.yaml --prune
```

Some optional analyses require additional tools or reference resources. See
`docs/optional_tool_environments.md` for optional method-specific guidance.

## Intake Table

Each row in the intake table represents one sequencing library. The example
templates are:

- `config/intake_rnaseq_project.example.tsv`
- `config/intake_smallrna_project.example.tsv`

Required columns are intentionally minimal:

- `library_id`: stable unique library identifier.
- `input_1`: local FASTQ path or public run accession.

Recommended columns for real projects:

- `biospecimen_id`: biological sample identifier used to match assays.
- `project`: project or biological system name.
- `assay`: `rnaseq` or `smallrna`.
- `input_2`: second FASTQ for paired-end libraries.
- `condition`: biological condition used in contrasts.
- `time_h`, `dose`, `unit`, `replicate`, `batch`: design metadata used for
  stratification, reporting, and model formulas.

Local files may be single-end or paired-end FASTQs. Public run accessions are
materialized with `sra-tools`; layout is inspected during materialization.

## Project Configuration

Start from one of the real-project templates:

```bash
cp config/aspis_rnaseq_project.example.yaml config/my_project.yaml
cp config/intake_rnaseq_project.example.tsv config/my_project_intake.tsv
```

or, for smallRNA-only projects:

```bash
cp config/aspis_smallrna_project.example.yaml config/my_project.yaml
cp config/intake_smallrna_project.example.tsv config/my_project_intake.tsv
```

Then edit `config/my_project.yaml` to set:

- `intake`: path to the project intake table.
- `paths`: a unique run namespace under `work/`, `meta/`, and `results/`.
- `resources`: genome, annotation, aligner indexes, miRNA references,
  feature-set resources, and target resources.
- `design`: condition column, control label, formula, covariates, and blocking
  factors.
- assay sections such as `rnaseq_alignment`, `rnaseq_quantification`,
  `rnaseq_differential`, `rnaseq_dtu`, and `smallrna`.

Use a new `paths.*` namespace for each independent run. This prevents results
from different configurations from being mixed.

## Reference Resources

ASPIS does not bundle biological databases or restricted resources. Configure
them explicitly in the project YAML.

Typical RNA-seq resources:

- genome FASTA.
- gene annotation GTF.
- STAR or HISAT2 index, or enough information to build one.
- optional feature-set GMT or table files for GO, Reactome, or other licensed
  resources.

Typical smallRNA resources:

- mature miRNA FASTA, for example from miRBase.
- contaminant FASTA for depletion.
- residual genome reference, when residual alignment is enabled.
- optional miRNA target tables and target feature-set tables.

Feature-set enrichment and miRNA target analysis run only when explicit
resources are configured. When resources are absent, reports mark the analysis
as `not_configured` instead of silently pretending it ran.

Resource licensing is the responsibility of the user. Open resources such as GO
and Reactome can be prepared for public workflows. Restricted resources such as
commercial databases, KEGG, or MSigDB collections should only be used when the
user has the appropriate rights.

See:

- `docs/resource_preparation.md`
- `docs/feature_set_resources.md`
- `docs/optional_tool_environments.md`

## Running Locally

Always start with a dry run:

```bash
snakemake --configfile config/my_project.yaml --cores 4 --dry-run
```

Run the workflow:

```bash
snakemake --configfile config/my_project.yaml --cores 8 --rerun-incomplete
```

To build a specific final dashboard, target the run index:

```bash
snakemake results/<run_id>/index.html \
  --configfile config/my_project.yaml \
  --cores 8 \
  --rerun-incomplete
```

Replace `<run_id>` with the namespace configured under `paths`.

## Running On SLURM

ASPIS includes a Snakemake SLURM profile in `profiles/slurm`.

```bash
snakemake results/<run_id>/index.html \
  --workflow-profile profiles/slurm \
  --configfile config/my_project.yaml \
  --rerun-incomplete
```

Cluster resource defaults should be adapted to the local site policy. Runtime,
memory, partition, account, and disk requests can be set through the profile,
the project configuration, or command-line `--set-resources` values.

## Outputs

ASPIS separates transient work, metadata, and final results:

```text
work/<run_id>/       Materialized inputs, scratch files, and temporary outputs
meta/<run_id>/       Manifests, analysis plan, environment report, audit data
results/<run_id>/    Reports, tables, plots, PDFs, and final branch outputs
```

Main report entry points:

- `results/<run_id>/index.html`: run-level dashboard.
- `results/<run_id>/projects/<project>/index.html`: integrated project report.
- `results/<run_id>/branches/rnaseq/<project>/report/index.html`: RNA-seq
  branch report.
- `results/<run_id>/branches/smallrna/<project>/report/index.html`:
  smallRNA branch report, when smallRNA is enabled.

The project report is organized around:

- project overview and contrast evidence matrix.
- RNA-seq differential expression.
- RNA-seq feature-set enrichment.
- DTU and splicing methods.
- isoform-switch candidates and DTU support.
- smallRNA differential expression.
- miRNA targets and matched miRNA-mRNA evidence.
- QC, design, provenance, and source artifacts.

HTML reports are intended for navigation and review. TSV manifests and result
tables are the source of truth for downstream analysis. PDFs are rendered as
portable summaries for sharing, not as replacements for the tabular outputs.

## Testing

Run the local smoke tests after editing workflow logic:

```bash
bash tests/run_local_smokes.sh
```

Selected contract tests can also be run directly, for example:

```bash
python3 tests/validate_project_report_contract.py
python3 tests/validate_biological_integration_contract.py
python3 tests/validate_resource_mapping_qa_contract.py
```

The test suite focuses on workflow contracts, report inventory, resource
handling, and representative smoke runs. Large biological validation runs are
expected to be executed by users on their own data and compute environment.

## Documentation

Detailed documentation lives in `docs/`:

- `docs/resource_preparation.md`: preparing open reference resources.
- `docs/feature_set_resources.md`: feature-set table and GMT contracts.
- `docs/rnaseq_real_project.md`: RNA-seq project setup.
- `docs/smallrna_real_project.md`: smallRNA project setup.
- `docs/optional_tool_environments.md`: optional DTU/splicing tools.
- `docs/real_data_readiness.md`: readiness checks before large runs.

Site-specific launcher scripts and examples are kept under `tests/` and
`docs/`. Treat them as examples to adapt, not as universal commands.

## Reproducibility Notes

- Keep raw sequencing data outside the Git repository.
- Use stable `library_id` and `biospecimen_id` values.
- Use one run namespace per project configuration.
- Record exact reference releases and checksums for public analyses.
- Configure licensed resources explicitly and document their provenance.
- Inspect `meta/<run_id>/environment_report.tsv` and
  `logs/execution/*.execution.tsv` when auditing a run.

## Scope

ASPIS automates sequencing analysis and evidence organization. It does not make
automated biological claims. Reports summarize statistical outputs, QC state,
resource availability, and cross-assay links so that domain experts can review
the evidence.

## License And Citation

No license or citation file is currently included in this repository. Add a
repository license and citation metadata before public release if redistribution
or formal citation is required.
