# Offline Resource Preparation

ASPIS should not contact biological databases or protein-domain services during
a production Snakemake run. Prepare open-license resources once, keep the frozen
files beside the project inputs, and point the project config to the prepared
TSVs.

This gives two useful properties:

- reports are reproducible because database versions are fixed;
- the resource provenance table records what was used;
- missing resources are reported as `not_configured`, not as silent biological
  zeros.

ASPIS defaults to open-license resources. Do not use restricted, license-gated,
academic-only, non-commercial-only, or non-redistributable resources in the
standard validation bundle. KEGG, MSigDB, SignalP, TMHMM/DeepTMHMM, and similar
resources can be useful in some labs, but they are not part of the open ASPIS
default because their terms are not equivalent to open data or open-source
software licenses. If a local project uses them anyway, it must be an explicit
project decision outside the default open-resource workflow.

The machine-readable source policy is
`config/aspis_open_resource_sources.example.yaml`. It documents the current
recommended open RNA-seq sources, the reviewed user-provided categories, and the
manual-only restricted resources. Keep this file and
`tests/validate_open_resource_policy.py` updated whenever the default resource
policy changes.

## Recommended Layout

Use one project-controlled resource directory on the machine that will run the
analysis:

```text
/path/to/aspis_resources/
  source/
    go-basic.obo
    goa_human.gaf.gz
    Ensembl2Reactome_All_Levels.txt
    wikipathways_human.gmt
    project_open_mirna_targets.tsv
  beas/
    feature_sets/
    smallrna_targets/
```

The `source/` files are frozen upstream exports. The `beas/` files are
ASPIS-ready normalized resources. The preparation scripts normalize files that
you already downloaded or exported; they do not silently download biological
databases during Snakemake execution.

## Open RNA-seq ORA/GSEA Feature Sets

Download frozen open source files once. The recommended open source plan is
recorded in `config/aspis_open_resource_sources.example.yaml`. Example source
endpoints:

```bash
mkdir -p /path/to/aspis_resources/source
cd /path/to/aspis_resources/source

curl -L -o go-basic.obo \
  https://current.geneontology.org/ontology/go-basic.obo

curl -L -o goa_human.gaf.gz \
  https://current.geneontology.org/annotations/goa_human.gaf.gz

curl -L -o Ensembl2Reactome_All_Levels.txt \
  https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt
```

Add open GMT files if desired. WikiPathways is the preferred pathway-like
addition because pathway content is open; export or prepare it as GMT and place
it in `source/`, for example:

```text
wikipathways_human.gmt
```

Then prepare ASPIS TSVs. Use the same GTF used for featureCounts/StringTie so
symbols and external IDs are mapped to the count-matrix gene IDs.

```bash
cd /path/to/aspis

python3 workflow/scripts/prepare_feature_set_resources.py \
  --gtf /path/to/reference/Homo_sapiens.GRCh38.112.chr.gtf \
  --outdir /path/to/aspis_resources/beas/feature_sets \
  --resource-version "GO_current_Reactome_current_WikiPathways_manual" \
  --go-gaf /path/to/aspis_resources/source/goa_human.gaf.gz \
  --go-obo /path/to/aspis_resources/source/go-basic.obo \
  --go-id-field symbol \
  --reactome /path/to/aspis_resources/source/Ensembl2Reactome_All_Levels.txt \
  --gmt /path/to/aspis_resources/source/wikipathways_human.gmt \
  --config-fragment /path/to/aspis_resources/beas/feature_sets/aspis_feature_sets.yaml
```

Important outputs:

```text
gene_id_map.tsv
gene_identifier_map.tsv
transcript_to_gene_map.tsv
go_bp.tsv
go_mf.tsv
go_cc.tsv
reactome.tsv
gmt_*.tsv
unmapped_features.tsv
resource_provenance.tsv
resource_summary.tsv
aspis_feature_sets.yaml
```

Review `gene_id_map.tsv`, `gene_identifier_map.tsv`, and
`transcript_to_gene_map.tsv` before trusting a resource bundle. These files show
which GTF gene IDs, gene symbols/synonyms, transcript IDs, and optional external
ID-map aliases were available for resource normalization.

Review `unmapped_features.tsv`. Some unmapped rows are expected when a source
database contains identifiers that are not present in the GTF, but a mostly
unmapped resource usually means the source ID namespace is wrong.

Review `resource_provenance.tsv` and `resource_summary.tsv` before submitting a
real run. If these files are declared under `resources.*.provenance` or
`resources.*.summary`, `validate_project_inputs.py` checks that they exist, are
non-empty, use controlled license status values, and still match the recorded
checksums for local source/output files.

Use the generated YAML fragment to set:

```yaml
resources:
  rnaseq_feature_sets:
    tables: /path/to/aspis_resources/beas/feature_sets/go_bp.tsv,/path/to/...
    provenance: /path/to/aspis_resources/beas/feature_sets/resource_provenance.tsv
    summary: /path/to/aspis_resources/beas/feature_sets/resource_summary.tsv

rnaseq_differential:
  report_feature_set_tables: /path/to/aspis_resources/beas/feature_sets/go_bp.tsv,/path/to/...
```

The same gene-level feature sets can also be used for smallRNA target-gene
enrichment after a miRNA target table has been configured.

## smallRNA Target Resources

smallRNA target enrichment needs a miRNA-to-target table. GO, Reactome, and
other feature-set resources are not target databases; they are used after
targets have been mapped to genes.

Prepare a frozen open-license or project-owned target export. Then normalize it:

```bash
python3 workflow/scripts/prepare_mirna_target_resources.py \
  --gtf /path/to/reference/Homo_sapiens.GRCh38.112.chr.gtf \
  --input /path/to/aspis_resources/source/project_open_mirna_targets.tsv \
  --outdir /path/to/aspis_resources/beas/smallrna_targets \
  --database project_open_targets \
  --evidence-type validated \
  --resource-version "manual_release_label" \
  --config-fragment /path/to/aspis_resources/beas/smallrna_targets/aspis_targets.yaml
```

On G100, the same preparation step is wrapped by a helper that mirrors the
RNA-seq feature-set helper:

```bash
MODE=dry-run bash tests/prepare_g100_smallrna_targets.sh "$ACCOUNT"
MODE=check bash tests/prepare_g100_smallrna_targets.sh "$ACCOUNT"
MODE=run bash tests/prepare_g100_smallrna_targets.sh "$ACCOUNT"
```

ASPIS intentionally does not bundle a universal default miRNA-target database.
miRBase is appropriate for miRNA sequences and names, but it is not a
miRNA-target interaction resource. Many target databases are free web resources
but have academic, non-commercial, citation, registration, or redistribution
terms that need a project-level review before they can be used. For the open
validation path, use a project-owned table or a target export whose license you
have reviewed and can redistribute with the analysis.

If the target export uses Entrez, UniProt, or another external identifier
namespace, provide an
`--id-map-table` with `source_id` and `target_id` columns, where `target_id`
resolves to the GTF gene ID.

The helper expects `project_open_mirna_targets.tsv` under
`/g100_work/$ACCOUNT/aspis_resources/source` by default. Override paths and
metadata with `ASPIS_MIRNA_TARGET_INPUT`, `ASPIS_MIRNA_TARGET_DATABASE`,
`ASPIS_MIRNA_TARGET_EVIDENCE_TYPE`, `ASPIS_MIRNA_TARGET_VERSION`,
`ASPIS_MIRNA_TARGET_LICENSE`, `ASPIS_MIRNA_TARGET_LICENSE_STATUS`,
`ASPIS_MIRNA_TARGET_ID_MAP_TABLES`, and column-specific overrides such as
`ASPIS_MIRNA_TARGET_MIRNA_COLUMN` or `ASPIS_MIRNA_TARGET_TARGET_COLUMN`.

Important outputs:

```text
<database>_targets.tsv
<database>_target_feature_sets.tsv
<database>_unmapped_targets.tsv
<database>_unmapped_mirnas.tsv
<database>_target_provenance.tsv
<database>_target_summary.tsv
aspis_targets.yaml
```

Review `<database>_unmapped_targets.tsv` for target IDs or symbols that could
not be resolved to the GTF gene namespace. Review
`<database>_unmapped_mirnas.tsv` for blank miRNA IDs and rows filtered by the
configured species. The generated `<database>_target_feature_sets.tsv` represents
each miRNA as a feature set of normalized target genes and can be used by the
smallRNA target-feature-set and miRNA-mRNA target-feature-set report layers.

Use the generated YAML fragment to set:

```yaml
resources:
  smallrna_targets:
    target_table: /path/to/aspis_resources/beas/smallrna_targets/project_open_targets_targets.tsv
    target_tables: [/path/to/aspis_resources/beas/smallrna_targets/project_open_targets_targets.tsv]
    target_feature_set_tables: [/path/to/aspis_resources/beas/smallrna_targets/project_open_targets_target_feature_sets.tsv]
    provenance: /path/to/aspis_resources/beas/smallrna_targets/project_open_targets_target_provenance.tsv
    summary: /path/to/aspis_resources/beas/smallrna_targets/project_open_targets_target_summary.tsv
smallrna:
  target_enrichment_mode: table
  target_table: /path/to/aspis_resources/beas/smallrna_targets/project_open_targets_targets.tsv
  target_tables: [/path/to/aspis_resources/beas/smallrna_targets/project_open_targets_targets.tsv]
  target_feature_set_tables: [/path/to/aspis_resources/beas/smallrna_targets/project_open_targets_target_feature_sets.tsv]
mirna_mrna_integration:
  run: true
```

## Isoform-Switch Functional Annotation

Core isoform-switch calling and exon diagrams do not require external domain
databases. More complete consequence annotation has two layers.

First, configure sequence extraction when an appropriate genome object is
available:

```yaml
rnaseq_differential:
  isoform_switch_genome_object: "BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38"
```

Use a genome object whose chromosome naming matches the GTF. If the GTF uses
Ensembl-style chromosomes (`1`, `2`, `X`) and the genome object uses UCSC-style
names (`chr1`, `chr2`, `chrX`), sequence extraction can fail or produce empty
FASTA files. Use `BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38`
only when the annotation uses UCSC-style `chr` sequence names.

Second, provide precomputed annotation tables or command templates. Precomputed
tables are the most reproducible interface:

```yaml
rnaseq_differential:
  isoform_switch_functional_annotation_tables:
    - /path/to/aspis_resources/beas/isoform_annotations/interproscan.tsv
    - /path/to/aspis_resources/beas/isoform_annotations/pfam.domtblout
    - /path/to/aspis_resources/beas/isoform_annotations/cpat.tsv
```

Alternatively, command templates can run site-installed tools from the report
step. Templates may use `{aa_fasta}`, `{nt_fasta}`, `{outdir}`, and
`{tool_name}`. Generated annotation files under `external_annotations/` are
parsed automatically.

```yaml
rnaseq_differential:
  isoform_switch_interproscan_command: >-
    interproscan.sh -i {aa_fasta} -f TSV -o {outdir}/interproscan.tsv
  isoform_switch_pfam_command: >-
    hmmscan --domtblout {outdir}/pfam.domtblout /path/to/Pfam-A.hmm {aa_fasta}
  isoform_switch_coding_potential_command: >-
    cpat.py -x /path/to/cpat_model/hexamer.tsv -d /path/to/cpat_model/model.dat
    -g {nt_fasta} -o {outdir}/cpat
```

InterPro/InterProScan and Pfam/HMMER are appropriate open defaults when
installed from open channels and paired with their open downloadable data. CPAT
is also suitable when its model files are kept under project provenance. Tools
or databases with academic-only, non-commercial, registration-only, or otherwise
restricted terms are not part of the ASPIS open-resource validation path.

## G100 BEAS Application

For the current BEAS validation, prepare resources under a shared path such as:

```text
/g100_work/<ACCOUNT>/aspis_resources/beas/
```

After running the preparation scripts, merge the generated YAML fragments into:

```text
config/aspis_g100_beas_full.yaml
```

For BEAS RNA-seq ORA/GSEA resources, the G100 helper wraps the approved GO and
Reactome preparation command and fails early when frozen source files are
missing:

```bash
MODE=dry-run bash tests/prepare_g100_beas_feature_sets.sh "$ACCOUNT"
MODE=check bash tests/prepare_g100_beas_feature_sets.sh "$ACCOUNT"
MODE=run bash tests/prepare_g100_beas_feature_sets.sh "$ACCOUNT"
```

The helper expects frozen source files under
`/g100_work/$ACCOUNT/aspis_resources/source` by default. Override paths with
`ASPIS_RESOURCE_SOURCE_DIR`, `ASPIS_RESOURCE_ROOT`, `ASPIS_RESOURCE_GTF`,
`ASPIS_GO_GAF`, `ASPIS_GO_OBO`, `ASPIS_REACTOME`, and optional
`ASPIS_OPEN_GMT` when needed. Inspect the generated `unmapped_features.tsv`,
`resource_summary.tsv`, and `aspis_feature_sets.yaml` before copying paths into
the real BEAS config.

The BEAS full config already enables the GRCh38 BSgenome object needed for
isoform-switch NT/AA sequence extraction. The default open-resource decision is
recorded in `config/aspis_open_resource_sources.example.yaml`: prepare GO
GAF/OBO and Reactome for RNA-seq enrichment, add only reviewed open/user-owned
GMTs when needed, and configure miRNA targets only from reviewed open or
project-owned local exports.

For BEAS smallRNA target enrichment resources, put the frozen reviewed target
export under `/g100_work/$ACCOUNT/aspis_resources/source` or point the helper at
the export explicitly:

```bash
export ASPIS_MIRNA_TARGET_INPUT=/g100_work/$ACCOUNT/aspis_resources/source/project_open_mirna_targets.tsv
export ASPIS_MIRNA_TARGET_DATABASE=project_open_targets
export ASPIS_MIRNA_TARGET_EVIDENCE_TYPE=user_provided
export ASPIS_MIRNA_TARGET_VERSION=manual_release_label
MODE=dry-run bash tests/prepare_g100_smallrna_targets.sh "$ACCOUNT"

MODE=check bash tests/prepare_g100_smallrna_targets.sh "$ACCOUNT"
MODE=run bash tests/prepare_g100_smallrna_targets.sh "$ACCOUNT"
```

ORA/GSEA and miRNA target enrichment remain
disabled until the prepared feature-set and target TSVs exist and their paths
are pasted into the config.

Then rerun:

```bash
MODE=dry-run bash tests/run_g100_full_project.sh "$ACCOUNT" config/aspis_g100_beas_full.yaml
MODE=run bash tests/run_g100_full_project.sh "$ACCOUNT" config/aspis_g100_beas_full.yaml
```
