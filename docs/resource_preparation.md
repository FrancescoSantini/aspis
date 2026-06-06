# Offline Resource Preparation

ASPIS should not contact GO, Reactome, KEGG, MSigDB, miRNA-target databases, or
protein-domain services during a production Snakemake run. Prepare resources
once, keep the frozen files beside the project inputs, and point the project
config to the prepared TSVs.

This gives two useful properties:

- reports are reproducible because database versions are fixed;
- missing resources are reported as `not_configured`, not as silent biological
  zeros.

## Recommended Layout

Use one project-controlled resource directory on the machine that will run the
analysis:

```text
/path/to/aspis_resources/
  source/
    go-basic.obo
    goa_human.gaf.gz
    Ensembl2Reactome_All_Levels.txt
    kegg_hsa_link_pathway.tsv
    kegg_hsa_pathway_names.tsv
    kegg_hsa_ensembl_conv.tsv
    msigdb_hallmark.gmt
    mirtarbase_hsa.csv
  beas/
    feature_sets/
    smallrna_targets/
```

The `source/` files are frozen upstream exports. The `beas/` files are
ASPIS-ready normalized resources.

## RNA-seq ORA/GSEA Feature Sets

Download frozen source files once. Example source endpoints:

```bash
mkdir -p /path/to/aspis_resources/source
cd /path/to/aspis_resources/source

curl -L -o go-basic.obo \
  https://current.geneontology.org/ontology/go-basic.obo

curl -L -o goa_human.gaf.gz \
  https://current.geneontology.org/annotations/goa_human.gaf.gz

curl -L -o Ensembl2Reactome_All_Levels.txt \
  https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt

curl -L -o kegg_hsa_link_pathway.tsv \
  https://rest.kegg.jp/link/pathway/hsa

curl -L -o kegg_hsa_pathway_names.tsv \
  https://rest.kegg.jp/list/pathway/hsa

curl -L -o kegg_hsa_ensembl_conv.tsv \
  https://rest.kegg.jp/conv/hsa/ensembl
```

MSigDB requires accepting its license terms, so download GMT files manually from
MSigDB and place them in `source/`, for example:

```text
msigdb_hallmark.gmt
msigdb_c2_cp_reactome.gmt
```

Then prepare ASPIS TSVs. Use the same GTF used for featureCounts/StringTie so
symbols and external IDs are mapped to the count-matrix gene IDs.

```bash
cd /path/to/aspis

python3 workflow/scripts/prepare_feature_set_resources.py \
  --gtf /path/to/reference/Homo_sapiens.GRCh38.112.chr.gtf \
  --outdir /path/to/aspis_resources/beas/feature_sets \
  --resource-version "GO_current_Reactome_current_KEGG_2026-06_MSigDB_manual" \
  --go-gaf /path/to/aspis_resources/source/goa_human.gaf.gz \
  --go-obo /path/to/aspis_resources/source/go-basic.obo \
  --go-id-field symbol \
  --reactome /path/to/aspis_resources/source/Ensembl2Reactome_All_Levels.txt \
  --kegg-link-table /path/to/aspis_resources/source/kegg_hsa_link_pathway.tsv \
  --kegg-name-table /path/to/aspis_resources/source/kegg_hsa_pathway_names.tsv \
  --kegg-conv-table /path/to/aspis_resources/source/kegg_hsa_ensembl_conv.tsv \
  --msigdb-gmt /path/to/aspis_resources/source/msigdb_hallmark.gmt \
  --msigdb-gmt /path/to/aspis_resources/source/msigdb_c2_cp_reactome.gmt \
  --config-fragment /path/to/aspis_resources/beas/feature_sets/aspis_feature_sets.yaml
```

Important outputs:

```text
gene_id_map.tsv
go_bp.tsv
go_mf.tsv
go_cc.tsv
reactome.tsv
kegg.tsv
msigdb_*.tsv
unmapped_features.tsv
resource_provenance.tsv
aspis_feature_sets.yaml
```

Review `unmapped_features.tsv`. Some unmapped rows are expected when a source
database contains identifiers that are not present in the GTF, but a mostly
unmapped resource usually means the source ID namespace is wrong.

Use the generated YAML fragment to set:

```yaml
resources:
  rnaseq_feature_sets:
    tables: /path/to/aspis_resources/beas/feature_sets/go_bp.tsv,/path/to/...

rnaseq_differential:
  report_feature_set_tables: /path/to/aspis_resources/beas/feature_sets/go_bp.tsv,/path/to/...
```

The same gene-level feature sets can also be used for smallRNA target-gene
enrichment after a miRNA target table has been configured.

## smallRNA Target Resources

smallRNA target enrichment needs a miRNA-to-target table. GO/Reactome/KEGG are
not target databases; they are used after targets have been mapped to genes.

Prepare a frozen validated or predicted target export, for example miRTarBase,
TarBase, TargetScan, or a project-curated target list. Then normalize it:

```bash
python3 workflow/scripts/prepare_mirna_target_resources.py \
  --gtf /path/to/reference/Homo_sapiens.GRCh38.112.chr.gtf \
  --input /path/to/aspis_resources/source/mirtarbase_hsa.csv \
  --outdir /path/to/aspis_resources/beas/smallrna_targets \
  --database miRTarBase \
  --evidence-type validated \
  --resource-version "manual_release_label" \
  --config-fragment /path/to/aspis_resources/beas/smallrna_targets/aspis_targets.yaml
```

If the target export uses Entrez/KEGG/UniProt identifiers, provide an
`--id-map-table` with `source_id` and `target_id` columns, where `target_id`
resolves to the GTF gene ID. The script also accepts `--kegg-conv-table` for
two-column KEGG-style conversion files.

Use the generated YAML fragment to set:

```yaml
resources:
  smallrna_targets:
    target_table: /path/to/aspis_resources/beas/smallrna_targets/mirtarbase_targets.tsv
smallrna:
  target_enrichment_mode: table
  target_table: /path/to/aspis_resources/beas/smallrna_targets/mirtarbase_targets.tsv
  target_feature_set_tables: /path/to/aspis_resources/beas/feature_sets/go_bp.tsv,/path/to/...
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
  isoform_switch_genome_object: "BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38"
```

Use a genome object whose chromosome naming matches the GTF. If the GTF uses
Ensembl-style chromosomes (`1`, `2`, `X`) and the genome object uses UCSC-style
names (`chr1`, `chr2`, `chrX`), sequence extraction can fail or produce empty
FASTA files.

Second, provide precomputed annotation tables or command templates. Precomputed
tables are the most reproducible interface:

```yaml
rnaseq_differential:
  isoform_switch_functional_annotation_tables:
    - /path/to/aspis_resources/beas/isoform_annotations/interproscan.tsv
    - /path/to/aspis_resources/beas/isoform_annotations/pfam.domtblout
    - /path/to/aspis_resources/beas/isoform_annotations/cpat.tsv
    - /path/to/aspis_resources/beas/isoform_annotations/signalp.tsv
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

InterProScan, Pfam HMM databases, SignalP, DeepTMHMM/TMHMM, DeepLoc2, and
IUPred/NetSurfP are not bundled with ASPIS because their databases, licenses, or
install procedures are site-specific. The stable ASPIS contract is the
normalized annotation table or the command-template output.

## G100 BEAS Application

For the current BEAS validation, prepare resources under a shared path such as:

```text
/g100_work/<ACCOUNT>/aspis_resources/beas/
```

After running the preparation scripts, merge the generated YAML fragments into:

```text
config/aspis_g100_beas_full.yaml
```

Then rerun:

```bash
MODE=dry-run bash tests/run_g100_full_project.sh "$ACCOUNT" config/aspis_g100_beas_full.yaml
MODE=run bash tests/run_g100_full_project.sh "$ACCOUNT" config/aspis_g100_beas_full.yaml
```
