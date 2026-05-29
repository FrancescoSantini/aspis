# Feature-Set Resource Configuration

ASPIS enrichment is deliberately offline and file-based. GO, Reactome, KEGG,
MSigDB, and project-specific resources should be exported before the workflow is
submitted, committed or archived with the project inputs, and referenced from
the project config.

This avoids hidden network/database state on clusters and makes the enrichment
universe auditable.

## Where Resources Are Used

RNA-seq differential reports use:

```yaml
rnaseq_differential:
  report_feature_sets: resources/feature_sets/msigdb_hallmark.gmt
  report_feature_set_tables: resources/feature_sets/go_bp.tsv,resources/feature_sets/reactome.tsv
```

smallRNA target-gene reports use:

```yaml
smallrna:
  target_feature_sets: resources/feature_sets/msigdb_hallmark.gmt
  target_feature_set_tables: resources/feature_sets/go_bp.tsv,resources/feature_sets/reactome.tsv
```

Use RNA-seq settings for gene/transcript DESeq2 reports. Use smallRNA settings
for target-gene enrichment after miRNA target mapping or miRNA-mRNA integration.
Use `smallrna.mirna_feature_sets` and `smallrna.mirna_feature_set_tables` only
for resources whose members are miRNA identifiers themselves, such as miRNA
families, seed groups, genomic clusters, or curated miRNA classes.

## Accepted GMT Format

GMT is accepted for RNA-seq and smallRNA target-gene resources:

```text
SET_ID<TAB>DESCRIPTION<TAB>FEATURE_ID<TAB>FEATURE_ID...
```

Example:

```text
HALLMARK_INTERFERON_GAMMA_RESPONSE	MSigDB hallmark interferon gamma response	STAT1	IRF1	CXCL10
```

ASPIS records the file stem as `feature_set_source` and uses `gmt` as
`feature_set_collection`.

## Accepted TSV Format

TSV resources are preferred when you want explicit source and collection labels.
Required columns:

```text
set_id	feature_id
```

Optional columns:

```text
source	collection	description	resource_version
```

Recommended convention:

```text
set_id	feature_id	source	collection	resource_version	description
GO:0006955	STAT1	GO	BP	2026-05	immune response
R-HSA-913531	STAT1	Reactome	Reactome	2026-05	interferon signaling
hsa04630	STAT1	KEGG	KEGG	2026-05	JAK-STAT signaling pathway
HALLMARK_INTERFERON_GAMMA_RESPONSE	STAT1	MSigDB	H	2026.1	hallmark interferon gamma response
```

ASPIS records resource versions in enrichment universe and result tables when
TSV resources provide one of these optional columns: `resource_version`,
`version`, `source_version`, `database_version`, `collection_version`, or
`release`. Plain GMT files do not carry structured version metadata, so they
are reported as `unknown` unless converted to TSV or encoded into the file stem.

## ID Matching Rules

Use the same identifier namespace as the DESeq2 result layer whenever possible.

For gene-level RNA-seq:

- `feature_id` should match the gene count/result ID, usually the `Geneid`
  column.
- If counts use Ensembl IDs, use Ensembl IDs in the resource.
- If counts use symbols, use symbols in the resource.

For transcript-level RNA-seq:

- pathway resources are usually gene-level, not transcript-level;
- ASPIS can map transcript result IDs back to parent genes when feature metadata
  contains `transcript_id` and `gene_id`;
- transcript-native feature sets are also allowed when `feature_id` matches
  transcript IDs directly.

For smallRNA target-gene enrichment:

- `feature_id` must match `target_id` after target-table normalization;
- use the same namespace as the miRNA target table, preferably stable gene IDs;
- target-symbol-only resources are fragile unless the target table also uses
  symbols as `target_id`.

For smallRNA miRNA-ID enrichment:

- `feature_id` must match the miRNA feature IDs in the DESeq2 result table,
  usually the `Geneid` values derived from the miRBase SAF/FASTA;
- this is separate from target-gene enrichment and must not use GO/Reactome/KEGG
  gene memberships unless those rows genuinely contain miRNA identifiers;
- appropriate resources include miRNA families, seed families, genomic clusters,
  conserved miRNA groups, or project-curated miRNA classes.

## miRNA Target Table Evidence Labels

Target tables may contain source-specific labels such as `source`,
`source_type`, `database`, and `evidence`. ASPIS preserves those raw fields, but
also emits a controlled `target_evidence_type` column so reports can be filtered
consistently across mixed resources.

Preferred target-table columns:

```text
mirna_id	target_id	target_symbol	database	source	source_type	target_evidence_type	resource_version	evidence
```

Accepted controlled `target_evidence_type` values:

```text
validated
predicted
conserved
user_provided
matched_expressed
inverse_integrated
unspecified
mixed
```

For database target tables, use:

- `validated` for experimentally supported target resources such as
  miRTarBase/TarBase-style exports;
- `predicted` for computational prediction resources such as
  TargetScan/miRanda/miRWalk-style exports;
- `conserved` when the source specifically represents conserved predictions or
  conserved target relationships;
- `user_provided` for project-curated target lists that do not come from a
  named target database.

When `target_evidence_type` is absent, ASPIS infers a controlled label from
`source_type`, `database`, `source`, and `evidence` where possible. Raw
`target_source_type` remains in the outputs, so source-specific detail is not
lost.

Matched RNA-seq integration uses additional controlled labels:

- `matched_expressed` means the database target is present in matched RNA-seq
  DESeq2/count outputs;
- `inverse_integrated` means the matched RNA-seq target also has the expected
  opposite miRNA/mRNA log2 fold-change direction. Anticorrelated subsets are
  reported separately when matched sample-level counts allow correlation
  testing.

Aggregate rows spanning multiple evidence classes use `mixed`.

## Recommended Resource Layout

Keep resources under a project-controlled directory:

```text
resources/
  feature_sets/
    go_bp.tsv
    go_mf.tsv
    go_cc.tsv
    reactome.tsv
    kegg.tsv
    msigdb_hallmark.gmt
    msigdb_c2_reactome.gmt
```

Then reference them explicitly:

```yaml
rnaseq_differential:
  report_feature_set_tables: >-
    resources/feature_sets/go_bp.tsv,
    resources/feature_sets/go_mf.tsv,
    resources/feature_sets/go_cc.tsv,
    resources/feature_sets/reactome.tsv,
    resources/feature_sets/kegg.tsv
  report_feature_sets: resources/feature_sets/msigdb_hallmark.gmt
```

For YAML values, comma-separated paths must not contain spaces unless quoted
carefully. Prefer path names without spaces.

## GO

Use one TSV row per gene-term membership.

Recommended columns:

```text
set_id	feature_id	source	collection	resource_version	description
GO:0006955	ENSG00000115415	GO	BP	2026-05	immune response
GO:0005515	ENSG00000115415	GO	MF	2026-05	protein binding
GO:0005634	ENSG00000115415	GO	CC	2026-05	nucleus
```

Keep BP, MF, and CC either as separate files or as one file with
`collection` set to `BP`, `MF`, or `CC`.

## Reactome

Reactome resources should be gene-to-pathway membership tables:

```text
set_id	feature_id	source	collection	resource_version	description
R-HSA-913531	ENSG00000115415	Reactome	Reactome	2026-05	interferon signaling
```

If you use MSigDB Reactome collections instead, GMT is fine, but record that in
the file name or source label so the report provenance is clear.

## KEGG

KEGG resources should be pathway membership tables:

```text
set_id	feature_id	source	collection	resource_version	description
hsa04630	ENSG00000115415	KEGG	KEGG	2026-05	JAK-STAT signaling pathway
```

Use organism-specific pathway IDs and an identifier namespace that matches your
counts or target table.

## MSigDB

MSigDB GMT files can be used directly:

```yaml
rnaseq_differential:
  report_feature_sets: resources/feature_sets/msigdb_hallmark.gmt,resources/feature_sets/msigdb_c2_cp_reactome.gmt
```

If you convert MSigDB to TSV, use `source` for the database family,
`collection` for the collection/subcollection, and `resource_version` for the
MSigDB release:

```text
set_id	feature_id	source	collection	resource_version	description
HALLMARK_INTERFERON_GAMMA_RESPONSE	STAT1	MSigDB	H	2026.1	hallmark interferon gamma response
REACTOME_INTERFERON_SIGNALING	STAT1	MSigDB	C2:CP:REACTOME	2026.1	interferon signaling
```

## Interpretation

RNA-seq enrichment is interpreted as feature-set enrichment among observed
DESeq2 genes/transcripts.

smallRNA target-gene enrichment is interpreted as potentially regulated target
processes. It is not direct evidence that a pathway is activated or repressed
unless matched RNA-seq integration supports the target expression direction.

smallRNA miRNA-ID feature-set enrichment is interpreted as enrichment among the
observed differentially expressed miRNAs themselves. It is useful for resources
such as seed families or genomic clusters and is deliberately reported
separately from target-gene enrichment.

When matched RNA-seq is available, ASPIS also writes ranked target feature-set
outputs for inverse miRNA-mRNA integrations. These outputs rank target genes by
the matched RNA-seq target DE statistic when available, falling back to signed
`-log10(pvalue/padj)` or target log2 fold change. They are GSEA-style
running-score summaries over matched target genes, not permutation-based fgsea
p-values. ORA-style target feature sets and ranked target feature sets are kept
as separate report assets.

## Minimum Practical Checks

Before a real run:

1. Confirm that `feature_id` overlaps the DESeq2 result IDs or target IDs.
2. Confirm that each source has enough mapped tested features.
3. Inspect `feature_set_universe.tsv` or `target_feature_set_universe.tsv` in
   the report outputs for mapping loss and final universe size.
4. Treat enrichments with tiny universes or high mapping loss as exploratory.
