#!/usr/bin/env python3
"""Prepare offline miRNA target resources for ASPIS smallRNA reports.

The output target table is intentionally simple and auditable:

    mirna_id  target_id  target_symbol  database  source  source_type
    target_evidence_type  resource_version  evidence

`target_id` is mapped through the same GTF gene map used by RNA-seq so matched
miRNA-mRNA integration can compare miRNA targets directly against gene-level
DESeq2 results.

Use open-license or project-owned target resources for the standard ASPIS
validation path. Restricted identifier conversion resources require explicit
opt-in after manual license review.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

from prepare_feature_set_resources import (
    GeneResolver,
    RESOURCE_SUMMARY_COLUMNS,
    add_id_map,
    add_kegg_conv,
    build_gene_resolver,
    provenance_row,
    sha256,
    sniff_delimiter,
)


TARGET_COLUMNS = [
    "mirna_id",
    "target_id",
    "target_symbol",
    "database",
    "source",
    "source_type",
    "target_evidence_type",
    "resource_version",
    "evidence",
]
TARGET_FEATURE_SET_COLUMNS = [
    "set_id",
    "feature_id",
    "source",
    "collection",
    "resource_version",
    "description",
    "mirna_id",
    "database",
    "target_evidence_type",
]
UNMAPPED_COLUMNS = [
    "mirna_id",
    "raw_target_id",
    "raw_target_symbol",
    "mapping_status",
    "evidence",
]
UNMAPPED_MIRNA_COLUMNS = [
    "mirna_id",
    "raw_target_id",
    "raw_target_symbol",
    "mapping_status",
    "reason",
    "species",
    "evidence",
]
PROVENANCE_COLUMNS = [
    "resource_id",
    "resource_kind",
    "path",
    "provider",
    "source",
    "collection",
    "release",
    "resource_version",
    "url",
    "source_path",
    "source_checksum_sha256",
    "checksum_sha256",
    "license",
    "license_status",
    "identifier_namespace",
    "prepared_by",
    "prepared_at",
    "notes",
]


MIRNA_CANDIDATES = [
    "mirna_id",
    "miRNA",
    "miRNA_ID",
    "miRNA ID",
    "miRNA name",
    "mirbase_id",
]
TARGET_ID_CANDIDATES = [
    "target_id",
    "target_gene_id",
    "Target Gene (Entrez ID)",
    "Target Entrez Gene ID",
    "Entrez Gene ID",
    "Gene ID",
    "gene_id",
    "ensembl_gene_id",
    "Ensembl Gene ID",
]
TARGET_SYMBOL_CANDIDATES = [
    "target_symbol",
    "target_gene_symbol",
    "Target Gene",
    "Target Gene Symbol",
    "Gene Symbol",
    "gene_name",
    "symbol",
]
EVIDENCE_CANDIDATES = [
    "evidence",
    "Experiments",
    "Support Type",
    "support_type",
    "method",
    "Reference",
    "References (PMID)",
]
SPECIES_CANDIDATES = [
    "species",
    "Species",
    "Species (Target Gene)",
    "target_species",
]
SUPPORTED_IDENTIFIER_NAMESPACES = {
    "gtf_gene_id",
    "gene_id",
    "ensembl_gene_id",
    "gene_symbol",
    "entrez_gene_id",
    "target_id",
    "toy_gene_id",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, help="Reference GTF used by RNA-seq quantification")
    parser.add_argument("--input", required=True, help="Frozen target database export TSV/CSV")
    parser.add_argument("--outdir", required=True, help="Output resource directory")
    parser.add_argument("--database", required=True, help="Database/resource label for provenance")
    parser.add_argument(
        "--evidence-type",
        required=True,
        choices=["validated", "predicted", "conserved", "user_provided", "unspecified"],
        help="Controlled ASPIS target evidence label",
    )
    parser.add_argument("--resource-version", default="", help="Database/release label")
    parser.add_argument("--prepared-by", default="", help="Provenance author label")
    parser.add_argument("--license", default="user_provided", help="License or usage label recorded in provenance")
    parser.add_argument(
        "--license-status",
        choices=["open", "user_provided", "restricted", "unknown"],
        default="user_provided",
        help="Controlled license status recorded in provenance",
    )
    parser.add_argument(
        "--identifier-namespace",
        default="gtf_gene_id",
        help="Identifier namespace used in normalized target_id values",
    )
    parser.add_argument("--mirna-column", default="", help="Override miRNA column")
    parser.add_argument("--target-column", default="", help="Override target ID column")
    parser.add_argument("--target-symbol-column", default="", help="Override target symbol column")
    parser.add_argument("--evidence-column", default="", help="Override evidence/details column")
    parser.add_argument("--species-column", default="", help="Optional species column")
    parser.add_argument("--species", default="Homo sapiens", help="Species filter when a species column is available")
    parser.add_argument("--id-map-table", action="append", default=[], help="Optional source_id/target_id ID map")
    parser.add_argument(
        "--kegg-conv-table",
        action="append",
        default=[],
        help="Restricted/non-default two-column KEGG conv-style ID map. Requires --allow-restricted-resources.",
    )
    parser.add_argument(
        "--allow-restricted-resources",
        action="store_true",
        help="Permit restricted ID-conversion inputs after the user has verified license terms.",
    )
    parser.add_argument("--unmapped-action", choices=["drop", "keep"], default="drop")
    parser.add_argument("--config-fragment", default="", help="Optional YAML fragment path to write")
    return parser.parse_args()


def choose_column(fieldnames: list[str], explicit: str, candidates: list[str], label: str, required: bool = True) -> str:
    if explicit:
        if explicit not in fieldnames:
            raise ValueError(f"Configured {label} column is absent: {explicit}")
        return explicit
    lowered = {field.lower(): field for field in fieldnames}
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    if required:
        raise ValueError(f"Could not infer {label} column from: {fieldnames}")
    return ""


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    delimiter = sniff_delimiter(path)
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Target table is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def resolve_target(resolver: GeneResolver, raw_id: str, raw_symbol: str) -> tuple[str, str]:
    for value in [raw_id, raw_symbol]:
        mapped, status = resolver.resolve(value)
        if status == "mapped":
            return mapped, "mapped"
        if status == "ambiguous":
            return "", "ambiguous"
    return "", "unmapped"


def identifier_namespace_supported(namespace: str) -> bool:
    return namespace in SUPPORTED_IDENTIFIER_NAMESPACES or namespace.startswith("custom:") or namespace.startswith("project:")


def validate_identifier_namespace(namespace: str) -> None:
    if not namespace:
        raise ValueError("--identifier-namespace must not be empty")
    if not identifier_namespace_supported(namespace):
        raise ValueError(
            f"Unsupported identifier namespace: {namespace!r}. "
            "Use a supported namespace or prefix project-specific namespaces with 'custom:'"
        )


def write_table(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def build_target_feature_sets(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    feature_rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        mirna_id = row.get("mirna_id", "")
        target_id = row.get("target_id", "")
        database = row.get("database", "")
        if not mirna_id or not target_id:
            continue
        key = (mirna_id, target_id, database)
        if key in seen:
            continue
        seen.add(key)
        evidence_type = row.get("target_evidence_type", "")
        feature_rows.append(
            {
                "set_id": mirna_id,
                "feature_id": target_id,
                "source": database,
                "collection": f"{evidence_type}_mirna_targets" if evidence_type else "mirna_targets",
                "resource_version": row.get("resource_version", ""),
                "description": f"Targets of {mirna_id} from {database}".strip(),
                "mirna_id": mirna_id,
                "database": database,
                "target_evidence_type": evidence_type,
            }
        )
    feature_rows.sort(key=lambda row: (row["source"], row["collection"], row["set_id"], row["feature_id"]))
    return feature_rows


def write_config_fragment(
    path: Path, target_table: Path, target_feature_set_table: Path, provenance_table: Path, summary_table: Path
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# Generated by workflow/scripts/prepare_mirna_target_resources.py\n"
        "# Paste or merge these keys into a project config after reviewing paths.\n\n"
        "resources:\n"
        "  smallrna_targets:\n"
        f"    target_table: {str(target_table)!r}\n"
        f"    target_tables: [{str(target_table)!r}]\n"
        f"    target_feature_set_tables: [{str(target_feature_set_table)!r}]\n"
        f"    provenance: {str(provenance_table)!r}\n"
        f"    summary: {str(summary_table)!r}\n"
        "smallrna:\n"
        "  target_enrichment_mode: table\n"
        f"  target_table: {str(target_table)!r}\n"
        f"  target_tables: [{str(target_table)!r}]\n"
        f"  target_feature_set_tables: [{str(target_feature_set_table)!r}]\n"
        "mirna_mrna_integration:\n"
        "  run: true\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    version = args.resource_version or "unknown"
    prepared_at = datetime.now(timezone.utc).isoformat()
    validate_identifier_namespace(args.identifier_namespace)
    if args.kegg_conv_table and not args.allow_restricted_resources:
        raise ValueError(
            "Restricted/non-open resource inputs require --allow-restricted-resources after "
            "manual license review: --kegg-conv-table"
        )

    resolver = build_gene_resolver(Path(args.gtf))
    for path in args.id_map_table:
        add_id_map(Path(path), resolver)
    for path in args.kegg_conv_table:
        add_kegg_conv(Path(path), resolver)

    fieldnames, raw_rows = read_rows(input_path)
    mirna_col = choose_column(fieldnames, args.mirna_column, MIRNA_CANDIDATES, "miRNA")
    target_col = choose_column(fieldnames, args.target_column, TARGET_ID_CANDIDATES, "target ID", required=False)
    symbol_col = choose_column(fieldnames, args.target_symbol_column, TARGET_SYMBOL_CANDIDATES, "target symbol", required=False)
    if not target_col and not symbol_col:
        raise ValueError("Target resource needs either a target ID column or a target symbol column")
    evidence_col = choose_column(fieldnames, args.evidence_column, EVIDENCE_CANDIDATES, "evidence", required=False)
    species_col = choose_column(fieldnames, args.species_column, SPECIES_CANDIDATES, "species", required=False)

    rows: list[dict[str, str]] = []
    unmapped: list[dict[str, str]] = []
    unmapped_mirnas: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for raw in raw_rows:
        raw_target = raw.get(target_col, "") if target_col else ""
        raw_symbol = raw.get(symbol_col, "") if symbol_col else ""
        evidence = raw.get(evidence_col, "") if evidence_col else ""
        species_value = raw.get(species_col, "") if species_col else ""
        mirna = raw.get(mirna_col, "").strip()
        if species_col and args.species and species_value and args.species.lower() not in species_value.lower():
            unmapped_mirnas.append(
                {
                    "mirna_id": mirna,
                    "raw_target_id": raw_target,
                    "raw_target_symbol": raw_symbol,
                    "mapping_status": "filtered_species",
                    "reason": f"species {species_value!r} does not match {args.species!r}",
                    "species": species_value,
                    "evidence": evidence,
                }
            )
            continue
        if not mirna:
            unmapped_mirnas.append(
                {
                    "mirna_id": mirna,
                    "raw_target_id": raw_target,
                    "raw_target_symbol": raw_symbol,
                    "mapping_status": "blank_mirna_id",
                    "reason": "miRNA identifier is empty",
                    "species": species_value,
                    "evidence": evidence,
                }
            )
            continue
        mapped, status = resolve_target(resolver, raw_target, raw_symbol)
        if status != "mapped":
            unmapped.append(
                {
                    "mirna_id": mirna,
                    "raw_target_id": raw_target,
                    "raw_target_symbol": raw_symbol,
                    "mapping_status": status,
                    "evidence": evidence,
                }
            )
            if args.unmapped_action == "drop":
                continue
            mapped = raw_target or raw_symbol
        key = (mirna, mapped, args.database)
        if not mirna or not mapped or key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "mirna_id": mirna,
                "target_id": mapped,
                "target_symbol": raw_symbol,
                "database": args.database,
                "source": args.database,
                "source_type": args.evidence_type,
                "target_evidence_type": args.evidence_type,
                "resource_version": version,
                "evidence": evidence,
            }
        )
    rows.sort(key=lambda row: (row["mirna_id"], row["target_id"], row["database"]))
    target_feature_sets = build_target_feature_sets(rows)
    target_table = outdir / f"{args.database.lower()}_targets.tsv"
    target_feature_set_table = outdir / f"{args.database.lower()}_target_feature_sets.tsv"
    unmapped_table = outdir / f"{args.database.lower()}_unmapped_targets.tsv"
    unmapped_mirna_table = outdir / f"{args.database.lower()}_unmapped_mirnas.tsv"
    provenance_table = outdir / f"{args.database.lower()}_target_provenance.tsv"
    summary_table = outdir / f"{args.database.lower()}_target_summary.tsv"
    write_table(target_table, TARGET_COLUMNS, rows)
    write_table(target_feature_set_table, TARGET_FEATURE_SET_COLUMNS, target_feature_sets)
    write_table(unmapped_table, UNMAPPED_COLUMNS, unmapped)
    write_table(unmapped_mirna_table, UNMAPPED_MIRNA_COLUMNS, unmapped_mirnas)
    provenance = [
        provenance_row(
            f"{args.database.lower()}_targets",
            "smallrna_target_table",
            target_table,
            args.database,
            args.database,
            args.evidence_type,
            version,
            args.prepared_by,
            prepared_at,
            f"{len(rows)} mapped miRNA-target pairs; {len(unmapped)} unmapped/ambiguous rows; source checksum {sha256(input_path)}",
            source_path=input_path,
            license_label=args.license,
            license_status=args.license_status,
            identifier_namespace=args.identifier_namespace,
        ),
        provenance_row(
            f"{args.database.lower()}_target_feature_sets",
            "smallrna_target_feature_set_table",
            target_feature_set_table,
            args.database,
            args.database,
            f"{args.evidence_type}_mirna_targets",
            version,
            args.prepared_by,
            prepared_at,
            f"{len(target_feature_sets)} miRNA target-set memberships from {len(rows)} mapped miRNA-target pairs",
            source_path=input_path,
            license_label=args.license,
            license_status=args.license_status,
            identifier_namespace=args.identifier_namespace,
        ),
    ]
    write_table(provenance_table, PROVENANCE_COLUMNS, provenance)
    summary = [
        {
            "resource_id": f"{args.database.lower()}_targets",
            "resource_kind": "smallrna_target_table",
            "path": str(target_table),
            "source_file": str(input_path),
            "source": args.database,
            "collection": args.evidence_type,
            "resource_version": version,
            "license_status": args.license_status,
            "identifier_namespace": args.identifier_namespace,
            "n_memberships": str(len(rows)),
            "n_sets": str(len({row["mirna_id"] for row in rows})),
            "n_features": str(len({row["target_id"] for row in rows})),
            "n_unmapped_or_ambiguous": str(len(unmapped) + len(unmapped_mirnas)),
            "mapping_status": "empty" if not rows else "ok",
        },
        {
            "resource_id": f"{args.database.lower()}_target_feature_sets",
            "resource_kind": "smallrna_target_feature_set_table",
            "path": str(target_feature_set_table),
            "source_file": str(input_path),
            "source": args.database,
            "collection": f"{args.evidence_type}_mirna_targets",
            "resource_version": version,
            "license_status": args.license_status,
            "identifier_namespace": args.identifier_namespace,
            "n_memberships": str(len(target_feature_sets)),
            "n_sets": str(len({row["set_id"] for row in target_feature_sets})),
            "n_features": str(len({row["feature_id"] for row in target_feature_sets})),
            "n_unmapped_or_ambiguous": str(len(unmapped) + len(unmapped_mirnas)),
            "mapping_status": "empty" if not target_feature_sets else "ok",
        },
    ]
    write_table(summary_table, RESOURCE_SUMMARY_COLUMNS, summary)
    if args.config_fragment:
        write_config_fragment(Path(args.config_fragment), target_table, target_feature_set_table, provenance_table, summary_table)
    print(f"Prepared {len(rows)} miRNA-target pairs: {target_table}")
    print(f"Target feature-set table: {target_feature_set_table}")
    print(f"Unmapped/ambiguous target rows: {len(unmapped)}")
    print(f"Unmapped/filtered miRNA rows: {len(unmapped_mirnas)}")
    print(f"Provenance: {provenance_table}")
    print(f"Resource summary: {summary_table}")
    if args.config_fragment:
        print(f"Config fragment: {args.config_fragment}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
