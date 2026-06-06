#!/usr/bin/env python3
"""Prepare offline feature-set resources for ASPIS ORA/GSEA reports.

The workflow itself should not query external biological databases at run time.
By default, this helper accepts open-license resources only, such as GO,
Reactome, WikiPathways/GMT, and project-owned custom tables. KEGG and MSigDB
inputs require an explicit opt-in because their redistribution/use terms are
not equivalent to open data licenses.

This helper converts frozen source exports into the ASPIS TSV contract:

    set_id  feature_id  source  collection  resource_version  description

`feature_id` is mapped to the gene identifier namespace used by the same GTF
that drives featureCounts/StringTie. For the current human RNA-seq configs this
is normally an Ensembl gene ID.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


OUTPUT_COLUMNS = [
    "set_id",
    "feature_id",
    "source",
    "collection",
    "resource_version",
    "description",
    "original_feature_id",
    "mapping_status",
]
UNMAPPED_COLUMNS = [
    "source_file",
    "resource_kind",
    "set_id",
    "original_feature_id",
    "mapping_status",
    "description",
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
    "checksum_sha256",
    "prepared_by",
    "prepared_at",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, help="Reference GTF used by the RNA-seq quantification layer")
    parser.add_argument("--outdir", required=True, help="Directory for ASPIS-ready resource outputs")
    parser.add_argument("--resource-version", default="", help="Version label to attach when source-specific versions are absent")
    parser.add_argument("--prepared-by", default=os.environ.get("USER", ""), help="Provenance author label")
    parser.add_argument("--unmapped-action", choices=["drop", "keep"], default="drop")
    parser.add_argument(
        "--id-map-table",
        action="append",
        default=[],
        help=(
            "Optional TSV/CSV mapping table with source_id and target_id columns. "
            "Use this for external identifiers that need to map to GTF gene_id."
        ),
    )
    parser.add_argument(
        "--kegg-conv-table",
        action="append",
        default=[],
        help=(
            "Restricted/non-default two-column KEGG conv output. Requires "
            "--allow-restricted-resources."
        ),
    )
    parser.add_argument("--go-gaf", default="", help="GO annotation GAF/GAF.GZ export")
    parser.add_argument("--go-obo", default="", help="Optional go-basic.obo/go.obo for GO term names")
    parser.add_argument(
        "--go-id-field",
        choices=["symbol", "db_object_id", "gene_product_form_id"],
        default="symbol",
        help="Which GAF identifier to map through the GTF gene map",
    )
    parser.add_argument("--reactome", default="", help="Reactome Ensembl2Reactome_All_Levels.txt-style export")
    parser.add_argument("--reactome-species", default="Homo sapiens", help="Species filter for Reactome rows")
    parser.add_argument(
        "--kegg-link-table",
        default="",
        help="Restricted/non-default KEGG REST link/pathway output. Requires --allow-restricted-resources.",
    )
    parser.add_argument(
        "--kegg-name-table",
        default="",
        help="Restricted/non-default KEGG pathway-name table. Requires --allow-restricted-resources.",
    )
    parser.add_argument(
        "--gmt",
        action="append",
        default=[],
        help="Open-license GMT-format gene-set file, e.g. WikiPathways",
    )
    parser.add_argument(
        "--msigdb-gmt",
        action="append",
        default=[],
        help="Restricted/non-default MSigDB GMT file. Requires --allow-restricted-resources.",
    )
    parser.add_argument(
        "--allow-restricted-resources",
        action="store_true",
        help="Permit KEGG/MSigDB inputs after the user has verified license terms.",
    )
    parser.add_argument(
        "--custom-table",
        action="append",
        default=[],
        help="ASPIS-style TSV with set_id and feature_id plus optional source/collection/version/description",
    )
    parser.add_argument("--config-fragment", default="", help="Optional YAML fragment path to write")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace")
    return path.open("r", encoding="utf-8-sig", errors="replace")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_gene_version(value: str) -> str:
    return re.sub(r"\.\d+$", "", value.strip())


def clean_id(value: str) -> str:
    value = value.strip()
    for prefix in ["path:", "hsa:", "ncbi-geneid:", "ensembl:", "uniprot:"]:
        if value.lower().startswith(prefix):
            return value.split(":", 1)[1]
    return value


def safe_stem(path: Path) -> str:
    name = path.name
    for suffix in [".gz", ".txt", ".tsv", ".csv", ".gaf", ".gmt"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "resource"


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


class GeneResolver:
    def __init__(self) -> None:
        self.alias_to_gene_ids: dict[str, set[str]] = defaultdict(set)
        self.gene_rows: dict[str, dict[str, str]] = {}

    def add_gene(self, gene_id: str, aliases: Iterable[str], row: dict[str, str]) -> None:
        gene_id = gene_id.strip()
        if not gene_id:
            return
        self.gene_rows.setdefault(gene_id, row)
        all_aliases = {gene_id, split_gene_version(gene_id), *[alias.strip() for alias in aliases if alias.strip()]}
        for alias in all_aliases:
            cleaned = clean_id(alias)
            if cleaned:
                self.alias_to_gene_ids[cleaned].add(gene_id)
                self.alias_to_gene_ids[split_gene_version(cleaned)].add(gene_id)

    def add_alias(self, source_id: str, target_id: str) -> None:
        source_id = clean_id(source_id)
        target_id = clean_id(target_id)
        resolved, status = self.resolve(target_id)
        if status == "mapped" and resolved:
            self.alias_to_gene_ids[source_id].add(resolved)
            self.alias_to_gene_ids[split_gene_version(source_id)].add(resolved)

    def resolve(self, value: str) -> tuple[str, str]:
        value = clean_id(value)
        if not value:
            return "", "empty"
        candidates = self.alias_to_gene_ids.get(value) or self.alias_to_gene_ids.get(split_gene_version(value)) or set()
        if len(candidates) == 1:
            return next(iter(candidates)), "mapped"
        if len(candidates) > 1:
            return "", "ambiguous"
        return "", "unmapped"


def build_gene_resolver(gtf: Path) -> GeneResolver:
    resolver = GeneResolver()
    with open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            gene_biotype = attrs.get("gene_biotype", attrs.get("gene_type", ""))
            synonyms = re.split(r"[,|]", attrs.get("gene_synonym", ""))
            resolver.add_gene(
                gene_id,
                [gene_name, *synonyms],
                {
                    "gene_id": gene_id,
                    "gene_id_stripped": split_gene_version(gene_id),
                    "gene_name": gene_name,
                    "gene_biotype": gene_biotype,
                    "source_gtf": str(gtf),
                },
            )
    return resolver


def sniff_delimiter(path: Path) -> str:
    with open_text(path) as handle:
        first = handle.readline()
    return "," if first.count(",") > first.count("\t") else "\t"


def add_id_map(path: Path, resolver: GeneResolver) -> int:
    delimiter = sniff_delimiter(path)
    added = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames and {"source_id", "target_id"} <= set(reader.fieldnames):
            for row in reader:
                resolver.add_alias(row.get("source_id", ""), row.get("target_id", ""))
                added += 1
        else:
            handle.seek(0)
            for line in handle:
                parts = line.rstrip("\n").split(delimiter)
                if len(parts) >= 2:
                    resolver.add_alias(parts[0], parts[1])
                    added += 1
    return added


def add_kegg_conv(path: Path, resolver: GeneResolver) -> int:
    added = 0
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            resolver.add_alias(parts[0], parts[1])
            resolver.add_alias(parts[1], parts[0])
            added += 1
    return added


def write_gene_map(path: Path, resolver: GeneResolver) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        columns = ["gene_id", "gene_id_stripped", "gene_name", "gene_biotype", "source_gtf"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene_id in sorted(resolver.gene_rows):
            writer.writerow(resolver.gene_rows[gene_id])


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in OUTPUT_COLUMNS})


def standardize(
    memberships: Iterable[dict[str, str]],
    resolver: GeneResolver,
    unmapped_action: str,
    source_file: Path,
    resource_kind: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    unmapped: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for member in memberships:
        original = member.get("feature_id", "").strip()
        mapped, status = resolver.resolve(original)
        if status != "mapped":
            unmapped.append(
                {
                    "source_file": str(source_file),
                    "resource_kind": resource_kind,
                    "set_id": member.get("set_id", ""),
                    "original_feature_id": original,
                    "mapping_status": status,
                    "description": member.get("description", ""),
                }
            )
            if unmapped_action == "drop":
                continue
            mapped = original
        key = (member.get("set_id", ""), mapped, member.get("source", ""), member.get("collection", ""))
        if key in seen:
            continue
        seen.add(key)
        row = dict(member)
        row["feature_id"] = mapped
        row["original_feature_id"] = original
        row["mapping_status"] = status
        rows.append(row)
    rows.sort(key=lambda row: (row.get("source", ""), row.get("collection", ""), row.get("set_id", ""), row.get("feature_id", "")))
    return rows, unmapped


def parse_obo(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    names: dict[str, str] = {}
    namespaces: dict[str, str] = {}
    current_id = ""
    current_name = ""
    current_namespace = ""
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "[Term]":
                if current_id:
                    names[current_id] = current_name
                    namespaces[current_id] = current_namespace
                current_id = current_name = current_namespace = ""
            elif line.startswith("id: "):
                current_id = line[4:].strip()
            elif line.startswith("name: "):
                current_name = line[6:].strip()
            elif line.startswith("namespace: "):
                current_namespace = line[11:].strip()
    if current_id:
        names[current_id] = current_name
        namespaces[current_id] = current_namespace
    return names, namespaces


def go_collection(aspect: str, namespace: str) -> str:
    aspect_map = {"P": "BP", "F": "MF", "C": "CC"}
    namespace_map = {
        "biological_process": "BP",
        "molecular_function": "MF",
        "cellular_component": "CC",
    }
    return aspect_map.get(aspect, namespace_map.get(namespace, aspect or "GO"))


def go_memberships(path: Path, obo: Path | None, go_id_field: str, version: str) -> list[dict[str, str]]:
    names, namespaces = parse_obo(obo) if obo else ({}, {})
    field_index = {"db_object_id": 1, "symbol": 2, "gene_product_form_id": 16}[go_id_field]
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            qualifier = parts[3]
            if "NOT" in {item.strip() for item in qualifier.split("|")}:
                continue
            go_id = parts[4]
            feature_id = parts[field_index] if len(parts) > field_index else ""
            if go_id == "" or feature_id == "":
                continue
            rows.append(
                {
                    "set_id": go_id,
                    "feature_id": feature_id,
                    "source": "GO",
                    "collection": go_collection(parts[8] if len(parts) > 8 else "", namespaces.get(go_id, "")),
                    "resource_version": version,
                    "description": names.get(go_id, ""),
                }
            )
    return rows


def reactome_memberships(path: Path, species: str, version: str) -> list[dict[str, str]]:
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            row_species = parts[5] if len(parts) > 5 else ""
            if species and row_species and species.lower() not in row_species.lower():
                continue
            rows.append(
                {
                    "set_id": parts[1],
                    "feature_id": parts[0],
                    "source": "Reactome",
                    "collection": "Reactome",
                    "resource_version": version,
                    "description": parts[3],
                }
            )
    return rows


def read_kegg_names(path: Path) -> dict[str, str]:
    names: dict[str, str] = {}
    if not path:
        return names
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                names[clean_id(parts[0])] = parts[1]
    return names


def kegg_memberships(path: Path, names_path: Path | None, version: str) -> list[dict[str, str]]:
    names = read_kegg_names(names_path) if names_path else {}
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene_id = clean_id(parts[0])
            pathway = clean_id(parts[1])
            rows.append(
                {
                    "set_id": pathway,
                    "feature_id": gene_id,
                    "source": "KEGG",
                    "collection": "KEGG",
                    "resource_version": version,
                    "description": names.get(pathway, ""),
                }
            )
    return rows


def gmt_memberships(path: Path, version: str, source: str = "GMT") -> list[dict[str, str]]:
    rows = []
    collection = safe_stem(path)
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            set_id, description, *members = parts
            for member in members:
                if member:
                    rows.append(
                        {
                            "set_id": set_id,
                            "feature_id": member,
                            "source": source,
                            "collection": collection,
                            "resource_version": version,
                            "description": description,
                        }
                    )
    return rows


def custom_memberships(path: Path, version: str) -> list[dict[str, str]]:
    delimiter = sniff_delimiter(path)
    rows = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            return rows
        required = {"set_id", "feature_id"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Custom table {path} is missing columns: {sorted(missing)}")
        for raw in reader:
            rows.append(
                {
                    "set_id": (raw.get("set_id") or "").strip(),
                    "feature_id": (raw.get("feature_id") or "").strip(),
                    "source": (raw.get("source") or safe_stem(path)).strip(),
                    "collection": (raw.get("collection") or "custom").strip(),
                    "resource_version": (
                        raw.get("resource_version")
                        or raw.get("version")
                        or raw.get("source_version")
                        or raw.get("database_version")
                        or version
                    ).strip(),
                    "description": (raw.get("description") or "").strip(),
                }
            )
    return rows


def provenance_row(
    resource_id: str,
    kind: str,
    path: Path,
    provider: str,
    source: str,
    collection: str,
    version: str,
    prepared_by: str,
    prepared_at: str,
    notes: str = "",
) -> dict[str, str]:
    return {
        "resource_id": resource_id,
        "resource_kind": kind,
        "path": str(path),
        "provider": provider,
        "source": source,
        "collection": collection,
        "release": version,
        "resource_version": version,
        "url": "",
        "checksum_sha256": sha256(path) if path.exists() and path.is_file() else "",
        "prepared_by": prepared_by,
        "prepared_at": prepared_at,
        "notes": notes,
    }


def write_unmapped(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=UNMAPPED_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in UNMAPPED_COLUMNS})


def write_provenance(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROVENANCE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in PROVENANCE_COLUMNS})


def yaml_list_value(paths: list[Path]) -> str:
    return ",".join(str(path) for path in paths)


def write_config_fragment(path: Path, tables: list[Path], provenance: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table_value = yaml_list_value(tables)
    path.write_text(
        "# Generated by workflow/scripts/prepare_feature_set_resources.py\n"
        "# Paste or merge these keys into a project config after reviewing paths.\n\n"
        "resources:\n"
        "  rnaseq_feature_sets:\n"
        "    gmt: ''\n"
        f"    tables: {table_value!r}\n"
        "rnaseq_differential:\n"
        "  report_feature_sets: ''\n"
        f"  report_feature_set_tables: {table_value!r}\n"
        "smallrna:\n"
        f"  target_feature_set_tables: {table_value!r}\n"
        f"# resource_provenance: {provenance}\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    version = args.resource_version or "unknown"
    prepared_at = datetime.now(timezone.utc).isoformat()
    restricted_inputs = []
    if args.kegg_link_table:
        restricted_inputs.append("--kegg-link-table")
    if args.kegg_name_table:
        restricted_inputs.append("--kegg-name-table")
    if args.kegg_conv_table:
        restricted_inputs.append("--kegg-conv-table")
    if args.msigdb_gmt:
        restricted_inputs.append("--msigdb-gmt")
    if restricted_inputs and not args.allow_restricted_resources:
        raise ValueError(
            "Restricted/non-open resource inputs require --allow-restricted-resources after "
            f"manual license review: {', '.join(restricted_inputs)}"
        )
    resolver = build_gene_resolver(Path(args.gtf))
    for path in args.id_map_table:
        add_id_map(Path(path), resolver)
    for path in args.kegg_conv_table:
        add_kegg_conv(Path(path), resolver)

    gene_map = outdir / "gene_id_map.tsv"
    write_gene_map(gene_map, resolver)
    output_tables: list[Path] = []
    unmapped_rows: list[dict[str, str]] = []
    provenance_rows = [
        provenance_row(
            "gene_id_map",
            "gene_id_map",
            gene_map,
            "ASPIS",
            "GTF",
            "gene_map",
            version,
            args.prepared_by,
            prepared_at,
            f"Prepared from {args.gtf}",
        )
    ]

    def emit(name: str, source_path: Path, kind: str, memberships: list[dict[str, str]], source: str, collection: str) -> None:
        rows, unmapped = standardize(memberships, resolver, args.unmapped_action, source_path, kind)
        output = outdir / f"{name}.tsv"
        write_rows(output, rows)
        output_tables.append(output)
        unmapped_rows.extend(unmapped)
        provenance_rows.append(
            provenance_row(
                name,
                kind,
                output,
                source,
                source,
                collection,
                version,
                args.prepared_by,
                prepared_at,
                f"{len(rows)} mapped memberships; {len(unmapped)} unmapped/ambiguous memberships from {source_path}",
            )
        )

    if args.go_gaf:
        go_path = Path(args.go_gaf)
        go_rows = go_memberships(go_path, Path(args.go_obo) if args.go_obo else None, args.go_id_field, version)
        for collection in ["BP", "MF", "CC"]:
            emit(f"go_{collection.lower()}", go_path, "go_feature_set_table", [row for row in go_rows if row["collection"] == collection], "GO", collection)
    if args.reactome:
        reactome_path = Path(args.reactome)
        emit(
            "reactome",
            reactome_path,
            "reactome_feature_set_table",
            reactome_memberships(reactome_path, args.reactome_species, version),
            "Reactome",
            "Reactome",
        )
    if args.kegg_link_table:
        kegg_path = Path(args.kegg_link_table)
        emit(
            "kegg",
            kegg_path,
            "kegg_feature_set_table",
            kegg_memberships(kegg_path, Path(args.kegg_name_table) if args.kegg_name_table else None, version),
            "KEGG",
            "KEGG",
        )
    for gmt in args.gmt:
        gmt_path = Path(gmt)
        emit(
            f"gmt_{safe_stem(gmt_path)}",
            gmt_path,
            "gmt_feature_set_table",
            gmt_memberships(gmt_path, version),
            "GMT",
            safe_stem(gmt_path),
        )
    for gmt in args.msigdb_gmt:
        gmt_path = Path(gmt)
        emit(f"msigdb_{safe_stem(gmt_path)}", gmt_path, "msigdb_feature_set_table", gmt_memberships(gmt_path, version, "MSigDB"), "MSigDB", safe_stem(gmt_path))
    for table in args.custom_table:
        table_path = Path(table)
        emit(f"custom_{safe_stem(table_path)}", table_path, "custom_feature_set_table", custom_memberships(table_path, version), safe_stem(table_path), "custom")

    unmapped_path = outdir / "unmapped_features.tsv"
    provenance_path = outdir / "resource_provenance.tsv"
    write_unmapped(unmapped_path, unmapped_rows)
    write_provenance(provenance_path, provenance_rows)
    if args.config_fragment:
        write_config_fragment(Path(args.config_fragment), output_tables, provenance_path)

    print(f"Prepared {len(output_tables)} feature-set table(s) in {outdir}")
    print(f"Unmapped/ambiguous memberships: {len(unmapped_rows)}")
    print(f"Provenance: {provenance_path}")
    if args.config_fragment:
        print(f"Config fragment: {args.config_fragment}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
