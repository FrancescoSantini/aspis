#!/usr/bin/env python3
"""Merge per-library RNA-seq featureCounts outputs into gene matrices."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from display_labels import gene_display_label


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project", "layout", "bam"}
REQUIRED_PLAN_COLUMNS = {"project", "assay", "status", "gene_counter", "annotation_gtf"}
METADATA_COLUMNS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
GENE_ANNOTATION_COLUMNS = ["gene_name", "gene_display", "gene_biotype"]
ATTR_RE = re.compile(r'([A-Za-z0-9_.-]+) "([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Aligned RNA-seq samples TSV")
    parser.add_argument("--plan", required=True, help="RNA-seq quantification plan TSV")
    parser.add_argument("--counts", required=True, help="Merged gene count matrix TSV")
    parser.add_argument("--metadata", required=True, help="Gene metadata TSV")
    parser.add_argument("--manifest", required=True, help="Combined featureCounts manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("manifests", nargs="+", help="Per-library featureCounts manifests")
    return parser.parse_args()


def read_table(path: Path, required: set[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def read_plan(path: Path) -> dict[str, str]:
    _, rows = read_table(path, REQUIRED_PLAN_COLUMNS)
    if len(rows) != 1:
        raise ValueError(f"Quantification plan must contain exactly one row: {path}")
    row = rows[0]
    if row.get("status") != "ready":
        raise ValueError("Quantification plan is not ready: " + row.get("reason", ""))
    if row.get("gene_counter") != "featurecounts":
        raise ValueError(f"featureCounts requested with gene_counter={row.get('gene_counter')!r}")
    if not Path(row["annotation_gtf"]).exists():
        raise FileNotFoundError(f"annotation_gtf does not exist: {row['annotation_gtf']}")
    return row


def parse_attrs(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(text)}


def read_gene_annotations(gtf_path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}
    biotype_keys = ("gene_biotype", "gene_type", "transcript_biotype", "transcript_type", "biotype")
    with gtf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get("gene_id", "")
            if not gene_id:
                continue
            current = annotations.setdefault(gene_id, {"gene_name": "", "gene_biotype": ""})
            if not current["gene_name"]:
                current["gene_name"] = attrs.get("gene_name", "") or attrs.get("gene", "")
            if not current["gene_biotype"]:
                current["gene_biotype"] = next((attrs[key] for key in biotype_keys if attrs.get(key)), "")
    return annotations


def read_featurecounts_counts(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, int]]:
    metadata: dict[str, dict[str, str]] = {}
    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        data_lines = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(data_lines, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"featureCounts output is empty: {path}")
        sample_column = reader.fieldnames[-1]
        for row in reader:
            gene_id = row["Geneid"]
            metadata[gene_id] = {column: row.get(column, "") for column in METADATA_COLUMNS}
            counts[gene_id] = int(float(row.get(sample_column, "0")))
    return metadata, counts


def ordered_manifest_rows(paths: list[str], expected_ids: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    by_library: dict[str, dict[str, str]] = {}
    columns: list[str] = []
    for text_path in paths:
        manifest_columns, rows = read_table(Path(text_path))
        columns.extend(column for column in manifest_columns if column not in columns)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one row in {text_path}, found {len(rows)}")
        row = rows[0]
        library_id = row.get("library_id", "")
        if not library_id:
            raise ValueError(f"{text_path}: library_id is empty")
        if library_id in by_library:
            raise ValueError(f"Duplicate featureCounts manifest for {library_id}")
        output = row.get("featurecounts_output", "")
        if not output or not Path(output).exists():
            raise ValueError(f"{library_id}: featurecounts_output does not exist: {output}")
        by_library[library_id] = row
    missing = [library_id for library_id in expected_ids if library_id not in by_library]
    extra = sorted(set(by_library) - set(expected_ids))
    if missing or extra:
        raise ValueError(f"featureCounts manifests do not match samples: missing={missing}, extra={extra}")
    return columns, [by_library[library_id] for library_id in expected_ids]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_matrices(
    counts_path: Path,
    metadata_path: Path,
    samples: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    annotations: dict[str, dict[str, str]],
) -> None:
    merged_metadata: dict[str, dict[str, str]] = {}
    counts_by_sample: dict[str, dict[str, int]] = {}
    for sample, manifest_row in zip(samples, manifest_rows, strict=True):
        metadata, counts = read_featurecounts_counts(Path(manifest_row["featurecounts_output"]))
        merged_metadata.update(metadata)
        counts_by_sample[sample["library_id"]] = counts
    gene_ids = sorted(merged_metadata)
    sample_ids = [row["library_id"] for row in samples]

    counts_path.parent.mkdir(parents=True, exist_ok=True)
    with counts_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_COLUMNS + sample_ids, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene_id in gene_ids:
            row = dict(merged_metadata[gene_id])
            for sample_id in sample_ids:
                row[sample_id] = str(counts_by_sample.get(sample_id, {}).get(gene_id, 0))
            writer.writerow(row)

    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=METADATA_COLUMNS + GENE_ANNOTATION_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for gene_id in gene_ids:
            row = dict(merged_metadata[gene_id])
            row.update(annotations.get(gene_id, {}))
            row["gene_display"] = gene_display_label(gene_id, row.get("gene_name", ""))
            writer.writerow({column: row.get(column, "") for column in METADATA_COLUMNS + GENE_ANNOTATION_COLUMNS})


def write_done(path: Path, samples: list[dict[str, str]]) -> None:
    paired = sum(1 for row in samples if row.get("layout") == "paired")
    single = sum(1 for row in samples if row.get("layout") == "single")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tsingle\tpaired\n")
        handle.write(f"ok\t{len(samples)}\t{single}\t{paired}\n")


def main() -> int:
    args = parse_args()
    _, samples = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    plan = read_plan(Path(args.plan))
    expected_ids = [row["library_id"] for row in samples]
    manifest_columns, manifest_rows = ordered_manifest_rows(args.manifests, expected_ids)
    write_matrices(
        Path(args.counts),
        Path(args.metadata),
        samples,
        manifest_rows,
        read_gene_annotations(Path(plan["annotation_gtf"])),
    )
    write_tsv(Path(args.manifest), manifest_columns, manifest_rows)
    write_done(Path(args.done), samples)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
