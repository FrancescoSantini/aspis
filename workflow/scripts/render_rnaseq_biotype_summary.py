#!/usr/bin/env python3
"""Summarize RNA-seq count and differential features by annotation biotype."""

from __future__ import annotations

import argparse
import csv
import html
import re
from collections import Counter, defaultdict
from pathlib import Path


BIOTYPE_KEYS = ("gene_biotype", "gene_type", "transcript_biotype", "transcript_type", "biotype")
STAT_COLUMNS = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
MANIFEST_COLUMNS = ["resource", "status", "path", "rows", "detail"]
COUNT_SUMMARY_COLUMNS = ["level", "biotype", "detected_features", "total_count"]
DIFF_SUMMARY_COLUMNS = ["level", "contrast_id", "biotype", "tested", "significant", "up", "down"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-gtf", default="")
    parser.add_argument("--gene-counts", default="")
    parser.add_argument("--gene-metadata", default="")
    parser.add_argument("--transcript-counts", default="")
    parser.add_argument("--transcript-metadata", default="")
    parser.add_argument("--gene-deseq2-manifest", default="")
    parser.add_argument("--transcript-deseq2-manifest", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--count-summary", required=True)
    parser.add_argument("--differential-summary", required=True)
    parser.add_argument("--html", required=True)
    parser.add_argument("--done", required=True)
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


def write_table(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


ATTR_RE = re.compile(r'([A-Za-z0-9_.-]+) "([^"]*)"')


def parse_attrs(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(text)}


def biotype_from_attrs(attrs: dict[str, str]) -> str:
    for key in BIOTYPE_KEYS:
        if attrs.get(key):
            return attrs[key]
    return "unclassified"


def read_gtf_biotypes(path_text: str) -> tuple[dict[str, str], dict[str, str]]:
    gene_biotypes: dict[str, str] = {}
    transcript_biotypes: dict[str, str] = {}
    if not path_text:
        return gene_biotypes, transcript_biotypes
    path = Path(path_text)
    if not path.exists():
        return gene_biotypes, transcript_biotypes
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            attrs = parse_attrs(fields[8])
            biotype = biotype_from_attrs(attrs)
            gene_id = attrs.get("gene_id") or attrs.get("ref_gene_id") or ""
            transcript_id = attrs.get("transcript_id") or attrs.get("ref_transcript_id") or ""
            if gene_id and feature_type in {"gene", "transcript", "exon"}:
                gene_biotypes.setdefault(gene_id, biotype)
            if transcript_id:
                transcript_biotypes.setdefault(transcript_id, biotype)
    return gene_biotypes, transcript_biotypes


def feature_id_column(columns: list[str], preferred: str) -> str:
    if preferred in columns:
        return preferred
    for column in columns:
        if column not in STAT_COLUMNS:
            return column
    return columns[0]


def metadata_biotypes(path_text: str, feature_column: str, gtf_biotypes: dict[str, str]) -> dict[str, str]:
    if not path_text:
        return {}
    columns, rows = read_table(Path(path_text))
    feature_column = feature_id_column(columns, feature_column)
    biotype_column = next((column for column in ["biotype", "gene_biotype", "gene_type", "feature_type"] if column in columns), "")
    gene_id_column = next((column for column in ["gene_id", "Geneid"] if column in columns), "")
    result = {}
    for row in rows:
        feature_id = row.get(feature_column, "")
        if not feature_id:
            continue
        biotype = row.get(biotype_column, "") if biotype_column else ""
        if not biotype and gene_id_column:
            biotype = gtf_biotypes.get(row.get(gene_id_column, ""), "")
        if not biotype:
            biotype = gtf_biotypes.get(feature_id, "unclassified")
        result[feature_id] = biotype or "unclassified"
    return result


def count_summary(path_text: str, level: str, feature_column: str, biotypes: dict[str, str]) -> list[dict[str, str]]:
    if not path_text:
        return []
    columns, rows = read_table(Path(path_text))
    feature_column = feature_id_column(columns, feature_column)
    sample_columns = [column for column in columns if column != feature_column and column not in {"Chr", "Start", "End", "Strand", "Length"}]
    detected: Counter[str] = Counter()
    total: Counter[str] = Counter()
    for row in rows:
        feature_id = row.get(feature_column, "")
        biotype = biotypes.get(feature_id, "unclassified")
        values = []
        for column in sample_columns:
            try:
                values.append(int(float(row.get(column, "0") or "0")))
            except ValueError:
                values.append(0)
        if any(value > 0 for value in values):
            detected[biotype] += 1
        total[biotype] += sum(values)
    return [
        {"level": level, "biotype": biotype, "detected_features": str(detected[biotype]), "total_count": str(total[biotype])}
        for biotype in sorted(set(detected) | set(total))
    ]


def direction(row: dict[str, str]) -> str:
    try:
        log2fc = float(row.get("log2FoldChange", "0") or "0")
    except ValueError:
        return "unchanged"
    if log2fc > 0:
        return "up"
    if log2fc < 0:
        return "down"
    return "unchanged"


def manifest_rows(path_text: str) -> list[dict[str, str]]:
    if not path_text:
        return []
    path = Path(path_text)
    if not path.exists():
        return []
    _columns, rows = read_table(path, {"contrast_id", "status", "results", "filtered"})
    return [row for row in rows if row.get("status") == "ok"]


def differential_summary(path_text: str, level: str, feature_column: str, biotypes: dict[str, str]) -> list[dict[str, str]]:
    rows = []
    for manifest_row in manifest_rows(path_text):
        contrast_id = manifest_row["contrast_id"]
        result_columns, result_rows = read_table(Path(manifest_row["results"]))
        filtered_columns, filtered_rows = read_table(Path(manifest_row["filtered"]))
        result_feature = feature_id_column(result_columns, feature_column)
        filtered_feature = feature_id_column(filtered_columns, feature_column)
        tested = Counter(biotypes.get(row.get(result_feature, ""), "unclassified") for row in result_rows)
        significant = Counter(biotypes.get(row.get(filtered_feature, ""), "unclassified") for row in filtered_rows)
        up = Counter(
            biotypes.get(row.get(filtered_feature, ""), "unclassified")
            for row in filtered_rows
            if direction(row) == "up"
        )
        down = Counter(
            biotypes.get(row.get(filtered_feature, ""), "unclassified")
            for row in filtered_rows
            if direction(row) == "down"
        )
        for biotype in sorted(set(tested) | set(significant) | set(up) | set(down)):
            rows.append(
                {
                    "level": level,
                    "contrast_id": contrast_id,
                    "biotype": biotype,
                    "tested": str(tested[biotype]),
                    "significant": str(significant[biotype]),
                    "up": str(up[biotype]),
                    "down": str(down[biotype]),
                }
            )
    return rows


def html_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    if not rows:
        return "<p>No rows.</p>"
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(row.get(column, ''))}</td>" for column in columns) + "</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def write_html(path: Path, count_rows: list[dict[str, str]], diff_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>RNA-seq biotype summary</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
    th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; }}
    th {{ background: #f2f2f2; }}
  </style>
</head>
<body>
  <h1>RNA-seq biotype summary</h1>
  <h2>Detected/count composition</h2>
  {html_table(count_rows, COUNT_SUMMARY_COLUMNS)}
  <h2>Differential composition</h2>
  {html_table(diff_rows, DIFF_SUMMARY_COLUMNS)}
</body>
</html>
"""
    path.write_text(content, encoding="utf-8")


def write_done(path: Path, count_rows: list[dict[str, str]], diff_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tcount_biotype_rows\tdifferential_biotype_rows\n")
        handle.write(f"ok\t{len(count_rows)}\t{len(diff_rows)}\n")


def main() -> int:
    args = parse_args()
    gene_gtf, transcript_gtf = read_gtf_biotypes(args.annotation_gtf)
    gene_biotypes = metadata_biotypes(args.gene_metadata, "Geneid", gene_gtf)
    transcript_biotypes = metadata_biotypes(args.transcript_metadata, "transcript_id", transcript_gtf or gene_gtf)
    count_rows = []
    count_rows.extend(count_summary(args.gene_counts, "gene", "Geneid", gene_biotypes))
    count_rows.extend(count_summary(args.transcript_counts, "transcript", "transcript_id", transcript_biotypes))
    diff_rows = []
    diff_rows.extend(differential_summary(args.gene_deseq2_manifest, "gene", "Geneid", gene_biotypes))
    diff_rows.extend(differential_summary(args.transcript_deseq2_manifest, "transcript", "transcript_id", transcript_biotypes))
    write_table(Path(args.count_summary), COUNT_SUMMARY_COLUMNS, count_rows)
    write_table(Path(args.differential_summary), DIFF_SUMMARY_COLUMNS, diff_rows)
    write_html(Path(args.html), count_rows, diff_rows)
    write_table(
        Path(args.manifest),
        MANIFEST_COLUMNS,
        [
            {"resource": "count_biotype_summary", "status": "ok", "path": args.count_summary, "rows": str(len(count_rows)), "detail": ""},
            {"resource": "differential_biotype_summary", "status": "ok", "path": args.differential_summary, "rows": str(len(diff_rows)), "detail": ""},
            {"resource": "biotype_summary_html", "status": "ok", "path": args.html, "rows": str(len(count_rows) + len(diff_rows)), "detail": ""},
        ],
    )
    write_done(Path(args.done), count_rows, diff_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
