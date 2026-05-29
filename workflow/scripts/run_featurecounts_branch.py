#!/usr/bin/env python3
"""Run featureCounts per RNA-seq sample and merge gene-level count matrices."""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project", "layout", "bam"}
REQUIRED_PLAN_COLUMNS = {"project", "assay", "status", "gene_counter", "annotation_gtf"}
METADATA_COLUMNS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
GENE_ANNOTATION_COLUMNS = ["gene_name", "gene_biotype"]
ATTR_RE = re.compile(r'([A-Za-z0-9_.-]+) "([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Aligned RNA-seq samples TSV")
    parser.add_argument("--plan", required=True, help="RNA-seq quantification plan TSV")
    parser.add_argument("--outdir", required=True, help="featureCounts per-sample output directory")
    parser.add_argument("--counts", required=True, help="Merged gene count matrix TSV")
    parser.add_argument("--metadata", required=True, help="Gene metadata TSV")
    parser.add_argument("--manifest", required=True, help="featureCounts output manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--featurecounts", default="featureCounts", help="featureCounts executable")
    parser.add_argument("--threads", type=int, default=1, help="Threads per sample")
    parser.add_argument("--single-extra-args", default="", help="Extra args for single-end BAMs")
    parser.add_argument("--paired-extra-args", default="-p --countReadPairs", help="Extra args for paired BAMs")
    parser.add_argument("--extra-args", default="", help="Extra args for every featureCounts run")
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
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


def validate_samples(rows: list[dict[str, str]], plan: dict[str, str]) -> None:
    errors = []
    if not rows:
        errors.append("aligned sample table has no rows")
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")
        if row.get("project") != plan.get("project"):
            errors.append(f"{library_id}: project does not match plan")
        if row.get("layout") not in {"single", "paired"}:
            errors.append(f"{library_id}: unsupported layout {row.get('layout')!r}")
        bam = row.get("bam", "")
        if not bam:
            errors.append(f"{library_id}: bam is empty")
        elif not Path(bam).exists():
            errors.append(f"{library_id}: bam does not exist: {bam}")
    if errors:
        raise ValueError("featureCounts cannot start:\n- " + "\n- ".join(errors))


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def parse_attrs(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(text)}


def read_gene_annotations(gtf_path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}
    if not gtf_path.exists():
        return annotations
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


def run_featurecounts(
    row: dict[str, str],
    plan: dict[str, str],
    outdir: Path,
    featurecounts: str,
    threads: int,
    single_extra_args: str,
    paired_extra_args: str,
    extra_args: str,
) -> dict[str, str]:
    sample_dir = outdir / row["library_id"]
    sample_dir.mkdir(parents=True, exist_ok=True)
    raw_output = sample_dir / "featurecounts.tsv"
    log_path = sample_dir / "featurecounts.log"
    summary_path = Path(str(raw_output) + ".summary")
    for path in (raw_output, log_path, summary_path):
        if path.exists():
            path.unlink()

    command = [
        featurecounts,
        "-T",
        str(max(1, threads)),
        "-a",
        plan["annotation_gtf"],
        "-o",
        str(raw_output),
    ]
    if row["layout"] == "paired":
        command.extend(shlex.split(paired_extra_args))
    else:
        command.extend(shlex.split(single_extra_args))
    command.extend(shlex.split(extra_args))
    command.append(row["bam"])

    print("[CMD] " + shlex.join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    with log_path.open("w", encoding="utf-8") as handle:
        if completed.stdout:
            handle.write(completed.stdout)
        if completed.stderr:
            handle.write(completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(f"{row['library_id']}: featureCounts exited with status {completed.returncode}")
    if not raw_output.exists():
        raise RuntimeError(f"{row['library_id']}: featureCounts did not produce {raw_output}")

    return {
        "library_id": row["library_id"],
        "layout": row["layout"],
        "bam": row["bam"],
        "featurecounts_output": str(raw_output),
        "featurecounts_summary": str(summary_path) if summary_path.exists() else "",
        "featurecounts_log": str(log_path),
        "status": "ok",
        "message": "",
    }


def read_featurecounts_counts(path: Path, library_id: str) -> tuple[dict[str, dict[str, str]], dict[str, int]]:
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
            raw_count = row.get(sample_column, "0")
            counts[gene_id] = int(float(raw_count))
    return metadata, counts


def write_merged_counts(
    counts_path: Path,
    metadata_path: Path,
    rows: list[dict[str, str]],
    fc_rows: list[dict[str, str]],
    gene_annotations: dict[str, dict[str, str]],
) -> None:
    merged_metadata: dict[str, dict[str, str]] = {}
    counts_by_sample: dict[str, dict[str, int]] = {}
    for sample, fc_row in zip(rows, fc_rows, strict=True):
        sample_metadata, sample_counts = read_featurecounts_counts(
            Path(fc_row["featurecounts_output"]),
            sample["library_id"],
        )
        merged_metadata.update(sample_metadata)
        counts_by_sample[sample["library_id"]] = sample_counts

    gene_ids = sorted(merged_metadata)
    sample_ids = [row["library_id"] for row in rows]

    counts_path.parent.mkdir(parents=True, exist_ok=True)
    with counts_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=METADATA_COLUMNS + sample_ids,
            delimiter="\t",
            lineterminator="\n",
        )
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
            output_row = dict(merged_metadata[gene_id])
            output_row.update(gene_annotations.get(gene_id, {}))
            writer.writerow({column: output_row.get(column, "") for column in METADATA_COLUMNS + GENE_ANNOTATION_COLUMNS})


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "library_id",
        "layout",
        "bam",
        "featurecounts_output",
        "featurecounts_summary",
        "featurecounts_log",
        "status",
        "message",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    paired = sum(1 for row in rows if row.get("layout") == "paired")
    single = sum(1 for row in rows if row.get("layout") == "single")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tsingle\tpaired\n")
        handle.write(f"ok\t{len(rows)}\t{single}\t{paired}\n")


def main() -> int:
    args = parse_args()
    _, rows = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    plan = read_plan(Path(args.plan))
    validate_samples(rows, plan)
    featurecounts = executable_path(args.featurecounts)

    outdir = Path(args.outdir)
    fc_rows = [
        run_featurecounts(
            row,
            plan,
            outdir,
            featurecounts,
            args.threads,
            args.single_extra_args,
            args.paired_extra_args,
            args.extra_args,
        )
        for row in rows
    ]
    write_merged_counts(
        Path(args.counts),
        Path(args.metadata),
        rows,
        fc_rows,
        read_gene_annotations(Path(plan["annotation_gtf"])),
    )
    write_manifest(Path(args.manifest), fc_rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
