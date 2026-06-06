#!/usr/bin/env python3
"""Align miRBase-unmapped smallRNA reads to a genome and classify biotypes."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shlex
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_COLUMNS = {"library_id", "assay", "project", "layout", "mirbase_unmapped_fastq_1"}
ADDED_COLUMNS = [
    "residual_pre_alignment_fastq_1",
    "residual_bam",
    "residual_unmapped_fastq_1",
    "residual_flagstat",
    "residual_alignment_log",
    "residual_alignment_tool",
]
MANIFEST_COLUMNS = [
    "library_id",
    "project",
    "assay",
    "input_fastq_1",
    "genome_unmapped_fastq_1",
    "bam",
    "flagstat",
    "alignment_log",
    "assignment_tsv",
    "status",
    "message",
    "input_reads",
    "genome_aligned_reads",
    "genome_unmapped_reads",
    "annotated_reads",
    "unassigned_reads",
]
REF_CIGAR_OPS = {"M", "D", "N", "=", "X"}
BIOTYPE_KEYS = ("gene_biotype", "gene_type", "transcript_biotype", "transcript_type", "biotype")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="miRBase-aligned smallRNA sample table TSV")
    parser.add_argument("--library-id", default="", help="Optional single library to process")
    parser.add_argument("--outdir", required=True, help="Residual genome output directory")
    parser.add_argument("--output", required=True, help="Residual-aligned sample table TSV")
    parser.add_argument("--manifest", required=True, help="Residual alignment manifest TSV")
    parser.add_argument("--biotype-counts", required=True, help="Output biotype count matrix TSV")
    parser.add_argument("--feature-counts", required=True, help="Output feature count matrix TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--index-prefix", required=True, help="Genome Bowtie index prefix")
    parser.add_argument("--annotation-gtf", default="", help="GTF used to classify residual genome alignments")
    parser.add_argument("--bowtie", default="bowtie", help="Bowtie executable")
    parser.add_argument("--samtools", default="samtools", help="samtools executable")
    parser.add_argument("--threads", type=int, default=1, help="Threads per sample")
    parser.add_argument("--mismatches", type=int, default=1, help="Bowtie -v mismatches")
    parser.add_argument("--multi-alignments", type=int, default=10, help="Bowtie -k alignments per read")
    parser.add_argument("--extra-args", default="--best --strata", help="Extra Bowtie arguments")
    return parser.parse_args()


def read_tsv(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def select_library(rows: list[dict[str, str]], library_id: str) -> list[dict[str, str]]:
    if not library_id:
        return rows
    matches = [row for row in rows if row.get("library_id") == library_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one row for {library_id!r}, found {len(matches)}")
    return matches


def validate_args(args: argparse.Namespace) -> list[str]:
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.mismatches < 0:
        raise ValueError("--mismatches cannot be negative")
    if args.multi_alignments < 1:
        raise ValueError("--multi-alignments must be >= 1")
    try:
        return shlex.split(args.extra_args)
    except ValueError as exc:
        raise ValueError(f"--extra-args is not valid shell-like syntax: {exc}") from exc


def validate_samples(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("smallRNA residual alignment received an empty sample table")
    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "smallrna":
            errors.append(f"{library_id}: expected assay='smallrna', got {row.get('assay')!r}")
        if row.get("layout") != "single":
            errors.append(f"{library_id}: residual smallRNA alignment expects single-end libraries")
        fastq = row.get("mirbase_unmapped_fastq_1", "")
        if not fastq:
            errors.append(f"{library_id}: mirbase_unmapped_fastq_1 is empty")
        elif not Path(fastq).exists():
            errors.append(f"{library_id}: miRBase-unmapped FASTQ does not exist: {fastq}")
    if errors:
        raise ValueError("smallRNA residual genome alignment cannot start:\n- " + "\n- ".join(errors))


def output_columns(input_columns: list[str]) -> list[str]:
    columns = list(input_columns)
    for column in ADDED_COLUMNS:
        if column not in columns:
            columns.append(column)
    return columns


def outputs_for(row: dict[str, str], outdir: Path) -> dict[str, Path]:
    library_dir = outdir / row["library_id"]
    return {
        "library_dir": library_dir,
        "sam": library_dir / "residual_genome.sam",
        "unsorted_bam": library_dir / "residual_genome.unsorted.bam",
        "bam": library_dir / "residual_genome.bam",
        "unmapped_fastq_1": library_dir / "genome_unmapped.fastq.gz",
        "tmp_unmapped_fastq_1": library_dir / "genome_unmapped.fastq",
        "flagstat": library_dir / "flagstat.txt",
        "alignment_log": library_dir / "bowtie.log",
        "assignment_tsv": library_dir / "residual_assignments.tsv",
    }


def remove_stale_outputs(outputs: dict[str, Path]) -> None:
    for key, path in outputs.items():
        if key == "library_dir":
            continue
        if path.exists():
            path.unlink()


def count_fastq_records(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(1 for _line in handle) // 4


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def log_tail(path: Path, max_lines: int = 20) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-max_lines:])


def run_command(command: list[str], *, stdout: Path | None = None, stderr: Path | None = None) -> None:
    print("[CMD] " + shlex.join(command))
    stdout_handle = stdout.open("w", encoding="utf-8") if stdout else None
    stderr_handle = stderr.open("w", encoding="utf-8") if stderr else None
    try:
        completed = subprocess.run(command, check=False, stdout=stdout_handle, stderr=stderr_handle, text=True)
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()
    if completed.returncode != 0:
        message = log_tail(stderr) if stderr else ""
        detail = f"\n{message}" if message else ""
        raise RuntimeError(f"{Path(command[0]).name} exited with status {completed.returncode}{detail}")


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in text.rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, value = part.split(" ", 1)
            attributes[key] = value.strip().strip('"')
        elif "=" in part:
            key, value = part.split("=", 1)
            attributes[key] = value.strip().strip('"')
    return attributes


def feature_biotype(attributes: dict[str, str]) -> str:
    for key in BIOTYPE_KEYS:
        value = attributes.get(key, "")
        if value:
            return value
    return "unannotated"


def read_gtf(path_text: str) -> dict[str, list[dict[str, object]]]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"Residual annotation GTF does not exist: {path}")
    features: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type not in {"gene", "transcript", "exon"}:
                continue
            attributes = parse_attributes(fields[8])
            feature_id = attributes.get("gene_id") or attributes.get("transcript_id") or fields[0]
            features[fields[0]].append(
                {
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "feature_id": feature_id,
                    "feature_name": attributes.get("gene_name", feature_id),
                    "biotype": feature_biotype(attributes),
                    "feature_type": feature_type,
                }
            )
    for contig in features:
        features[contig].sort(key=lambda item: (item["feature_type"] != "gene", item["start"], item["end"]))
    return features


def cigar_ref_length(cigar: str, sequence: str) -> int:
    if not cigar or cigar == "*":
        return max(1, len(sequence))
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in REF_CIGAR_OPS:
            total += int(length)
    return max(1, total)


def assign_feature(
    contig: str,
    start: int,
    end: int,
    features: dict[str, list[dict[str, object]]],
) -> tuple[str, str, str]:
    for feature in features.get(contig, []):
        if int(feature["end"]) < start or int(feature["start"]) > end:
            continue
        return str(feature["feature_id"]), str(feature["feature_name"]), str(feature["biotype"])
    return "", "", "unassigned"


def parse_sam_assignments(
    sam: Path,
    gtf_features: dict[str, list[dict[str, object]]],
) -> tuple[list[dict[str, str]], Counter[str], Counter[str]]:
    seen_reads = set()
    assignments: list[dict[str, str]] = []
    biotype_counts: Counter[str] = Counter()
    feature_counts: Counter[str] = Counter()
    with sam.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            read_id = fields[0]
            if read_id in seen_reads:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            seen_reads.add(read_id)
            contig = fields[2]
            start = int(fields[3])
            end = start + cigar_ref_length(fields[5], fields[9]) - 1
            feature_id, feature_name, biotype = assign_feature(contig, start, end, gtf_features)
            assignments.append(
                {
                    "read_id": read_id,
                    "contig": contig,
                    "start": str(start),
                    "end": str(end),
                    "feature_id": feature_id,
                    "feature_name": feature_name,
                    "biotype": biotype,
                }
            )
            biotype_counts[biotype] += 1
            if feature_id:
                feature_counts["\t".join([feature_id, feature_name, biotype])] += 1
    return assignments, biotype_counts, feature_counts


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_empty_fastq(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb"):
        pass


def gzip_fastq(source: Path, target: Path) -> None:
    with source.open("rb") as input_handle, gzip.open(target, "wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle)
    source.unlink(missing_ok=True)


def run_alignment(
    row: dict[str, str],
    outputs: dict[str, Path],
    args: argparse.Namespace,
    extra_args: list[str],
    gtf_features: dict[str, list[dict[str, object]]],
) -> tuple[dict[str, str], Counter[str], Counter[str]]:
    input_fastq = Path(row["mirbase_unmapped_fastq_1"])
    input_reads = count_fastq_records(input_fastq)
    outputs["library_dir"].mkdir(parents=True, exist_ok=True)
    remove_stale_outputs(outputs)

    if input_reads == 0:
        outputs["sam"].write_text("", encoding="utf-8")
        outputs["flagstat"].write_text("0 residual reads\n", encoding="utf-8")
        outputs["alignment_log"].write_text("No miRBase-unmapped reads to align\n", encoding="utf-8")
        write_empty_fastq(outputs["unmapped_fastq_1"])
        assignments: list[dict[str, str]] = []
        biotype_counts: Counter[str] = Counter()
        feature_counts: Counter[str] = Counter()
        bam_path = ""
    else:
        bowtie = executable_path(args.bowtie)
        samtools = executable_path(args.samtools)
        command = [
            bowtie,
            "-v",
            str(args.mismatches),
            "-k",
            str(args.multi_alignments),
            "-p",
            str(args.threads),
            "--un",
            str(outputs["tmp_unmapped_fastq_1"]),
            "-S",
            *extra_args,
            args.index_prefix,
            str(input_fastq),
        ]
        run_command(command, stdout=outputs["sam"], stderr=outputs["alignment_log"])
        run_command([samtools, "view", "-bS", "-o", str(outputs["unsorted_bam"]), str(outputs["sam"])])
        run_command(
            [
                samtools,
                "sort",
                "-@",
                str(args.threads),
                "-o",
                str(outputs["bam"]),
                str(outputs["unsorted_bam"]),
            ]
        )
        run_command([samtools, "flagstat", str(outputs["bam"])], stdout=outputs["flagstat"])
        if outputs["tmp_unmapped_fastq_1"].exists():
            gzip_fastq(outputs["tmp_unmapped_fastq_1"], outputs["unmapped_fastq_1"])
        else:
            write_empty_fastq(outputs["unmapped_fastq_1"])
        assignments, biotype_counts, feature_counts = parse_sam_assignments(outputs["sam"], gtf_features)
        outputs["unsorted_bam"].unlink(missing_ok=True)
        bam_path = str(outputs["bam"])

    write_tsv(
        outputs["assignment_tsv"],
        ["read_id", "contig", "start", "end", "feature_id", "feature_name", "biotype"],
        assignments,
    )
    genome_unmapped_reads = count_fastq_records(outputs["unmapped_fastq_1"])
    aligned_reads = max(0, input_reads - genome_unmapped_reads)
    unassigned_reads = biotype_counts.get("unassigned", 0)
    return (
        {
            "library_id": row["library_id"],
            "project": row.get("project", ""),
            "assay": row.get("assay", ""),
            "input_fastq_1": str(input_fastq),
            "genome_unmapped_fastq_1": str(outputs["unmapped_fastq_1"]),
            "bam": bam_path,
            "flagstat": str(outputs["flagstat"]),
            "alignment_log": str(outputs["alignment_log"]),
            "assignment_tsv": str(outputs["assignment_tsv"]),
            "status": "ok",
            "message": "",
            "input_reads": str(input_reads),
            "genome_aligned_reads": str(aligned_reads),
            "genome_unmapped_reads": str(genome_unmapped_reads),
            "annotated_reads": str(aligned_reads - unassigned_reads),
            "unassigned_reads": str(unassigned_reads),
        },
        biotype_counts,
        feature_counts,
    )


def write_count_matrix(path: Path, first_columns: list[str], rows: list[dict[str, str]]) -> None:
    write_tsv(path, first_columns, rows)


def write_done(path: Path, manifest_rows: list[dict[str, str]]) -> None:
    total_input = sum(int(row["input_reads"]) for row in manifest_rows)
    total_aligned = sum(int(row["genome_aligned_reads"]) for row in manifest_rows)
    total_unmapped = sum(int(row["genome_unmapped_reads"]) for row in manifest_rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tinput_reads\tgenome_aligned_reads\tgenome_unmapped_reads\n")
        handle.write(f"ok\t{len(manifest_rows)}\t{total_input}\t{total_aligned}\t{total_unmapped}\n")


def main() -> int:
    args = parse_args()
    extra_args = validate_args(args)
    input_columns, rows = read_tsv(Path(args.samples), REQUIRED_COLUMNS)
    rows = select_library(rows, args.library_id)
    validate_samples(rows)
    gtf_features = read_gtf(args.annotation_gtf)

    outdir = Path(args.outdir)
    output_rows = []
    manifest_rows = []
    biotypes_by_sample: dict[str, Counter[str]] = {}
    features_by_sample: dict[str, Counter[str]] = {}
    for row in rows:
        outputs = outputs_for(row, outdir)
        manifest_row, biotype_counts, feature_counts = run_alignment(row, outputs, args, extra_args, gtf_features)
        sample_id = row["library_id"]
        output_rows.append(
            {
                **row,
                "residual_pre_alignment_fastq_1": row["mirbase_unmapped_fastq_1"],
                "residual_bam": manifest_row["bam"],
                "residual_unmapped_fastq_1": manifest_row["genome_unmapped_fastq_1"],
                "residual_flagstat": manifest_row["flagstat"],
                "residual_alignment_log": manifest_row["alignment_log"],
                "residual_alignment_tool": "bowtie",
            }
        )
        manifest_rows.append(manifest_row)
        biotypes_by_sample[sample_id] = biotype_counts
        features_by_sample[sample_id] = feature_counts

    sample_ids = [row["library_id"] for row in rows]
    biotypes = sorted({key for counts in biotypes_by_sample.values() for key in counts} | {"genome_unmapped"})
    biotype_rows = []
    for biotype in biotypes:
        output_row = {"biotype": biotype}
        for sample_id in sample_ids:
            if biotype == "genome_unmapped":
                manifest_row = next(row for row in manifest_rows if row["library_id"] == sample_id)
                output_row[sample_id] = manifest_row["genome_unmapped_reads"]
            else:
                output_row[sample_id] = str(biotypes_by_sample[sample_id].get(biotype, 0))
        biotype_rows.append(output_row)

    feature_keys = sorted({key for counts in features_by_sample.values() for key in counts})
    feature_rows = []
    for key in feature_keys:
        feature_id, feature_name, biotype = key.split("\t", 2)
        output_row = {"feature_id": feature_id, "feature_name": feature_name, "biotype": biotype}
        for sample_id in sample_ids:
            output_row[sample_id] = str(features_by_sample[sample_id].get(key, 0))
        feature_rows.append(output_row)

    write_tsv(Path(args.output), output_columns(input_columns), output_rows)
    write_tsv(Path(args.manifest), MANIFEST_COLUMNS, manifest_rows)
    write_count_matrix(Path(args.biotype_counts), ["biotype", *sample_ids], biotype_rows)
    write_count_matrix(Path(args.feature_counts), ["feature_id", "feature_name", "biotype", *sample_ids], feature_rows)
    write_done(Path(args.done), manifest_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
