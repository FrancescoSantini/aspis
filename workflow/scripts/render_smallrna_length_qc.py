#!/usr/bin/env python3
"""Render smallRNA read-length, arm, and mapped-length QC summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


MANIFEST_COLUMNS = ["resource", "status", "path", "rows", "detail"]
LENGTH_COLUMNS = ["stage", "library_id", "length", "reads", "fraction"]
STAGE_COLUMNS = [
    "stage",
    "library_id",
    "total_reads",
    "modal_length",
    "mean_length",
    "min_length",
    "max_length",
    "n_lengths",
]
ARM_COLUMNS = ["arm", "detected_mirnas", "total_count", "fraction"]
ISOMIR_LENGTH_COLUMNS = ["length", "estimated_mirbase_mapped_reads", "fraction"]
STAGE_PATH_COLUMNS = {
    "raw": "fastq_1",
    "trimmed": "fastq_1",
    "depleted": "fastq_1",
    "mirbase_unmapped": "mirbase_unmapped_fastq_1",
}
ARM_RE = re.compile(r"(?i)(?:^|[-_])([53]p)(?:$|[-_])")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-samples", required=True)
    parser.add_argument("--trimmed-samples", required=True)
    parser.add_argument("--depleted-samples", required=True)
    parser.add_argument("--aligned-samples", required=True)
    parser.add_argument("--mirna-counts", required=True)
    parser.add_argument("--mirna-metadata", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--length-distribution", required=True)
    parser.add_argument("--stage-summary", required=True)
    parser.add_argument("--arm-summary", required=True)
    parser.add_argument("--isomir-length-summary", required=True)
    parser.add_argument("--length-plot", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--max-reads", type=int, default=200000)
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


def fastq_lengths(path: Path, max_reads: int) -> Counter[int]:
    counts: Counter[int] = Counter()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            if index % 4 == 1:
                counts[len(line.strip())] += 1
                if sum(counts.values()) >= max_reads:
                    break
    return counts


def stage_counts(stage: str, rows: list[dict[str, str]], path_column: str, max_reads: int) -> dict[str, Counter[int]]:
    counts: dict[str, Counter[int]] = {}
    for row in rows:
        library_id = row.get("library_id", "")
        fastq = row.get(path_column, "")
        if not library_id or not fastq:
            continue
        path = Path(fastq)
        if not path.exists():
            raise FileNotFoundError(f"{stage} FASTQ for {library_id} does not exist: {fastq}")
        counts[library_id] = fastq_lengths(path, max_reads)
    return counts


def length_rows(all_counts: dict[str, dict[str, Counter[int]]]) -> list[dict[str, str]]:
    rows = []
    for stage, by_library in all_counts.items():
        for library_id, counts in by_library.items():
            total = sum(counts.values()) or 1
            for length in sorted(counts):
                rows.append(
                    {
                        "stage": stage,
                        "library_id": library_id,
                        "length": str(length),
                        "reads": str(counts[length]),
                        "fraction": f"{counts[length] / total:.8g}",
                    }
                )
    return rows


def stage_summary_rows(all_counts: dict[str, dict[str, Counter[int]]]) -> list[dict[str, str]]:
    rows = []
    for stage, by_library in all_counts.items():
        for library_id, counts in by_library.items():
            total = sum(counts.values())
            if total:
                modal_length = max(counts, key=lambda length: (counts[length], -length))
                mean_length = sum(length * value for length, value in counts.items()) / total
                min_length = min(counts)
                max_length = max(counts)
            else:
                modal_length = mean_length = min_length = max_length = 0
            rows.append(
                {
                    "stage": stage,
                    "library_id": library_id,
                    "total_reads": str(total),
                    "modal_length": str(modal_length),
                    "mean_length": f"{float(mean_length):.6g}",
                    "min_length": str(min_length),
                    "max_length": str(max_length),
                    "n_lengths": str(len(counts)),
                }
            )
    return rows


def feature_id_column(columns: list[str]) -> str:
    for column in ["Geneid", "mirna_id", "feature_id", "id"]:
        if column in columns:
            return column
    return columns[0]


def feature_arm(feature_id: str) -> str:
    match = ARM_RE.search(feature_id)
    if not match:
        return "unannotated_arm"
    return match.group(1).lower()


def parse_count(value: str) -> int:
    try:
        return int(float(value))
    except ValueError:
        return 0


def arm_summary_rows(counts_path: Path) -> list[dict[str, str]]:
    columns, rows = read_table(counts_path)
    feature_column = feature_id_column(columns)
    sample_columns = [
        column
        for column in columns
        if column
        not in {
            feature_column,
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
            "feature_type",
        }
    ]
    detected: Counter[str] = Counter()
    total: Counter[str] = Counter()
    for row in rows:
        feature_id = row.get(feature_column, "")
        if not feature_id:
            continue
        arm = feature_arm(feature_id)
        feature_total = sum(parse_count(row.get(sample, "")) for sample in sample_columns)
        if feature_total > 0:
            detected[arm] += 1
            total[arm] += feature_total
    grand_total = sum(total.values()) or 1
    return [
        {
            "arm": arm,
            "detected_mirnas": str(detected[arm]),
            "total_count": str(total[arm]),
            "fraction": f"{total[arm] / grand_total:.8g}",
        }
        for arm in sorted(set(detected) | set(total))
    ]


def aggregate_by_length(by_library: dict[str, Counter[int]]) -> Counter[int]:
    total: Counter[int] = Counter()
    for counts in by_library.values():
        total.update(counts)
    return total


def isomir_length_rows(all_counts: dict[str, dict[str, Counter[int]]]) -> list[dict[str, str]]:
    depleted = aggregate_by_length(all_counts.get("depleted", {}))
    unmapped = aggregate_by_length(all_counts.get("mirbase_unmapped", {}))
    estimated = Counter({length: max(0, depleted[length] - unmapped[length]) for length in set(depleted) | set(unmapped)})
    total = sum(estimated.values()) or 1
    return [
        {
            "length": str(length),
            "estimated_mirbase_mapped_reads": str(estimated[length]),
            "fraction": f"{estimated[length] / total:.8g}",
        }
        for length in sorted(estimated)
        if estimated[length] > 0
    ]


def write_length_svg(path: Path, all_counts: dict[str, dict[str, Counter[int]]]) -> None:
    aggregate = {stage: aggregate_by_length(counts) for stage, counts in all_counts.items()}
    lengths = sorted({length for counts in aggregate.values() for length in counts})
    width = 980
    height = 420
    left = 70
    right = 30
    top = 42
    bottom = 72
    plot_width = width - left - right
    plot_height = height - top - bottom
    colors = {
        "raw": "#4d4d4d",
        "trimmed": "#2166ac",
        "depleted": "#1b7837",
        "mirbase_unmapped": "#b2182b",
    }
    max_reads = max((max(counts.values()) for counts in aggregate.values() if counts), default=1)
    min_length = min(lengths) if lengths else 0
    max_length = max(lengths) if lengths else 1

    def scale_x(length: int) -> float:
        if max_length == min_length:
            return left + plot_width / 2
        return left + (length - min_length) * plot_width / (max_length - min_length)

    def scale_y(reads: int) -> float:
        return top + plot_height - (reads / max_reads) * plot_height

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="28" y="26" font-family="sans-serif" font-size="18" font-weight="700">smallRNA read-length distribution</text>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#777"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#777"/>',
    ]
    for index, (stage, counts) in enumerate(aggregate.items()):
        if not counts:
            continue
        points = " ".join(f"{scale_x(length):.1f},{scale_y(counts[length]):.1f}" for length in sorted(counts))
        color = colors.get(stage, "#4d4d4d")
        elements.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        elements.append(
            f'<text x="{left + index * 170}" y="{height - 26}" font-family="sans-serif" font-size="12" fill="{color}">{html.escape(stage)}</text>'
        )
    elements.extend(
        [
            f'<text x="{left}" y="{height - 50}" font-family="sans-serif" font-size="12">{min_length}</text>',
            f'<text x="{left + plot_width - 20}" y="{height - 50}" font-family="sans-serif" font-size="12">{max_length}</text>',
            f'<text x="{width / 2 - 45:.1f}" y="{height - 16}" font-family="sans-serif" font-size="12">read length</text>',
            '<text x="8" y="210" transform="rotate(-90 8 210)" font-family="sans-serif" font-size="12">reads</text>',
            "</svg>",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    write_table(path, MANIFEST_COLUMNS, rows)


def write_done(path: Path, length_rows_count: int, arm_rows_count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlength_rows\tarm_rows\n")
        handle.write(f"ok\t{length_rows_count}\t{arm_rows_count}\n")


def main() -> int:
    args = parse_args()
    if args.max_reads < 1:
        raise ValueError("--max-reads must be >= 1")
    stage_specs = {
        "raw": (Path(args.raw_samples), STAGE_PATH_COLUMNS["raw"]),
        "trimmed": (Path(args.trimmed_samples), STAGE_PATH_COLUMNS["trimmed"]),
        "depleted": (Path(args.depleted_samples), STAGE_PATH_COLUMNS["depleted"]),
        "mirbase_unmapped": (Path(args.aligned_samples), STAGE_PATH_COLUMNS["mirbase_unmapped"]),
    }
    all_counts: dict[str, dict[str, Counter[int]]] = {}
    for stage, (path, column) in stage_specs.items():
        _columns, rows = read_table(path, {"library_id", column})
        all_counts[stage] = stage_counts(stage, rows, column, args.max_reads)

    length = length_rows(all_counts)
    stage_summary = stage_summary_rows(all_counts)
    arm = arm_summary_rows(Path(args.mirna_counts))
    isomir_lengths = isomir_length_rows(all_counts)
    write_table(Path(args.length_distribution), LENGTH_COLUMNS, length)
    write_table(Path(args.stage_summary), STAGE_COLUMNS, stage_summary)
    write_table(Path(args.arm_summary), ARM_COLUMNS, arm)
    write_table(Path(args.isomir_length_summary), ISOMIR_LENGTH_COLUMNS, isomir_lengths)
    write_length_svg(Path(args.length_plot), all_counts)
    manifest_rows = [
        {"resource": "length_distribution", "status": "ok", "path": args.length_distribution, "rows": str(len(length)), "detail": "per-sample read-length histograms by processing stage"},
        {"resource": "stage_summary", "status": "ok", "path": args.stage_summary, "rows": str(len(stage_summary)), "detail": "modal and mean length by sample and stage"},
        {"resource": "arm_summary", "status": "ok", "path": args.arm_summary, "rows": str(len(arm)), "detail": "5p/3p abundance inferred from miRNA feature names"},
        {"resource": "isomir_length_summary", "status": "ok", "path": args.isomir_length_summary, "rows": str(len(isomir_lengths)), "detail": "estimated miRBase-mapped read-length spectrum from depleted minus unmapped reads"},
        {"resource": "length_plot", "status": "ok", "path": args.length_plot, "rows": str(len(length)), "detail": "aggregate read-length SVG"},
    ]
    write_manifest(Path(args.manifest), manifest_rows)
    write_done(Path(args.done), len(length), len(arm))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
