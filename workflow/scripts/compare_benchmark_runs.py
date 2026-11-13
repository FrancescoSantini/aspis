#!/usr/bin/env python3
"""Create publication-oriented tables and SVG figures from benchmark bundles."""

from __future__ import annotations

import argparse
import csv
import html
import re
from pathlib import Path


RUN_FIELDS = [
    "label", "benchmark_dir", "benchmark_mode", "submitted_jobs", "completed_jobs",
    "failed_jobs", "workflow_wall_seconds", "total_cpu_seconds", "cpu_hours",
    "allocated_cpu_seconds", "cpu_efficiency", "maximum_job_rss_mb",
    "disk_read_sum_mb", "disk_write_sum_mb", "queue_median_seconds",
    "queue_max_seconds", "disk_io_observed",
]
LABEL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def read_key_values(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["field"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def svg_bar_chart(
    path: Path,
    title: str,
    labels: list[str],
    values: list[float],
    unit: str,
    color: str = "#2f6f9f",
) -> None:
    width = 1000
    row_height = 34
    top = 75
    left = 330
    right = 150
    height = top + max(len(labels), 1) * row_height + 55
    plot_width = width - left - right
    maximum = max(values, default=0.0) or 1.0
    pieces = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:24px;font-weight:bold}.label{font-size:13px}.value{font-size:12px}</style>',
        f'<text x="30" y="38" class="title">{html.escape(title)}</text>',
    ]
    for index, (label, value) in enumerate(zip(labels, values)):
        y = top + index * row_height
        bar_width = plot_width * value / maximum
        pieces.append(f'<text x="{left - 12}" y="{y + 17}" text-anchor="end" class="label">{html.escape(label)}</text>')
        pieces.append(f'<rect x="{left}" y="{y}" width="{bar_width:.2f}" height="22" fill="{color}"/>')
        pieces.append(f'<text x="{left + bar_width + 8:.2f}" y="{y + 16}" class="value">{value:.3g} {html.escape(unit)}</text>')
    pieces.append("</svg>\n")
    path.write_text("\n".join(pieces), encoding="utf-8")


def parse_run(specification: str) -> tuple[str, Path]:
    if "=" not in specification:
        raise argparse.ArgumentTypeError("--run must use LABEL=BENCHMARK_DIRECTORY")
    label, directory = specification.split("=", 1)
    if not LABEL.fullmatch(label) or not directory:
        raise argparse.ArgumentTypeError("--run must use LABEL=BENCHMARK_DIRECTORY")
    return label, Path(directory)


def rule_category(rule: str) -> str:
    groups = [
        ("input", ("materialize", "inspect_fastq", "build_analysis_plan", "build_branch_design")),
        ("QC and preprocessing", ("fastqc", "multiqc", "preprocess", "deplete", "qc_", "length_qc")),
        ("alignment", ("align_", "alignment", "star_index", "hisat2_index", "residual_genome_index")),
        ("quantification and assembly", ("featurecounts", "stringtie", "gffcompare", "quantification", "transcript_matrix")),
        ("differential analysis", ("deseq2", "differential")),
        ("DTU and isoform analysis", ("dtu", "dexseq", "drimseq", "suppa", "rmats", "isoform_switch")),
        ("integration and enrichment", ("mirna_mrna", "target", "enrichment", "feature_set")),
        ("reporting and provenance", ("render_", "report", "dashboard", "pdf", "provenance")),
    ]
    lowered = rule.lower()
    return next((category for category, tokens in groups if any(token in lowered for token in tokens)), "other")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True, type=parse_run)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top-rules", type=int, default=20)
    args = parser.parse_args()
    if args.top_rules < 1:
        parser.error("--top-rules must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    comparison: list[dict[str, object]] = []
    combined_rules: list[dict[str, object]] = []
    for label, directory in args.run:
        required = [directory / "slurm_run_summary.tsv", directory / "slurm_rule_summary.tsv"]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            parser.error("missing benchmark summaries: " + ", ".join(missing))
        run = read_key_values(required[0])
        metadata = read_key_values(directory / "metadata.tsv") if (directory / "metadata.tsv").is_file() else {}
        cpu_seconds = number(run.get("total_cpu_seconds", ""))
        row: dict[str, object] = {field: "" for field in RUN_FIELDS}
        row.update({
            "label": label,
            "benchmark_dir": str(directory),
            "benchmark_mode": metadata.get("benchmark_mode", ""),
            "submitted_jobs": run.get("submitted_jobs", ""),
            "completed_jobs": run.get("completed_jobs", ""),
            "failed_jobs": run.get("failed_jobs", ""),
            "workflow_wall_seconds": run.get("workflow_wall_seconds", ""),
            "total_cpu_seconds": run.get("total_cpu_seconds", ""),
            "cpu_hours": "" if cpu_seconds is None else f"{cpu_seconds / 3600:.3f}",
            "allocated_cpu_seconds": run.get("allocated_cpu_seconds", ""),
            "cpu_efficiency": run.get("cpu_efficiency", ""),
            "maximum_job_rss_mb": run.get("maximum_job_rss_mb", ""),
            "disk_read_sum_mb": run.get("disk_read_sum_mb", ""),
            "disk_write_sum_mb": run.get("disk_write_sum_mb", ""),
            "queue_median_seconds": run.get("queue_median_seconds", ""),
            "queue_max_seconds": run.get("queue_max_seconds", ""),
            "disk_io_observed": run.get("disk_io_observed", ""),
        })
        comparison.append(row)
        for rule in read_rows(required[1]):
            combined_rules.append({"run_label": label, "category": rule_category(rule["rule"]), **rule})

    write_tsv(args.output_dir / "benchmark_comparison.tsv", RUN_FIELDS, comparison)
    rule_fields = list(combined_rules[0]) if combined_rules else ["run_label", "rule"]
    write_tsv(args.output_dir / "benchmark_rule_comparison.tsv", rule_fields, combined_rules)

    run_labels = [str(row["label"]) for row in comparison]
    metrics = [
        ("workflow_wall_seconds", "Workflow wall time", "seconds", "run_wall_time.svg"),
        ("cpu_hours", "Consumed CPU time", "CPU-hours", "run_cpu_hours.svg"),
        ("maximum_job_rss_mb", "Maximum job resident memory", "MiB", "run_peak_memory.svg"),
    ]
    for field, title, unit, filename in metrics:
        values = [number(str(row[field])) for row in comparison]
        if all(value is not None for value in values):
            svg_bar_chart(args.output_dir / filename, title, run_labels, [float(value) for value in values], unit)
    efficiencies = [number(str(row["cpu_efficiency"])) for row in comparison]
    if all(value is not None for value in efficiencies):
        svg_bar_chart(
            args.output_dir / "run_cpu_efficiency.svg",
            "CPU efficiency",
            run_labels,
            [100 * float(value) for value in efficiencies],
            "%",
            "#3d8b5f",
        )
    for field, title, filename in (
        ("disk_read_sum_mb", "Aggregate disk reads", "run_disk_read.svg"),
        ("disk_write_sum_mb", "Aggregate disk writes", "run_disk_write.svg"),
    ):
        values = [number(str(row[field])) for row in comparison]
        if all(value is not None for value in values) and any(float(value) > 0 for value in values):
            svg_bar_chart(args.output_dir / filename, title, run_labels, [float(value) for value in values], "MiB")

    for label, _ in args.run:
        rows = [row for row in combined_rules if row["run_label"] == label]
        rule_metrics = [
            ("elapsed_sum_seconds", "Rule aggregate wall time", "seconds", "rule_runtime"),
            ("cpu_sum_seconds", "Rule CPU consumption", "CPU-seconds", "rule_cpu"),
            ("cpu_efficiency", "Rule CPU efficiency", "fraction", "rule_cpu_efficiency"),
            ("rss_max_mb", "Rule peak resident memory", "MiB", "rule_memory"),
            ("disk_read_sum_mb", "Rule disk reads", "MiB", "rule_disk_read"),
            ("disk_write_sum_mb", "Rule disk writes", "MiB", "rule_disk_write"),
        ]
        for field, title, unit, stem in rule_metrics:
            ranked = [(str(row["rule"]), number(str(row.get(field, "")))) for row in rows]
            ranked = [(name, value) for name, value in ranked if value is not None and value > 0]
            ranked.sort(key=lambda item: item[1], reverse=True)
            ranked = ranked[: args.top_rules]
            if ranked:
                svg_bar_chart(
                    args.output_dir / f"{label}_{stem}.svg",
                    f"{title}: {label}",
                    [name for name, _ in ranked],
                    [float(value) for _, value in ranked],
                    unit,
                )

        categories: dict[str, float] = {}
        for row in rows:
            elapsed = number(str(row.get("elapsed_sum_seconds", "")))
            if elapsed is not None:
                category = str(row["category"])
                categories[category] = categories.get(category, 0.0) + elapsed
        ranked_categories = sorted(categories.items(), key=lambda item: item[1], reverse=True)
        if ranked_categories:
            svg_bar_chart(
                args.output_dir / f"{label}_workflow_composition.svg",
                f"Workflow rule-category composition: {label}",
                [category for category, _ in ranked_categories],
                [value for _, value in ranked_categories],
                "aggregate job-seconds",
                "#7a5aa6",
            )

    (args.output_dir / "comparison.done").touch()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
