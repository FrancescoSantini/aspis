#!/usr/bin/env python3
"""Plan transcript-level DESeq2 contrasts from RNA-seq transcript counts."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from plan_gene_differential import (
    REQUIRED_SAMPLE_COLUMNS,
    contrast_id,
    grouped_rows,
    read_table,
    validate_samples,
    write_plan,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--transcript-counts", required=True, help="StringTie transcript count matrix")
    parser.add_argument("--differential-plan", required=True, help="Differential layer plan TSV")
    parser.add_argument("--output", required=True, help="Transcript contrast plan TSV")
    parser.add_argument("--outdir", required=True, help="Transcript DESeq2 output directory")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--condition-col", default="condition", help="Condition column")
    parser.add_argument("--control-label", default="control", help="Control condition label")
    parser.add_argument("--contrast-by", nargs="*", default=[], help="Optional stratifying columns")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum samples per group")
    return parser.parse_args()


def count_sample_columns(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Transcript count matrix is empty: {path}") from exc
    if "transcript_id" not in header:
        raise ValueError(f"Transcript count matrix lacks transcript_id column: {path}")
    return [column for column in header if column != "transcript_id"]


def validate_differential_plan(path: Path) -> list[str]:
    _, rows = read_table(path, {"level", "method", "status"})
    transcript_rows = [
        row for row in rows
        if row.get("level") == "transcript" and row.get("method") == "deseq2"
    ]
    if not transcript_rows:
        return [f"differential plan has no transcript/deseq2 row: {path}"]
    if transcript_rows[0].get("status") != "ready":
        return ["transcript/deseq2 differential layer is not ready: " + transcript_rows[0].get("reason", "")]
    return []


def build_plan_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    sample_columns, samples = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    if not samples:
        raise ValueError("Branch samples table has no rows")
    count_columns = count_sample_columns(Path(args.transcript_counts))
    contrast_by = [column for column in args.contrast_by if column]
    errors = validate_samples(
        sample_columns,
        samples,
        args.project,
        args.condition_col,
        contrast_by,
        count_columns,
    )
    errors.extend(validate_differential_plan(Path(args.differential_plan)))
    if args.min_replicates < 1:
        errors.append("--min-replicates must be >= 1")
    if errors:
        raise ValueError("Transcript differential plan cannot be built:\n- " + "\n- ".join(errors))

    rows = []
    outdir = Path(args.outdir)
    for values, by_condition in sorted(grouped_rows(samples, args.condition_col, contrast_by).items()):
        controls = by_condition.get(args.control_label, [])
        for condition in sorted(condition for condition in by_condition if condition != args.control_label):
            tested = by_condition[condition]
            cid = contrast_id(condition, args.control_label, contrast_by, values)
            selected = sorted(row["library_id"] for row in controls + tested)
            reasons = []
            if len(controls) < args.min_replicates:
                reasons.append(
                    f"control group has {len(controls)} sample(s); {args.min_replicates} required"
                )
            if len(tested) < args.min_replicates:
                reasons.append(
                    f"{condition!r} group has {len(tested)} sample(s); {args.min_replicates} required"
                )
            status = "blocked" if reasons else "ready"
            contrast_dir = outdir / "contrasts" / cid
            rows.append(
                {
                    "project": args.project,
                    "assay": "rnaseq",
                    "level": "transcript",
                    "method": "deseq2",
                    "contrast_id": cid,
                    "status": status,
                    "reason": "; ".join(reasons),
                    "condition_col": args.condition_col,
                    "control_label": args.control_label,
                    "test_label": condition,
                    "contrast_by": ",".join(contrast_by),
                    "contrast_values": ",".join(values),
                    "n_control": str(len(controls)),
                    "n_test": str(len(tested)),
                    "samples": ",".join(selected),
                    "contrast_dir": str(contrast_dir),
                    "counts": str(contrast_dir / "counts.tsv"),
                    "coldata": str(contrast_dir / "coldata.tsv"),
                    "results": str(contrast_dir / "deseq2_results.tsv"),
                    "filtered": str(contrast_dir / "deseq2_significant.tsv"),
                    "normalized_counts": str(contrast_dir / "normalized_counts.tsv"),
                    "summary": str(contrast_dir / "summary.tsv"),
                    "log": str(contrast_dir / "deseq2.log"),
                }
            )

    if not rows:
        rows.append(
            {
                "project": args.project,
                "assay": "rnaseq",
                "level": "transcript",
                "method": "deseq2",
                "contrast_id": "no_contrasts",
                "status": "blocked",
                "reason": f"no non-control conditions found for control label {args.control_label!r}",
                "condition_col": args.condition_col,
                "control_label": args.control_label,
                "test_label": "",
                "contrast_by": ",".join(contrast_by),
                "contrast_values": "",
                "n_control": "0",
                "n_test": "0",
                "samples": "",
                "contrast_dir": "",
                "counts": "",
                "coldata": "",
                "results": "",
                "filtered": "",
                "normalized_counts": "",
                "summary": "",
                "log": "",
            }
        )
    return rows


def main() -> int:
    args = parse_args()
    write_plan(Path(args.output), build_plan_rows(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
