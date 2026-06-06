#!/usr/bin/env python3
"""Merge per-item status manifests and write a compatible done sentinel."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

DONE_COLUMNS = {
    "deseq2": ("contrasts", ["ready"]),
    "isoform_switch": ("contrasts", []),
    "plots": ("plots", []),
    "enrichment": ("enrichment", []),
    "summaries": ("summaries", []),
}

REQUIRED_INPUT_COLUMNS = {
    "deseq2": {"contrast_id", "status", "reason", "n_significant"},
    "isoform_switch": {"contrast_id", "status", "reason", "n_transcripts", "n_genes"},
    "plots": {"project", "level", "contrast_id", "status", "n_features", "n_significant"},
    "enrichment": {"project", "level", "contrast_id", "status", "n_ranked", "n_feature_sets"},
    "summaries": {"project", "level", "contrast_id", "status", "n_features", "n_significant", "n_up", "n_down"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", required=True, choices=sorted(DONE_COLUMNS), help="Done sentinel schema")
    parser.add_argument("--manifest", required=True, help="Merged manifest TSV")
    parser.add_argument("--done", required=True, help="Merged done sentinel")
    parser.add_argument("inputs", nargs="+", help="Per-item manifest TSV files")
    parser.add_argument("--summary-top-n", type=int, default=20, help="Top features for summary fallback rendering")
    return parser.parse_args()


def read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Manifest is empty: {path}")
        return list(reader.fieldnames), [{key: (value or "") for key, value in row.items()} for row in reader]


def sort_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (row.get("level", ""), row.get("contrast_id", ""), row.get("resource", ""))


def validate_columns(kind: str, path: Path, columns: list[str]) -> None:
    missing = REQUIRED_INPUT_COLUMNS[kind] - set(columns)
    if missing:
        raise ValueError(
            f"{kind} merge input {path} is missing columns: {sorted(missing)}. "
            "Expected a per-item manifest, not a planning table."
        )


def write_single_plan_row(path: Path, done: Path, columns: list[str], row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({column: row.get(column, "") for column in columns})
    with done.open("w", encoding="utf-8") as handle:
        handle.write("status\tlevel\tcontrast_id\n")
        handle.write(f"ok\t{row.get('level', '')}\t{row.get('contrast_id', '')}\n")


def render_missing_summary_items(
    output_manifest: Path,
    columns: list[str],
    rows: list[dict[str, str]],
    summary_top_n: int,
) -> None:
    reports_dir = output_manifest.parent.parent
    summary_script = Path(__file__).with_name("render_rnaseq_differential_summary.py")
    for row in rows:
        level = row.get("level", "")
        contrast_id = row.get("contrast_id", "")
        if not level or not contrast_id:
            continue
        item_plan = reports_dir / "items" / level / contrast_id / "report_plan.tsv"
        item_plan_done = reports_dir / "items" / level / contrast_id / "report_plan.done"
        item_manifest = reports_dir / "summaries" / "items" / level / contrast_id / "summaries_manifest.tsv"
        item_done = reports_dir / "summaries" / "items" / level / contrast_id / "summaries.done"
        if item_manifest.exists() and item_done.exists():
            continue
        write_single_plan_row(item_plan, item_plan_done, columns, row)
        subprocess.run(
            [
                sys.executable,
                str(summary_script),
                "--plan",
                str(item_plan),
                "--manifest",
                str(item_manifest),
                "--done",
                str(item_done),
                "--top-n",
                str(summary_top_n),
            ],
            check=True,
        )


def expand_plan_placeholder(
    kind: str,
    output_manifest: Path,
    input_paths: list[Path],
    summary_top_n: int,
) -> list[Path]:
    if len(input_paths) != 1:
        return input_paths

    input_path = input_paths[0]
    columns, rows = read_manifest(input_path)
    if not rows:
        return input_paths
    if input_path.name not in {"report_plan.tsv", "contrast_plan.tsv"}:
        return input_paths
    if not (REQUIRED_INPUT_COLUMNS[kind] - set(columns)):
        return input_paths

    if kind in {"plots", "enrichment", "summaries"}:
        required = {"level", "contrast_id"}
        if required - set(columns):
            return input_paths
        section_dir = output_manifest.parent
        section = section_dir.name
        expanded = [
            section_dir / "items" / row.get("level", "") / row.get("contrast_id", "") / f"{section}_manifest.tsv"
            for row in rows
            if row.get("level", "") and row.get("contrast_id", "")
        ]
    elif kind in {"deseq2", "isoform_switch"}:
        if "contrast_id" not in columns:
            return input_paths
        expanded = [
            output_manifest.parent / "contrast_manifests" / f"{row.get('contrast_id', '')}.manifest.tsv"
            for row in rows
            if row.get("contrast_id", "")
        ]
    else:
        return input_paths

    if not expanded:
        return input_paths
    missing_paths = [path for path in expanded if not path.exists()]
    if missing_paths and kind == "summaries":
        render_missing_summary_items(output_manifest, columns, rows, summary_top_n)
        missing_paths = [path for path in expanded if not path.exists()]
    if missing_paths:
        preview = ", ".join(str(path) for path in missing_paths[:5])
        if len(missing_paths) > 5:
            preview += f", ... ({len(missing_paths)} missing total)"
        raise FileNotFoundError(
            f"{kind} merge received planning table {input_path}, but expected per-item manifest(s) "
            f"are missing: {preview}"
        )
    return expanded


def write_manifest(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, kind: str, rows: list[dict[str, str]]) -> None:
    prefix, extra_statuses = DONE_COLUMNS[kind]
    ok = sum(1 for row in rows if row.get("status") == "ok")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    extra_counts = {status: sum(1 for row in rows if row.get("status") == status) for status in extra_statuses}
    status = "failed" if failed else "ok" if ok and not blocked and not any(extra_counts.values()) else "blocked"

    columns = ["status", f"{prefix}_ok"]
    values = [status, str(ok)]
    for extra_status in extra_statuses:
        columns.append(f"{prefix}_{extra_status}")
        values.append(str(extra_counts[extra_status]))
    columns.extend([f"{prefix}_blocked", f"{prefix}_failed"])
    values.extend([str(blocked), str(failed)])
    if kind != "deseq2":
        columns.append(f"{prefix}_total")
        values.append(str(len(rows)))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.write("\t".join(values) + "\n")

    if failed:
        failed_ids = ", ".join(row.get("contrast_id", "") for row in rows if row.get("status") == "failed")
        raise RuntimeError(f"Merged {kind} manifest contains failed item(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if not args.inputs:
        raise ValueError("At least one input manifest is required")

    input_paths = expand_plan_placeholder(args.kind, Path(args.manifest), [Path(path) for path in args.inputs], args.summary_top_n)

    columns: list[str] | None = None
    rows: list[dict[str, str]] = []
    for input_path in input_paths:
        current_columns, current_rows = read_manifest(input_path)
        validate_columns(args.kind, input_path, current_columns)
        if columns is None:
            columns = current_columns
        elif current_columns != columns:
            raise ValueError(f"Manifest columns do not match first input: {input_path}")
        rows.extend(current_rows)

    if columns is None or not rows:
        raise ValueError("Input manifests did not contain rows")
    rows.sort(key=sort_key)
    write_manifest(Path(args.manifest), columns, rows)
    write_done(Path(args.done), args.kind, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
