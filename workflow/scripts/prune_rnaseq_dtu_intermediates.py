#!/usr/bin/env python3
"""Prune re-creatable RNA-seq DTU intermediate files after reporting."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PRUNE_COLUMNS = [
    "project",
    "method",
    "contrast_id",
    "file_kind",
    "path",
    "status",
    "reason",
    "bytes",
]

DONE_COLUMNS = [
    "status",
    "removed",
    "missing",
    "skipped",
    "failed",
    "bytes_removed",
    "total",
    "reason",
]

INPUT_SLICE_FILES = {
    "dtu_counts.tsv": "dtu_counts",
    "dtu_coldata.tsv": "dtu_coldata",
}

LOG_FILES = {
    "stdout.log": "stdout_log",
    "stderr.log": "stderr_log",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method-manifest", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument(
        "--delete-logs",
        action="store_true",
        help="Also delete per-method stdout/stderr logs. By default, logs are preserved.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write the prune manifest without removing files.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, rows: list[dict[str, str]], dry_run: bool) -> None:
    removed = [row for row in rows if row.get("status") == "removed"]
    missing = [row for row in rows if row.get("status") == "missing"]
    skipped = [row for row in rows if row.get("status") == "skipped"]
    failed = [row for row in rows if row.get("status") == "failed"]
    bytes_removed = sum(int(row.get("bytes") or 0) for row in removed)
    if failed:
        status = "failed"
        reason = f"{len(failed)} DTU intermediate file(s) could not be removed"
    elif dry_run:
        status = "dry_run"
        reason = f"{len(rows)} DTU intermediate prune candidate(s) inspected"
    else:
        status = "ok"
        reason = f"{len(removed)} DTU intermediate file(s) removed"
    write_tsv(
        path,
        [
            {
                "status": status,
                "removed": str(len(removed)),
                "missing": str(len(missing)),
                "skipped": str(len(skipped)),
                "failed": str(len(failed)),
                "bytes_removed": str(bytes_removed),
                "total": str(len(rows)),
                "reason": reason,
            }
        ],
        DONE_COLUMNS,
    )


def candidate_files(row: dict[str, str], delete_logs: bool) -> list[tuple[Path, str]]:
    outdir = Path(row.get("output_dir", ""))
    candidates = [(outdir / filename, kind) for filename, kind in INPUT_SLICE_FILES.items()]
    if delete_logs:
        candidates.extend((outdir / filename, kind) for filename, kind in LOG_FILES.items())
    return candidates


def prune_file(path: Path, dry_run: bool) -> tuple[str, str, int]:
    if not path.exists():
        return "missing", "file already absent", 0
    if not path.is_file():
        return "skipped", "not a regular file", 0
    size = path.stat().st_size
    if dry_run:
        return "skipped", "dry run", size
    try:
        path.unlink()
    except OSError as exc:
        return "failed", str(exc), size
    return "removed", "", size


def main() -> int:
    args = parse_args()
    method_manifest = Path(args.method_manifest)
    rows = read_tsv(method_manifest)
    prune_rows: list[dict[str, str]] = []

    for row in rows:
        if row.get("status") != "completed":
            continue
        project = row.get("project", "")
        method = row.get("method", "")
        contrast_id = row.get("contrast_id", "")
        if not row.get("output_dir", ""):
            prune_rows.append(
                {
                    "project": project,
                    "method": method,
                    "contrast_id": contrast_id,
                    "file_kind": "output_dir",
                    "path": "",
                    "status": "skipped",
                    "reason": "missing output_dir",
                    "bytes": "0",
                }
            )
            continue
        for path, file_kind in candidate_files(row, args.delete_logs):
            status, reason, size = prune_file(path, args.dry_run)
            prune_rows.append(
                {
                    "project": project,
                    "method": method,
                    "contrast_id": contrast_id,
                    "file_kind": file_kind,
                    "path": str(path),
                    "status": status,
                    "reason": reason,
                    "bytes": str(size),
                }
            )

    write_tsv(Path(args.manifest), prune_rows, PRUNE_COLUMNS)
    write_done(Path(args.done), prune_rows, args.dry_run)
    return 1 if any(row.get("status") == "failed" for row in prune_rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
