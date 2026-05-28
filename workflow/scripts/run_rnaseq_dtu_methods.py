#!/usr/bin/env python3
"""Run or manifest optional RNA-seq event-level DTU/splicing engines."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


METHOD_ALIASES = {
    "drimseq": "DRIMSeq",
    "dexseq": "DEXSeq",
    "suppa2": "SUPPA2",
    "suppa": "SUPPA2",
    "rmats": "rMATS",
    "rmats-turbo": "rMATS",
    "rmats_turbo": "rMATS",
}

MANIFEST_COLUMNS = [
    "project",
    "assay",
    "level",
    "method",
    "status",
    "reason",
    "command",
    "output_dir",
    "stdout",
    "stderr",
    "plan",
    "samples",
    "aligned_samples",
    "transcript_counts",
    "transcript_metadata",
    "annotation_gtf",
]


class SafeFormat(dict):
    def __missing__(self, key: str) -> str:
        return "{" + key + "}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--aligned-samples", required=True)
    parser.add_argument("--transcript-counts", required=True)
    parser.add_argument("--transcript-metadata", required=True)
    parser.add_argument("--annotation-gtf", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--method", default="planned")
    parser.add_argument("--methods", default="DRIMSeq,DEXSeq,SUPPA2,rMATS")
    parser.add_argument("--drimseq-command", default="")
    parser.add_argument("--dexseq-command", default="")
    parser.add_argument("--suppa2-command", default="")
    parser.add_argument("--rmats-command", default="")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in MANIFEST_COLUMNS})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    failures = [row for row in rows if row["status"] == "failed"]
    completed = [row for row in rows if row["status"] == "completed"]
    planned = [row for row in rows if row["status"] == "planned"]
    if failures:
        status = "failed"
        reason = ",".join(row["method"] for row in failures)
    elif completed:
        status = "ok"
        reason = f"{len(completed)} method(s) completed"
    else:
        status = "planned"
        reason = f"{len(planned)} method(s) have no configured command template"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tcompleted_methods\tplanned_methods\treason\n")
        handle.write(f"{status}\t{len(completed)}\t{len(planned)}\t{reason}\n")


def normalize_methods(method: str, methods: str) -> list[str]:
    selected = (method or "").strip()
    if selected and selected.lower() not in {"planned", "all", "auto", "none"}:
        method_tokens = [selected]
    else:
        method_tokens = [item.strip() for item in methods.split(",") if item.strip()]
    normalized: list[str] = []
    for token in method_tokens:
        canonical = METHOD_ALIASES.get(token.strip().lower(), token.strip())
        if canonical and canonical not in normalized:
            normalized.append(canonical)
    return normalized


def command_for_method(args: argparse.Namespace, method: str) -> str:
    commands = {
        "DRIMSeq": args.drimseq_command,
        "DEXSeq": args.dexseq_command,
        "SUPPA2": args.suppa2_command,
        "rMATS": args.rmats_command,
    }
    return commands.get(method, "")


def context_for(args: argparse.Namespace, method: str, method_dir: Path) -> SafeFormat:
    context = SafeFormat()
    context.update(
        {
            "project": args.project,
            "method": method,
            "outdir": str(method_dir.resolve()),
            "plan": str(Path(args.plan).resolve()),
            "samples": str(Path(args.samples).resolve()),
            "aligned_samples": str(Path(args.aligned_samples).resolve()),
            "transcript_counts": str(Path(args.transcript_counts).resolve()),
            "transcript_metadata": str(Path(args.transcript_metadata).resolve()),
            "annotation_gtf": str(Path(args.annotation_gtf).resolve()),
        }
    )
    return context


def run_method(args: argparse.Namespace, method: str, plan: Path) -> dict[str, str]:
    method_dir = Path(args.outdir) / method.lower().replace("-", "_")
    method_dir.mkdir(parents=True, exist_ok=True)
    stdout = method_dir / "stdout.log"
    stderr = method_dir / "stderr.log"
    command_template = command_for_method(args, method).strip()
    row = {
        "project": args.project,
        "assay": "rnaseq",
        "level": "differential_transcript_usage",
        "method": method,
        "command": command_template,
        "output_dir": str(method_dir),
        "stdout": "",
        "stderr": "",
        "plan": str(plan),
        "samples": args.samples,
        "aligned_samples": args.aligned_samples,
        "transcript_counts": args.transcript_counts,
        "transcript_metadata": args.transcript_metadata,
        "annotation_gtf": args.annotation_gtf,
    }
    if not command_template:
        row.update(
            {
                "status": "planned",
                "reason": "no command template configured for this optional method",
            }
        )
        return row
    command = command_template.format_map(context_for(args, method, method_dir))
    completed = subprocess.run(command, shell=True, text=True, capture_output=True, check=False)
    stdout.write_text(completed.stdout, encoding="utf-8")
    stderr.write_text(completed.stderr, encoding="utf-8")
    row["command"] = command
    row["stdout"] = str(stdout)
    row["stderr"] = str(stderr)
    if completed.returncode == 0:
        row["status"] = "completed"
        row["reason"] = ""
    else:
        row["status"] = "failed"
        row["reason"] = f"command exited with status {completed.returncode}"
    return row


def main() -> int:
    args = parse_args()
    required = [
        args.plan,
        args.samples,
        args.aligned_samples,
        args.transcript_counts,
        args.transcript_metadata,
        args.annotation_gtf,
    ]
    for path_text in required:
        if not Path(path_text).exists():
            raise FileNotFoundError(path_text)
    read_tsv(Path(args.plan))
    methods = normalize_methods(args.method, args.methods)
    if not methods:
        raise ValueError("no DTU methods selected")
    rows = [run_method(args, method, Path(args.plan)) for method in methods]
    write_tsv(Path(args.manifest), rows)
    write_done(Path(args.done), rows)
    failed = [row["method"] for row in rows if row["status"] == "failed"]
    if failed:
        raise RuntimeError(f"DTU method command(s) failed: {','.join(failed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
