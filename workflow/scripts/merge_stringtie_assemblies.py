#!/usr/bin/env python3
"""Merge per-sample StringTie assemblies into a project transcriptome."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_ASSEMBLY_COLUMNS = {"library_id", "assembly_gtf", "status"}
REQUIRED_PLAN_COLUMNS = {"status", "annotation_gtf", "transcriptome_mode"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-manifest", required=True, help="StringTie assembly manifest TSV")
    parser.add_argument("--plan", required=True, help="RNA-seq quantification plan TSV")
    parser.add_argument("--assemblies-list", required=True, help="Output list of assembly GTFs")
    parser.add_argument("--merged-gtf", required=True, help="Merged StringTie GTF")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--stringtie", default="stringtie", help="StringTie executable")
    parser.add_argument("--extra-args", default="", help="Extra StringTie --merge args")
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def read_plan(path: Path) -> dict[str, str]:
    rows = read_table(path, REQUIRED_PLAN_COLUMNS)
    if len(rows) != 1:
        raise ValueError(f"Quantification plan must contain exactly one row: {path}")
    row = rows[0]
    if row.get("status") != "ready":
        raise ValueError("Quantification plan is not ready: " + row.get("reason", ""))
    if row.get("transcriptome_mode") != "reference_guided_novel":
        raise ValueError("StringTie merge requires transcriptome_mode='reference_guided_novel'")
    if not Path(row["annotation_gtf"]).exists():
        raise FileNotFoundError(f"annotation_gtf does not exist: {row['annotation_gtf']}")
    return row


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def assembly_gtfs(rows: list[dict[str, str]]) -> list[Path]:
    gtfs = []
    errors = []
    for row in rows:
        if row.get("status") != "ok":
            errors.append(f"{row.get('library_id', '<unknown>')}: assembly status is not ok")
            continue
        gtf = Path(row.get("assembly_gtf", ""))
        if not gtf.exists():
            errors.append(f"{row.get('library_id', '<unknown>')}: assembly_gtf does not exist: {gtf}")
        else:
            gtfs.append(gtf)
    if errors:
        raise ValueError("StringTie merge cannot start:\n- " + "\n- ".join(errors))
    if not gtfs:
        raise ValueError("StringTie merge received no assembly GTFs")
    return gtfs


def main() -> int:
    args = parse_args()
    plan = read_plan(Path(args.plan))
    rows = read_table(Path(args.assembly_manifest), REQUIRED_ASSEMBLY_COLUMNS)
    gtfs = assembly_gtfs(rows)
    stringtie = executable_path(args.stringtie)

    assemblies_list = Path(args.assemblies_list)
    merged_gtf = Path(args.merged_gtf)
    done = Path(args.done)
    assemblies_list.parent.mkdir(parents=True, exist_ok=True)
    merged_gtf.parent.mkdir(parents=True, exist_ok=True)
    assemblies_list.write_text("\n".join(str(gtf) for gtf in gtfs) + "\n", encoding="utf-8")

    command = [
        stringtie,
        "--merge",
        "-G",
        plan["annotation_gtf"],
        "-o",
        str(merged_gtf),
    ]
    command.extend(shlex.split(args.extra_args))
    command.append(str(assemblies_list))

    print("[CMD] " + shlex.join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")
    if completed.returncode != 0:
        raise RuntimeError(f"StringTie --merge failed with status {completed.returncode}")
    if not merged_gtf.exists() or merged_gtf.stat().st_size == 0:
        raise RuntimeError(f"StringTie --merge produced an empty GTF: {merged_gtf}")

    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(f"status\tassemblies\tmerged_gtf\nok\t{len(gtfs)}\t{merged_gtf}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
