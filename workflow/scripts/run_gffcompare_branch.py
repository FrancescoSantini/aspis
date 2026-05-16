#!/usr/bin/env python3
"""Annotate a merged transcriptome with gffcompare and normalize output names."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


REQUIRED_PLAN_COLUMNS = {"status", "annotation_gtf", "transcriptome_mode"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--merged-gtf", required=True, help="Merged StringTie GTF")
    parser.add_argument("--plan", required=True, help="RNA-seq quantification plan TSV")
    parser.add_argument("--outdir", required=True, help="gffcompare output directory")
    parser.add_argument("--prefix", required=True, help="gffcompare output prefix basename")
    parser.add_argument("--annotated-gtf", required=True, help="Stable annotated GTF output")
    parser.add_argument("--tracking", required=True, help="Stable tracking output")
    parser.add_argument("--tmap", required=True, help="Stable tmap output")
    parser.add_argument("--refmap", required=True, help="Stable refmap output")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--gffcompare", default="gffcompare", help="gffcompare executable")
    return parser.parse_args()


def read_plan(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Quantification plan is empty: {path}")
        missing = REQUIRED_PLAN_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Quantification plan is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if len(rows) != 1:
        raise ValueError(f"Quantification plan must contain exactly one row: {path}")
    row = rows[0]
    if row.get("status") != "ready":
        raise ValueError("Quantification plan is not ready: " + row.get("reason", ""))
    if row.get("transcriptome_mode") != "reference_guided_novel":
        raise ValueError("gffcompare requires transcriptome_mode='reference_guided_novel'")
    if not Path(row["annotation_gtf"]).exists():
        raise FileNotFoundError(f"annotation_gtf does not exist: {row['annotation_gtf']}")
    return row


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def require_nonempty(path: Path, label: str) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"gffcompare did not produce {label}: {path}")


def log_tail(path: Path, max_lines: int = 20) -> str:
    if not path.exists() or path.stat().st_size == 0:
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-max_lines:])


def copy_first(paths: list[Path], output: Path, label: str, required: bool = True) -> None:
    for path in paths:
        if path.exists() and path.stat().st_size > 0:
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, output)
            return
    if required:
        raise RuntimeError(f"Could not find gffcompare {label}: {', '.join(str(path) for path in paths)}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("", encoding="utf-8")


def gffcompare_map_candidates(outdir: Path, merged_gtf: Path, prefix: str, suffix: str) -> list[Path]:
    """Return possible gffcompare per-query map files.

    gffcompare writes .tmap/.refmap per input query, and those files can appear
    next to the query GTF rather than next to the main output prefix.
    """
    candidates = [
        outdir / f"{prefix}.{merged_gtf.name}.{suffix}",
        merged_gtf.parent / f"{prefix}.{merged_gtf.name}.{suffix}",
        Path(f"{prefix}.{merged_gtf.name}.{suffix}"),
    ]
    candidates.extend(sorted(outdir.glob(f"{prefix}.*.{suffix}")))
    candidates.extend(sorted(merged_gtf.parent.glob(f"{prefix}.*.{suffix}")))

    unique = []
    seen = set()
    for path in candidates:
        key = str(path)
        if key not in seen:
            unique.append(path)
            seen.add(key)
    return unique


def main() -> int:
    args = parse_args()
    plan = read_plan(Path(args.plan))
    merged_gtf = Path(args.merged_gtf)
    if not merged_gtf.exists():
        raise FileNotFoundError(f"merged_gtf does not exist: {merged_gtf}")

    gffcompare = executable_path(args.gffcompare)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix_path = outdir / args.prefix
    stderr_log = outdir / f"{args.prefix}.stderr.log"
    command = [
        gffcompare,
        "-r",
        plan["annotation_gtf"],
        "-o",
        str(prefix_path),
        str(merged_gtf),
    ]

    print("[CMD] " + " ".join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    stderr_log.write_text(
        (completed.stdout or "") + (completed.stderr or ""),
        encoding="utf-8",
    )
    if completed.returncode != 0:
        tail = log_tail(stderr_log)
        detail = f"\nLast lines from {stderr_log}:\n{tail}" if tail else f"\nLog file: {stderr_log}"
        raise RuntimeError(f"gffcompare failed with status {completed.returncode}{detail}")

    raw_annotated = prefix_path.with_suffix(".annotated.gtf")
    raw_tracking = prefix_path.with_suffix(".tracking")
    require_nonempty(raw_annotated, "annotated GTF")
    require_nonempty(raw_tracking, "tracking file")

    annotated_gtf = Path(args.annotated_gtf)
    tracking = Path(args.tracking)
    tmap = Path(args.tmap)
    refmap = Path(args.refmap)
    annotated_gtf.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(raw_annotated, annotated_gtf)
    shutil.copyfile(raw_tracking, tracking)

    tmap_candidates = gffcompare_map_candidates(outdir, merged_gtf, args.prefix, "tmap")
    refmap_candidates = gffcompare_map_candidates(outdir, merged_gtf, args.prefix, "refmap")
    copy_first(tmap_candidates, tmap, "tmap", required=True)
    copy_first(refmap_candidates, refmap, "refmap", required=False)

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(
        "status\tannotated_gtf\ttracking\ttmap\n"
        f"ok\t{annotated_gtf}\t{tracking}\t{tmap}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
