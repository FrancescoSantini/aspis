#!/usr/bin/env python3
"""Materialize one ASPIS intake row into canonical FASTQ files."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


INSDC_RUN_RE = re.compile(r"^[SED]RR\d+$")
ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_READ_RUN_FIELDS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "secondary_study_accession",
    "tax_id",
    "scientific_name",
    "instrument_platform",
    "instrument_model",
    "library_name",
    "library_layout",
    "library_strategy",
    "library_source",
    "library_selection",
    "nominal_length",
    "read_count",
    "base_count",
    "experiment_title",
    "study_title",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]
ASSAY_ALIASES = {
    "rnaseq": "rnaseq",
    "rna-seq": "rnaseq",
    "rna_seq": "rnaseq",
    "rna seq": "rnaseq",
    "mrna": "rnaseq",
    "mrnaseq": "rnaseq",
    "mrna-seq": "rnaseq",
    "mrna_seq": "rnaseq",
    "mrna seq": "rnaseq",
    "longrna": "rnaseq",
    "longrna-seq": "rnaseq",
    "longrna seq": "rnaseq",
    "smallrna": "smallrna",
    "smallrna-seq": "smallrna",
    "small-rna": "smallrna",
    "small-rna-seq": "smallrna",
    "small rna": "smallrna",
    "small rna seq": "smallrna",
    "mirna": "smallrna",
    "mirnaseq": "smallrna",
    "mirna-seq": "smallrna",
    "mirna seq": "smallrna",
}
NORMALIZED_ASSAY_ALIASES = {
    re.sub(r"[^a-z0-9]+", "", key.lower()): value
    for key, value in ASSAY_ALIASES.items()
}
ASSAY_DECLARATION_FIELDS = ("assay_hint", "assay")
ASSAY_METADATA_FIELDS = ("library_strategy", "library_selection")
LIBRARY_STRATEGY_ASSAYS = {
    "rnaseq": "rnaseq",
    "mrnaseq": "rnaseq",
    "mirnaseq": "smallrna",
    "smallrnaseq": "smallrna",
}
LIBRARY_SELECTION_ASSAYS = {
    "mirna": "smallrna",
    "mirnasizefractionation": "smallrna",
}
RESERVED_METADATA_KEYS = {
    "library_id",
    "input_1",
    "input_2",
    "assay_hint",
    "assay",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intake", required=True, help="TSV intake sheet")
    parser.add_argument("--library-id", required=True, help="Library ID to materialize")
    parser.add_argument("--outdir", required=True, help="Canonical raw output directory")
    parser.add_argument("--metadata", required=True, help="Per-library metadata JSON")
    parser.add_argument("--sra-cache-dir", default="cache/sra", help="SRA cache directory")
    parser.add_argument("--scratch-dir", default="work/tmp", help="Temporary scratch directory")
    parser.add_argument(
        "--local-link-mode",
        choices=("symlink", "copy"),
        default="symlink",
        help="How to materialize local gzipped FASTQ inputs",
    )
    parser.add_argument("--sra-max-size", default="40G", help="prefetch -X value")
    parser.add_argument(
        "--sra-spot-limit",
        type=int,
        default=0,
        help="If >0, run a partial SRA smoke extraction with fastq-dump -X",
    )
    parser.add_argument(
        "--no-validate-sra",
        action="store_true",
        help="Skip vdb-validate even if available",
    )
    parser.add_argument(
        "--public-metadata-mode",
        choices=("auto", "off", "required"),
        default="auto",
        help="Resolve INSDC run metadata through ENA before materialization",
    )
    parser.add_argument(
        "--ena-api-url",
        default=ENA_API_URL,
        help="ENA filereport API endpoint",
    )
    parser.add_argument(
        "--public-metadata-timeout",
        type=int,
        default=60,
        help="Seconds to wait for public metadata resolution",
    )
    return parser.parse_args()


def read_intake_row(path: Path, library_id: str) -> dict[str, str]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Intake sheet is empty: {path}")
        required = {"library_id", "input_1"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Intake sheet is missing required columns: {sorted(missing)}")
        for row in reader:
            rows.append({key: (value or "").strip() for key, value in row.items()})

    matches = [row for row in rows if row.get("library_id") == library_id]
    if not matches:
        raise ValueError(f"Library {library_id!r} not found in {path}")
    if len(matches) > 1:
        raise ValueError(f"Library {library_id!r} appears more than once in {path}")
    return matches[0]


def resolve_path(raw_path: str, intake_path: Path) -> Path:
    path = Path(raw_path).expanduser()
    if path.is_absolute() and path.exists():
        return path

    cwd_path = Path.cwd() / path
    if cwd_path.exists():
        return cwd_path

    intake_relative = intake_path.parent / path
    if intake_relative.exists():
        return intake_relative

    return cwd_path


def clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def is_gzipped_fastq(path: Path) -> bool:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    return suffixes[-2:] in ([".fastq", ".gz"], [".fq", ".gz"])


def is_plain_fastq(path: Path) -> bool:
    return path.suffix.lower() in {".fastq", ".fq"}


def validate_fastq_path(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Input FASTQ does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Input FASTQ is not a file: {path}")
    if not (is_gzipped_fastq(path) or is_plain_fastq(path)):
        raise ValueError(
            "Unsupported local input extension for "
            f"{path}. Expected .fastq, .fq, .fastq.gz, or .fq.gz."
        )


def link_or_copy(src: Path, dest: Path, mode: str) -> str:
    ensure_parent(dest)
    if dest.exists() or dest.is_symlink():
        dest.unlink()

    if mode == "copy":
        shutil.copy2(src, dest)
        return "copy"

    try:
        os.symlink(src.resolve(), dest)
        return "symlink"
    except OSError:
        shutil.copy2(src, dest)
        return "copy"


def gzip_copy(src: Path, dest: Path) -> None:
    ensure_parent(dest)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    with src.open("rb") as in_handle, gzip.open(dest, "wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)


def materialize_local_fastq(
    row: dict[str, str],
    intake_path: Path,
    outdir: Path,
    link_mode: str,
) -> tuple[str, str | None, str, list[str]]:
    input_1 = resolve_path(row["input_1"], intake_path)
    input_2_raw = row.get("input_2", "")
    input_2 = resolve_path(input_2_raw, intake_path) if input_2_raw else None

    validate_fastq_path(input_1)
    if input_2 is not None:
        validate_fastq_path(input_2)

    operations: list[str] = []
    r1 = outdir / "R1.fastq.gz"
    if is_gzipped_fastq(input_1):
        op = link_or_copy(input_1, r1, link_mode)
        operations.append(f"{op}:{input_1}->{r1}")
    else:
        gzip_copy(input_1, r1)
        operations.append(f"gzip:{input_1}->{r1}")

    r2_value = None
    layout = "single"
    if input_2 is not None:
        r2 = outdir / "R2.fastq.gz"
        if is_gzipped_fastq(input_2):
            op = link_or_copy(input_2, r2, link_mode)
            operations.append(f"{op}:{input_2}->{r2}")
        else:
            gzip_copy(input_2, r2)
            operations.append(f"gzip:{input_2}->{r2}")
        r2_value = str(r2)
        layout = "paired"

    return str(r1), r2_value, layout, operations


def require_command(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Required command not found on PATH: {name}")
    return executable


def run_command(command: list[str]) -> None:
    print("[CMD]", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def read_url_text(url: str, timeout: int) -> str:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "aspis-materializer/0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            charset = response.headers.get_content_charset() or "utf-8"
            return response.read().decode(charset)
    except urllib.error.URLError as exc:
        raise RuntimeError(f"could not query public metadata endpoint: {exc}") from exc


def fetch_ena_read_run_metadata(
    accession: str,
    api_url: str,
    timeout: int,
) -> dict[str, str]:
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(ENA_READ_RUN_FIELDS),
            "format": "tsv",
            "download": "false",
        }
    )
    text = read_url_text(f"{api_url}?{query}", timeout=timeout)
    if not text.strip():
        raise RuntimeError(f"ENA returned no metadata for {accession}")

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise RuntimeError(f"ENA returned an empty metadata table for {accession}")

    rows = [
        {key: (value or "").strip() for key, value in row.items()}
        for row in reader
    ]
    if not rows:
        raise RuntimeError(f"ENA returned no read_run rows for {accession}")

    exact = [row for row in rows if row.get("run_accession") == accession]
    if len(exact) == 1:
        selected = exact[0]
    elif len(rows) == 1:
        selected = rows[0]
    else:
        observed = sorted(row.get("run_accession", "") for row in rows)
        raise RuntimeError(
            f"ENA returned multiple read_run rows for {accession}: {observed}"
        )

    return {key: value for key, value in selected.items() if value}


def merge_missing_values(
    row: dict[str, str],
    values: dict[str, str],
) -> tuple[dict[str, str], list[str]]:
    merged = dict(row)
    added = []
    for key, value in values.items():
        if value and not merged.get(key, ""):
            merged[key] = value
            added.append(key)
    return merged, sorted(added)


def resolve_public_metadata(
    row: dict[str, str],
    mode: str,
    api_url: str,
    timeout: int,
) -> tuple[dict[str, str], dict[str, str], list[str]]:
    input_1 = row.get("input_1", "")
    status = {
        "public_metadata_source": "",
        "public_metadata_status": "not_applicable",
        "public_metadata_error": "",
    }
    operations: list[str] = []

    if mode == "off" or not INSDC_RUN_RE.match(input_1):
        return row, status, operations

    try:
        metadata = fetch_ena_read_run_metadata(input_1, api_url=api_url, timeout=timeout)
    except Exception as exc:  # noqa: BLE001 - preserve user hint fallback if present.
        if mode == "required":
            raise RuntimeError(f"Public metadata resolution failed for {input_1}: {exc}") from exc
        status.update(
            {
                "public_metadata_source": "ena",
                "public_metadata_status": "failed",
                "public_metadata_error": str(exc),
            }
        )
        operations.append(f"resolve-metadata:ena:{input_1}:failed")
        return row, status, operations

    merged, added_keys = merge_missing_values(row, metadata)
    status.update(
        {
            "public_metadata_source": "ena",
            "public_metadata_status": "resolved",
            "public_metadata_error": "",
        }
    )
    operations.append(
        f"resolve-metadata:ena:{input_1}:added={','.join(added_keys) or 'none'}"
    )
    return merged, status, operations


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def gzip_sra_fastq(src: Path, dest: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"Expected SRA conversion output does not exist: {src}")
    gzip_copy(src, dest)


def canonicalize_dumped_fastqs(
    accession: str,
    outdir: Path,
    dump_dir: Path,
    operations: list[str],
) -> tuple[str, str | None, str, list[str]]:
    r1_source = first_existing(
        [
            dump_dir / f"{accession}_1.fastq",
            dump_dir / f"{accession}.fastq",
        ]
    )
    r2_source = first_existing([dump_dir / f"{accession}_2.fastq"])

    if r1_source is None:
        observed = sorted(path.name for path in dump_dir.glob("*.fastq"))
        raise RuntimeError(
            f"Could not find FASTQ output for {accession}. Observed: {observed}"
        )

    r1 = outdir / "R1.fastq.gz"
    gzip_sra_fastq(r1_source, r1)
    operations.append(f"gzip:{r1_source}->{r1}")

    r2_value = None
    layout = "single"
    if r2_source is not None:
        r2 = outdir / "R2.fastq.gz"
        gzip_sra_fastq(r2_source, r2)
        operations.append(f"gzip:{r2_source}->{r2}")
        r2_value = str(r2)
        layout = "paired"

    return str(r1), r2_value, layout, operations


def materialize_limited_insdc_run(
    accession: str,
    outdir: Path,
    scratch_dir: Path,
    spot_limit: int,
) -> tuple[str, str | None, str, list[str]]:
    require_command("fastq-dump")

    scratch_dir.mkdir(parents=True, exist_ok=True)
    dump_dir = scratch_dir / "fastq-dump" / accession
    clean_dir(dump_dir)

    run_command(
        [
            "fastq-dump",
            "--split-files",
            "-X",
            str(spot_limit),
            "-O",
            str(dump_dir),
            accession,
        ]
    )
    operations = [f"fastq-dump:{accession}:max_spots={spot_limit}"]
    result = canonicalize_dumped_fastqs(accession, outdir, dump_dir, operations)
    shutil.rmtree(dump_dir, ignore_errors=True)
    return result


def materialize_insdc_run(
    accession: str,
    outdir: Path,
    cache_dir: Path,
    scratch_dir: Path,
    max_size: str,
    validate_sra: bool,
    spot_limit: int,
) -> tuple[str, str | None, str, list[str]]:
    if spot_limit < 0:
        raise ValueError("--sra-spot-limit cannot be negative")
    if spot_limit > 0:
        return materialize_limited_insdc_run(
            accession=accession,
            outdir=outdir,
            scratch_dir=scratch_dir,
            spot_limit=spot_limit,
        )

    require_command("prefetch")
    require_command("fasterq-dump")

    cache_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    dump_dir = scratch_dir / "fasterq" / accession
    clean_dir(dump_dir)

    run_command(["prefetch", "-X", max_size, "-O", str(cache_dir), accession])
    accession_dir = cache_dir / accession

    if validate_sra and shutil.which("vdb-validate") is not None:
        run_command(["vdb-validate", str(accession_dir)])

    run_command(
        [
            "fasterq-dump",
            "--split-files",
            "-O",
            str(dump_dir),
            "-t",
            str(scratch_dir),
            str(accession_dir),
        ]
    )

    operations = [f"prefetch:{accession}", f"fasterq-dump:{accession}"]
    result = canonicalize_dumped_fastqs(accession, outdir, dump_dir, operations)
    shutil.rmtree(dump_dir, ignore_errors=True)
    return result


def normalize_metadata_value(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.strip().lower())


def canonical_assay_value(value: str) -> str | None:
    normalized = normalize_metadata_value(value)
    if not normalized:
        return None
    return NORMALIZED_ASSAY_ALIASES.get(normalized)


def assay_from_metadata_field(field: str, value: str) -> str | None:
    normalized = normalize_metadata_value(value)
    if field == "library_strategy":
        return LIBRARY_STRATEGY_ASSAYS.get(normalized)
    if field == "library_selection":
        return LIBRARY_SELECTION_ASSAYS.get(normalized)
    return None


def assay_from_row(row: dict[str, str]) -> tuple[str, str, str]:
    declared = []
    for field in ASSAY_DECLARATION_FIELDS:
        value = row.get(field, "").strip()
        if not value:
            continue
        assay = canonical_assay_value(value)
        if assay is None:
            allowed = ", ".join(sorted(set(ASSAY_ALIASES.values())))
            raise ValueError(
                f"Unsupported {field} {value!r}. Use one of: {allowed}; or leave blank."
            )
        declared.append((field, value, assay))

    if declared:
        assays = {assay for _, _, assay in declared}
        if len(assays) > 1:
            details = ", ".join(f"{field}={value!r}" for field, value, _ in declared)
            raise ValueError(f"Conflicting assay declarations: {details}")
        field, value, assay = declared[0]
        return assay, "user_hint", f"{field}={value}"

    metadata_hits = []
    for field in ASSAY_METADATA_FIELDS:
        value = row.get(field, "").strip()
        if not value:
            continue
        assay = assay_from_metadata_field(field, value)
        if assay is not None:
            metadata_hits.append((field, value, assay))

    if metadata_hits:
        assays = {assay for _, _, assay in metadata_hits}
        details = ", ".join(f"{field}={value!r}" for field, value, _ in metadata_hits)
        if len(assays) == 1:
            return metadata_hits[0][2], "metadata", details
        return "unknown", "ambiguous_metadata", f"conflicting metadata: {details}"

    return (
        "unknown",
        "unclassified",
        "no assay_hint, assay, or recognized library_strategy/library_selection metadata",
    )


def write_metadata(path: Path, payload: dict[str, object]) -> None:
    ensure_parent(path)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    args = parse_args()

    intake_path = Path(args.intake)
    row = read_intake_row(intake_path, args.library_id)
    if not row.get("input_1"):
        raise ValueError(f"Library {args.library_id!r} has empty input_1")

    outdir = Path(args.outdir)
    metadata_path = Path(args.metadata)

    row, public_metadata_status, metadata_operations = resolve_public_metadata(
        row,
        mode=args.public_metadata_mode,
        api_url=args.ena_api_url,
        timeout=args.public_metadata_timeout,
    )
    assay, assay_confidence, assay_reason = assay_from_row(row)
    input_1 = row["input_1"]
    input_2 = row.get("input_2", "")

    if INSDC_RUN_RE.match(input_1) and assay == "unknown":
        detail = assay_reason
        if public_metadata_status["public_metadata_error"]:
            detail = public_metadata_status["public_metadata_error"]
        raise RuntimeError(
            f"Public run assay could not be classified before downloading {input_1}. "
            "Provide assay_hint/assay, or provide/resolve recognized library metadata. "
            f"Reason: {detail}"
        )

    clean_dir(outdir)

    if resolve_path(input_1, intake_path).exists():
        source_type = "local_fastq"
        source_id = input_1
        archive = ""
        fastq_1, fastq_2, layout, operations = materialize_local_fastq(
            row,
            intake_path,
            outdir,
            args.local_link_mode,
        )
    elif INSDC_RUN_RE.match(input_1):
        if input_2:
            raise ValueError("input_2 is not supported for public run accessions")
        source_type = "insdc_run"
        source_id = input_1
        archive = "INSDC"
        fastq_1, fastq_2, layout, operations = materialize_insdc_run(
            accession=input_1,
            outdir=outdir,
            cache_dir=Path(args.sra_cache_dir),
            scratch_dir=Path(args.scratch_dir),
            max_size=args.sra_max_size,
            validate_sra=not args.no_validate_sra,
            spot_limit=args.sra_spot_limit,
        )
    else:
        raise ValueError(
            f"input_1 is neither an existing FASTQ path nor a supported INSDC run accession: {input_1}"
        )

    metadata = {
        key: value
        for key, value in row.items()
        if key not in RESERVED_METADATA_KEYS and value != ""
    }
    operations = metadata_operations + operations
    payload: dict[str, object] = {
        "library_id": args.library_id,
        "source_type": source_type,
        "source_id": source_id,
        "archive": archive,
        **public_metadata_status,
        "assay": assay,
        "assay_confidence": assay_confidence,
        "assay_reason": assay_reason,
        "layout": layout,
        "fastq_1": fastq_1,
        "fastq_2": fastq_2 or "",
        "input_1": input_1,
        "input_2": input_2,
        "metadata": metadata,
        "operations": operations,
        "materialized_at": datetime.now(timezone.utc).isoformat(),
    }
    write_metadata(metadata_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
