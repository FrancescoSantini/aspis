#!/usr/bin/env python3
"""Prepare a local smallRNA reference FASTA and SAF annotation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


MANIFEST_COLUMNS = [
    "source_fasta",
    "output_fasta",
    "saf",
    "species_prefix",
    "replace_u_with_t",
    "n_records",
    "total_bases",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="Input miRBase-like FASTA")
    parser.add_argument("--output-fasta", required=True, help="Normalized output FASTA")
    parser.add_argument("--saf", required=True, help="Output SAF annotation")
    parser.add_argument("--manifest", required=True, help="Reference manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument(
        "--species-prefix",
        default="",
        help="Optional FASTA record ID prefix to keep, e.g. hsa for human miRBase records",
    )
    parser.add_argument(
        "--replace-u-with-t",
        action="store_true",
        help="Convert RNA U bases to DNA T bases in output FASTA",
    )
    return parser.parse_args()


def fasta_records(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header = ""
    seq_chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append(record_from_parts(header, seq_chunks))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header:
        records.append(record_from_parts(header, seq_chunks))
    return records


def record_from_parts(header: str, seq_chunks: list[str]) -> tuple[str, str, str]:
    record_id = header.split()[0]
    if not record_id:
        raise ValueError("FASTA record has an empty identifier")
    sequence = "".join(seq_chunks).replace(" ", "").upper()
    if not sequence:
        raise ValueError(f"FASTA record {record_id!r} has an empty sequence")
    return record_id, header, sequence


def filter_records(
    records: list[tuple[str, str, str]],
    species_prefix: str,
    replace_u_with_t: bool,
) -> list[tuple[str, str, str]]:
    kept = []
    seen = set()
    for record_id, header, sequence in records:
        if species_prefix and not record_id.startswith(species_prefix):
            continue
        if record_id in seen:
            raise ValueError(f"Duplicate FASTA record ID after filtering: {record_id}")
        seen.add(record_id)
        if replace_u_with_t:
            sequence = sequence.replace("U", "T")
        invalid = sorted(set(sequence) - set("ACGTNRYKMSWBDHVU"))
        if invalid:
            raise ValueError(f"FASTA record {record_id!r} has unsupported base(s): {invalid}")
        kept.append((record_id, header, sequence))
    if not kept:
        prefix = f" with prefix {species_prefix!r}" if species_prefix else ""
        raise ValueError(f"No FASTA records kept{prefix}")
    return kept


def write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record_id, header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def write_saf(path: Path, records: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["GeneID", "Chr", "Start", "End", "Strand"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record_id, _header, sequence in records:
            writer.writerow(
                {
                    "GeneID": record_id,
                    "Chr": record_id,
                    "Start": "1",
                    "End": str(len(sequence)),
                    "Strand": "+",
                }
            )


def write_manifest(path: Path, args: argparse.Namespace, records: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "source_fasta": args.fasta,
        "output_fasta": args.output_fasta,
        "saf": args.saf,
        "species_prefix": args.species_prefix,
        "replace_u_with_t": str(args.replace_u_with_t).lower(),
        "n_records": str(len(records)),
        "total_bases": str(sum(len(sequence) for _record_id, _header, sequence in records)),
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def write_done(path: Path, records: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\trecords\n")
        handle.write(f"ok\t{len(records)}\n")


def main() -> int:
    args = parse_args()
    source = Path(args.fasta)
    if not source.exists():
        raise FileNotFoundError(f"Input FASTA does not exist: {source}")
    records = filter_records(
        fasta_records(source),
        args.species_prefix,
        args.replace_u_with_t,
    )
    write_fasta(Path(args.output_fasta), records)
    write_saf(Path(args.saf), records)
    write_manifest(Path(args.manifest), args, records)
    write_done(Path(args.done), records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
