#!/usr/bin/env python3
"""Generate the deterministic, redistribution-safe ASPIS synthetic project."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
READ_LENGTH = 75
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCC"
SEED = 20261022


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def deterministic_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fastq_gz(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0, compresslevel=9) as compressed:
        for name, sequence in records:
            quality = "I" * len(sequence)
            compressed.write(f"@{name}\n{sequence}\n+\n{quality}\n".encode())
    path.write_bytes(buffer.getvalue())


def wrap_fasta(name: str, sequence: str) -> str:
    wrapped = "\n".join(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    return f">{name}\n{wrapped}\n"


def rnaseq_records(
    genome: str, sample: str, counts: dict[str, int], replicate: int
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    intervals = {"syn_gene_a": (200, 1200), "syn_gene_b": (1600, 2600), "syn_gene_c": (3000, 4200)}
    first: list[tuple[str, str]] = []
    second: list[tuple[str, str]] = []
    for gene, count in counts.items():
        start, end = intervals[gene]
        span = end - start
        for index in range(count):
            fragment_start = start + ((index * 17 + replicate * 11) % (span - 230))
            fragment = genome[fragment_start : fragment_start + 225]
            name = f"{sample}:{gene}:{index + 1}"
            first.append((f"{name}/1", fragment[:READ_LENGTH]))
            second.append((f"{name}/2", reverse_complement(fragment[-READ_LENGTH:])))
    return first, second


def smallrna_records(sample: str, counts: dict[str, int], sequences: dict[str, str]) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for feature, count in counts.items():
        insert = sequences[feature]
        for index in range(count):
            records.append((f"{sample}:{feature}:{index + 1}", insert + ADAPTER))
    return records


def main() -> int:
    genome = deterministic_dna(5000, SEED)
    mirnas = {
        "syn-miR-1-5p": "TGGAATGTAAAGAAGTATGTAT",
        "syn-miR-2-3p": "TAAAGTGCTTATAGTGCAGGTA",
    }
    contaminant = "GGTTCGATCCCGGCTCAG"
    residual = genome[4500:4524]
    small_sequences = {**mirnas, "contaminant_tRNA": contaminant, "residual_fragment": residual}

    write_text(DATA / "reference" / "genome.fa", wrap_fasta("syn_chr1", genome))
    gtf_rows = []
    for gene, start, end, biotype in (
        ("syn_gene_a", 201, 1200, "protein_coding"),
        ("syn_gene_b", 1601, 2600, "protein_coding"),
        ("syn_gene_c", 3001, 4200, "lncRNA"),
        ("syn_snoRNA", 4501, 4540, "snoRNA"),
    ):
        tx = f"{gene}_tx"
        attrs = f'gene_id "{gene}"; gene_name "{gene}"; gene_biotype "{biotype}";'
        tx_attrs = attrs + f' transcript_id "{tx}"; transcript_biotype "{biotype}";'
        gtf_rows.extend(
            [
                f"syn_chr1\tASPIS\tgene\t{start}\t{end}\t.\t+\t.\t{attrs}",
                f"syn_chr1\tASPIS\ttranscript\t{start}\t{end}\t.\t+\t.\t{tx_attrs}",
                f"syn_chr1\tASPIS\texon\t{start}\t{end}\t.\t+\t.\t{tx_attrs} exon_number \"1\";",
            ]
        )
    write_text(DATA / "reference" / "annotation.gtf", "\n".join(gtf_rows) + "\n")
    write_text(
        DATA / "reference" / "mature_mirna.fa",
        "".join(wrap_fasta(name, sequence) for name, sequence in mirnas.items()),
    )
    write_text(DATA / "reference" / "contaminants.fa", wrap_fasta("synthetic_tRNA", contaminant))
    write_tsv(
        DATA / "reference" / "mirna_targets.tsv",
        ["mirna_id", "target_id", "target_symbol", "database", "source", "source_type", "target_evidence_type", "resource_version", "evidence"],
        [
            dict(mirna_id="syn-miR-1-5p", target_id="syn_gene_b", target_symbol="syn_gene_b", database="ASPIS_SYNTHETIC", source="deterministic_fixture", source_type="user_provided", target_evidence_type="user_provided", resource_version="1", evidence="inverse_input_signal"),
            dict(mirna_id="syn-miR-2-3p", target_id="syn_gene_a", target_symbol="syn_gene_a", database="ASPIS_SYNTHETIC", source="deterministic_fixture", source_type="user_provided", target_evidence_type="user_provided", resource_version="1", evidence="inverse_input_signal"),
        ],
    )

    rnaseq_truth = {
        "control": {"syn_gene_a": 120, "syn_gene_b": 30, "syn_gene_c": 60},
        "treated": {"syn_gene_a": 30, "syn_gene_b": 120, "syn_gene_c": 60},
    }
    smallrna_truth = {
        "control": {"syn-miR-1-5p": 120, "syn-miR-2-3p": 30, "contaminant_tRNA": 20, "residual_fragment": 20},
        "treated": {"syn-miR-1-5p": 30, "syn-miR-2-3p": 120, "contaminant_tRNA": 20, "residual_fragment": 20},
    }
    intake_rows: list[dict[str, object]] = []
    for condition in ("control", "treated"):
        for replicate in (1, 2):
            biospecimen = f"synthetic_{condition}_r{replicate}"
            rna_id = f"{biospecimen}_rnaseq"
            r1, r2 = rnaseq_records(genome, rna_id, rnaseq_truth[condition], replicate)
            r1_path = DATA / "fastq" / f"{rna_id}_R1.fastq.gz"
            r2_path = DATA / "fastq" / f"{rna_id}_R2.fastq.gz"
            write_fastq_gz(r1_path, r1)
            write_fastq_gz(r2_path, r2)
            intake_rows.append(
                dict(library_id=rna_id, biospecimen_id=biospecimen, project="ASPIS_SYNTHETIC", input_1=r1_path.relative_to(ROOT.parent.parent).as_posix(), input_2=r2_path.relative_to(ROOT.parent.parent).as_posix(), assay_hint="rnaseq", condition=condition, treatment="synthetic_signal" if condition == "treated" else "none", dose=1 if condition == "treated" else 0, dose_unit="arbitrary_unit", time_h=24, replicate=replicate, batch="synthetic_batch")
            )
            small_id = f"{biospecimen}_smallrna"
            small_path = DATA / "fastq" / f"{small_id}.fastq.gz"
            write_fastq_gz(small_path, smallrna_records(small_id, smallrna_truth[condition], small_sequences))
            intake_rows.append(
                dict(library_id=small_id, biospecimen_id=biospecimen, project="ASPIS_SYNTHETIC", input_1=small_path.relative_to(ROOT.parent.parent).as_posix(), input_2="", assay_hint="smallrna", condition=condition, treatment="synthetic_signal" if condition == "treated" else "none", dose=1 if condition == "treated" else 0, dose_unit="arbitrary_unit", time_h=24, replicate=replicate, batch="synthetic_batch")
            )

    columns = ["library_id", "biospecimen_id", "project", "input_1", "input_2", "assay_hint", "condition", "treatment", "dose", "dose_unit", "time_h", "replicate", "batch"]
    write_tsv(ROOT / "intake.tsv", columns, intake_rows)
    write_tsv(
        ROOT / "truth_rnaseq_counts.tsv",
        ["condition", "feature_id", "fragments_per_replicate", "expected_direction"],
        [dict(condition=condition, feature_id=feature, fragments_per_replicate=count, expected_direction="down" if feature == "syn_gene_a" and condition == "treated" else "up" if feature == "syn_gene_b" and condition == "treated" else "unchanged") for condition, counts in rnaseq_truth.items() for feature, count in counts.items()],
    )
    write_tsv(
        ROOT / "truth_smallrna_counts.tsv",
        ["condition", "feature_id", "reads_per_replicate", "expected_direction"],
        [dict(condition=condition, feature_id=feature, reads_per_replicate=count, expected_direction="down" if feature == "syn-miR-1-5p" and condition == "treated" else "up" if feature == "syn-miR-2-3p" and condition == "treated" else "unchanged") for condition, counts in smallrna_truth.items() for feature, count in counts.items()],
    )
    checksums = []
    for path in sorted(DATA.rglob("*")):
        if path.is_file():
            checksums.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT).as_posix()}")
    write_text(ROOT / "checksums.sha256", "\n".join(checksums) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
