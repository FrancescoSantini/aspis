#!/usr/bin/env python3
"""Build transcript count and metadata tables from StringTie quantified GTFs."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_QUANT_COLUMNS = {"library_id", "quant_gtf", "status"}
REQUIRED_PLAN_COLUMNS = {"status", "read_length"}
ATTR_RE = re.compile(r'([A-Za-z0-9_.-]+) "([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quant-manifest", required=True, help="StringTie quantification manifest TSV")
    parser.add_argument("--plan", required=True, help="RNA-seq quantification plan TSV")
    parser.add_argument("--tmap", default="", help="Optional gffcompare tmap")
    parser.add_argument("--counts", required=True, help="Transcript count matrix TSV")
    parser.add_argument("--metadata", required=True, help="Transcript metadata TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--known-codes-strict", default="=")
    parser.add_argument("--known-codes-lenient", default="=,c,k,m,n,y")
    parser.add_argument("--gene-type-view", default="strict", choices=("strict", "lenient"))
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def read_plan(path: Path) -> dict[str, str]:
    _, rows = read_table(path, REQUIRED_PLAN_COLUMNS)
    if len(rows) != 1:
        raise ValueError(f"Quantification plan must contain exactly one row: {path}")
    row = rows[0]
    if row.get("status") != "ready":
        raise ValueError("Quantification plan is not ready: " + row.get("reason", ""))
    try:
        read_length = int(row["read_length"])
    except ValueError as exc:
        raise ValueError(f"read_length must be an integer: {row['read_length']!r}") from exc
    if read_length <= 0:
        raise ValueError(f"read_length must be > 0: {read_length}")
    row["read_length"] = str(read_length)
    return row


def parse_attrs(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(text)}


def mode(values: list[str]) -> str:
    clean = [value for value in values if value]
    if not clean:
        return ""
    return Counter(clean).most_common(1)[0][0]


def class_code_rank(code: str) -> int:
    order = ["=", "j", "c", "k", "m", "n", "y", "o", "e", "i", "u", "x", "s", "p", "r"]
    return order.index(code) if code in order else len(order)


def read_tmap(path: Path) -> dict[str, str]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        if reader.fieldnames is None:
            return {}
        query_column = "q_id" if "q_id" in reader.fieldnames else "qry_id"
        if query_column not in reader.fieldnames or "class_code" not in reader.fieldnames:
            return {}
        by_tx: dict[str, list[str]] = defaultdict(list)
        for row in reader:
            tx = (row.get(query_column) or "").strip()
            code = (row.get("class_code") or "").strip()
            if tx and code:
                by_tx[tx].append(code)
    return {
        tx: sorted(set(codes), key=class_code_rank)[0]
        for tx, codes in by_tx.items()
        if codes
    }


def transcript_length(feature: dict[str, object]) -> int:
    exons = feature["exons"]
    if exons:
        return sum(end - start + 1 for start, end in exons)
    return int(feature["end"]) - int(feature["start"]) + 1


def parse_gtf(path: Path, sample_id: str, read_length: int) -> tuple[dict[str, int], dict[str, dict[str, str]]]:
    transcripts: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature_type, start, end, _score, strand, _frame, attr_text = fields
            attrs = parse_attrs(attr_text)
            tx_id = attrs.get("transcript_id", "")
            gene_id = attrs.get("ref_gene_id") or attrs.get("gene_id", "")
            if not tx_id or not gene_id:
                continue
            feature = transcripts.setdefault(
                tx_id,
                {
                    "gene_id": gene_id,
                    "gene_name": attrs.get("ref_gene_name") or attrs.get("gene_name", ""),
                    "chr": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "cov": "",
                    "exons": [],
                },
            )
            feature["gene_id"] = gene_id or feature["gene_id"]
            if attrs.get("ref_gene_name") or attrs.get("gene_name"):
                feature["gene_name"] = attrs.get("ref_gene_name") or attrs.get("gene_name", "")
            if feature_type == "transcript":
                feature["chr"] = chrom
                feature["start"] = int(start)
                feature["end"] = int(end)
                feature["strand"] = strand
                feature["cov"] = attrs.get("cov", feature.get("cov", ""))
            elif feature_type == "exon":
                feature["exons"].append((int(start), int(end)))

    counts: dict[str, int] = {}
    metadata: dict[str, dict[str, str]] = {}
    for tx_id, feature in transcripts.items():
        length = max(0, transcript_length(feature))
        try:
            cov = float(feature["cov"] or 0)
        except ValueError:
            cov = 0.0
        count = int(math.ceil(cov * length / read_length))
        counts[tx_id] = max(0, count)
        metadata[tx_id] = {
            "transcript_id": tx_id,
            "gene_id": str(feature["gene_id"]),
            "gene_name": str(feature["gene_name"]),
            "Chr": str(feature["chr"]),
            "Start": str(feature["start"]),
            "End": str(feature["end"]),
            "Strand": str(feature["strand"]),
            "Length": str(length),
        }
    if not metadata:
        raise ValueError(f"No transcripts found in quantified GTF for {sample_id}: {path}")
    return counts, metadata


def classify_gene_type(
    gene_id: str,
    class_code: str,
    strict: set[str],
    lenient: set[str],
    view: str,
) -> tuple[str, str, str]:
    looks_known = bool(gene_id) and not gene_id.startswith(("MSTRG", "STRG"))
    strict_known = (class_code in strict) if class_code else looks_known
    lenient_known = (class_code in lenient) if class_code else looks_known
    strict_value = "Known" if strict_known else "Novel"
    lenient_value = "Known" if lenient_known else "Novel"
    return strict_value, lenient_value, strict_value if view == "strict" else lenient_value


def classify_transcript_discovery(gene_id: str, class_code: str) -> tuple[str, str, str, str, str, str]:
    """Map gffcompare class codes to an explicit transcript discovery class."""
    if class_code == "=":
        return (
            "known_transcript",
            "known",
            "no",
            "known_compatible",
            "Known/reference-compatible",
            "exact reference transcript match",
        )
    if class_code == "j":
        return (
            "novel_isoform_known_gene",
            "novel_isoform",
            "yes",
            "novel_isoform",
            "Novel isoform",
            "novel splice junction compatible with a known gene",
        )
    if class_code == "u":
        return (
            "intergenic_novel_locus",
            "novel_locus",
            "yes",
            "novel_locus",
            "Novel locus",
            "intergenic transcript with no reference overlap",
        )
    if class_code == "i":
        return (
            "intronic_novel_candidate",
            "ambiguous_overlap",
            "no",
            "ambiguous",
            "Ambiguous overlap",
            "transcript contained within a reference intron",
        )
    if class_code in {"x", "s"}:
        return (
            "antisense_novel_candidate",
            "ambiguous_overlap",
            "no",
            "ambiguous",
            "Ambiguous overlap",
            "antisense or opposite-strand overlap with reference annotation",
        )
    if class_code in {"c", "k", "m", "n", "y"}:
        return (
            "reference_contained_or_containing",
            "reference_overlap",
            "no",
            "known_compatible",
            "Known/reference-compatible",
            "contained, containing, retained-intron, or reference-compatible overlap",
        )
    if class_code in {"o", "e"}:
        return (
            "ambiguous_reference_overlap",
            "ambiguous_overlap",
            "no",
            "ambiguous",
            "Ambiguous overlap",
            "generic or exonic overlap with reference annotation",
        )
    if class_code in {"p", "r"}:
        return (
            "likely_artifact_or_repeat",
            "low_confidence",
            "no",
            "artifact",
            "Artifact/repeat",
            "possible polymerase run-on, pre-mRNA, or repeat-associated transcript",
        )
    if class_code:
        return (
            "unclassified_gffcompare_code",
            "unclassified",
            "no",
            "ambiguous",
            "Ambiguous overlap",
            f"unmapped gffcompare class code {class_code}",
        )
    if gene_id.startswith(("MSTRG", "STRG")):
        return (
            "unclassified_novel_candidate",
            "novel_locus",
            "yes",
            "novel_locus",
            "Novel locus",
            "StringTie novel gene without gffcompare class code",
        )
    return (
        "unclassified_reference_compatible",
        "reference_overlap",
        "no",
        "known_compatible",
        "Known/reference-compatible",
        "no gffcompare class code for reference-like gene id",
    )


def main() -> int:
    args = parse_args()
    plan = read_plan(Path(args.plan))
    read_length = int(plan["read_length"])
    _, quant_rows = read_table(Path(args.quant_manifest), REQUIRED_QUANT_COLUMNS)
    if not quant_rows:
        raise ValueError("StringTie quantification manifest has no rows")
    tmap_codes = read_tmap(Path(args.tmap)) if args.tmap else {}
    strict_codes = set(code.strip() for code in args.known_codes_strict.split(",") if code.strip())
    lenient_codes = set(code.strip() for code in args.known_codes_lenient.split(",") if code.strip())

    sample_ids = []
    counts_by_sample: dict[str, dict[str, int]] = {}
    metadata_by_tx: dict[str, dict[str, str]] = {}
    for row in quant_rows:
        if row.get("status") != "ok":
            raise ValueError(f"{row.get('library_id', '<unknown>')}: quantification status is not ok")
        sample_id = row["library_id"]
        quant_gtf = Path(row["quant_gtf"])
        if not quant_gtf.exists():
            raise FileNotFoundError(f"{sample_id}: quant_gtf does not exist: {quant_gtf}")
        sample_counts, sample_metadata = parse_gtf(quant_gtf, sample_id, read_length)
        sample_ids.append(sample_id)
        counts_by_sample[sample_id] = sample_counts
        metadata_by_tx.update(sample_metadata)

    tx_ids = sorted(metadata_by_tx)
    counts_path = Path(args.counts)
    counts_path.parent.mkdir(parents=True, exist_ok=True)
    with counts_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["transcript_id"] + sample_ids,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for tx_id in tx_ids:
            row = {"transcript_id": tx_id}
            for sample_id in sample_ids:
                row[sample_id] = str(counts_by_sample.get(sample_id, {}).get(tx_id, 0))
            writer.writerow(row)

    metadata_columns = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "class_code",
        "transcript_discovery_class",
        "transcript_novelty",
        "true_novel_candidate",
        "transcript_plot_group",
        "transcript_plot_label",
        "gffcompare_description",
        "gene_type_strict",
        "gene_type_lenient",
        "gene_type",
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length",
    ]
    metadata_path = Path(args.metadata)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for tx_id in tx_ids:
            row = dict(metadata_by_tx[tx_id])
            class_code = tmap_codes.get(tx_id, "")
            strict_value, lenient_value, selected = classify_gene_type(
                row["gene_id"],
                class_code,
                strict_codes,
                lenient_codes,
                args.gene_type_view,
            )
            discovery_class, novelty, true_novel, plot_group, plot_label, description = classify_transcript_discovery(row["gene_id"], class_code)
            row.update(
                {
                    "class_code": class_code,
                    "transcript_discovery_class": discovery_class,
                    "transcript_novelty": novelty,
                    "true_novel_candidate": true_novel,
                    "transcript_plot_group": plot_group,
                    "transcript_plot_label": plot_label,
                    "gffcompare_description": description,
                    "gene_type_strict": strict_value,
                    "gene_type_lenient": lenient_value,
                    "gene_type": selected,
                }
            )
            writer.writerow({column: row.get(column, "") for column in metadata_columns})

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(f"status\ttranscripts\tlibraries\nok\t{len(tx_ids)}\t{len(sample_ids)}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
