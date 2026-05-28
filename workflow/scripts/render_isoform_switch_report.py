#!/usr/bin/env python3
"""Render event-level isoform-switch report assets.

This script is intentionally downstream of IsoformSwitchAnalyzeR. It does not
re-score switches; it turns the per-contrast switch outputs into a stable,
auditable report with candidate tables, sequence summaries, optional functional
annotation imports, and one HTML/SVG page per switch event.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


MANIFEST_REQUIRED = {
    "contrast_id",
    "status",
    "reason",
    "detailed",
    "consequences",
    "nt_fasta",
    "aa_fasta",
}
TRANSCRIPT_METADATA_REQUIRED = {"transcript_id", "gene_id"}
CANDIDATE_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "isoform_id",
    "switch_role",
    "dIF",
    "switch_statistic",
    "switch_statistic_name",
    "candidate_status",
    "transcript_discovery_class",
    "transcript_novelty",
    "transcript_plot_group",
    "gffcompare_class_code",
    "consequence_summary",
    "source_detailed",
]
EVENT_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "status",
    "reason",
    "switch_in_isoform",
    "switch_out_isoform",
    "switch_in_dIF",
    "switch_out_dIF",
    "max_abs_dIF",
    "best_switch_statistic",
    "best_switch_statistic_name",
    "n_isoforms_in_gene",
    "n_candidate_isoforms",
    "n_switch_consequences",
    "n_functional_annotations",
    "plot_svg",
    "event_html",
]
SEQUENCE_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "isoform_id",
    "switch_role",
    "nt_length",
    "aa_length",
    "nt_sequence",
    "aa_sequence",
    "sequence_status",
]
ANNOTATION_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "isoform_id",
    "source",
    "feature_type",
    "feature_id",
    "feature_name",
    "start_aa",
    "end_aa",
    "score",
    "description",
    "status",
]
PLOT_MANIFEST_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "status",
    "reason",
    "plot_svg",
    "event_html",
    "n_isoforms",
    "n_candidate_isoforms",
    "nt_fasta",
    "aa_fasta",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Isoform-switch manifest TSV")
    parser.add_argument("--transcript-metadata", required=True, help="Transcript metadata TSV")
    parser.add_argument("--annotated-gtf", required=True, help="Annotated transcriptome GTF")
    parser.add_argument("--outdir", required=True, help="Report output directory")
    parser.add_argument("--candidate-table", required=True, help="Output candidate isoform TSV")
    parser.add_argument("--event-summary", required=True, help="Output event summary TSV")
    parser.add_argument("--sequence-table", required=True, help="Output event sequence TSV")
    parser.add_argument("--functional-annotation-table", required=True, help="Output normalized annotation TSV")
    parser.add_argument("--plot-manifest", required=True, help="Output plot manifest TSV")
    parser.add_argument("--html", required=True, help="Output project-level HTML report")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--padj", type=float, default=0.1, help="Maximum switch q-value/p-value")
    parser.add_argument("--dif", type=float, default=0.1, help="Minimum absolute dIF")
    parser.add_argument("--top-n", type=int, default=30, help="Maximum switch events to render")
    parser.add_argument(
        "--functional-annotation-tables",
        default="",
        help=(
            "Optional comma-separated TSVs with protein/domain annotations. "
            "Supported identifiers: isoform_id, transcript_id, or protein_id."
        ),
    )
    return parser.parse_args()


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    cleaned = cleaned.strip("._")
    return cleaned or "unknown"


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def read_table(path: Path, required: Optional[set[str]] = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Input TSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        fieldnames = list(reader.fieldnames)
        missing = (required or set()) - set(fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    return fieldnames, rows


def write_table(path: Path, columns: list[str], rows: Iterable[dict[str, object]]) -> None:
    ensure_parent(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def first_existing(row: dict[str, str], names: Iterable[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value != "":
            return value
    return ""


def first_column(fieldnames: Iterable[str], candidates: Iterable[str]) -> str:
    fields = list(fieldnames)
    lowered = {field.lower(): field for field in fields}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return ""


def find_statistic_column(fieldnames: Iterable[str]) -> str:
    preferred = [
        "isoform_switch_q_value",
        "switch_q_value",
        "q_value",
        "qvalue",
        "adj_pvalue",
        "adj_p_value",
        "padj",
        "FDR",
        "p_value",
        "pvalue",
        "pval",
        "p_val",
    ]
    return first_column(fieldnames, preferred)


def to_float(value: str) -> Optional[float]:
    if value is None or str(value).strip() == "":
        return None
    try:
        parsed = float(str(value))
    except ValueError:
        return None
    if math.isnan(parsed):
        return None
    return parsed


def parse_fasta(path: Path) -> dict[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    sequences: dict[str, list[str]] = {}
    current = ""
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                sequences.setdefault(current, [])
            elif current:
                sequences[current].append(line)
    joined = {key: "".join(value) for key, value in sequences.items()}
    aliases = dict(joined)
    for key, sequence in joined.items():
        for alias in {key.split("|")[0], key.split(".")[0], key.rsplit("|", 1)[-1]}:
            if alias and alias not in aliases:
                aliases[alias] = sequence
    return aliases


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]


def parse_gtf(path: Path) -> dict[str, TranscriptModel]:
    models: dict[str, TranscriptModel] = {}
    if not path.exists() or path.stat().st_size == 0:
        return models
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start_text, end_text, _, strand, _, attributes = parts
            if feature != "exon":
                continue
            attrs = parse_gtf_attributes(attributes)
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id)
            model = models.setdefault(
                transcript_id,
                TranscriptModel(transcript_id, gene_id, gene_name, chrom, strand, []),
            )
            model.exons.append((int(start_text), int(end_text)))
    for model in models.values():
        model.exons.sort()
    return models


def metadata_by_transcript(path: Path) -> dict[str, dict[str, str]]:
    _, rows = read_table(path, TRANSCRIPT_METADATA_REQUIRED)
    indexed = {}
    for row in rows:
        transcript_id = row.get("transcript_id", "")
        if transcript_id:
            indexed[transcript_id] = row
    return indexed


def consequence_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, str], list[str]]:
    lookup: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in rows:
        isoform_id = first_existing(row, ["isoform_id", "transcript_id", "isoform"])
        gene_id = first_existing(row, ["gene_id", "geneID", "gene"])
        text_parts = []
        for key in ["consequence", "switch_consequence", "category", "featureCompared", "condition_1", "condition_2"]:
            if row.get(key, ""):
                text_parts.append(f"{key}={row[key]}")
        if not text_parts:
            text_parts = [f"{key}={value}" for key, value in row.items() if value]
        text = "; ".join(text_parts)
        if isoform_id:
            lookup[(gene_id, isoform_id)].append(text)
            lookup[("", isoform_id)].append(text)
        if gene_id:
            lookup[(gene_id, "")].append(text)
    return lookup


def normalized_annotations(
    paths_text: str,
    candidate_isoforms: set[str],
    event_by_isoform: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for item in [part.strip() for part in paths_text.split(",") if part.strip()]:
        path = Path(item)
        if not path.exists():
            raise FileNotFoundError(f"Functional annotation table does not exist: {path}")
        fieldnames, raw_rows = read_table(path)
        iso_col = first_column(fieldnames, ["isoform_id", "transcript_id", "protein_id", "query", "query_id"])
        if not iso_col:
            raise ValueError(f"Annotation table lacks isoform/transcript/protein identifier column: {path}")
        for raw in raw_rows:
            isoform_id = raw.get(iso_col, "")
            if not isoform_id:
                continue
            isoform_aliases = {isoform_id, isoform_id.split("|")[0], isoform_id.split(".")[0]}
            matched = sorted(alias for alias in isoform_aliases if alias in candidate_isoforms)
            if not matched:
                continue
            source = first_existing(raw, ["source", "analysis", "database", "signature_desc"])
            feature_id = first_existing(raw, ["feature_id", "accession", "signature_accession", "interpro_accession", "id"])
            feature_name = first_existing(raw, ["feature_name", "name", "signature_description", "interpro_description"])
            feature_type = first_existing(raw, ["feature_type", "type", "analysis", "database"]) or "protein_feature"
            start_aa = first_existing(raw, ["start_aa", "start", "from", "ali_start", "hmm_start"])
            end_aa = first_existing(raw, ["end_aa", "end", "to", "ali_end", "hmm_end"])
            score = first_existing(raw, ["score", "evalue", "e_value", "bitscore"])
            description = first_existing(raw, ["description", "desc", "interpro_description", "signature_description"])
            for matched_id in matched:
                for event in event_by_isoform.get(matched_id, []):
                    rows.append(
                        {
                            "event_id": event["event_id"],
                            "contrast_id": event["contrast_id"],
                            "gene_id": event["gene_id"],
                            "gene_name": event["gene_name"],
                            "isoform_id": matched_id,
                            "source": source or path.name,
                            "feature_type": feature_type,
                            "feature_id": feature_id,
                            "feature_name": feature_name,
                            "start_aa": start_aa,
                            "end_aa": end_aa,
                            "score": score,
                            "description": description,
                            "status": "ok",
                        }
                    )
    rows.sort(key=lambda row: (row["event_id"], row["isoform_id"], row["source"], row["start_aa"]))
    return rows


def role_for_isoform(isoform_id: str, switch_in: str, switch_out: str) -> str:
    if isoform_id == switch_in:
        return "switch_in"
    if isoform_id == switch_out:
        return "switch_out"
    return "same_gene"


def build_events(
    manifest_rows: list[dict[str, str]],
    metadata: dict[str, dict[str, str]],
    padj: float,
    dif: float,
    top_n: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, dict[str, str]], dict[str, str], dict[str, str]]:
    candidate_rows: list[dict[str, str]] = []
    event_rows: list[dict[str, str]] = []
    event_context: dict[str, dict[str, str]] = {}
    nt_sequences: dict[str, str] = {}
    aa_sequences: dict[str, str] = {}
    event_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    all_gene_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)

    for manifest in manifest_rows:
        if manifest.get("status") != "ok":
            continue
        detailed_path = Path(manifest.get("detailed", ""))
        if not detailed_path.exists() or detailed_path.stat().st_size == 0:
            continue
        fieldnames, detailed_rows = read_table(detailed_path)
        d_if_col = first_column(fieldnames, ["dIF", "d_if", "delta_isoform_fraction"])
        statistic_col = find_statistic_column(fieldnames)
        isoform_col = first_column(fieldnames, ["isoform_id", "transcript_id"])
        gene_col = first_column(fieldnames, ["gene_id", "gene_id.x", "gene_id.y", "geneID"])
        gene_name_col = first_column(fieldnames, ["gene_name", "gene_name.x", "gene_name.y", "gene_symbol"])
        if not d_if_col or not statistic_col or not isoform_col or not gene_col:
            continue

        nt_sequences.update(parse_fasta(Path(manifest.get("nt_fasta", ""))))
        aa_sequences.update(parse_fasta(Path(manifest.get("aa_fasta", ""))))

        for row in detailed_rows:
            isoform_id = row.get(isoform_col, "")
            gene_id = row.get(gene_col, "")
            if not isoform_id or not gene_id:
                continue
            d_if = to_float(row.get(d_if_col, ""))
            statistic = to_float(row.get(statistic_col, ""))
            if d_if is None:
                continue
            clean_row = {
                **row,
                "_contrast_id": manifest["contrast_id"],
                "_dIF": str(d_if),
                "_statistic": "" if statistic is None else str(statistic),
                "_statistic_name": statistic_col,
                "_isoform_id": isoform_id,
                "_gene_id": gene_id,
                "_gene_name": row.get(gene_name_col, "") if gene_name_col else "",
                "_detailed": str(detailed_path),
            }
            all_gene_rows[(manifest["contrast_id"], gene_id)].append(clean_row)
            if statistic is not None and statistic <= padj and abs(d_if) >= dif:
                event_groups[(manifest["contrast_id"], gene_id)].append(clean_row)

    ranked_events = sorted(
        event_groups.items(),
        key=lambda item: max(abs(float(row["_dIF"])) for row in item[1]),
        reverse=True,
    )
    for (contrast_id, gene_id), candidates in ranked_events[:top_n]:
        gene_rows = all_gene_rows[(contrast_id, gene_id)]
        gene_rows = sorted(gene_rows, key=lambda row: float(row["_dIF"]))
        switch_out = gene_rows[0]
        switch_in = gene_rows[-1]
        gene_name = switch_in.get("_gene_name") or switch_out.get("_gene_name") or gene_id
        event_id = f"{safe_id(contrast_id)}__{safe_id(gene_id)}"
        candidate_isoforms = {row["_isoform_id"] for row in candidates}
        best_stat = min(
            [to_float(row.get("_statistic", "")) for row in candidates if to_float(row.get("_statistic", "")) is not None],
            default=None,
        )
        event = {
            "event_id": event_id,
            "contrast_id": contrast_id,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "status": "ok",
            "reason": "",
            "switch_in_isoform": switch_in["_isoform_id"],
            "switch_out_isoform": switch_out["_isoform_id"],
            "switch_in_dIF": switch_in["_dIF"],
            "switch_out_dIF": switch_out["_dIF"],
            "max_abs_dIF": str(max(abs(float(row["_dIF"])) for row in candidates)),
            "best_switch_statistic": "" if best_stat is None else str(best_stat),
            "best_switch_statistic_name": candidates[0].get("_statistic_name", ""),
            "n_isoforms_in_gene": str(len({row["_isoform_id"] for row in gene_rows})),
            "n_candidate_isoforms": str(len(candidate_isoforms)),
            "n_switch_consequences": "0",
            "n_functional_annotations": "0",
            "plot_svg": "",
            "event_html": "",
        }
        event_rows.append(event)
        event_context[event_id] = {
            **event,
            "_gene_rows": gene_rows,  # type: ignore[dict-item]
            "_candidate_isoforms": candidate_isoforms,  # type: ignore[dict-item]
        }
        for row in gene_rows:
            isoform_id = row["_isoform_id"]
            meta = metadata.get(isoform_id, {})
            status = "candidate" if isoform_id in candidate_isoforms else "same_gene_context"
            candidate_rows.append(
                {
                    "event_id": event_id,
                    "contrast_id": contrast_id,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "isoform_id": isoform_id,
                    "switch_role": role_for_isoform(isoform_id, switch_in["_isoform_id"], switch_out["_isoform_id"]),
                    "dIF": row["_dIF"],
                    "switch_statistic": row["_statistic"],
                    "switch_statistic_name": row["_statistic_name"],
                    "candidate_status": status,
                    "transcript_discovery_class": meta.get("transcript_discovery_class", ""),
                    "transcript_novelty": meta.get("transcript_novelty", ""),
                    "transcript_plot_group": meta.get("transcript_plot_group", ""),
                    "gffcompare_class_code": meta.get("gffcompare_class_code", meta.get("class_code", "")),
                    "consequence_summary": "",
                    "source_detailed": row["_detailed"],
                }
            )
    return candidate_rows, event_rows, event_context, nt_sequences, aa_sequences


def attach_consequences(
    candidate_rows: list[dict[str, str]],
    event_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
) -> None:
    all_consequences: list[dict[str, str]] = []
    for manifest in manifest_rows:
        if manifest.get("status") == "ok" and manifest.get("consequences", ""):
            path = Path(manifest["consequences"])
            if path.exists() and path.stat().st_size > 0:
                _, rows = read_table(path)
                all_consequences.extend(rows)
    lookup = consequence_lookup(all_consequences)
    event_counts: dict[str, int] = defaultdict(int)
    for row in candidate_rows:
        gene_id = row["gene_id"]
        isoform_id = row["isoform_id"]
        pieces = []
        for key in [(gene_id, isoform_id), ("", isoform_id), (gene_id, "")]:
            pieces.extend(lookup.get(key, []))
        unique = []
        for item in pieces:
            if item and item not in unique:
                unique.append(item)
        row["consequence_summary"] = " | ".join(unique[:8])
        event_counts[row["event_id"]] += len(unique)
    for row in event_rows:
        row["n_switch_consequences"] = str(event_counts.get(row["event_id"], 0))


def build_sequence_rows(
    candidate_rows: list[dict[str, str]],
    nt_sequences: dict[str, str],
    aa_sequences: dict[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen = set()
    for row in candidate_rows:
        key = (row["event_id"], row["isoform_id"])
        if key in seen:
            continue
        seen.add(key)
        isoform_id = row["isoform_id"]
        nt_seq = nt_sequences.get(isoform_id, "")
        aa_seq = aa_sequences.get(isoform_id, "")
        if nt_seq and aa_seq:
            status = "nt_and_aa"
        elif nt_seq:
            status = "nt_only"
        elif aa_seq:
            status = "aa_only"
        else:
            status = "missing"
        rows.append(
            {
                "event_id": row["event_id"],
                "contrast_id": row["contrast_id"],
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "isoform_id": isoform_id,
                "switch_role": row["switch_role"],
                "nt_length": str(len(nt_seq)) if nt_seq else "0",
                "aa_length": str(len(aa_seq)) if aa_seq else "0",
                "nt_sequence": nt_seq,
                "aa_sequence": aa_seq,
                "sequence_status": status,
            }
        )
    return rows


def relative(path: str, base: Path) -> str:
    return os.path.relpath(path, start=base.parent)


def svg_text(x: float, y: float, text: str, size: int = 12, weight: str = "400") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-weight="{weight}" font-family="Arial, sans-serif">{html.escape(text)}</text>'
    )


def render_event_svg(
    event: dict[str, str],
    event_context: dict[str, object],
    models: dict[str, TranscriptModel],
    annotations: list[dict[str, str]],
    path: Path,
) -> None:
    gene_rows: list[dict[str, str]] = event_context["_gene_rows"]  # type: ignore[assignment]
    isoforms = [row["_isoform_id"] for row in gene_rows]
    isoforms = list(dict.fromkeys(isoforms))
    relevant_models = [models[isoform] for isoform in isoforms if isoform in models and models[isoform].exons]
    width = 1180
    row_height = 58
    margin_left = 210
    margin_right = 40
    top = 72
    height = top + max(1, len(isoforms)) * row_height + 70
    path.parent.mkdir(parents=True, exist_ok=True)
    pieces = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(24, 30, f"{event['gene_name']} ({event['gene_id']})", 18, "700"),
        svg_text(24, 52, f"{event['contrast_id']} | max abs dIF {event['max_abs_dIF']}", 12),
    ]
    if not relevant_models:
        pieces.append(svg_text(24, 100, "No exon models were available in the annotated GTF.", 14, "700"))
        pieces.append("</svg>")
        path.write_text("\n".join(pieces), encoding="utf-8")
        return
    start = min(start for model in relevant_models for start, _ in model.exons)
    end = max(end for model in relevant_models for _, end in model.exons)
    span = max(1, end - start + 1)
    scale_width = width - margin_left - margin_right

    def x_for(position: int) -> float:
        return margin_left + ((position - start) / span) * scale_width

    annotations_by_isoform: dict[str, list[dict[str, str]]] = defaultdict(list)
    for annotation in annotations:
        annotations_by_isoform[annotation["isoform_id"]].append(annotation)

    d_if_by_isoform = {row["_isoform_id"]: row["_dIF"] for row in gene_rows}
    for index, isoform_id in enumerate(isoforms):
        y = top + index * row_height
        role = role_for_isoform(isoform_id, event["switch_in_isoform"], event["switch_out_isoform"])
        color = {"switch_in": "#1f77b4", "switch_out": "#d62728"}.get(role, "#57606a")
        model = models.get(isoform_id)
        label = f"{isoform_id} ({role}, dIF {d_if_by_isoform.get(isoform_id, '')})"
        pieces.append(svg_text(24, y + 14, label, 12, "700" if role != "same_gene" else "400"))
        pieces.append(f'<line x1="{margin_left}" x2="{width - margin_right}" y1="{y}" y2="{y}" stroke="#8c959f" stroke-width="1"/>')
        if model:
            for exon_start, exon_end in model.exons:
                exon_x = x_for(exon_start)
                exon_w = max(3, x_for(exon_end) - exon_x)
                pieces.append(
                    f'<rect x="{exon_x:.1f}" y="{y - 9:.1f}" width="{exon_w:.1f}" height="18" '
                    f'fill="{color}" fill-opacity="0.82" stroke="#24292f" stroke-width="0.5"/>'
                )
        anno_y = y + 22
        anno_rows = annotations_by_isoform.get(isoform_id, [])
        if anno_rows:
            pieces.append(f'<line x1="{margin_left}" x2="{width - margin_right}" y1="{anno_y}" y2="{anno_y}" stroke="#d0d7de"/>')
            aa_lengths = [
                to_float(a.get("end_aa", "")) or to_float(a.get("start_aa", "")) or 0
                for a in anno_rows
            ]
            aa_span = max(aa_lengths + [1])
            for annotation in anno_rows[:12]:
                start_aa = to_float(annotation.get("start_aa", "")) or 1
                end_aa = to_float(annotation.get("end_aa", "")) or start_aa
                feature_x = margin_left + ((start_aa - 1) / aa_span) * scale_width
                feature_w = max(8, ((end_aa - start_aa + 1) / aa_span) * scale_width)
                label_text = annotation.get("feature_name") or annotation.get("feature_id") or annotation.get("source")
                pieces.append(
                    f'<rect x="{feature_x:.1f}" y="{anno_y - 7:.1f}" width="{feature_w:.1f}" height="14" '
                    'fill="#8250df" fill-opacity="0.7" stroke="#6639ba" stroke-width="0.5">'
                    f'<title>{html.escape(label_text)}</title></rect>'
                )
        else:
            pieces.append(svg_text(margin_left, anno_y + 4, "No imported protein/domain annotations", 10))
    pieces.append(svg_text(margin_left, height - 22, f"Genomic span: {start}-{end}", 11))
    pieces.append("</svg>")
    path.write_text("\n".join(pieces), encoding="utf-8")


def render_event_html(
    event: dict[str, str],
    candidate_rows: list[dict[str, str]],
    sequence_rows: list[dict[str, str]],
    annotation_rows: list[dict[str, str]],
    out_path: Path,
    svg_path: Path,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    event_candidates = [row for row in candidate_rows if row["event_id"] == event["event_id"]]
    event_sequences = [row for row in sequence_rows if row["event_id"] == event["event_id"]]
    event_annotations = [row for row in annotation_rows if row["event_id"] == event["event_id"]]

    def table(headers: list[str], rows: list[dict[str, str]]) -> str:
        if not rows:
            return "<p>No rows.</p>"
        body = []
        for row in rows:
            body.append("<tr>" + "".join(f"<td>{html.escape(row.get(header, ''))}</td>" for header in headers) + "</tr>")
        return (
            "<table><thead><tr>"
            + "".join(f"<th>{html.escape(header)}</th>" for header in headers)
            + "</tr></thead><tbody>"
            + "\n".join(body)
            + "</tbody></table>"
        )

    sequence_blocks = []
    for row in event_sequences:
        sequence_blocks.append(
            "<details>"
            f"<summary>{html.escape(row['isoform_id'])} | {html.escape(row['switch_role'])} | "
            f"nt {html.escape(row['nt_length'])}, aa {html.escape(row['aa_length'])}</summary>"
            "<h3>Nucleotide sequence</h3>"
            f"<pre>{html.escape(row['nt_sequence'] or 'not available')}</pre>"
            "<h3>Amino-acid sequence</h3>"
            f"<pre>{html.escape(row['aa_sequence'] or 'not available')}</pre>"
            "</details>"
        )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(event['event_id'])}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1280px; }}
    h1 {{ margin-bottom: 4px; }}
    .muted {{ color: #57606a; }}
    table {{ border-collapse: collapse; width: 100%; margin: 12px 0 24px; }}
    th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; }}
    pre {{ white-space: pre-wrap; word-break: break-word; background: #f6f8fa; padding: 12px; }}
    img {{ max-width: 100%; border: 1px solid #d0d7de; }}
  </style>
</head>
<body>
  <h1>{html.escape(event['gene_name'])} isoform switch</h1>
  <div class="muted">{html.escape(event['event_id'])} | contrast {html.escape(event['contrast_id'])}</div>
  <p>Switch-in isoform: <strong>{html.escape(event['switch_in_isoform'])}</strong>;
     switch-out isoform: <strong>{html.escape(event['switch_out_isoform'])}</strong>.</p>
  <img src="{html.escape(relative(str(svg_path), out_path))}" alt="Isoform switch plot">
  <h2>Candidate Isoforms</h2>
  {table(['isoform_id', 'switch_role', 'dIF', 'switch_statistic', 'switch_statistic_name', 'candidate_status', 'transcript_plot_group', 'consequence_summary'], event_candidates)}
  <h2>Functional Annotations</h2>
  {table(['isoform_id', 'source', 'feature_type', 'feature_id', 'feature_name', 'start_aa', 'end_aa', 'score', 'description'], event_annotations)}
  <h2>Sequences</h2>
  {''.join(sequence_blocks) if sequence_blocks else '<p>No sequence rows available.</p>'}
</body>
</html>
"""
    out_path.write_text(html_text, encoding="utf-8")


def render_project_html(event_rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    body = []
    for row in event_rows:
        link = ""
        if row.get("event_html") and Path(row["event_html"]).exists():
            link = f'<a href="{html.escape(relative(row["event_html"], output))}">{html.escape(row["event_id"])}</a>'
        else:
            link = html.escape(row["event_id"])
        body.append(
            "<tr>"
            f"<td>{link}</td>"
            f"<td>{html.escape(row['contrast_id'])}</td>"
            f"<td>{html.escape(row['gene_name'])}</td>"
            f"<td>{html.escape(row['gene_id'])}</td>"
            f"<td>{html.escape(row['switch_in_isoform'])}</td>"
            f"<td>{html.escape(row['switch_out_isoform'])}</td>"
            f"<td>{html.escape(row['max_abs_dIF'])}</td>"
            f"<td>{html.escape(row['best_switch_statistic'])}</td>"
            f"<td>{html.escape(row['n_functional_annotations'])}</td>"
            "</tr>"
        )
    if not body:
        body.append('<tr><td colspan="9">No significant isoform-switch events passed the configured thresholds.</td></tr>')
    output.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Isoform-switch report</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1440px; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
  </style>
</head>
<body>
  <h1>Isoform-switch report</h1>
  <p>Events are ranked by absolute isoform fraction change and linked to exon/annotation/sequence pages.</p>
  <table>
    <thead>
      <tr>
        <th>event</th><th>contrast</th><th>gene</th><th>gene_id</th>
        <th>switch-in</th><th>switch-out</th><th>max abs dIF</th>
        <th>best statistic</th><th>annotations</th>
      </tr>
    </thead>
    <tbody>
      {''.join(body)}
    </tbody>
  </table>
</body>
</html>
""",
        encoding="utf-8",
    )


def write_done(path: Path, events: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tevents_total\tevents_ok\n")
        handle.write(f"ok\t{len(events)}\t{sum(1 for row in events if row.get('status') == 'ok')}\n")


def main() -> int:
    args = parse_args()
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")
    _, manifest_rows = read_table(Path(args.manifest), MANIFEST_REQUIRED)
    metadata = metadata_by_transcript(Path(args.transcript_metadata))
    gtf_models = parse_gtf(Path(args.annotated_gtf))

    candidate_rows, event_rows, event_context, nt_sequences, aa_sequences = build_events(
        manifest_rows,
        metadata,
        args.padj,
        args.dif,
        args.top_n,
    )
    attach_consequences(candidate_rows, event_rows, manifest_rows)
    sequence_rows = build_sequence_rows(candidate_rows, nt_sequences, aa_sequences)

    event_by_isoform: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidate_rows:
        event_by_isoform[row["isoform_id"]].append(
            {
                "event_id": row["event_id"],
                "contrast_id": row["contrast_id"],
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
            }
        )
    annotation_rows = normalized_annotations(
        args.functional_annotation_tables,
        {row["isoform_id"] for row in candidate_rows},
        event_by_isoform,
    )
    annotation_counts = defaultdict(int)
    for row in annotation_rows:
        annotation_counts[row["event_id"]] += 1

    outdir = Path(args.outdir)
    plot_rows = []
    for event in event_rows:
        event_dir = outdir / "events" / safe_id(event["event_id"])
        svg_path = event_dir / "switch.svg"
        html_path = event_dir / "index.html"
        event_annotations = [row for row in annotation_rows if row["event_id"] == event["event_id"]]
        render_event_svg(event, event_context[event["event_id"]], gtf_models, event_annotations, svg_path)
        render_event_html(event, candidate_rows, sequence_rows, annotation_rows, html_path, svg_path)
        event["plot_svg"] = str(svg_path)
        event["event_html"] = str(html_path)
        event["n_functional_annotations"] = str(annotation_counts[event["event_id"]])
        plot_rows.append(
            {
                "event_id": event["event_id"],
                "contrast_id": event["contrast_id"],
                "gene_id": event["gene_id"],
                "gene_name": event["gene_name"],
                "status": event["status"],
                "reason": event["reason"],
                "plot_svg": str(svg_path),
                "event_html": str(html_path),
                "n_isoforms": event["n_isoforms_in_gene"],
                "n_candidate_isoforms": event["n_candidate_isoforms"],
                "nt_fasta": "",
                "aa_fasta": "",
            }
        )

    write_table(Path(args.candidate_table), CANDIDATE_COLUMNS, candidate_rows)
    write_table(Path(args.event_summary), EVENT_COLUMNS, event_rows)
    write_table(Path(args.sequence_table), SEQUENCE_COLUMNS, sequence_rows)
    write_table(Path(args.functional_annotation_table), ANNOTATION_COLUMNS, annotation_rows)
    write_table(Path(args.plot_manifest), PLOT_MANIFEST_COLUMNS, plot_rows)
    render_project_html(event_rows, Path(args.html))
    write_done(Path(args.done), event_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
