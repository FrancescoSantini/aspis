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
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from report_navigation import report_map_css, report_map_item, report_shell_close, report_shell_open


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
    "gene_display",
    "gene_biotype",
    "transcript_biotype",
    "switch_biotype_class",
    "switch_rank",
    "isoform_id",
    "transcript_display",
    "switch_role",
    "dIF",
    "padj_qvalue",
    "switch_statistic",
    "switch_statistic_name",
    "isoform_fraction_control",
    "isoform_fraction_test",
    "switch_direction",
    "novelty_group",
    "reason_selected",
    "candidate_status",
    "transcript_discovery_class",
    "transcript_novelty",
    "transcript_plot_group",
    "gffcompare_class_code",
    "is_stringtie_assembly",
    "assembly_evidence_class",
    "assembly_evidence_label",
    "assembly_evidence_note",
    "consequence_summary",
    "source_detailed",
]
EVENT_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "gene_display",
    "gene_biotype",
    "switch_biotype_class",
    "switch_interpretation_label",
    "coding_priority_rank",
    "coding_priority_score",
    "coding_priority_tier",
    "coding_priority_reasons",
    "switch_rank",
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
    "event_nt_fasta",
    "event_aa_fasta",
    "plot_svg",
    "event_html",
]
NCRNA_SWITCH_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "transcript_biotype",
    "switch_biotype_class",
    "isoform_id",
    "paired_isoform_id",
    "switch_role",
    "dIF",
    "padj_qvalue",
    "transcript_length_change",
    "exon_gain_loss",
    "intron_retention_change",
    "gained_splice_junctions",
    "lost_splice_junctions",
    "n_gained_splice_junctions",
    "n_lost_splice_junctions",
    "TSS_change",
    "TES_change",
    "promoter_context_change",
    "isoform_proximal_gene_context",
    "paired_isoform_proximal_gene_context",
    "antisense_overlap",
    "conserved_exon_change",
    "motif_change",
    "host_smallrna_change",
    "resource_antisense_overlap",
    "ncrna_resource_annotations",
    "pseudogene_caution",
    "coding_potential_change",
    "interpretation_label",
]
CODING_SWITCH_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "switch_rank",
    "coding_priority_rank",
    "coding_priority_score",
    "coding_priority_tier",
    "coding_priority_reasons",
    "functional_consequence_summary",
    "switch_in_isoform",
    "switch_out_isoform",
    "switch_in_dIF",
    "switch_out_dIF",
    "max_abs_dIF",
    "best_switch_statistic",
    "nmd_change",
    "coding_potential_change",
    "orf_length_change_aa",
    "gained_domain",
    "lost_domain",
    "gained_signal_peptide",
    "lost_signal_peptide",
    "gained_transmembrane_region",
    "lost_transmembrane_region",
    "gained_disorder_region",
    "lost_disorder_region",
    "localization_change",
    "n_functional_annotations",
    "event_html",
]
SEQUENCE_COLUMNS = [
    "switch_pair_id",
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "isoform_id",
    "paired_isoform_id",
    "switch_role",
    "pair_role",
    "orf_length_aa",
    "orf_start",
    "orf_end",
    "cds_coordinates",
    "premature_stop",
    "nmd_status",
    "coding_potential",
    "gained_exon_coordinates",
    "lost_exon_coordinates",
    "gained_aa_interval",
    "lost_aa_interval",
    "nt_length",
    "aa_length",
    "nt_sequence",
    "aa_sequence",
    "affected_nt_sequence",
    "affected_aa_sequence",
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
    "feature_change",
    "status",
    "match_type",
    "reason",
    "source_path",
]
ANNOTATION_QA_COLUMNS = [
    "source_path",
    "source_name",
    "source_kind",
    "status",
    "reason",
    "input_rows",
    "parsed_rows",
    "rows_with_identifier",
    "matched_rows",
    "unmatched_rows",
    "duplicate_identifier_rows",
    "unsupported_column_count",
    "unsupported_columns",
    "mapping_loss_fraction",
    "output_rows",
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
    "plots_pdf",
]
EXTERNAL_TOOL_COLUMNS = [
    "tool_group",
    "tool_name",
    "status",
    "returncode",
    "command",
    "produced_files",
    "parser_status",
    "parser_detail",
    "stdout_log",
    "stderr_log",
    "detail",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Isoform-switch manifest TSV")
    parser.add_argument("--transcript-metadata", required=True, help="Transcript metadata TSV")
    parser.add_argument("--annotated-gtf", required=True, help="Annotated transcriptome GTF")
    parser.add_argument("--outdir", required=True, help="Report output directory")
    parser.add_argument("--candidate-table", required=True, help="Output candidate isoform TSV")
    parser.add_argument("--event-summary", required=True, help="Output event summary TSV")
    parser.add_argument("--ncrna-switch-table", required=True, help="Output ncRNA-aware switch interpretation TSV")
    parser.add_argument("--coding-switch-summary", required=True, help="Output coding-switch consequence-priority TSV")
    parser.add_argument("--sequence-table", required=True, help="Output event sequence TSV")
    parser.add_argument("--functional-annotation-table", required=True, help="Output normalized annotation TSV")
    parser.add_argument("--functional-annotation-qa", required=True, help="Output per-source annotation import QA TSV")
    parser.add_argument("--plot-manifest", required=True, help="Output plot manifest TSV")
    parser.add_argument("--external-tool-manifest", required=True, help="Output optional external-tool manifest TSV")
    parser.add_argument("--plots-pdf", required=True, help="Output multi-page switch plot summary PDF")
    parser.add_argument("--html", required=True, help="Output project-level HTML report")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--padj", type=float, default=0.1, help="Maximum switch q-value/p-value")
    parser.add_argument("--dif", type=float, default=0.1, help="Minimum absolute dIF")
    parser.add_argument("--top-n", type=int, default=30, help="Maximum switch events to render")
    parser.add_argument(
        "--functional-annotation-tables",
        default="",
        help=(
            "Optional comma-separated annotation files. Supports generic TSVs "
            "with isoform_id/transcript_id/protein_id plus native InterProScan "
            "TSV, HMMER/Pfam domtblout, CPAT/CPC2 TSV, SignalP summary, "
            "DeepTMHMM/TMHMM GFF-like output, DeepLoc2 TSV, and IUPred2A table/raw output."
        ),
    )
    parser.add_argument(
        "--ncrna-annotation-tables",
        default="",
        help=(
            "Optional comma-separated TSVs with resource-backed ncRNA annotations. "
            "Rows can match by transcript_id/isoform_id, gene_id, or chrom/start/end overlap "
            "with gained/lost switch intervals. Recognized feature_type values include "
            "conserved_exon, rbp_motif, mirna_motif, host_small_rna, and antisense_overlap."
        ),
    )
    parser.add_argument("--interproscan-command", default="", help="Optional command template for InterProScan")
    parser.add_argument("--pfam-command", default="", help="Optional command template for Pfam/hmmscan")
    parser.add_argument("--coding-potential-command", default="", help="Optional command template for CPAT/CPC2")
    parser.add_argument("--signalp-command", default="", help="Optional command template for SignalP")
    parser.add_argument("--tm-command", default="", help="Optional command template for DeepTMHMM/TMHMM")
    parser.add_argument("--localization-command", default="", help="Optional command template for DeepLoc2/localization")
    parser.add_argument("--disorder-command", default="", help="Optional command template for IUPred2A/NetSurfP")
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


def normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower().lstrip("#")).strip("_")


def value_by_normalized_key(row: dict[str, str], names: Iterable[str]) -> str:
    lookup = {normalized_key(key): value for key, value in row.items()}
    for name in names:
        value = lookup.get(normalized_key(name), "")
        if value != "":
            return value
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


def annotation_row(
    *,
    isoform_id: str,
    source: str,
    feature_type: str,
    feature_id: str = "",
    feature_name: str = "",
    start_aa: str = "",
    end_aa: str = "",
    score: str = "",
    description: str = "",
) -> dict[str, str]:
    return {
        "isoform_id": isoform_id,
        "source": source,
        "feature_type": feature_type,
        "feature_id": feature_id,
        "feature_name": feature_name,
        "start_aa": start_aa,
        "end_aa": end_aa,
        "score": score,
        "description": description,
    }


def parse_signalp_cleavage(value: str) -> tuple[str, str]:
    if not value:
        return "", ""
    match = re.search(r"CS\s+pos:\s*(\d+)\s*-\s*(\d+)", value, flags=re.IGNORECASE)
    if match:
        return "1", match.group(1)
    match = re.search(r"\b(\d+)\s*-\s*(\d+)\b", value)
    if match:
        return "1", match.group(1)
    return "", ""


def contiguous_intervals(positions: list[int]) -> list[tuple[int, int]]:
    if not positions:
        return []
    sorted_positions = sorted(set(positions))
    intervals = []
    start = previous = sorted_positions[0]
    for position in sorted_positions[1:]:
        if position == previous + 1:
            previous = position
            continue
        intervals.append((start, previous))
        start = previous = position
    intervals.append((start, previous))
    return intervals


def parse_interproscan_headerless(path: Path, lines: list[str]) -> list[dict[str, str]]:
    rows = []
    parsed_any = False
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 11:
            return []
        if not (parts[6].isdigit() and parts[7].isdigit()):
            return []
        parsed_any = True
        interpro_accession = parts[11] if len(parts) > 11 and parts[11] != "-" else ""
        interpro_description = parts[12] if len(parts) > 12 and parts[12] != "-" else ""
        signature_accession = parts[4]
        signature_description = parts[5]
        analysis = parts[3]
        rows.append(
            annotation_row(
                isoform_id=parts[0],
                source=f"interproscan:{analysis}",
                feature_type="protein_domain",
                feature_id=interpro_accession or signature_accession,
                feature_name=interpro_description or signature_description,
                start_aa=parts[6],
                end_aa=parts[7],
                score=parts[8],
                description=f"signature={signature_accession}; signature_description={signature_description}",
            )
        )
    return rows if parsed_any else []


def parse_hmmer_domtblout(path: Path, lines: list[str]) -> list[dict[str, str]]:
    rows = []
    parsed_any = False
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.split(maxsplit=22)
        if len(parts) < 22:
            continue
        if not (parts[17].isdigit() and parts[18].isdigit()):
            continue
        parsed_any = True
        target_name = parts[0]
        target_accession = "" if parts[1] == "-" else parts[1]
        description = parts[22] if len(parts) > 22 else target_name
        rows.append(
            annotation_row(
                isoform_id=parts[3],
                source="hmmer_domtblout",
                feature_type="protein_domain",
                feature_id=target_accession or target_name,
                feature_name=target_name,
                start_aa=parts[17],
                end_aa=parts[18],
                score=parts[12],
                description=description,
            )
        )
    return rows if parsed_any else []


def parse_deeptmhmm_gff(path: Path, lines: list[str]) -> list[dict[str, str]]:
    rows = []
    parsed_any = False
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        feature = parts[2]
        if feature.lower() not in {"tmhelix", "transmembrane", "inside", "outside", "signal", "signal_peptide"}:
            continue
        parsed_any = True
        feature_type = "transmembrane" if feature.lower() in {"tmhelix", "transmembrane"} else "topology"
        rows.append(
            annotation_row(
                isoform_id=parts[0],
                source=parts[1] or "deeptmhmm",
                feature_type=feature_type,
                feature_id=feature,
                feature_name=feature,
                start_aa=parts[3],
                end_aa=parts[4],
                score="" if parts[5] == "." else parts[5],
                description=parts[8],
            )
        )
    return rows if parsed_any else []


def parse_iupred_raw(path: Path, lines: list[str]) -> list[dict[str, str]]:
    rows = []
    current_id = ""
    positions_by_id: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for line in lines:
        if line.startswith("#"):
            match = re.search(r"['\"]([^'\"]+)['\"]", line)
            if match:
                current_id = match.group(1)
            continue
        parts = line.split()
        if len(parts) < 3 or not parts[0].isdigit() or not current_id:
            continue
        score = to_float(parts[2])
        if score is not None:
            positions_by_id[current_id].append((int(parts[0]), score))
    for isoform_id, values in positions_by_id.items():
        disordered = [position for position, score in values if score >= 0.5]
        for start, end in contiguous_intervals(disordered):
            max_score = max(score for position, score in values if start <= position <= end)
            rows.append(
                annotation_row(
                    isoform_id=isoform_id,
                    source="iupred2a",
                    feature_type="disorder",
                    feature_id="iupred2a_disorder",
                    feature_name="predicted_disordered_region",
                    start_aa=str(start),
                    end_aa=str(end),
                    score=f"{max_score:.6g}",
                    description="IUPred2A score >= 0.5",
                )
            )
    return rows


def parse_structured_annotation_table(path: Path) -> list[dict[str, str]]:
    fieldnames, raw_rows = read_table(path)
    normalized_fields = {normalized_key(field) for field in fieldnames}
    rows = []
    is_interpro = {"protein_accession", "analysis", "signature_accession", "start", "stop"} <= normalized_fields
    is_cpat = bool({"coding_prob", "coding_probability"} & normalized_fields and {"coding_label", "label"} & normalized_fields)
    is_signalp = "prediction" in normalized_fields and bool({"cs_position", "cs_pos"} & normalized_fields)
    is_deeploc = bool({"localizations", "localization"} & normalized_fields and {"protein_id", "protein_id_seq"} & normalized_fields)
    is_iupred = "position" in normalized_fields and bool({"iupred2", "iupred_score", "score"} & normalized_fields)
    if is_iupred:
        grouped: dict[str, list[tuple[int, float]]] = defaultdict(list)
        for raw in raw_rows:
            isoform_id = value_by_normalized_key(raw, ["protein_id", "isoform_id", "transcript_id", "id", "query"])
            position = value_by_normalized_key(raw, ["position", "pos", "aa_position"])
            score = to_float(value_by_normalized_key(raw, ["iupred2", "iupred_score", "disorder_score", "score"]))
            if isoform_id and position.isdigit() and score is not None:
                grouped[isoform_id].append((int(position), score))
        for isoform_id, values in grouped.items():
            disordered = [position for position, score in values if score >= 0.5]
            for start, end in contiguous_intervals(disordered):
                max_score = max(score for position, score in values if start <= position <= end)
                rows.append(
                    annotation_row(
                        isoform_id=isoform_id,
                        source="iupred2a",
                        feature_type="disorder",
                        feature_id="iupred2a_disorder",
                        feature_name="predicted_disordered_region",
                        start_aa=str(start),
                        end_aa=str(end),
                        score=f"{max_score:.6g}",
                        description="IUPred2A score >= 0.5",
                    )
                )
        return rows
    for raw in raw_rows:
        if is_interpro:
            analysis = value_by_normalized_key(raw, ["analysis"])
            interpro_accession = value_by_normalized_key(raw, ["interpro_accession"])
            signature_accession = value_by_normalized_key(raw, ["signature_accession"])
            interpro_description = value_by_normalized_key(raw, ["interpro_description"])
            signature_description = value_by_normalized_key(raw, ["signature_description"])
            rows.append(
                annotation_row(
                    isoform_id=value_by_normalized_key(raw, ["protein_accession", "isoform_id", "query"]),
                    source=f"interproscan:{analysis}",
                    feature_type="protein_domain",
                    feature_id=interpro_accession or signature_accession,
                    feature_name=interpro_description or signature_description,
                    start_aa=value_by_normalized_key(raw, ["start"]),
                    end_aa=value_by_normalized_key(raw, ["stop", "end"]),
                    score=value_by_normalized_key(raw, ["score", "evalue"]),
                    description=f"signature={signature_accession}; signature_description={signature_description}",
                )
            )
        elif is_cpat:
            label = value_by_normalized_key(raw, ["coding_label", "label", "prediction"])
            rows.append(
                annotation_row(
                    isoform_id=value_by_normalized_key(raw, ["id", "mrna", "seq_id", "transcript_id", "isoform_id", "query"]),
                    source="coding_potential",
                    feature_type="coding_potential",
                    feature_id=label,
                    feature_name=label,
                    score=value_by_normalized_key(raw, ["coding_prob", "coding_probability"]),
                    description="CPAT/CPC2 coding-potential prediction",
                )
            )
        elif is_signalp:
            prediction = value_by_normalized_key(raw, ["prediction"])
            cleavage = value_by_normalized_key(raw, ["cs_position", "cs_pos", "cleavage_site"])
            start, end = parse_signalp_cleavage(cleavage) if prediction and prediction.upper() != "OTHER" else ("", "")
            rows.append(
                annotation_row(
                    isoform_id=value_by_normalized_key(raw, ["id", "protein_id", "isoform_id", "transcript_id"]),
                    source="signalp",
                    feature_type="signal_peptide",
                    feature_id=prediction,
                    feature_name=prediction,
                    start_aa=start,
                    end_aa=end,
                    score=value_by_normalized_key(raw, ["sp_sec_spi", "signal_peptide_probability", "probability", "score"]),
                    description=cleavage,
                )
            )
        elif is_deeploc:
            localization = value_by_normalized_key(raw, ["localizations", "localization"])
            rows.append(
                annotation_row(
                    isoform_id=value_by_normalized_key(raw, ["protein_id", "id", "isoform_id", "transcript_id"]),
                    source="deeploc2",
                    feature_type="localization",
                    feature_id=localization,
                    feature_name=localization,
                    score=value_by_normalized_key(raw, ["score", "probability"]),
                    description=value_by_normalized_key(raw, ["signals", "signal"]),
                )
            )
        else:
            rows.append(raw)
    return rows


def native_annotation_rows(path: Path) -> list[dict[str, str]]:
    lines = [line.rstrip("\n") for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    if not lines:
        return []
    for parser in [parse_deeptmhmm_gff, parse_hmmer_domtblout, parse_interproscan_headerless, parse_iupred_raw]:
        parsed = parser(path, lines)
        if parsed:
            return parsed
    return parse_structured_annotation_table(path)


ANNOTATION_ID_FIELDS = ["isoform_id", "transcript_id", "protein_id", "query", "query_id", "id"]
ANNOTATION_GENE_FIELDS = ["gene_id", "gene", "geneid", "gene_name", "gene_symbol"]
SUPPORTED_GENERIC_ANNOTATION_FIELDS = {
    "analysis",
    "accession",
    "aa_position",
    "ali_end",
    "ali_start",
    "coding_label",
    "coding_potential",
    "coding_prob",
    "coding_probability",
    "cs_pos",
    "cs_position",
    "database",
    "description",
    "desc",
    "disorder_score",
    "e_value",
    "end",
    "end_aa",
    "feature_id",
    "feature_name",
    "feature_type",
    "from",
    "gene",
    "gene_id",
    "gene_name",
    "gene_symbol",
    "geneid",
    "id",
    "interpro_accession",
    "interpro_description",
    "isoform_id",
    "iupred2",
    "iupred_score",
    "label",
    "localization",
    "localizations",
    "name",
    "prediction",
    "probability",
    "protein_accession",
    "protein_id",
    "query",
    "query_id",
    "score",
    "seq_id",
    "signature_accession",
    "signature_desc",
    "signature_description",
    "signals",
    "source",
    "sp_sec_spi",
    "start",
    "start_aa",
    "stop",
    "to",
    "transcript_id",
    "type",
}


def nonempty_source_lines(path: Path) -> list[str]:
    return [
        line.rstrip("\n")
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    ]


def looks_like_annotation_header(line: str) -> bool:
    fields = [normalized_key(part) for part in line.lstrip("#").split("\t")]
    if len(fields) < 2:
        return False
    recognizable = set(fields) & SUPPORTED_GENERIC_ANNOTATION_FIELDS
    return bool(recognizable & {normalized_key(field) for field in ANNOTATION_ID_FIELDS + ANNOTATION_GENE_FIELDS})


def source_input_row_count(path: Path) -> int:
    lines = [line for line in nonempty_source_lines(path) if not line.startswith("#") or looks_like_annotation_header(line)]
    if not lines:
        return 0
    if looks_like_annotation_header(lines[0]):
        return max(0, len(lines) - 1)
    return len([line for line in lines if not line.startswith("#")])


def source_header_columns(path: Path) -> list[str]:
    for line in nonempty_source_lines(path):
        if looks_like_annotation_header(line):
            return [part.strip().lstrip("#").strip() for part in line.split("\t")]
        if line.startswith("#"):
            continue
        break
    return []


def unsupported_annotation_columns(path: Path) -> list[str]:
    unsupported = []
    for column in source_header_columns(path):
        if normalized_key(column) not in SUPPORTED_GENERIC_ANNOTATION_FIELDS:
            unsupported.append(column)
    return unsupported


def annotation_source_kind(path: Path, rows: list[dict[str, str]]) -> str:
    sources = {first_existing(row, ["source", "analysis", "database"]).lower() for row in rows[:20]}
    feature_types = {first_existing(row, ["feature_type", "type"]).lower() for row in rows[:20]}
    if any(source.startswith("interproscan:") for source in sources):
        return "interproscan_tsv"
    if "hmmer_domtblout" in sources or path.suffix.lower() == ".domtblout":
        return "hmmer_pfam_domtblout"
    if "coding_potential" in sources or "coding_potential" in feature_types:
        return "coding_potential_tsv"
    if "signalp" in sources or "signal_peptide" in feature_types:
        return "signalp_table"
    if any(source in {"deeptmhmm", "tmhmm"} for source in sources) or "transmembrane" in feature_types:
        return "tm_topology_table"
    if "deeploc2" in sources or "localization" in feature_types:
        return "localization_table"
    if "iupred2a" in sources or "disorder" in feature_types:
        return "disorder_table"
    return "generic_tsv"


def annotation_aliases(value: str) -> set[str]:
    if not value:
        return set()
    aliases = {value, value.split("|")[0], value.rsplit("|", 1)[-1]}
    aliases.add(value.split(".")[0])
    return {alias for alias in aliases if alias}


def annotation_qa_row(
    *,
    path: Path | None,
    source_name: str,
    source_kind: str,
    status: str,
    reason: str,
    input_rows: int = 0,
    parsed_rows: int = 0,
    rows_with_identifier: int = 0,
    matched_rows: int = 0,
    duplicate_identifier_rows: int = 0,
    unsupported_columns: list[str] | None = None,
    output_rows: int = 0,
) -> dict[str, str]:
    unmatched_rows = max(0, parsed_rows - matched_rows)
    mapping_loss = (unmatched_rows / parsed_rows) if parsed_rows else 0.0
    unsupported = unsupported_columns or []
    return {
        "source_path": "" if path is None else str(path),
        "source_name": source_name,
        "source_kind": source_kind,
        "status": status,
        "reason": reason,
        "input_rows": str(input_rows),
        "parsed_rows": str(parsed_rows),
        "rows_with_identifier": str(rows_with_identifier),
        "matched_rows": str(matched_rows),
        "unmatched_rows": str(unmatched_rows),
        "duplicate_identifier_rows": str(duplicate_identifier_rows),
        "unsupported_column_count": str(len(unsupported)),
        "unsupported_columns": ";".join(unsupported),
        "mapping_loss_fraction": f"{mapping_loss:.6g}",
        "output_rows": str(output_rows),
    }


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
    gene_biotype: str
    transcript_biotype: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]
    cds: list[tuple[int, int]]


@dataclass
class GeneLocus:
    gene_id: str
    gene_name: str
    gene_biotype: str
    chrom: str
    strand: str
    start: int
    end: int


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
            if feature not in {"exon", "CDS"}:
                continue
            attrs = parse_gtf_attributes(attributes)
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id)
            gene_biotype = first_existing(
                attrs,
                ["gene_biotype", "gene_type", "gene_type_name", "biotype"],
            )
            transcript_biotype = first_existing(
                attrs,
                ["transcript_biotype", "transcript_type", "transcript_type_name", "biotype"],
            )
            model = models.setdefault(
                transcript_id,
                TranscriptModel(transcript_id, gene_id, gene_name, gene_biotype, transcript_biotype, chrom, strand, [], []),
            )
            if gene_biotype and not model.gene_biotype:
                model.gene_biotype = gene_biotype
            if transcript_biotype and not model.transcript_biotype:
                model.transcript_biotype = transcript_biotype
            if feature == "exon":
                model.exons.append((int(start_text), int(end_text)))
            else:
                model.cds.append((int(start_text), int(end_text)))
    for model in models.values():
        model.exons.sort()
        model.cds.sort()
    return models


def metadata_by_transcript(path: Path) -> dict[str, dict[str, str]]:
    _, rows = read_table(path, TRANSCRIPT_METADATA_REQUIRED)
    indexed = {}
    for row in rows:
        transcript_id = row.get("transcript_id", "")
        if transcript_id:
            indexed[transcript_id] = row
    return indexed


def looks_like_stringtie_id(value: str) -> bool:
    return bool(value) and value.startswith(("MSTRG", "STRG"))


def display_gene(gene_id: str, gene_name: str) -> str:
    gene_id = (gene_id or "").strip()
    gene_name = (gene_name or "").strip()
    if gene_name and gene_id and gene_name != gene_id:
        return f"{gene_name} ({gene_id})"
    return gene_name or gene_id


def display_transcript(transcript_id: str, gene_id: str, gene_name: str) -> str:
    transcript_id = (transcript_id or "").strip()
    gene_label = display_gene(gene_id, gene_name)
    if gene_label and transcript_id:
        return f"{gene_label} | {transcript_id}"
    return transcript_id or gene_label


def fallback_assembly_evidence(gene_id: str, transcript_id: str, discovery_class: str) -> tuple[str, str, str, str]:
    is_stringtie = "yes" if looks_like_stringtie_id(gene_id) or looks_like_stringtie_id(transcript_id) else "no"
    if discovery_class in {"known_transcript", "reference_contained_or_containing", "unclassified_reference_compatible"}:
        return (
            is_stringtie,
            "reference_compatible",
            "Reference-compatible transcript",
            "Annotated or reference-compatible transcript model; not a novel-assembly claim.",
        )
    if discovery_class == "novel_isoform_known_gene":
        return (
            is_stringtie,
            "candidate_novel_isoform",
            "Candidate novel isoform",
            "RNA-seq assembly supports a candidate novel isoform of a known gene; independent validation is recommended before treating it as a validated transcript.",
        )
    if discovery_class in {"intergenic_novel_locus", "unclassified_novel_candidate"}:
        return (
            is_stringtie,
            "candidate_novel_locus",
            "Candidate novel locus",
            "RNA-seq assembly supports an intergenic transcript model; independent validation is required before naming a new gene or transcript.",
        )
    if discovery_class == "likely_artifact_or_repeat":
        return (
            is_stringtie,
            "low_confidence_assembly",
            "Low-confidence assembly model",
            "Possible artifact, repeat, pre-mRNA, or run-on signal; do not treat as a validated transcript without external evidence.",
        )
    return (
        is_stringtie,
        "ambiguous_assembly_model",
        "Ambiguous assembled model",
        "Assembly overlaps annotation ambiguously; review splice structure, strand, expression, and independent evidence before interpretation.",
    )


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
    event_by_gene: dict[str, list[dict[str, str]]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    qa_rows: list[dict[str, str]] = []
    paths = [part.strip() for part in paths_text.split(",") if part.strip()]
    if not paths:
        qa_rows.append(
            annotation_qa_row(
                path=None,
                source_name="functional_annotation_tables",
                source_kind="precomputed_tables",
                status="not_configured",
                reason="No functional annotation table or parseable command output configured",
            )
        )
        return rows, qa_rows

    for item in paths:
        path = Path(item)
        if not path.exists():
            qa_rows.append(
                annotation_qa_row(
                    path=path,
                    source_name=path.name,
                    source_kind="unknown",
                    status="blocked",
                    reason="Functional annotation table does not exist",
                )
            )
            continue
        input_rows = source_input_row_count(path)
        unsupported_columns = unsupported_annotation_columns(path)
        try:
            raw_rows = native_annotation_rows(path)
        except Exception as exc:
            qa_rows.append(
                annotation_qa_row(
                    path=path,
                    source_name=path.name,
                    source_kind="unknown",
                    status="failed",
                    reason=f"Could not parse annotation table: {exc}",
                    input_rows=input_rows,
                    unsupported_columns=unsupported_columns,
                )
            )
            continue
        source_kind = annotation_source_kind(path, raw_rows)
        if not raw_rows:
            qa_rows.append(
                annotation_qa_row(
                    path=path,
                    source_name=path.name,
                    source_kind=source_kind,
                    status="blocked",
                    reason="No parseable annotation rows found",
                    input_rows=input_rows,
                    unsupported_columns=unsupported_columns,
                )
            )
            continue
        identifiers = [
            first_existing(raw, ANNOTATION_ID_FIELDS) or first_existing(raw, ANNOTATION_GENE_FIELDS)
            for raw in raw_rows
        ]
        rows_with_identifier = sum(1 for value in identifiers if value)
        duplicate_identifier_rows = sum(count - 1 for count in Counter(value for value in identifiers if value).values() if count > 1)
        if not rows_with_identifier:
            qa_rows.append(
                annotation_qa_row(
                    path=path,
                    source_name=path.name,
                    source_kind=source_kind,
                    status="blocked",
                    reason="Annotation table lacks a recognized isoform/transcript/protein/gene identifier",
                    input_rows=input_rows,
                    parsed_rows=len(raw_rows),
                    rows_with_identifier=rows_with_identifier,
                    duplicate_identifier_rows=duplicate_identifier_rows,
                    unsupported_columns=unsupported_columns,
                )
            )
            continue
        source_output_rows = 0
        source_matched_rows = 0
        for raw in raw_rows:
            isoform_id = first_existing(raw, ANNOTATION_ID_FIELDS)
            gene_id = first_existing(raw, ANNOTATION_GENE_FIELDS)
            if not isoform_id and not gene_id:
                continue
            matched_isoforms = sorted(alias for alias in annotation_aliases(isoform_id) if alias in candidate_isoforms)
            matched_events = []
            match_type = ""
            if matched_isoforms:
                source_matched_rows += 1
                match_type = "isoform_id"
                for matched_id in matched_isoforms:
                    for event in event_by_isoform.get(matched_id, []):
                        matched_events.append((matched_id, event, "isoform_id"))
            elif gene_id and gene_id in event_by_gene:
                source_matched_rows += 1
                match_type = "gene_id"
                for event in event_by_gene.get(gene_id, []):
                    matched_events.append((isoform_id, event, "gene_id"))
            if not matched_events:
                continue
            source = first_existing(raw, ["source", "analysis", "database", "signature_desc"])
            feature_id = first_existing(raw, ["feature_id", "accession", "signature_accession", "interpro_accession", "id"])
            feature_name = first_existing(raw, ["feature_name", "name", "signature_description", "interpro_description"])
            feature_type = first_existing(raw, ["feature_type", "type", "analysis", "database"]) or "protein_feature"
            start_aa = first_existing(raw, ["start_aa", "start", "from", "ali_start", "hmm_start"])
            end_aa = first_existing(raw, ["end_aa", "end", "to", "ali_end", "hmm_end"])
            score = first_existing(raw, ["score", "evalue", "e_value", "bitscore"])
            description = first_existing(raw, ["description", "desc", "interpro_description", "signature_description"])
            for matched_id, event, row_match_type in matched_events:
                source_output_rows += 1
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
                        "feature_change": "",
                        "status": "ok",
                        "match_type": row_match_type,
                        "reason": f"matched by {row_match_type}",
                        "source_path": str(path),
                    }
                )
        if source_matched_rows:
            status = "ok"
            reason = "One or more annotation rows matched rendered isoform-switch candidates"
        else:
            status = "ok_no_matches"
            reason = "Annotation table parsed, but no rows matched rendered isoform-switch candidates"
        qa_rows.append(
            annotation_qa_row(
                path=path,
                source_name=path.name,
                source_kind=source_kind,
                status=status,
                reason=reason,
                input_rows=input_rows,
                parsed_rows=len(raw_rows),
                rows_with_identifier=rows_with_identifier,
                matched_rows=source_matched_rows,
                duplicate_identifier_rows=duplicate_identifier_rows,
                unsupported_columns=unsupported_columns,
                output_rows=source_output_rows,
            )
        )
    rows.sort(key=lambda row: (row["event_id"], row["isoform_id"], row["source"], row["start_aa"]))
    return rows, qa_rows


def intervals_overlap_text(feature_start: str, feature_end: str, interval_text_value: str) -> bool:
    start = to_float(feature_start)
    end = to_float(feature_end)
    if start is None or end is None or not interval_text_value:
        return False
    for part in interval_text_value.split(","):
        if "-" not in part:
            continue
        left, right = part.split("-", 1)
        left_value = to_float(left)
        right_value = to_float(right)
        if left_value is None or right_value is None:
            continue
        if not (end < left_value or start > right_value):
            return True
    return False


def annotate_feature_changes(annotation_rows: list[dict[str, str]], sequence_rows: list[dict[str, str]]) -> None:
    sequence_by_key = {(row["event_id"], row["isoform_id"]): row for row in sequence_rows}
    for row in annotation_rows:
        sequence = sequence_by_key.get((row["event_id"], row["isoform_id"]), {})
        if not sequence:
            row["feature_change"] = "unclassified"
            continue
        role = sequence.get("pair_role", "")
        if role == "gained_isoform" and intervals_overlap_text(
            row.get("start_aa", ""),
            row.get("end_aa", ""),
            sequence.get("gained_aa_interval", ""),
        ):
            row["feature_change"] = "gained_domain_overlap"
        elif role == "lost_isoform" and intervals_overlap_text(
            row.get("start_aa", ""),
            row.get("end_aa", ""),
            sequence.get("lost_aa_interval", ""),
        ):
            row["feature_change"] = "lost_domain_overlap"
        else:
            row["feature_change"] = "retained_or_unmapped"


def role_for_isoform(isoform_id: str, switch_in: str, switch_out: str) -> str:
    if isoform_id == switch_in:
        return "switch_in"
    if isoform_id == switch_out:
        return "switch_out"
    return "same_gene"


def switch_direction(d_if: str) -> str:
    value = to_float(d_if)
    if value is None:
        return "unknown"
    if value > 0:
        return "test_gain"
    if value < 0:
        return "test_loss"
    return "unchanged"


def fraction_value(row: dict[str, str], role: str) -> str:
    if role == "control":
        names = [
            "isoform_fraction_control",
            "isoformFraction_control",
            "isoform_fraction_1",
            "isoformFraction_1",
            "condition_1_isoform_fraction",
            "IF1",
        ]
    else:
        names = [
            "isoform_fraction_test",
            "isoformFraction_test",
            "isoform_fraction_2",
            "isoformFraction_2",
            "condition_2_isoform_fraction",
            "IF2",
        ]
    value = first_existing(row, names)
    parsed = to_float(value)
    if parsed is None:
        return ""
    return f"{parsed:.6g}"


def interval_text(chrom: str, intervals: list[tuple[int, int]]) -> str:
    return ",".join(f"{chrom}:{start}-{end}" for start, end in intervals)


def interval_span_text(intervals: list[tuple[int, int]]) -> str:
    if not intervals:
        return ""
    return ",".join(f"{start}-{end}" for start, end in intervals)


def subtract_intervals(source: list[tuple[int, int]], blockers: list[tuple[int, int]]) -> list[tuple[int, int]]:
    remaining: list[tuple[int, int]] = []
    for start, end in source:
        pieces = [(start, end)]
        for block_start, block_end in blockers:
            next_pieces = []
            for piece_start, piece_end in pieces:
                if block_end < piece_start or block_start > piece_end:
                    next_pieces.append((piece_start, piece_end))
                    continue
                if block_start > piece_start:
                    next_pieces.append((piece_start, block_start - 1))
                if block_end < piece_end:
                    next_pieces.append((block_end + 1, piece_end))
            pieces = next_pieces
        remaining.extend(piece for piece in pieces if piece[0] <= piece[1])
    return remaining


def ordered_intervals(model: TranscriptModel, feature: str = "exon") -> list[tuple[int, int]]:
    intervals = model.cds if feature == "cds" else model.exons
    if model.strand == "-":
        return sorted(intervals, reverse=True)
    return sorted(intervals)


def genomic_to_transcript_intervals(
    model: TranscriptModel,
    genomic_intervals: list[tuple[int, int]],
    feature: str = "exon",
) -> list[tuple[int, int]]:
    projected: list[tuple[int, int]] = []
    cursor = 1
    intervals = ordered_intervals(model, feature)
    for feat_start, feat_end in intervals:
        feat_len = feat_end - feat_start + 1
        for interval_start, interval_end in genomic_intervals:
            overlap_start = max(feat_start, interval_start)
            overlap_end = min(feat_end, interval_end)
            if overlap_start > overlap_end:
                continue
            if model.strand == "-":
                local_start = cursor + (feat_end - overlap_end)
                local_end = cursor + (feat_end - overlap_start)
            else:
                local_start = cursor + (overlap_start - feat_start)
                local_end = cursor + (overlap_end - feat_start)
            projected.append((local_start, local_end))
        cursor += feat_len
    projected.sort()
    return projected


def nt_subsequence(sequence: str, intervals: list[tuple[int, int]]) -> str:
    if not sequence or not intervals:
        return ""
    pieces = []
    for start, end in intervals:
        start_idx = max(0, start - 1)
        end_idx = min(len(sequence), end)
        if start_idx < end_idx:
            pieces.append(sequence[start_idx:end_idx])
    return "".join(pieces)


def aa_intervals_from_cds(model: Optional[TranscriptModel], intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if model is None or not intervals or not model.cds:
        return []
    cds_nt = genomic_to_transcript_intervals(model, intervals, "cds")
    aa_intervals = []
    for start, end in cds_nt:
        aa_start = max(1, math.floor((start - 1) / 3) + 1)
        aa_end = max(aa_start, math.ceil(end / 3))
        aa_intervals.append((aa_start, aa_end))
    return aa_intervals


def aa_subsequence(sequence: str, intervals: list[tuple[int, int]]) -> str:
    if not sequence or not intervals:
        return ""
    return "".join(sequence[max(0, start - 1): min(len(sequence), end)] for start, end in intervals)


def model_unique_context(
    switch_in: str,
    switch_out: str,
    models: dict[str, TranscriptModel],
) -> dict[str, object]:
    in_model = models.get(switch_in)
    out_model = models.get(switch_out)
    if in_model is None or out_model is None:
        return {
            "gained_exons": [],
            "lost_exons": [],
            "gained_aa": [],
            "lost_aa": [],
            "chrom": "",
        }
    gained = subtract_intervals(in_model.exons, out_model.exons)
    lost = subtract_intervals(out_model.exons, in_model.exons)
    return {
        "gained_exons": gained,
        "lost_exons": lost,
        "gained_aa": aa_intervals_from_cds(in_model, gained),
        "lost_aa": aa_intervals_from_cds(out_model, lost),
        "chrom": in_model.chrom or out_model.chrom,
    }


def first_available_feature(row: dict[str, str], names: Iterable[str]) -> str:
    value = first_existing(row, names)
    return value if value else ""


CODING_BIOTYPES = {"protein_coding"}
ARTIFACT_NOVELTY_GROUPS = {"artifact", "ambiguous"}
NONCODING_HINTS = [
    "lncrna",
    "lincrna",
    "antisense",
    "pseudogene",
    "snorna",
    "snrna",
    "rrna",
    "mirna",
    "ncrna",
    "scrna",
    "scarna",
    "ribozyme",
    "retained_intron",
    "processed_transcript",
    "sense_intronic",
    "sense_overlapping",
]


def normalized_biotype(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", (value or "").strip().lower()).strip("_")


def metadata_biotype(meta: dict[str, str], model: Optional[TranscriptModel]) -> tuple[str, str]:
    gene_biotype = first_existing(
        meta,
        ["gene_biotype", "gene_type", "gene_type_name", "biotype"],
    )
    transcript_biotype = first_existing(
        meta,
        ["transcript_biotype", "transcript_type", "transcript_type_name", "biotype"],
    )
    if model is not None:
        gene_biotype = gene_biotype or model.gene_biotype
        transcript_biotype = transcript_biotype or model.transcript_biotype
    return gene_biotype, transcript_biotype


def is_noncoding_biotype(value: str) -> bool:
    normalized = normalized_biotype(value)
    if not normalized or normalized in CODING_BIOTYPES:
        return False
    return any(hint in normalized for hint in NONCODING_HINTS) or normalized.endswith("_rna")


def is_pseudogene_biotype(value: str) -> bool:
    return "pseudogene" in normalized_biotype(value)


def is_pseudogene_switch(row: dict[str, str]) -> bool:
    return is_pseudogene_biotype(row.get("gene_biotype", "")) or is_pseudogene_biotype(row.get("transcript_biotype", ""))


def is_artifact_or_ambiguous(row: dict[str, str]) -> bool:
    values = {
        normalized_biotype(row.get("transcript_plot_group", "")),
        normalized_biotype(row.get("transcript_discovery_class", "")),
        normalized_biotype(row.get("transcript_novelty", "")),
    }
    return bool(values & ARTIFACT_NOVELTY_GROUPS)


def classify_switch_biotype(rows: list[dict[str, str]]) -> tuple[str, str, str]:
    if not rows:
        return "unclassified", "", "unclassified_switch"
    gene_biotype = first_existing(rows[0], ["gene_biotype"])
    if any(is_artifact_or_ambiguous(row) for row in rows):
        return "ambiguous_artifact", gene_biotype, "ambiguous_or_artifact"
    all_biotypes = [
        row.get("gene_biotype", "") for row in rows
    ] + [
        row.get("transcript_biotype", "") for row in rows
    ]
    has_coding = any(normalized_biotype(value) in CODING_BIOTYPES for value in all_biotypes)
    has_noncoding = any(is_noncoding_biotype(value) for value in all_biotypes)
    has_pseudogene = any(is_pseudogene_biotype(value) for value in all_biotypes)
    if has_coding and has_noncoding:
        return "mixed_coding_noncoding", gene_biotype, "coding_potential_transition"
    if has_pseudogene:
        return "noncoding", gene_biotype, "pseudogene_transcript_architecture_change"
    if has_noncoding:
        return "noncoding", gene_biotype, "noncoding_structure_change"
    if has_coding:
        return "coding", gene_biotype, "coding_consequence_candidate"
    return "unclassified", gene_biotype, "unclassified_switch"


def build_events(
    manifest_rows: list[dict[str, str]],
    metadata: dict[str, dict[str, str]],
    models: dict[str, TranscriptModel],
    padj: float,
    dif: float,
    top_n: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, dict[str, object]], dict[str, str], dict[str, str]]:
    candidate_rows: list[dict[str, str]] = []
    event_rows: list[dict[str, str]] = []
    event_context: dict[str, dict[str, object]] = {}
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
    for switch_rank, ((contrast_id, gene_id), candidates) in enumerate(ranked_events[:top_n], start=1):
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
            "gene_biotype": "",
            "switch_biotype_class": "",
            "switch_interpretation_label": "",
            "switch_rank": str(switch_rank),
            "coding_priority_rank": "",
            "coding_priority_score": "0",
            "coding_priority_tier": "",
            "coding_priority_reasons": "",
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
            "event_nt_fasta": "",
            "event_aa_fasta": "",
            "plot_svg": "",
            "event_html": "",
        }
        event_rows.append(event)
        event_context[event_id] = {
            **event,
            "_gene_rows": gene_rows,  # type: ignore[dict-item]
            "_candidate_isoforms": candidate_isoforms,  # type: ignore[dict-item]
            "_unique_context": model_unique_context(switch_in["_isoform_id"], switch_out["_isoform_id"], models),
        }
        current_candidate_rows = []
        for row in gene_rows:
            isoform_id = row["_isoform_id"]
            meta = metadata.get(isoform_id, {})
            model = models.get(isoform_id)
            gene_biotype, transcript_biotype = metadata_biotype(meta, model)
            gene_display = meta.get("gene_display") or display_gene(gene_id, gene_name)
            transcript_display = meta.get("transcript_display") or display_transcript(isoform_id, gene_id, gene_name)
            discovery_class = meta.get("transcript_discovery_class", "")
            (
                is_stringtie_assembly,
                assembly_evidence_class,
                assembly_evidence_label,
                assembly_evidence_note,
            ) = fallback_assembly_evidence(gene_id, isoform_id, discovery_class)
            status = "candidate" if isoform_id in candidate_isoforms else "same_gene_context"
            candidate = (
                {
                    "event_id": event_id,
                    "contrast_id": contrast_id,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_display": gene_display,
                    "gene_biotype": gene_biotype,
                    "transcript_biotype": transcript_biotype,
                    "switch_biotype_class": "",
                    "switch_rank": str(switch_rank),
                    "isoform_id": isoform_id,
                    "transcript_display": transcript_display,
                    "switch_role": role_for_isoform(isoform_id, switch_in["_isoform_id"], switch_out["_isoform_id"]),
                    "dIF": row["_dIF"],
                    "padj_qvalue": row["_statistic"],
                    "switch_statistic": row["_statistic"],
                    "switch_statistic_name": row["_statistic_name"],
                    "isoform_fraction_control": fraction_value(row, "control"),
                    "isoform_fraction_test": fraction_value(row, "test"),
                    "switch_direction": switch_direction(row["_dIF"]),
                    "novelty_group": meta.get("transcript_plot_group", ""),
                    "reason_selected": f"abs(dIF)>={dif}; {row['_statistic_name']}<={padj}"
                    if isoform_id in candidate_isoforms
                    else "same-gene context for selected switch",
                    "candidate_status": status,
                    "transcript_discovery_class": discovery_class,
                    "transcript_novelty": meta.get("transcript_novelty", ""),
                    "transcript_plot_group": meta.get("transcript_plot_group", ""),
                    "gffcompare_class_code": meta.get("gffcompare_class_code", meta.get("class_code", "")),
                    "is_stringtie_assembly": meta.get("is_stringtie_assembly", is_stringtie_assembly),
                    "assembly_evidence_class": meta.get("assembly_evidence_class", assembly_evidence_class),
                    "assembly_evidence_label": meta.get("assembly_evidence_label", assembly_evidence_label),
                    "assembly_evidence_note": meta.get("assembly_evidence_note", assembly_evidence_note),
                    "consequence_summary": "",
                    "source_detailed": row["_detailed"],
                    "ORF_length": first_available_feature(row, ["ORF_length", "orf_length", "ORF_length_aa", "orf_length_aa"]),
                    "ORF_start": first_available_feature(row, ["ORF_start", "orf_start"]),
                    "ORF_end": first_available_feature(row, ["ORF_end", "orf_end"]),
                    "premature_stop": first_available_feature(row, ["premature_stop", "prematureStop", "PTC", "has_premature_stop"]),
                    "NMD_status": first_available_feature(row, ["NMD_status", "nmd_status", "NMD"]),
                    "coding_potential": first_available_feature(row, ["coding_potential", "codingPotential", "CPAT_prediction", "CPC2_prediction"]),
                }
            )
            candidate_rows.append(candidate)
            current_candidate_rows.append(candidate)
        classified_rows = [
            row for row in current_candidate_rows if row["switch_role"] in {"switch_in", "switch_out"}
        ] or current_candidate_rows
        switch_class, gene_biotype, interpretation_label = classify_switch_biotype(classified_rows)
        event["gene_display"] = first_existing(current_candidate_rows[0], ["gene_display"]) if current_candidate_rows else display_gene(gene_id, gene_name)
        event["gene_biotype"] = gene_biotype
        event["switch_biotype_class"] = switch_class
        event["switch_interpretation_label"] = interpretation_label
        for candidate in current_candidate_rows:
            candidate["switch_biotype_class"] = switch_class
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
    event_rows: list[dict[str, str]],
    event_context: dict[str, dict[str, object]],
    models: dict[str, TranscriptModel],
    nt_sequences: dict[str, str],
    aa_sequences: dict[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen = set()
    event_by_id = {row["event_id"]: row for row in event_rows}
    for row in candidate_rows:
        key = (row["event_id"], row["isoform_id"])
        if key in seen:
            continue
        seen.add(key)
        isoform_id = row["isoform_id"]
        event = event_by_id[row["event_id"]]
        paired = event["switch_out_isoform"] if isoform_id == event["switch_in_isoform"] else event["switch_in_isoform"]
        pair_role = "gained_isoform" if isoform_id == event["switch_in_isoform"] else "lost_isoform" if isoform_id == event["switch_out_isoform"] else "context_isoform"
        unique = event_context.get(row["event_id"], {}).get("_unique_context", {})
        if not isinstance(unique, dict):
            unique = {}
        model = models.get(isoform_id)
        nt_seq = nt_sequences.get(isoform_id, "")
        aa_seq = aa_sequences.get(isoform_id, "")
        if pair_role == "gained_isoform":
            affected_exons = unique.get("gained_exons", [])
            affected_aa = unique.get("gained_aa", [])
            opposite_exons = unique.get("lost_exons", [])
            opposite_aa = unique.get("lost_aa", [])
        elif pair_role == "lost_isoform":
            affected_exons = unique.get("lost_exons", [])
            affected_aa = unique.get("lost_aa", [])
            opposite_exons = unique.get("gained_exons", [])
            opposite_aa = unique.get("gained_aa", [])
        else:
            affected_exons = []
            affected_aa = []
            opposite_exons = []
            opposite_aa = []
        affected_exons = list(affected_exons) if isinstance(affected_exons, list) else []
        affected_aa = list(affected_aa) if isinstance(affected_aa, list) else []
        opposite_exons = list(opposite_exons) if isinstance(opposite_exons, list) else []
        opposite_aa = list(opposite_aa) if isinstance(opposite_aa, list) else []
        affected_nt_intervals = genomic_to_transcript_intervals(model, affected_exons, "exon") if model else []
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
                "switch_pair_id": row["event_id"],
                "paired_isoform_id": paired,
                "switch_role": row["switch_role"],
                "pair_role": pair_role,
                "orf_length_aa": first_available_feature(row, ["ORF_length", "orf_length", "ORF_length_aa", "orf_length_aa"])
                or (str(len(aa_seq)) if aa_seq else ""),
                "orf_start": first_available_feature(row, ["ORF_start", "orf_start"]),
                "orf_end": first_available_feature(row, ["ORF_end", "orf_end"]),
                "cds_coordinates": interval_text(model.chrom, model.cds) if model and model.cds else "",
                "premature_stop": first_available_feature(row, ["premature_stop", "prematureStop", "PTC", "has_premature_stop"]),
                "nmd_status": first_available_feature(row, ["NMD_status", "nmd_status", "NMD"]),
                "coding_potential": first_available_feature(row, ["coding_potential", "codingPotential", "CPAT_prediction", "CPC2_prediction"]),
                "gained_exon_coordinates": interval_text(str(unique.get("chrom", "")), list(unique.get("gained_exons", [])))
                if isinstance(unique.get("gained_exons", []), list)
                else "",
                "lost_exon_coordinates": interval_text(str(unique.get("chrom", "")), list(unique.get("lost_exons", [])))
                if isinstance(unique.get("lost_exons", []), list)
                else "",
                "gained_aa_interval": interval_span_text(list(unique.get("gained_aa", [])))
                if isinstance(unique.get("gained_aa", []), list)
                else "",
                "lost_aa_interval": interval_span_text(list(unique.get("lost_aa", [])))
                if isinstance(unique.get("lost_aa", []), list)
                else "",
                "nt_length": str(len(nt_seq)) if nt_seq else "0",
                "aa_length": str(len(aa_seq)) if aa_seq else "0",
                "nt_sequence": nt_seq,
                "aa_sequence": aa_seq,
                "affected_nt_sequence": nt_subsequence(nt_seq, affected_nt_intervals),
                "affected_aa_sequence": aa_subsequence(aa_seq, affected_aa),
                "sequence_status": status,
            }
        )
    return rows


def transcript_model_length(model: Optional[TranscriptModel]) -> Optional[int]:
    if model is None or not model.exons:
        return None
    return sum(end - start + 1 for start, end in model.exons)


def row_transcript_length(row: dict[str, str], models: dict[str, TranscriptModel]) -> Optional[int]:
    value = to_float(row.get("nt_length", ""))
    if value is not None and value > 0:
        return int(value)
    return transcript_model_length(models.get(row.get("isoform_id", "")))


def tss_tes(model: Optional[TranscriptModel]) -> tuple[Optional[int], Optional[int]]:
    if model is None or not model.exons:
        return None, None
    starts = [start for start, _ in model.exons]
    ends = [end for _, end in model.exons]
    if model.strand == "-":
        return max(ends), min(starts)
    return min(starts), max(ends)


def signed_change(value: Optional[int]) -> str:
    if value is None:
        return ""
    return f"{value:+d}"


def coordinate_change(current: Optional[int], paired: Optional[int]) -> str:
    if current is None or paired is None:
        return "not_available"
    delta = current - paired
    return "unchanged" if delta == 0 else signed_change(delta)


def model_splice_junctions(model: Optional[TranscriptModel]) -> set[tuple[str, int, int, str]]:
    if model is None or len(model.exons) < 2:
        return set()
    junctions: set[tuple[str, int, int, str]] = set()
    for left, right in zip(sorted(model.exons), sorted(model.exons)[1:]):
        intron_start = left[1] + 1
        intron_end = right[0] - 1
        if intron_start <= intron_end:
            junctions.add((model.chrom, intron_start, intron_end, model.strand))
    return junctions


def splice_junction_text(junctions: Iterable[tuple[str, int, int, str]]) -> str:
    return "; ".join(f"{chrom}:{start}-{end}:{strand}" for chrom, start, end, strand in sorted(junctions))


def splice_junction_gain_loss(
    current: Optional[TranscriptModel],
    paired: Optional[TranscriptModel],
) -> tuple[set[tuple[str, int, int, str]], set[tuple[str, int, int, str]]]:
    current_junctions = model_splice_junctions(current)
    paired_junctions = model_splice_junctions(paired)
    return current_junctions - paired_junctions, paired_junctions - current_junctions


def gene_loci_from_models(models: dict[str, TranscriptModel]) -> list[GeneLocus]:
    loci: dict[str, GeneLocus] = {}
    for model in models.values():
        if not model.gene_id or not model.exons:
            continue
        start = min(exon_start for exon_start, _ in model.exons)
        end = max(exon_end for _, exon_end in model.exons)
        locus = loci.get(model.gene_id)
        if locus is None:
            loci[model.gene_id] = GeneLocus(
                gene_id=model.gene_id,
                gene_name=model.gene_name,
                gene_biotype=model.gene_biotype,
                chrom=model.chrom,
                strand=model.strand,
                start=start,
                end=end,
            )
            continue
        locus.start = min(locus.start, start)
        locus.end = max(locus.end, end)
        if model.gene_name and not locus.gene_name:
            locus.gene_name = model.gene_name
        if model.gene_biotype and not locus.gene_biotype:
            locus.gene_biotype = model.gene_biotype
    return list(loci.values())


def promoter_context_change(current: Optional[TranscriptModel], paired: Optional[TranscriptModel]) -> str:
    if current is None or paired is None:
        return "not_available"
    current_tss, _ = tss_tes(current)
    paired_tss, _ = tss_tes(paired)
    if current_tss is None or paired_tss is None:
        return "not_available"
    if current.chrom != paired.chrom or current.strand != paired.strand:
        return "different_locus_or_strand"
    delta = current_tss - paired_tss
    if delta == 0:
        return "unchanged"
    if current.strand == "-":
        direction = "upstream" if delta > 0 else "downstream"
    else:
        direction = "upstream" if delta < 0 else "downstream"
    return f"{direction}_TSS_shift:{abs(delta)}bp"


def distance_to_locus(point: int, locus: GeneLocus) -> int:
    if locus.start <= point <= locus.end:
        return 0
    return min(abs(point - locus.start), abs(point - locus.end))


def proximal_gene_context(
    model: Optional[TranscriptModel],
    gene_loci: list[GeneLocus],
    radius_bp: int = 1000,
) -> str:
    if model is None:
        return "not_available"
    tss, _ = tss_tes(model)
    if tss is None:
        return "not_available"
    matches = []
    for locus in gene_loci:
        if locus.gene_id == model.gene_id or locus.chrom != model.chrom:
            continue
        distance = distance_to_locus(tss, locus)
        if distance > radius_bp:
            continue
        relation = "overlap" if distance == 0 else "proximal"
        strand_relation = "same_strand" if locus.strand == model.strand else "opposite_strand"
        label = (
            f"{relation}:{locus.gene_id}"
            f"({locus.gene_name or locus.gene_id};"
            f"{locus.gene_biotype or 'biotype_unknown'};"
            f"{strand_relation};distance={distance}bp)"
        )
        matches.append((distance, label))
    if not matches:
        return f"no_other_gene_within_{radius_bp}bp"
    return "; ".join(label for _, label in sorted(matches)[:5])


def retained_intron_label(row: dict[str, str]) -> str:
    text = " ".join(
        [
            row.get("transcript_biotype", ""),
            row.get("transcript_discovery_class", ""),
            row.get("gffcompare_class_code", ""),
            row.get("consequence_summary", ""),
        ]
    ).lower()
    if "retained_intron" in text or "retained intron" in text:
        return "retained_intron_annotation"
    return "not_assessed"


def antisense_label(row: dict[str, str]) -> str:
    text = " ".join([row.get("gene_biotype", ""), row.get("transcript_biotype", "")]).lower()
    if "antisense" in text:
        return "annotated_antisense"
    return "not_assessed"


def coding_potential_change(row: dict[str, str], paired: dict[str, str]) -> str:
    current = row.get("coding_potential", "")
    other = paired.get("coding_potential", "")
    if not current and not other:
        return "not_available"
    if current == other:
        return "unchanged"
    return f"{other or 'unknown'}->{current or 'unknown'}"


def parse_coordinate_int(value: str) -> Optional[int]:
    cleaned = (value or "").replace(",", "").strip()
    if not cleaned:
        return None
    try:
        return int(float(cleaned))
    except ValueError:
        return None


def parse_genomic_intervals(text: str) -> list[tuple[str, int, int]]:
    intervals: list[tuple[str, int, int]] = []
    for match in re.finditer(r"([^:,\s;]+):(\d[\d,]*)-(\d[\d,]*)", text or ""):
        start = parse_coordinate_int(match.group(2))
        end = parse_coordinate_int(match.group(3))
        if start is None or end is None:
            continue
        intervals.append((match.group(1), min(start, end), max(start, end)))
    return intervals


def row_genomic_intervals(row: dict[str, str]) -> list[tuple[str, int, int]]:
    chrom = value_by_normalized_key(row, ["chrom", "chr", "seqname", "seqnames", "sequence_name"])
    start = parse_coordinate_int(value_by_normalized_key(row, ["start", "start_position", "genomic_start"]))
    end = parse_coordinate_int(value_by_normalized_key(row, ["end", "stop", "end_position", "genomic_end"]))
    intervals: list[tuple[str, int, int]] = []
    if chrom and start is not None and end is not None:
        intervals.append((chrom, min(start, end), max(start, end)))
    for field in ["coordinates", "coordinate", "genomic_coordinates", "interval", "locus", "region"]:
        intervals.extend(parse_genomic_intervals(value_by_normalized_key(row, [field])))
    return list(dict.fromkeys(intervals))


def ncrna_resource_category(feature_type: str) -> str:
    normalized = normalized_key(feature_type)
    if "antisense" in normalized:
        return "antisense"
    if any(token in normalized for token in ["host_small_rna", "host_smallrna", "embedded_small_rna", "embedded_snorna"]):
        return "host_smallrna"
    if any(token in normalized for token in ["snorna", "scrna", "trna_fragment", "small_rna_locus"]):
        return "host_smallrna"
    if any(token in normalized for token in ["conserved", "phastcons", "phylop", "gerp"]):
        return "conserved"
    if any(token in normalized for token in ["motif", "rbp", "mirna_binding", "seed_match"]):
        return "motif"
    return "other"


def split_optional_paths(paths_text: str) -> list[Path]:
    return [Path(part.strip()) for part in (paths_text or "").split(",") if part.strip()]


def read_ncrna_annotation_rows(paths_text: str) -> list[dict[str, object]]:
    annotations: list[dict[str, object]] = []
    for path in split_optional_paths(paths_text):
        _, rows = read_table(path)
        for row in rows:
            transcript_id = value_by_normalized_key(
                row,
                ["transcript_id", "isoform_id", "target_id", "query_id"],
            )
            gene_id = value_by_normalized_key(row, ["gene_id", "gene", "target_gene", "geneid"])
            intervals = row_genomic_intervals(row)
            if not transcript_id and not gene_id and not intervals:
                continue
            feature_type = value_by_normalized_key(
                row,
                ["feature_type", "type", "category", "annotation_type", "resource_type"],
            )
            source = value_by_normalized_key(row, ["source", "database", "resource", "resource_name"]) or path.stem
            feature_id = value_by_normalized_key(row, ["feature_id", "id", "accession", "motif_id", "target_id"])
            feature_name = value_by_normalized_key(row, ["feature_name", "name", "label", "motif_name", "description"])
            description = value_by_normalized_key(row, ["description", "note", "details", "evidence"])
            annotations.append(
                {
                    "source": source,
                    "feature_type": feature_type or "annotation",
                    "category": ncrna_resource_category(feature_type or feature_name or description),
                    "feature_id": feature_id,
                    "feature_name": feature_name,
                    "description": description,
                    "transcript_id": transcript_id,
                    "gene_id": gene_id,
                    "intervals": intervals,
                }
            )
    return annotations


def intervals_overlap(left: tuple[str, int, int], right: tuple[str, int, int]) -> bool:
    return left[0] == right[0] and left[1] <= right[2] and right[1] <= left[2]


def ncrna_annotation_matches(
    row: dict[str, str],
    sequence: dict[str, str],
    annotation: dict[str, object],
) -> bool:
    transcript_id = str(annotation.get("transcript_id", ""))
    gene_id = str(annotation.get("gene_id", ""))
    id_match = (transcript_id and transcript_id == row.get("isoform_id", "")) or (
        gene_id and gene_id == row.get("gene_id", "")
    )
    annotation_intervals = annotation.get("intervals", [])
    if not isinstance(annotation_intervals, list):
        annotation_intervals = []
    affected_intervals = parse_genomic_intervals(
        ",".join(
            [
                sequence.get("gained_exon_coordinates", ""),
                sequence.get("lost_exon_coordinates", ""),
            ]
        )
    )
    if annotation_intervals and affected_intervals:
        interval_match = any(
            intervals_overlap(left, right)
            for left in annotation_intervals
            for right in affected_intervals
        )
        return interval_match and (id_match or not transcript_id and not gene_id)
    if annotation_intervals and not affected_intervals:
        return bool(id_match)
    return bool(id_match)


def ncrna_annotation_label(annotation: dict[str, object]) -> str:
    source = str(annotation.get("source", "resource"))
    feature_type = str(annotation.get("feature_type", "annotation"))
    feature_id = str(annotation.get("feature_id", ""))
    feature_name = str(annotation.get("feature_name", ""))
    identity = feature_id or feature_name or str(annotation.get("description", "")) or "annotation"
    if feature_name and feature_name != identity:
        identity = f"{identity} ({feature_name})"
    return f"{source}:{feature_type}:{identity}"


def summarize_ncrna_resource_annotations(matches: list[dict[str, object]]) -> dict[str, str]:
    grouped: dict[str, list[str]] = defaultdict(list)
    for annotation in matches:
        grouped[str(annotation.get("category", "other"))].append(ncrna_annotation_label(annotation))

    def joined(category: str) -> str:
        return "; ".join(dict.fromkeys(grouped.get(category, [])))

    return {
        "conserved_exon_change": joined("conserved"),
        "motif_change": joined("motif"),
        "host_smallrna_change": joined("host_smallrna"),
        "resource_antisense_overlap": joined("antisense"),
        "ncrna_resource_annotations": "; ".join(dict.fromkeys(ncrna_annotation_label(row) for row in matches)),
    }


def ncrna_interpretation_label(
    row: dict[str, str],
    exon_gain_loss: str,
    tss_change: str,
    tes_change: str,
    gained_junctions: set[tuple[str, int, int, str]],
    lost_junctions: set[tuple[str, int, int, str]],
    promoter_change: str,
) -> str:
    switch_class = row.get("switch_biotype_class", "")
    if is_pseudogene_switch(row):
        return "pseudogene_transcript_architecture_change"
    if switch_class == "ambiguous_artifact":
        return "ambiguous_or_artifact"
    if switch_class == "mixed_coding_noncoding":
        return "coding_potential_transition"
    if (
        exon_gain_loss
        or gained_junctions
        or lost_junctions
        or tss_change not in {"", "unchanged", "not_available"}
        or tes_change not in {"", "unchanged", "not_available"}
        or promoter_change not in {"", "unchanged", "not_available"}
    ):
        return "noncoding_structure_change"
    return "noncoding_isoform_fraction_change"


def pseudogene_caution_label(row: dict[str, str]) -> str:
    if not is_pseudogene_switch(row):
        return ""
    return "interpret_as_transcript_architecture_not_protein_consequence"


def build_ncrna_switch_rows(
    candidate_rows: list[dict[str, str]],
    sequence_rows: list[dict[str, str]],
    models: dict[str, TranscriptModel],
    ncrna_annotations: list[dict[str, object]],
) -> list[dict[str, str]]:
    sequence_by_key = {(row["event_id"], row["isoform_id"]): row for row in sequence_rows}
    candidate_by_key = {(row["event_id"], row["isoform_id"]): row for row in candidate_rows}
    gene_loci = gene_loci_from_models(models)
    rows: list[dict[str, str]] = []
    for row in candidate_rows:
        if row.get("switch_role") not in {"switch_in", "switch_out"}:
            continue
        if row.get("switch_biotype_class") == "coding":
            continue
        sequence = sequence_by_key.get((row["event_id"], row["isoform_id"]), {})
        paired_isoform = sequence.get("paired_isoform_id", "")
        paired_candidate = candidate_by_key.get((row["event_id"], paired_isoform), {})
        paired_sequence = sequence_by_key.get((row["event_id"], paired_isoform), {})
        current_length = row_transcript_length(sequence, models)
        paired_length = row_transcript_length(paired_sequence, models)
        length_change = current_length - paired_length if current_length is not None and paired_length is not None else None
        current_model = models.get(row["isoform_id"])
        paired_model = models.get(paired_isoform)
        current_tss, current_tes = tss_tes(current_model)
        paired_tss, paired_tes = tss_tes(paired_model)
        tss_change = coordinate_change(current_tss, paired_tss)
        tes_change = coordinate_change(current_tes, paired_tes)
        gained_junctions, lost_junctions = splice_junction_gain_loss(current_model, paired_model)
        promoter_change = promoter_context_change(current_model, paired_model)
        exon_gain_loss = "; ".join(
            part
            for part in [
                f"gained={sequence.get('gained_exon_coordinates', '')}" if sequence.get("gained_exon_coordinates", "") else "",
                f"lost={sequence.get('lost_exon_coordinates', '')}" if sequence.get("lost_exon_coordinates", "") else "",
            ]
            if part
        )
        resource_summary = summarize_ncrna_resource_annotations(
            [
                annotation
                for annotation in ncrna_annotations
                if ncrna_annotation_matches(row, sequence, annotation)
            ]
        )
        rows.append(
            {
                "event_id": row["event_id"],
                "contrast_id": row["contrast_id"],
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "gene_biotype": row.get("gene_biotype", ""),
                "transcript_biotype": row.get("transcript_biotype", ""),
                "switch_biotype_class": row.get("switch_biotype_class", ""),
                "isoform_id": row["isoform_id"],
                "paired_isoform_id": paired_isoform,
                "switch_role": row["switch_role"],
                "dIF": row.get("dIF", ""),
                "padj_qvalue": row.get("padj_qvalue", ""),
                "transcript_length_change": signed_change(length_change),
                "exon_gain_loss": exon_gain_loss,
                "intron_retention_change": retained_intron_label(row),
                "gained_splice_junctions": splice_junction_text(gained_junctions),
                "lost_splice_junctions": splice_junction_text(lost_junctions),
                "n_gained_splice_junctions": str(len(gained_junctions)),
                "n_lost_splice_junctions": str(len(lost_junctions)),
                "TSS_change": tss_change,
                "TES_change": tes_change,
                "promoter_context_change": promoter_change,
                "isoform_proximal_gene_context": proximal_gene_context(current_model, gene_loci),
                "paired_isoform_proximal_gene_context": proximal_gene_context(paired_model, gene_loci),
                "antisense_overlap": antisense_label(row),
                "conserved_exon_change": resource_summary["conserved_exon_change"],
                "motif_change": resource_summary["motif_change"],
                "host_smallrna_change": resource_summary["host_smallrna_change"],
                "resource_antisense_overlap": resource_summary["resource_antisense_overlap"],
                "ncrna_resource_annotations": resource_summary["ncrna_resource_annotations"],
                "pseudogene_caution": pseudogene_caution_label(row),
                "coding_potential_change": coding_potential_change(row, paired_candidate),
                "interpretation_label": ncrna_interpretation_label(
                    row,
                    exon_gain_loss,
                    tss_change,
                    tes_change,
                    gained_junctions,
                    lost_junctions,
                    promoter_change,
                ),
            }
        )
    rows.sort(key=lambda item: (item["switch_biotype_class"], item["event_id"], item["switch_role"], item["isoform_id"]))
    return rows


def unique_join(values: Iterable[str]) -> str:
    return "; ".join(dict.fromkeys(value for value in values if value))


def annotation_display_label(row: dict[str, str]) -> str:
    feature_id = row.get("feature_id", "")
    feature_name = row.get("feature_name", "")
    source = row.get("source", "")
    label = feature_id or feature_name or row.get("description", "") or "annotation"
    if feature_name and feature_name != label:
        label = f"{label} ({feature_name})"
    return f"{source}:{label}" if source else label


def annotation_category(row: dict[str, str]) -> str:
    feature_type = normalized_key(row.get("feature_type", ""))
    feature_name = normalized_key(row.get("feature_name", ""))
    feature_id = normalized_key(row.get("feature_id", ""))
    text = "_".join([feature_type, feature_name, feature_id])
    if "signal_peptide" in text or text == "signal" or "signalp" in normalized_key(row.get("source", "")):
        return "signal_peptide"
    if "transmembrane" in text or "tmhelix" in text or "tm_region" in text:
        return "transmembrane"
    if "localization" in text or "deeploc" in normalized_key(row.get("source", "")):
        return "localization"
    if "disorder" in text or "iupred" in normalized_key(row.get("source", "")):
        return "disorder"
    if "domain" in text or "pfam" in text or "interpro" in normalized_key(row.get("source", "")):
        return "domain"
    return feature_type or "other"


def annotations_by_change(
    rows: list[dict[str, str]],
    category: str,
    change_prefix: str,
) -> list[str]:
    labels = []
    for row in rows:
        if annotation_category(row) != category:
            continue
        if row.get("feature_change", "").startswith(change_prefix):
            labels.append(annotation_display_label(row))
    return list(dict.fromkeys(labels))


def sequence_for_isoform(sequence_rows: list[dict[str, str]], event_id: str, isoform_id: str) -> dict[str, str]:
    for row in sequence_rows:
        if row.get("event_id") == event_id and row.get("isoform_id") == isoform_id:
            return row
    return {}


def positive_nmd_status(value: str) -> bool:
    normalized = normalized_key(value)
    if not normalized:
        return False
    if normalized.startswith("not_nmd") or normalized in {"no_nmd", "non_nmd", "nmd_negative"}:
        return False
    return "nmd" in normalized


def nmd_change_label(switch_out: dict[str, str], switch_in: dict[str, str]) -> str:
    before = switch_out.get("nmd_status", "")
    after = switch_in.get("nmd_status", "")
    if not before and not after:
        return ""
    if positive_nmd_status(before) == positive_nmd_status(after) and normalized_key(before) == normalized_key(after):
        return ""
    return f"{before or 'unknown'}->{after or 'unknown'}"


def coding_potential_class(value: str) -> str:
    normalized = normalized_key(value)
    if not normalized:
        return ""
    if any(token in normalized for token in ["noncoding", "non_coding", "low_potential", "no_orf", "no_coding"]):
        return "noncoding_or_low_potential"
    if "coding" in normalized or normalized in {"cp", "protein_coding"}:
        return "coding"
    return normalized


def coding_potential_transition_label(switch_out: dict[str, str], switch_in: dict[str, str]) -> str:
    before = coding_potential_class(switch_out.get("coding_potential", ""))
    after = coding_potential_class(switch_in.get("coding_potential", ""))
    if not before and not after:
        return ""
    if before == after:
        return ""
    return f"{before or 'unknown'}->{after or 'unknown'}"


def orf_length_change_label(switch_out: dict[str, str], switch_in: dict[str, str]) -> tuple[str, bool]:
    before = to_float(switch_out.get("orf_length_aa", ""))
    after = to_float(switch_in.get("orf_length_aa", ""))
    if before is None or after is None:
        return "", False
    delta = int(after - before)
    if delta == 0:
        return "0", False
    denominator = max(abs(before), 1.0)
    large = abs(delta) >= 50 or abs(delta) / denominator >= 0.25
    return signed_change(delta), large


def localization_change_label(
    event_annotations: list[dict[str, str]],
    switch_out_isoform: str,
    switch_in_isoform: str,
) -> str:
    before = {
        annotation_display_label(row)
        for row in event_annotations
        if row.get("isoform_id") == switch_out_isoform and annotation_category(row) == "localization"
    }
    after = {
        annotation_display_label(row)
        for row in event_annotations
        if row.get("isoform_id") == switch_in_isoform and annotation_category(row) == "localization"
    }
    if before == after:
        return ""
    gained = sorted(after - before)
    lost = sorted(before - after)
    pieces = []
    if gained:
        pieces.append("gained=" + "; ".join(gained))
    if lost:
        pieces.append("lost=" + "; ".join(lost))
    return " | ".join(pieces)


def priority_tier(score: int) -> str:
    if score >= 8:
        return "high"
    if score >= 4:
        return "moderate"
    if score > 0:
        return "low"
    return "annotation_limited"


def build_coding_switch_rows(
    event_rows: list[dict[str, str]],
    candidate_rows: list[dict[str, str]],
    sequence_rows: list[dict[str, str]],
    annotation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    candidate_by_event = defaultdict(list)
    for row in candidate_rows:
        candidate_by_event[row["event_id"]].append(row)
    annotation_by_event = defaultdict(list)
    for row in annotation_rows:
        annotation_by_event[row["event_id"]].append(row)

    rows = []
    for event in event_rows:
        if event.get("switch_biotype_class") != "coding":
            continue
        event_id = event["event_id"]
        switch_in = sequence_for_isoform(sequence_rows, event_id, event["switch_in_isoform"])
        switch_out = sequence_for_isoform(sequence_rows, event_id, event["switch_out_isoform"])
        event_annotations = annotation_by_event[event_id]
        gained_domain = annotations_by_change(event_annotations, "domain", "gained")
        lost_domain = annotations_by_change(event_annotations, "domain", "lost")
        gained_signal = annotations_by_change(event_annotations, "signal_peptide", "gained")
        lost_signal = annotations_by_change(event_annotations, "signal_peptide", "lost")
        gained_tm = annotations_by_change(event_annotations, "transmembrane", "gained")
        lost_tm = annotations_by_change(event_annotations, "transmembrane", "lost")
        gained_disorder = annotations_by_change(event_annotations, "disorder", "gained")
        lost_disorder = annotations_by_change(event_annotations, "disorder", "lost")
        localization_change = localization_change_label(event_annotations, event["switch_out_isoform"], event["switch_in_isoform"])
        nmd_change = nmd_change_label(switch_out, switch_in)
        coding_change = coding_potential_transition_label(switch_out, switch_in)
        orf_change, large_orf_change = orf_length_change_label(switch_out, switch_in)

        score = 0
        reasons = []
        if gained_domain or lost_domain:
            score += 3
            reasons.append("protein_domain_gain_loss")
        if gained_signal or lost_signal:
            score += 3
            reasons.append("signal_peptide_gain_loss")
        if gained_tm or lost_tm:
            score += 3
            reasons.append("transmembrane_region_gain_loss")
        if nmd_change:
            score += 4
            reasons.append("NMD_status_change")
        if coding_change:
            score += 4
            reasons.append("coding_potential_transition")
        if large_orf_change:
            score += 2
            reasons.append("large_ORF_length_change")
        if gained_disorder or lost_disorder:
            score += 1
            reasons.append("disorder_region_gain_loss")
        if localization_change:
            score += 2
            reasons.append("localization_change")

        consequence_summary = unique_join(
            row.get("consequence_summary", "")
            for row in candidate_by_event[event_id]
            if row.get("switch_role") in {"switch_in", "switch_out"}
        )
        rows.append(
            {
                "event_id": event_id,
                "contrast_id": event["contrast_id"],
                "gene_id": event["gene_id"],
                "gene_name": event["gene_name"],
                "gene_biotype": event.get("gene_biotype", ""),
                "switch_rank": event.get("switch_rank", ""),
                "coding_priority_rank": "",
                "coding_priority_score": str(score),
                "coding_priority_tier": priority_tier(score),
                "coding_priority_reasons": unique_join(reasons),
                "functional_consequence_summary": consequence_summary,
                "switch_in_isoform": event["switch_in_isoform"],
                "switch_out_isoform": event["switch_out_isoform"],
                "switch_in_dIF": event["switch_in_dIF"],
                "switch_out_dIF": event["switch_out_dIF"],
                "max_abs_dIF": event["max_abs_dIF"],
                "best_switch_statistic": event["best_switch_statistic"],
                "nmd_change": nmd_change,
                "coding_potential_change": coding_change,
                "orf_length_change_aa": orf_change,
                "gained_domain": unique_join(gained_domain),
                "lost_domain": unique_join(lost_domain),
                "gained_signal_peptide": unique_join(gained_signal),
                "lost_signal_peptide": unique_join(lost_signal),
                "gained_transmembrane_region": unique_join(gained_tm),
                "lost_transmembrane_region": unique_join(lost_tm),
                "gained_disorder_region": unique_join(gained_disorder),
                "lost_disorder_region": unique_join(lost_disorder),
                "localization_change": localization_change,
                "n_functional_annotations": event.get("n_functional_annotations", "0"),
                "event_html": event.get("event_html", ""),
            }
        )

    rows.sort(
        key=lambda item: (
            -int(item["coding_priority_score"]),
            -float(item.get("max_abs_dIF") or "0"),
            item["event_id"],
        )
    )
    event_by_id = {row["event_id"]: row for row in event_rows}
    for rank, row in enumerate(rows, start=1):
        row["coding_priority_rank"] = str(rank)
        event = event_by_id.get(row["event_id"], {})
        event["coding_priority_rank"] = str(rank)
        event["coding_priority_score"] = row["coding_priority_score"]
        event["coding_priority_tier"] = row["coding_priority_tier"]
        event["coding_priority_reasons"] = row["coding_priority_reasons"]
    return rows


def relative(path: str, base: Path) -> str:
    return os.path.relpath(path, start=base.parent)


def path_has_content(path_text: str) -> bool:
    if not path_text:
        return False
    path = Path(path_text)
    return path.is_file() and path.stat().st_size > 0


def file_link(label: str, path_text: str, html_path: Path) -> str:
    if not path_text:
        return f'<span class="muted">{html.escape(label)} not generated</span>'
    path = Path(path_text)
    if not path.exists():
        return f'<span class="muted">{html.escape(label)} missing</span>'
    if path.is_file() and path.stat().st_size == 0:
        return f'<span class="muted">{html.escape(label)} empty</span>'
    return f'<a href="{html.escape(relative(path_text, html_path))}">{html.escape(label)}</a>'


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
    row_height = 82
    margin_left = 210
    margin_right = 260
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
    fraction_by_isoform = {
        row["_isoform_id"]: (fraction_value(row, "control"), fraction_value(row, "test"))
        for row in gene_rows
    }
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
            for cds_start, cds_end in model.cds:
                cds_x = x_for(cds_start)
                cds_w = max(3, x_for(cds_end) - cds_x)
                pieces.append(
                    f'<rect x="{cds_x:.1f}" y="{y - 3:.1f}" width="{cds_w:.1f}" height="6" '
                    'fill="#24292f" fill-opacity="0.85"><title>CDS/ORF segment</title></rect>'
                )
        control_fraction, test_fraction = fraction_by_isoform.get(isoform_id, ("", ""))
        bar_x = width - 225
        bar_y = y - 15
        pieces.append(svg_text(bar_x, bar_y - 5, "IF control/test", 10))
        for offset, value, fill in [(0, control_fraction, "#8c959f"), (14, test_fraction, color)]:
            parsed = to_float(value) or 0
            parsed = parsed / 100 if parsed > 1 else parsed
            pieces.append(f'<rect x="{bar_x}" y="{bar_y + offset}" width="100" height="10" fill="#f6f8fa" stroke="#d0d7de"/>')
            pieces.append(f'<rect x="{bar_x}" y="{bar_y + offset}" width="{max(0, min(100, parsed * 100)):.1f}" height="10" fill="{fill}"/>')
            pieces.append(svg_text(bar_x + 108, bar_y + offset + 9, value or "NA", 10))
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
                if annotation.get("feature_change", "") in {"gained_domain_overlap", "lost_domain_overlap"}:
                    label_text = f"{annotation['feature_change'].replace('_domain_overlap', '')}: {label_text}"
                pieces.append(
                    f'<rect x="{feature_x:.1f}" y="{anno_y - 7:.1f}" width="{feature_w:.1f}" height="14" '
                    'fill="#8250df" fill-opacity="0.7" stroke="#6639ba" stroke-width="0.5">'
                    f'<title>{html.escape(label_text)}</title></rect>'
                )
                if annotation.get("feature_change", "") in {"gained_domain_overlap", "lost_domain_overlap"}:
                    pieces.append(svg_text(feature_x, anno_y + 22, label_text[:42], 9, "700"))
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
    ncrna_rows: list[dict[str, str]],
    coding_switch_rows: list[dict[str, str]],
    out_path: Path,
    svg_path: Path,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    event_candidates = [row for row in candidate_rows if row["event_id"] == event["event_id"]]
    event_sequences = [row for row in sequence_rows if row["event_id"] == event["event_id"]]
    event_annotations = [row for row in annotation_rows if row["event_id"] == event["event_id"]]
    event_ncrna_rows = [row for row in ncrna_rows if row["event_id"] == event["event_id"]]
    event_coding_rows = [row for row in coding_switch_rows if row["event_id"] == event["event_id"]]

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
            f"<p>ORF length: {html.escape(row.get('orf_length_aa', ''))}; "
            f"CDS: {html.escape(row.get('cds_coordinates', ''))}; "
            f"NMD: {html.escape(row.get('nmd_status', ''))}; "
            f"coding potential: {html.escape(row.get('coding_potential', ''))}</p>"
            f"<p>Gained exons: {html.escape(row.get('gained_exon_coordinates', ''))}; "
            f"lost exons: {html.escape(row.get('lost_exon_coordinates', ''))}; "
            f"gained AA: {html.escape(row.get('gained_aa_interval', ''))}; "
            f"lost AA: {html.escape(row.get('lost_aa_interval', ''))}</p>"
            "<h3>Nucleotide sequence</h3>"
            f"<pre>{html.escape(row['nt_sequence'] or 'not available')}</pre>"
            "<h3>Amino-acid sequence</h3>"
            f"<pre>{html.escape(row['aa_sequence'] or 'not available')}</pre>"
            "<h3>Affected nucleotide sequence</h3>"
            f"<pre>{html.escape(row.get('affected_nt_sequence', '') or 'not available')}</pre>"
            "<h3>Affected amino-acid sequence</h3>"
            f"<pre>{html.escape(row.get('affected_aa_sequence', '') or 'not available')}</pre>"
            "</details>"
        )
    if not event_sequences:
        sequence_note = "No sequence rows were produced for this event."
    elif not any(row.get("nt_sequence", "") or row.get("aa_sequence", "") for row in event_sequences):
        sequence_note = (
            "Nucleotide and amino-acid sequences were not generated. Configure "
            "rnaseq_differential.isoform_switch_genome_object with a BSgenome object "
            "matching the reference annotation, then rerun the isoform-switch report "
            "targets to enable sequence links and sequence-dependent consequence tools."
        )
    else:
        sequence_note = (
            "Sequence blocks expose the nucleotide, amino-acid, and affected-region "
            "sequences used for manual inspection or external validation of the switch event."
        )
    report_index = out_path.parent.parent.parent / "index.html"
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
    th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; vertical-align: top; overflow-wrap: anywhere; }}
    th {{ background: #f6f8fa; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    pre {{ white-space: pre-wrap; word-break: break-word; background: #f6f8fa; padding: 12px; }}
    img {{ max-width: 100%; border: 1px solid #d0d7de; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    table {{ display: block; overflow-x: auto; }}
    .breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    .toc {{ display: flex; flex-wrap: wrap; gap: 0.5rem 0.85rem; margin: 1rem 0 1.25rem; }}
    .toc a {{ border: 1px solid #d0d7de; border-radius: 999px; padding: 0.25rem 0.65rem; }}
  </style>
</head>
<body>
  <nav class="breadcrumbs"><a href="{html.escape(relative(str(report_index), out_path))}">Isoform-switch overview</a> / {html.escape(event['contrast_id'])} / {html.escape(event['event_id'])}</nav>
  <h1>{html.escape(event.get('gene_display', '') or event['gene_name'])} isoform switch</h1>
  <div class="muted">{html.escape(event['event_id'])} | contrast {html.escape(event['contrast_id'])}</div>
  <nav class="toc" aria-label="Page sections">
    <a href="#diagram">Diagram</a>
    <a href="#candidates">Candidates</a>
    <a href="#coding">Coding evidence</a>
    <a href="#ncrna">ncRNA evidence</a>
    <a href="#annotations">Functional annotations</a>
    <a href="#sequences">Sequences</a>
  </nav>
  <p>Class: <strong>{html.escape(event.get('switch_biotype_class', ''))}</strong>;
     gene biotype: <strong>{html.escape(event.get('gene_biotype', '') or 'not annotated')}</strong>;
     interpretation: <strong>{html.escape(event.get('switch_interpretation_label', ''))}</strong>.</p>
  <p>Switch-in isoform: <strong>{html.escape(event['switch_in_isoform'])}</strong>;
     switch-out isoform: <strong>{html.escape(event['switch_out_isoform'])}</strong>.</p>
  <p class="note">The exon diagram compares the switch-in and switch-out isoforms in genomic coordinates. It is a structural view of which transcript model gains relative usage and which loses usage in the test group.</p>
  <img id="diagram" src="{html.escape(relative(str(svg_path), out_path))}" alt="Isoform switch plot">
  <h2 id="candidates">Candidate Isoforms</h2>
  <p class="note">Candidate isoforms are the transcripts selected from the IsoformSwitchAnalyzeR output for this event. dIF is the change in isoform fraction between groups; positive and negative roles identify switch-in and switch-out transcripts. StringTie/MSTRG identifiers are RNA-seq assembled transcript models; evidence labels are conservative review classes, not independent validation.</p>
  {table(['switch_rank', 'transcript_display', 'isoform_id', 'switch_role', 'assembly_evidence_label', 'gffcompare_class_code', 'gene_biotype', 'transcript_biotype', 'switch_biotype_class', 'dIF', 'padj_qvalue', 'isoform_fraction_control', 'isoform_fraction_test', 'switch_direction', 'novelty_group', 'reason_selected', 'assembly_evidence_note', 'consequence_summary'], event_candidates)}
  <h2 id="coding">Coding Switch Prioritization</h2>
  <p class="note">This table ranks coding switches by predicted consequence evidence such as NMD, coding-potential changes, ORF length changes, domains, signal peptides, transmembrane regions, or localization annotations when those resources are available.</p>
  {table(['coding_priority_rank', 'coding_priority_score', 'coding_priority_tier', 'coding_priority_reasons', 'nmd_change', 'coding_potential_change', 'orf_length_change_aa', 'gained_domain', 'lost_domain', 'gained_signal_peptide', 'lost_signal_peptide', 'gained_transmembrane_region', 'lost_transmembrane_region', 'localization_change'], event_coding_rows)}
  <h2 id="ncrna">ncRNA Switch Interpretation</h2>
  <p class="note">For noncoding or mixed coding-potential events, the report emphasizes transcript architecture: length, exon structure, junction changes, TSS/TES shifts, overlap context, host-smallRNA context, and coding-potential cautions.</p>
  {table(['isoform_id', 'paired_isoform_id', 'switch_role', 'gene_biotype', 'transcript_biotype', 'switch_biotype_class', 'transcript_length_change', 'exon_gain_loss', 'intron_retention_change', 'gained_splice_junctions', 'lost_splice_junctions', 'TSS_change', 'TES_change', 'promoter_context_change', 'isoform_proximal_gene_context', 'paired_isoform_proximal_gene_context', 'antisense_overlap', 'conserved_exon_change', 'motif_change', 'host_smallrna_change', 'resource_antisense_overlap', 'pseudogene_caution', 'coding_potential_change', 'interpretation_label'], event_ncrna_rows)}
  <h2 id="annotations">Functional Annotations</h2>
  <p class="note">Functional annotations are optional external evidence imported from tools or tables such as domain, disorder, localization, signal-peptide, or transmembrane predictions. Empty rows mean the event lacks configured external annotation evidence, not necessarily that no function exists.</p>
  {table(['isoform_id', 'source', 'feature_type', 'feature_id', 'feature_name', 'start_aa', 'end_aa', 'score', 'feature_change', 'description'], event_annotations)}
  <h2 id="sequences">Sequences</h2>
  <p class="note">{html.escape(sequence_note)}</p>
  {''.join(sequence_blocks) if sequence_blocks else '<p>No sequence rows available.</p>'}
</body>
</html>
"""
    out_path.write_text(html_text, encoding="utf-8")


def write_fasta(path: Path, records: list[tuple[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            if not sequence:
                continue
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")
            written += 1
    return written


def write_event_fastas(
    event: dict[str, str],
    sequence_rows: list[dict[str, str]],
    event_dir: Path,
) -> tuple[str, str]:
    event_sequences = [row for row in sequence_rows if row["event_id"] == event["event_id"]]
    nt_path = event_dir / "switch_isoforms.nt.fa"
    aa_path = event_dir / "switch_isoforms.aa.fa"
    nt_records = write_fasta(
        nt_path,
        [
            (f"{row['isoform_id']}|{row['switch_role']}|{event['event_id']}", row["nt_sequence"])
            for row in event_sequences
        ],
    )
    aa_records = write_fasta(
        aa_path,
        [
            (f"{row['isoform_id']}|{row['switch_role']}|{event['event_id']}", row["aa_sequence"])
            for row in event_sequences
        ],
    )
    return str(nt_path) if nt_records else "", str(aa_path) if aa_records else ""


def write_selected_fastas(
    outdir: Path,
    sequence_rows: list[dict[str, str]],
) -> tuple[str, str]:
    nt_path = outdir / "switch_selected.nt.fa"
    aa_path = outdir / "switch_selected.aa.fa"
    nt_records = write_fasta(
        nt_path,
        [
            (f"{row['isoform_id']}|{row['switch_role']}|{row['event_id']}", row["nt_sequence"])
            for row in sequence_rows
        ],
    )
    aa_records = write_fasta(
        aa_path,
        [
            (f"{row['isoform_id']}|{row['switch_role']}|{row['event_id']}", row["aa_sequence"])
            for row in sequence_rows
        ],
    )
    return str(nt_path) if nt_records else "", str(aa_path) if aa_records else ""


def pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def write_simple_pdf(path: Path, event_rows: list[dict[str, str]]) -> None:
    """Write a lightweight multi-page PDF index for switch plots.

    This is intentionally dependency-free. The authoritative event plots remain
    the SVG files; the PDF gives reviewers one page per event with the same key
    identifiers and file pointers.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    events = event_rows or [
        {
            "event_id": "no_events",
            "gene_name": "No significant isoform switches",
            "gene_id": "",
            "contrast_id": "",
            "switch_in_isoform": "",
            "switch_out_isoform": "",
            "max_abs_dIF": "",
            "plot_svg": "",
            "event_html": "",
        }
    ]
    objects: list[bytes] = []
    page_refs = []
    for event in events:
        lines = [
            f"Isoform switch: {event.get('gene_name', '')} ({event.get('gene_id', '')})",
            f"Event: {event.get('event_id', '')}",
            f"Contrast: {event.get('contrast_id', '')}",
            f"Switch-in: {event.get('switch_in_isoform', '')}",
            f"Switch-out: {event.get('switch_out_isoform', '')}",
            f"Max abs dIF: {event.get('max_abs_dIF', '')}",
            f"SVG: {event.get('plot_svg', '')}",
            f"HTML: {event.get('event_html', '')}",
        ]
        content_lines = ["BT", "/F1 12 Tf", "50 770 Td"]
        first = True
        for line in lines:
            if not first:
                content_lines.append("0 -22 Td")
            first = False
            content_lines.append(f"({pdf_escape(line)}) Tj")
        content_lines.append("ET")
        stream = "\n".join(content_lines).encode("utf-8")
        content_obj = len(objects) + 1
        objects.append(b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"\nendstream")
        page_obj = len(objects) + 1
        page_refs.append(page_obj)
        objects.append(
            f"<< /Type /Page /Parent 0 0 R /MediaBox [0 0 612 792] "
            f"/Resources << /Font << /F1 0 0 R >> >> /Contents {content_obj} 0 R >>".encode("ascii")
        )
    font_obj = len(objects) + 1
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    pages_obj = len(objects) + 1
    kids = " ".join(f"{page} 0 R" for page in page_refs)
    objects.append(f"<< /Type /Pages /Kids [{kids}] /Count {len(page_refs)} >>".encode("ascii"))
    catalog_obj = len(objects) + 1
    objects.append(f"<< /Type /Catalog /Pages {pages_obj} 0 R >>".encode("ascii"))
    patched = []
    for index, obj in enumerate(objects, start=1):
        text = obj.replace(b"/Parent 0 0 R", f"/Parent {pages_obj} 0 R".encode("ascii"))
        text = text.replace(b"/F1 0 0 R", f"/F1 {font_obj} 0 R".encode("ascii"))
        patched.append(text)
    offsets = []
    output = bytearray(b"%PDF-1.4\n")
    for index, obj in enumerate(patched, start=1):
        offsets.append(len(output))
        output.extend(f"{index} 0 obj\n".encode("ascii"))
        output.extend(obj)
        output.extend(b"\nendobj\n")
    xref_offset = len(output)
    output.extend(f"xref\n0 {len(patched) + 1}\n0000000000 65535 f \n".encode("ascii"))
    for offset in offsets:
        output.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    output.extend(
        f"trailer\n<< /Size {len(patched) + 1} /Root {catalog_obj} 0 R >>\n"
        f"startxref\n{xref_offset}\n%%EOF\n".encode("ascii")
    )
    path.write_bytes(bytes(output))


def external_annotation_candidates(external_dir: Path) -> list[Path]:
    candidates = []
    supported_suffixes = {
        ".tsv",
        ".txt",
        ".domtblout",
        ".out",
        ".gff",
        ".gff3",
    }
    for path in sorted(external_dir.iterdir()):
        if not path.is_file():
            continue
        if path.name.endswith((".stdout.log", ".stderr.log", ".log")):
            continue
        if path.suffix.lower() in supported_suffixes:
            candidates.append(path)
    return candidates


def run_external_tool_commands(
    args: argparse.Namespace,
    nt_fasta: str,
    aa_fasta: str,
    outdir: Path,
) -> tuple[list[dict[str, str]], list[Path]]:
    external_dir = outdir / "external_annotations"
    external_dir.mkdir(parents=True, exist_ok=True)
    commands = [
        ("protein_domain", "interproscan", args.interproscan_command, "aa"),
        ("protein_domain", "pfam", args.pfam_command, "aa"),
        ("coding_potential", "coding_potential", args.coding_potential_command, "nt"),
        ("signal_peptide", "signalp", args.signalp_command, "aa"),
        ("transmembrane", "tm_topology", args.tm_command, "aa"),
        ("localization", "localization", args.localization_command, "aa"),
        ("disorder", "disorder", args.disorder_command, "aa"),
    ]
    rows = []
    fasta_by_kind = {"nt": nt_fasta, "aa": aa_fasta}
    for group, name, template, required_sequence in commands:
        if not template:
            rows.append(
                {
                    "tool_group": group,
                    "tool_name": name,
                    "status": "not_configured",
                    "returncode": "",
                    "command": "",
                    "produced_files": "",
                    "parser_status": "not_run",
                    "parser_detail": "No command output to parse",
                    "stdout_log": "",
                    "stderr_log": "",
                    "detail": "No command template configured",
                }
            )
            continue
        required_fasta = fasta_by_kind[required_sequence]
        if not path_has_content(required_fasta):
            rows.append(
                {
                    "tool_group": group,
                    "tool_name": name,
                    "status": "blocked",
                    "returncode": "",
                    "command": "",
                    "produced_files": "",
                    "parser_status": "not_run",
                    "parser_detail": "Required FASTA missing or empty",
                    "stdout_log": "",
                    "stderr_log": "",
                    "detail": (
                        f"Required {required_sequence.upper()} FASTA is empty or missing. "
                        "Configure rnaseq_differential.isoform_switch_genome_object with "
                        "a genome object matching the reference annotation, then rerun the "
                        "isoform-switch report targets."
                    ),
                }
            )
            continue
        stdout_log = external_dir / f"{safe_id(name)}.stdout.log"
        stderr_log = external_dir / f"{safe_id(name)}.stderr.log"
        command = template.format(
            nt_fasta=nt_fasta,
            aa_fasta=aa_fasta,
            outdir=str(external_dir),
            tool_name=safe_id(name),
        )
        before_files = {path.resolve() for path in external_annotation_candidates(external_dir)}
        completed = subprocess.run(command, shell=True, text=True, capture_output=True, check=False)
        stdout_log.write_text(completed.stdout or "", encoding="utf-8")
        stderr_log.write_text(completed.stderr or "", encoding="utf-8")
        after_files = {path.resolve() for path in external_annotation_candidates(external_dir)}
        produced_files = sorted(str(path) for path in after_files - before_files)
        if completed.returncode == 0 and not produced_files:
            status = "blocked"
            detail = "Command completed but produced no supported annotation output file"
            parser_status = "not_run"
            parser_detail = detail
        elif completed.returncode == 0:
            status = "ok"
            detail = ""
            parser_status = "pending"
            parser_detail = "Parser status is recorded in annotation QA"
        else:
            status = "failed"
            detail = "Command exited non-zero"
            parser_status = "not_run"
            parser_detail = "Command failed before parsing"
        rows.append(
            {
                "tool_group": group,
                "tool_name": name,
                "status": status,
                "returncode": str(completed.returncode),
                "command": command,
                "produced_files": ",".join(produced_files),
                "parser_status": parser_status,
                "parser_detail": parser_detail,
                "stdout_log": str(stdout_log),
                "stderr_log": str(stderr_log),
                "detail": detail,
            }
        )
    return rows, external_annotation_candidates(external_dir)


def update_external_tool_parser_status(
    external_tool_rows: list[dict[str, str]],
    annotation_qa_rows: list[dict[str, str]],
) -> None:
    qa_by_path = {str(Path(row.get("source_path", "")).resolve()): row for row in annotation_qa_rows if row.get("source_path")}
    for row in external_tool_rows:
        produced = [part.strip() for part in row.get("produced_files", "").split(",") if part.strip()]
        if row.get("status") != "ok" or not produced:
            continue
        matched_qa = [qa_by_path.get(str(Path(path).resolve())) for path in produced]
        matched_qa = [qa for qa in matched_qa if qa]
        if not matched_qa:
            row["parser_status"] = "not_assessed"
            row["parser_detail"] = "No annotation QA row matched produced file(s)"
            continue
        statuses = [qa.get("status", "") for qa in matched_qa]
        if any(status == "ok" for status in statuses):
            row["parser_status"] = "ok"
        elif any(status == "ok_no_matches" for status in statuses):
            row["parser_status"] = "ok_no_matches"
        elif any(status == "blocked" for status in statuses):
            row["parser_status"] = "blocked"
        elif any(status == "failed" for status in statuses):
            row["parser_status"] = "failed"
        else:
            row["parser_status"] = "unknown"
        row["parser_detail"] = "; ".join(
            f"{Path(qa.get('source_path', '')).name}:{qa.get('status', '')}:{qa.get('reason', '')}"
            for qa in matched_qa
        )


def render_project_html(
    event_rows: list[dict[str, str]],
    coding_switch_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    annotation_qa_rows: list[dict[str, str]],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    def status_badge(status: str) -> str:
        value = status or "unknown"
        return f'<span class="status {html.escape(value)}">{html.escape(value)}</span>'

    def metric_card(label: str, value: object) -> str:
        return f'<div class="metric"><strong>{html.escape(label)}</strong><span>{html.escape(str(value))}</span></div>'

    def table_source_links() -> str:
        links = [
            file_link("candidate table", str(output.parent / "switch_candidates.tsv"), output),
            file_link("event summary", str(output.parent / "switch_event_summary.tsv"), output),
            file_link("coding switch summary", str(output.parent / "coding_switch_summary.tsv"), output),
            file_link("ncRNA switch interpretation", str(output.parent / "ncrna_switch_interpretation.tsv"), output),
            file_link("sequence summary", str(output.parent / "switch_sequence_summary.tsv"), output),
            file_link("functional annotations", str(output.parent / "functional_annotation_summary.tsv"), output),
            file_link("annotation QA", str(output.parent / "functional_annotation_qa.tsv"), output),
            file_link("external tool manifest", str(output.parent / "external_tool_manifest.tsv"), output),
            file_link("plot manifest", str(output.parent / "switch_plot_manifest.tsv"), output),
            file_link("plot PDF", str(output.parent / "switch_plots.pdf"), output),
        ]
        return " | ".join(link for link in links if link)

    def event_sort_key(row: dict[str, str]) -> tuple[int, float, str]:
        rank_text = row.get("switch_rank", "")
        try:
            rank = int(float(rank_text))
        except (TypeError, ValueError):
            rank = 999999
        effect = abs(to_float(row.get("max_abs_dIF", "")) or 0.0)
        return (rank, -effect, row.get("event_id", ""))

    def top_event_cards(rows: list[dict[str, str]]) -> str:
        cards = []
        for row in sorted(rows, key=event_sort_key)[:12]:
            links = " ".join(
                link
                for link in [
                    file_link("event page", row.get("event_html", ""), output),
                    file_link("diagram", row.get("plot_svg", ""), output),
                    file_link("nt FASTA", row.get("event_nt_fasta", ""), output),
                    file_link("aa FASTA", row.get("event_aa_fasta", ""), output),
                ]
                if link
            )
            cards.append(
                '<article class="event-card">'
                f'<h3>{file_link(row.get("event_id", "event"), row.get("event_html", ""), output)}</h3>'
                f'<p><strong>{html.escape(row.get("gene_display") or row.get("gene_name") or row.get("gene_id", ""))}</strong> '
                f'<code>{html.escape(row.get("gene_id", ""))}</code></p>'
                f'<p>{html.escape(row.get("contrast_id", ""))}</p>'
                f'<p>{status_badge(row.get("status", ""))} '
                f'{html.escape(row.get("switch_biotype_class", ""))}; '
                f'max abs dIF {html.escape(row.get("max_abs_dIF", ""))}; '
                f'priority {html.escape(row.get("coding_priority_tier", ""))}</p>'
                f'<p class="asset-links">{links}</p>'
                "</article>"
            )
        if not cards:
            return '<p class="muted">No event cards are available because no switch events passed the configured thresholds.</p>'
        return '<div class="event-grid">' + "".join(cards) + "</div>"

    def contrast_status_table() -> str:
        events_by_contrast = Counter(row.get("contrast_id", "") for row in event_rows)
        body = []
        for row in sorted(manifest_rows, key=lambda item: item.get("contrast_id", "")):
            contrast_id = row.get("contrast_id", "")
            body.append(
                "<tr>"
                f"<td><code>{html.escape(contrast_id)}</code></td>"
                f"<td>{status_badge(row.get('status', ''))}</td>"
                f"<td>{html.escape(str(events_by_contrast.get(contrast_id, 0)))}</td>"
                f"<td>{html.escape(row.get('reason', ''))}</td>"
                f"<td>{file_link('detailed', row.get('detailed', ''), output)}</td>"
                f"<td>{file_link('consequences', row.get('consequences', ''), output)}</td>"
                "</tr>"
            )
        if not body:
            return '<p class="muted">No isoform-switch manifest rows are available.</p>'
        return (
            "<table><thead><tr><th>contrast</th><th>status</th><th>rendered events</th>"
            "<th>reason</th><th>detailed table</th><th>consequence table</th></tr></thead><tbody>"
            + "".join(body)
            + "</tbody></table>"
        )

    def annotation_status_table() -> str:
        if not annotation_qa_rows:
            return '<p class="muted">No optional annotation QA rows were written.</p>'
        body = []
        for row in annotation_qa_rows:
            body.append(
                "<tr>"
                f"<td>{html.escape(row.get('source_name', ''))}</td>"
                f"<td>{html.escape(row.get('source_kind', ''))}</td>"
                f"<td>{status_badge(row.get('status', ''))}</td>"
                f"<td>{html.escape(row.get('input_rows', ''))}</td>"
                f"<td>{html.escape(row.get('parsed_rows', ''))}</td>"
                f"<td>{html.escape(row.get('matched_rows', ''))}</td>"
                f"<td>{html.escape(row.get('unmatched_rows', ''))}</td>"
                f"<td>{html.escape(row.get('mapping_loss_fraction', ''))}</td>"
                f"<td>{html.escape(row.get('unsupported_column_count', ''))}</td>"
                f"<td>{html.escape(row.get('reason', ''))}</td>"
                "</tr>"
            )
        return (
            "<table><thead><tr><th>source</th><th>kind</th><th>status</th>"
            "<th>input</th><th>parsed</th><th>matched</th><th>unmatched</th>"
            "<th>mapping loss</th><th>unsupported columns</th><th>reason</th></tr></thead><tbody>"
            + "".join(body)
            + "</tbody></table>"
        )

    def overview_html() -> str:
        class_counts = Counter(row.get("switch_biotype_class", "unknown") or "unknown" for row in event_rows)
        status_counts = Counter(row.get("status", "unknown") or "unknown" for row in manifest_rows)
        annotation_status_counts = Counter(row.get("status", "unknown") or "unknown" for row in annotation_qa_rows)
        annotation_attention_count = annotation_status_counts.get("blocked", 0) + annotation_status_counts.get("failed", 0)
        source_links = table_source_links()
        metric_html = "".join(
            [
                metric_card("events", len(event_rows)),
                metric_card("coding", class_counts.get("coding", 0)),
                metric_card("noncoding/mixed", class_counts.get("noncoding", 0) + class_counts.get("mixed_coding_noncoding", 0)),
                metric_card("ambiguous", sum(count for key, count in class_counts.items() if key not in {"coding", "noncoding", "mixed_coding_noncoding"})),
                metric_card("contrasts ok", status_counts.get("ok", 0)),
                metric_card("contrasts blocked/failed", status_counts.get("blocked", 0) + status_counts.get("failed", 0)),
                metric_card("annotation sources ok", annotation_status_counts.get("ok", 0)),
                metric_card("annotation not configured", annotation_status_counts.get("not_configured", 0)),
                metric_card("annotation no matches", annotation_status_counts.get("ok_no_matches", 0)),
                metric_card("annotation blocked/failed", annotation_attention_count),
            ]
        )
        source_block = f'<p class="asset-links">{source_links}</p>' if source_links else ""
        return (
            '<section class="overview">'
            '<h2 id="isoform-overview">Isoform-Switch Overview</h2>'
            '<p class="note">This overview is the first page to inspect. It summarizes which contrasts ran, how many candidate switch events were rendered, and which events have direct diagram and sequence links. The detailed tables below remain available for exhaustive review.</p>'
            f'<div class="metrics">{metric_html}</div>'
            f"{source_block}"
            '<h2 id="contrast-status">Contrast Status</h2>'
            f"{contrast_status_table()}"
            '<h2 id="annotation-status">Optional Consequence Annotation Status</h2>'
            '<p class="note">Optional consequence annotations are supportive evidence. A status of not_configured or ok_no_matches does not invalidate the isoform-switch calls; blocked or failed rows identify resource or parser issues to fix before interpreting annotation-derived consequences.</p>'
            f"{annotation_status_table()}"
            '<h2 id="top-candidate-events">Top Candidate Events</h2>'
            f"{top_event_cards(event_rows)}"
            "</section>"
        )

    def coding_priority_rows(rows: list[dict[str, str]]) -> str:
        body = []
        for row in rows:
            link = ""
            if row.get("event_html") and Path(row["event_html"]).exists():
                link = f'<a href="{html.escape(relative(row["event_html"], output))}">{html.escape(row["event_id"])}</a>'
            else:
                link = html.escape(row["event_id"])
            body.append(
                "<tr>"
                f"<td>{link}</td>"
                f"<td>{html.escape(row['contrast_id'])}</td>"
                f"<td>{html.escape(row.get('coding_priority_rank', ''))}</td>"
                f"<td>{html.escape(row.get('coding_priority_score', ''))}</td>"
                f"<td>{html.escape(row.get('coding_priority_tier', ''))}</td>"
                f"<td>{html.escape(row.get('gene_display') or row['gene_name'])}</td>"
                f"<td>{html.escape(row['gene_id'])}</td>"
                f"<td>{html.escape(row.get('coding_priority_reasons', ''))}</td>"
                f"<td>{html.escape(row.get('gained_domain', ''))}</td>"
                f"<td>{html.escape(row.get('lost_domain', ''))}</td>"
                f"<td>{html.escape(row.get('nmd_change', ''))}</td>"
                f"<td>{html.escape(row.get('coding_potential_change', ''))}</td>"
                f"<td>{html.escape(row.get('orf_length_change_aa', ''))}</td>"
                "</tr>"
            )
        if not body:
            return '<tr><td colspan="13">No coding switches with consequence-priority rows.</td></tr>'
        return "".join(body)

    def section_rows(rows: list[dict[str, str]]) -> str:
        body = []
        for row in rows:
            link = ""
            if row.get("event_html") and Path(row["event_html"]).exists():
                link = f'<a href="{html.escape(relative(row["event_html"], output))}">{html.escape(row["event_id"])}</a>'
            else:
                link = html.escape(row["event_id"])
            body.append(
                "<tr>"
                f"<td>{link}</td>"
                f"<td>{html.escape(row['contrast_id'])}</td>"
                f"<td>{html.escape(row.get('switch_rank', ''))}</td>"
                f"<td>{html.escape(row.get('gene_display') or row['gene_name'])}</td>"
                f"<td>{html.escape(row['gene_id'])}</td>"
                f"<td>{html.escape(row.get('gene_biotype', ''))}</td>"
                f"<td>{html.escape(row.get('switch_biotype_class', ''))}</td>"
                f"<td>{html.escape(row.get('switch_interpretation_label', ''))}</td>"
                f"<td>{html.escape(row.get('coding_priority_tier', ''))}</td>"
                f"<td>{html.escape(row.get('coding_priority_reasons', ''))}</td>"
                f"<td>{html.escape(row['switch_in_isoform'])}</td>"
                f"<td>{html.escape(row['switch_out_isoform'])}</td>"
                f"<td>{html.escape(row['max_abs_dIF'])}</td>"
                f"<td>{html.escape(row['best_switch_statistic'])}</td>"
                f"<td>{html.escape(row['n_functional_annotations'])}</td>"
                f"<td>{file_link('nt', row.get('event_nt_fasta', ''), output)} {file_link('aa', row.get('event_aa_fasta', ''), output)}</td>"
                "</tr>"
            )
        if not body:
            return '<tr><td colspan="16">No events in this section.</td></tr>'
        return "".join(body)

    coding_events = [row for row in event_rows if row.get("switch_biotype_class") == "coding"]
    noncoding_events = [row for row in event_rows if row.get("switch_biotype_class") in {"noncoding", "mixed_coding_noncoding"}]
    ambiguous_events = [row for row in event_rows if row.get("switch_biotype_class") not in {"coding", "noncoding", "mixed_coding_noncoding"}]
    if not event_rows:
        ambiguous_events = []
    run_root = output.parents[6] if len(output.parents) > 6 else output.parent
    project = output.parents[3].name if len(output.parents) > 3 else ""
    sidebar = report_shell_open(
        "Report Map",
        [
            report_map_item("Run dashboard", run_root / "index.html"),
            report_map_item(
                "Project",
                children=[
                    report_map_item("Integrated project report", run_root / "projects" / project / "index.html"),
                    report_map_item("Combined project PDF", run_root / "projects" / project / "technical_report.pdf"),
                ],
            ),
            report_map_item(
                "RNA-seq",
                children=[
                    report_map_item("Differential index", output.parents[2] / "reports/index.html"),
                    report_map_item("DTU methods", output.parents[2] / "dtu/dtu_method_manifest.tsv"),
                    report_map_item("GO/Reactome", output.parents[2] / "reports/enrichment/index.html"),
                    report_map_item("RNA-seq technical PDF", output.parents[2] / "reports/technical_report.pdf"),
                ],
            ),
            report_map_item(
                "Isoform switch",
                children=[
                    report_map_item("Overview", "#isoform-overview"),
                    report_map_item("Contrast status", "#contrast-status"),
                    report_map_item("Annotation status", "#annotation-status"),
                    report_map_item("Top candidate events", "#top-candidate-events"),
                    report_map_item("Coding switches", "#coding-switches"),
                    report_map_item("Noncoding/mixed switches", "#noncoding-and-mixed-coding-potential-switches"),
                    report_map_item("Candidate table", output.parent / "switch_candidates.tsv"),
                    report_map_item("Interpretation consensus", output.parent / "isoform_interpretation_consensus.tsv"),
                ],
            ),
        ],
        output.parent,
    )
    sections = [
        ("Coding Switches", coding_events),
        ("Noncoding And Mixed Coding-Potential Switches", noncoding_events),
        ("Ambiguous, Artifact, Or Unclassified Switches", ambiguous_events),
    ]
    tables = []
    for title, rows in sections:
        tables.append(
            f"""
  {f'''
  <h2>Coding Switch Priority Summary</h2>
  <p class="note">This table is a compact ranking of coding-potential events by predicted functional impact. It is intended for prioritization before manual biological interpretation.</p>
  <table>
    <thead>
      <tr>
        <th>event</th><th>contrast</th><th>priority rank</th><th>score</th><th>tier</th>
        <th>gene display</th><th>gene_id</th><th>priority reasons</th>
        <th>gained domains</th><th>lost domains</th><th>NMD change</th>
        <th>coding-potential change</th><th>ORF length change</th>
      </tr>
    </thead>
    <tbody>
      {coding_priority_rows(coding_switch_rows)}
    </tbody>
  </table>
  ''' if title == 'Coding Switches' else ''}
  <h2 id="{html.escape(title.lower().replace(' ', '-'))}">{html.escape(title)}</h2>
  <p class="note">Each row is one isoform-switch event. The table gives the contrast, gene, switch class, leading switch-in/switch-out isoforms, effect size, annotation count, and links to event-level FASTA files when available.</p>
  <table>
    <thead>
      <tr>
        <th>event</th><th>contrast</th><th>rank</th><th>gene display</th><th>gene_id</th>
        <th>gene biotype</th><th>switch class</th><th>interpretation</th>
        <th>coding priority</th><th>priority reasons</th>
        <th>switch-in</th><th>switch-out</th><th>max abs dIF</th>
        <th>best statistic</th><th>annotations</th><th>FASTA</th>
      </tr>
    </thead>
    <tbody>
      {section_rows(rows)}
    </tbody>
  </table>
"""
        )
    if not event_rows:
        blocked_rows = [row for row in manifest_rows if row.get("status") == "blocked"]
        failed_rows = [row for row in manifest_rows if row.get("status") == "failed"]
        if blocked_rows:
            reasons = sorted({row.get("reason", "") for row in blocked_rows if row.get("reason", "")})
            reason_items = "".join(f"<li>{html.escape(reason)}</li>" for reason in reasons) or "<li>No reason recorded.</li>"
            tables = [
                '<p class="note"><strong>Isoform-switch analysis was blocked.</strong> '
                'No switch events were rendered because one or more planned contrasts did not run.</p>'
                f'<ul>{reason_items}</ul>'
                '<p>Install or configure the missing dependency/resource, then rerun the isoform-switch targets.</p>'
            ]
        elif failed_rows:
            reasons = sorted({row.get("reason", "") for row in failed_rows if row.get("reason", "")})
            reason_items = "".join(f"<li>{html.escape(reason)}</li>" for reason in reasons) or "<li>No reason recorded.</li>"
            tables = [
                '<p class="note"><strong>Isoform-switch analysis failed.</strong> '
                'No switch events were rendered because one or more planned contrasts failed.</p>'
                f'<ul>{reason_items}</ul>'
            ]
        else:
            tables = ['<p>No significant isoform-switch events passed the configured thresholds.</p>']
    output.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Isoform-switch report</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1680px; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; vertical-align: top; overflow-wrap: anywhere; }}
    th {{ background: #f6f8fa; }}
    code {{ white-space: normal; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 0.75rem; margin: 1rem 0; }}
    .metric {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    .metric strong {{ color: #57606a; display: block; font-size: 0.85rem; }}
    .metric span {{ display: block; font-size: 1.25rem; font-weight: 700; margin-top: 0.25rem; }}
    .event-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 0.85rem; }}
    .event-card {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.8rem; }}
    .event-card h3 {{ margin: 0 0 0.35rem; font-size: 1rem; }}
    .event-card p {{ margin: 0.25rem 0; }}
    .asset-links {{ color: #57606a; display: flex; flex-wrap: wrap; gap: 0.35rem 0.45rem; font-size: 0.92rem; }}
    .asset-links a {{ background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 4px; display: inline-block; line-height: 1.25; padding: 0.16rem 0.42rem; white-space: nowrap; }}
    .muted {{ color: #57606a; }}
    .status {{ font-weight: 700; }}
    .status.ok {{ color: #1a7f37; }}
    .status.not_configured {{ color: #57606a; }}
    .status.ok_no_matches {{ color: #57606a; }}
    .status.blocked {{ color: #9a6700; }}
    .status.failed {{ color: #cf222e; }}
    table {{ display: block; overflow-x: auto; }}
    nav.breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    {report_map_css()}
  </style>
</head>
<body>
  {sidebar}
  <nav class="breadcrumbs">ASPIS / RNA-seq / Isoform-switch overview</nav>
  <h1>Isoform-switch report</h1>
  <p>Events are ranked by absolute isoform fraction change and split into coding, noncoding/mixed, and ambiguous sections. Coding switches are additionally prioritized by predicted functional consequences. Noncoding switches are interpreted through transcript architecture rather than requiring ORF/domain evidence.</p>
  <p class="note">Use this page as an index of candidate events. The event links open per-gene pages with exon diagrams, switch-in/switch-out isoforms, coding or noncoding interpretation tables, optional functional annotations, and sequence extracts.</p>
  <p class="note">Isoform-switch calls are candidate transcript-usage changes. Prioritize events by effect size, replicate support, annotation quality, and biological plausibility; optional consequence annotation is supportive evidence, not a required condition for a valid switch.</p>
  {overview_html()}
  {''.join(tables)}
  {report_shell_close()}
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


def clean_stale_report_outputs(outdir: Path) -> None:
    events_dir = outdir / "events"
    if events_dir.exists():
        shutil.rmtree(events_dir)
    external_dir = outdir / "external_annotations"
    if external_dir.exists():
        shutil.rmtree(external_dir)
    for stale_name in ["switch_selected.nt.fa", "switch_selected.aa.fa"]:
        stale_path = outdir / stale_name
        if stale_path.exists():
            stale_path.unlink()


def main() -> int:
    args = parse_args()
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")
    outdir = Path(args.outdir)
    clean_stale_report_outputs(outdir)
    _, manifest_rows = read_table(Path(args.manifest), MANIFEST_REQUIRED)
    metadata = metadata_by_transcript(Path(args.transcript_metadata))
    gtf_models = parse_gtf(Path(args.annotated_gtf))

    candidate_rows, event_rows, event_context, nt_sequences, aa_sequences = build_events(
        manifest_rows,
        metadata,
        gtf_models,
        args.padj,
        args.dif,
        args.top_n,
    )
    attach_consequences(candidate_rows, event_rows, manifest_rows)
    sequence_rows = build_sequence_rows(candidate_rows, event_rows, event_context, gtf_models, nt_sequences, aa_sequences)
    ncrna_annotations = read_ncrna_annotation_rows(args.ncrna_annotation_tables)
    ncrna_switch_rows = build_ncrna_switch_rows(candidate_rows, sequence_rows, gtf_models, ncrna_annotations)

    event_by_isoform: dict[str, list[dict[str, str]]] = defaultdict(list)
    event_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidate_rows:
        event_stub = {
            "event_id": row["event_id"],
            "contrast_id": row["contrast_id"],
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
        }
        event_by_isoform[row["isoform_id"]].append(event_stub)
        event_by_gene[row["gene_id"]].append(event_stub)
    for gene_id, events in list(event_by_gene.items()):
        deduped = []
        seen_event_ids = set()
        for event in events:
            if event["event_id"] in seen_event_ids:
                continue
            seen_event_ids.add(event["event_id"])
            deduped.append(event)
        event_by_gene[gene_id] = deduped
    selected_nt_fasta, selected_aa_fasta = write_selected_fastas(outdir, sequence_rows)
    external_tool_rows, generated_annotation_paths = run_external_tool_commands(
        args,
        selected_nt_fasta,
        selected_aa_fasta,
        outdir,
    )
    annotation_path_text = ",".join(
        [args.functional_annotation_tables]
        + [str(path) for path in generated_annotation_paths]
    )

    annotation_rows, annotation_qa_rows = normalized_annotations(
        annotation_path_text,
        {row["isoform_id"] for row in candidate_rows},
        event_by_isoform,
        event_by_gene,
    )
    update_external_tool_parser_status(external_tool_rows, annotation_qa_rows)
    annotate_feature_changes(annotation_rows, sequence_rows)
    annotation_counts = defaultdict(int)
    for row in annotation_rows:
        annotation_counts[row["event_id"]] += 1
    for event in event_rows:
        event["n_functional_annotations"] = str(annotation_counts[event["event_id"]])

    coding_switch_rows = build_coding_switch_rows(event_rows, candidate_rows, sequence_rows, annotation_rows)

    plot_rows = []
    for event in event_rows:
        event_dir = outdir / "events" / safe_id(event["event_id"])
        svg_path = event_dir / "switch.svg"
        html_path = event_dir / "index.html"
        nt_fasta, aa_fasta = write_event_fastas(event, sequence_rows, event_dir)
        event["event_nt_fasta"] = nt_fasta
        event["event_aa_fasta"] = aa_fasta
        event_annotations = [row for row in annotation_rows if row["event_id"] == event["event_id"]]
        render_event_svg(event, event_context[event["event_id"]], gtf_models, event_annotations, svg_path)
        render_event_html(
            event,
            candidate_rows,
            sequence_rows,
            annotation_rows,
            ncrna_switch_rows,
            coding_switch_rows,
            html_path,
            svg_path,
        )
        event["plot_svg"] = str(svg_path)
        event["event_html"] = str(html_path)
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
                "nt_fasta": nt_fasta,
                "aa_fasta": aa_fasta,
                "plots_pdf": args.plots_pdf,
            }
        )

    event_by_id = {row["event_id"]: row for row in event_rows}
    for row in coding_switch_rows:
        event = event_by_id.get(row["event_id"], {})
        row["event_html"] = event.get("event_html", "")
        row["n_functional_annotations"] = event.get("n_functional_annotations", row.get("n_functional_annotations", "0"))

    write_simple_pdf(Path(args.plots_pdf), event_rows)
    write_table(Path(args.candidate_table), CANDIDATE_COLUMNS, candidate_rows)
    write_table(Path(args.event_summary), EVENT_COLUMNS, event_rows)
    write_table(Path(args.ncrna_switch_table), NCRNA_SWITCH_COLUMNS, ncrna_switch_rows)
    write_table(Path(args.coding_switch_summary), CODING_SWITCH_COLUMNS, coding_switch_rows)
    write_table(Path(args.sequence_table), SEQUENCE_COLUMNS, sequence_rows)
    write_table(Path(args.functional_annotation_table), ANNOTATION_COLUMNS, annotation_rows)
    write_table(Path(args.functional_annotation_qa), ANNOTATION_QA_COLUMNS, annotation_qa_rows)
    write_table(Path(args.plot_manifest), PLOT_MANIFEST_COLUMNS, plot_rows)
    write_table(Path(args.external_tool_manifest), EXTERNAL_TOOL_COLUMNS, external_tool_rows)
    render_project_html(event_rows, coding_switch_rows, manifest_rows, annotation_qa_rows, Path(args.html))
    write_done(Path(args.done), event_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
