"""Shared display-label helpers for report-facing feature identifiers."""

from __future__ import annotations

import re


MISSING_DISPLAY_VALUES = {"", "na", "n/a", "none", "null", "nan", "."}


def clean_display_value(value: object) -> str:
    cleaned = str(value or "").strip()
    if cleaned.lower() in MISSING_DISPLAY_VALUES:
        return ""
    missing_wrapped = re.fullmatch(r"(?i)(?:na|n/a|none|null|nan)\s*\(([^)]+)\)", cleaned)
    if missing_wrapped:
        return clean_display_value(missing_wrapped.group(1))
    return cleaned


def first_clean(row: dict[str, str], columns: list[str]) -> str:
    for column in columns:
        cleaned = clean_display_value(row.get(column, ""))
        if cleaned:
            return cleaned
    return ""


def gene_display_label(gene_id: str, gene_name: str) -> str:
    gene_id = clean_display_value(gene_id)
    gene_name = clean_display_value(gene_name)
    if gene_name and gene_id and gene_name != gene_id:
        return f"{gene_name} ({gene_id})"
    return gene_name or gene_id


def transcript_display_label(transcript_id: str, gene_id: str, gene_name: str) -> str:
    transcript_id = clean_display_value(transcript_id)
    gene_label = gene_display_label(gene_id, gene_name)
    if gene_label and transcript_id:
        return f"{gene_label} | {transcript_id}"
    return transcript_id or gene_label


def feature_display_label(row: dict[str, str], feature_id_column: str = "") -> str:
    existing = first_clean(row, ["feature_display", "transcript_display", "gene_display"])
    if existing:
        return existing

    feature_id = clean_display_value(row.get(feature_id_column, "")) if feature_id_column else ""
    transcript_id = first_clean(row, ["transcript_id", "isoform_id"])
    if feature_id_column in {"transcript_id", "isoform_id"} and feature_id:
        transcript_id = transcript_id or feature_id

    gene_id = first_clean(row, ["gene_id", "Geneid", "gene"])
    if feature_id_column == "Geneid" and feature_id:
        gene_id = gene_id or feature_id
    gene_name = first_clean(row, ["gene_name", "GeneName", "gene_symbol", "symbol"])

    if transcript_id:
        return transcript_display_label(transcript_id, gene_id, gene_name)

    gene_label = gene_display_label(gene_id, gene_name)
    if gene_label:
        return gene_label
    return feature_id or first_clean(row, ["feature_id", "event_id", "id"])
