#!/usr/bin/env python3
"""Build a lightweight ASPIS run provenance bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import yaml


MANIFEST_COLUMNS = [
    "label",
    "path",
    "status",
    "size_bytes",
    "sha256",
    "rows",
    "columns",
    "detail",
]
CONTEXT_COLUMNS = ["section", "metric", "value", "source"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--config-snapshot", required=True)
    parser.add_argument("--intake-snapshot", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--assay", required=True, choices=["rnaseq", "smallrna"])
    parser.add_argument("--project", required=True)
    parser.add_argument("--configfile", default="")
    parser.add_argument("--intake", default="")
    parser.add_argument("--preflight-report", default="")
    parser.add_argument("--artifacts", nargs="*", default=[])
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def table_shape(path: Path) -> tuple[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return "", ""
    try:
        with path.open(encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            rows = sum(1 for _ in reader)
    except UnicodeDecodeError:
        return "", ""
    except OSError:
        return "", ""
    return str(rows), str(len(header)) if header else ""


def manifest_row(label: str, path: Path, detail: str = "") -> dict[str, str]:
    if not path:
        return {
            "label": label,
            "path": "",
            "status": "not_configured",
            "size_bytes": "",
            "sha256": "",
            "rows": "",
            "columns": "",
            "detail": detail,
        }
    if not path.is_file():
        return {
            "label": label,
            "path": str(path),
            "status": "missing",
            "size_bytes": "",
            "sha256": "",
            "rows": "",
            "columns": "",
            "detail": detail,
        }
    rows, columns = table_shape(path) if path.suffix.lower() in {".tsv", ".txt", ".csv"} else ("", "")
    return {
        "label": label,
        "path": str(path),
        "status": "present",
        "size_bytes": str(path.stat().st_size),
        "sha256": sha256_file(path),
        "rows": rows,
        "columns": columns,
        "detail": detail,
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def copy_snapshot(source: str, dest: Path) -> Path | None:
    if not source:
        return None
    source_path = Path(source)
    if not source_path.is_file():
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_path, dest)
    return dest


def load_yaml(path: Path | None) -> dict:
    if path is None or not path.is_file():
        return {}
    with path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    return loaded if isinstance(loaded, dict) else {}


def add_context(rows: list[dict[str, str]], section: str, metric: str, value: object, source: str) -> None:
    if isinstance(value, (list, tuple)):
        value = ",".join(str(item) for item in value)
    rows.append(
        {
            "section": section,
            "metric": metric,
            "value": str(value if value is not None else ""),
            "source": source,
        }
    )


def add_file_context(rows: list[dict[str, str]], section: str, metric: str, path_text: str, source: str) -> None:
    if not path_text:
        add_context(rows, section, f"{metric}_status", "not_configured", source)
        return
    path = Path(path_text)
    add_context(rows, section, f"{metric}_path", path_text, source)
    add_context(rows, section, f"{metric}_status", "present" if path.exists() else "missing", source)
    if path.exists() and path.is_file():
        add_context(rows, section, f"{metric}_size_bytes", path.stat().st_size, source)


def table_by_name(artifacts: list[Path], suffix: str) -> Path | None:
    matches = [path for path in artifacts if str(path).endswith(suffix)]
    return matches[0] if matches else None


def summarize_design(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    samples = read_tsv(table_by_name(artifacts, "/samples.tsv") or Path(""))
    design = read_tsv(table_by_name(artifacts, "/design.tsv") or Path(""))
    if samples:
        add_context(rows, "design", "library_count", len(samples), "samples.tsv")
        if "layout" in samples[0]:
            add_context(rows, "design", "layouts", dict_count(samples, "layout"), "samples.tsv")
        if "condition" in samples[0]:
            add_context(rows, "design", "condition_counts", dict_count(samples, "condition"), "samples.tsv")
    if design:
        add_context(rows, "design", "design_rows", len(design), "design.tsv")
        add_context(rows, "design", "design_columns", list(design[0].keys()), "design.tsv")
        if "condition" in design[0]:
            add_context(rows, "design", "design_condition_counts", dict_count(design, "condition"), "design.tsv")
        for key in ("model_formula", "covariates", "batch_factors", "blocking_factors", "interaction_terms", "contrast_by"):
            if key in design[0]:
                add_context(rows, "design", key, design[0].get(key, ""), "design.tsv")


def dict_count(records: list[dict[str, str]], key: str) -> str:
    counts = Counter(record.get(key, "") or "missing" for record in records)
    return ",".join(f"{name}:{count}" for name, count in sorted(counts.items()))


def summarize_plan(rows: list[dict[str, str]], artifacts: list[Path], project: str, assay: str) -> None:
    plan = read_tsv(table_by_name(artifacts, "analysis_plan.tsv") or Path(""))
    for record in plan:
        if record.get("project") == project and record.get("assay") == assay:
            for key in ("status", "n_libraries", "libraries", "condition_col", "control_label"):
                if key in record:
                    add_context(rows, "analysis_plan", key, record.get(key, ""), "analysis_plan.tsv")
            break
    for suffix, section in [
        ("alignment_plan.tsv", "alignment"),
        ("quantification_plan.tsv", "quantification"),
        ("differential_plan.tsv", "differential"),
        ("smallrna_plan.tsv", "smallrna"),
        ("report_plan.tsv", "report"),
    ]:
        path = table_by_name(artifacts, suffix)
        records = read_tsv(path or Path(""))
        if records:
            add_context(rows, section, "plan_rows", len(records), str(path))
            if "status" in records[0]:
                add_context(rows, section, "status_counts", dict_count(records, "status"), str(path))
            if "level" in records[0]:
                add_context(rows, section, "levels", sorted({record.get("level", "") for record in records}), str(path))


def summarize_counts(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    for suffix, label in [
        ("gene_counts.tsv", "gene_counts"),
        ("transcript_counts.tsv", "transcript_counts"),
        ("mirna_counts.tsv", "mirna_counts"),
    ]:
        path = table_by_name(artifacts, suffix)
        if not path:
            continue
        n_rows, n_cols = table_shape(path)
        add_context(rows, "quantification", f"{label}_features", n_rows, str(path))
        add_context(rows, "quantification", f"{label}_columns", n_cols, str(path))
    for path in artifacts:
        if not str(path).endswith("metadata.tsv"):
            continue
        records = read_tsv(path)
        if not records:
            continue
        add_context(rows, "annotation", f"{path.name}_rows", len(records), str(path))
        for key in ("feature_type", "gene_type", "class_code"):
            if key in records[0]:
                add_context(rows, "annotation", f"{path.name}_{key}_counts", dict_count(records, key), str(path))


def summarize_sample_qc(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    metrics = read_tsv(table_by_name(artifacts, "sample_qc_metrics.tsv") or Path(""))
    if metrics:
        add_context(rows, "sample_qc", "sample_count", len(metrics), "sample_qc_metrics.tsv")
        add_context(rows, "sample_qc", "condition_counts", dict_count(metrics, "condition"), "sample_qc_metrics.tsv")
        total_counts = 0
        detected = []
        for record in metrics:
            try:
                total_counts += int(float(record.get("library_size", "0") or "0"))
            except ValueError:
                pass
            try:
                detected.append(int(float(record.get("detected_features", "0") or "0")))
            except ValueError:
                pass
        add_context(rows, "sample_qc", "total_counts", total_counts, "sample_qc_metrics.tsv")
        if detected:
            add_context(rows, "sample_qc", "detected_features_range", f"{min(detected)}-{max(detected)}", "sample_qc_metrics.tsv")
    correlations = read_tsv(table_by_name(artifacts, "sample_correlations.tsv") or Path(""))
    if correlations:
        values = []
        for record in correlations:
            if record.get("sample_a") == record.get("sample_b"):
                continue
            try:
                values.append(float(record.get("pearson_log_cpm", "")))
            except ValueError:
                pass
        if values:
            add_context(rows, "sample_qc", "correlation_range", f"{min(values):.4g}-{max(values):.4g}", "sample_correlations.tsv")


def summarize_differential(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    for path in artifacts:
        if not path.name.endswith("manifest.tsv"):
            continue
        if "deseq2" not in str(path) and "isoform_switch" not in str(path):
            continue
        records = read_tsv(path)
        if not records:
            continue
        section = "isoform_switch" if "isoform_switch" in str(path) else "deseq2"
        add_context(rows, section, f"{path.parent.name}_contrast_count", len(records), str(path))
        if "status" in records[0]:
            add_context(rows, section, f"{path.parent.name}_status_counts", dict_count(records, "status"), str(path))


def summarize_smallrna_residual(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    biotypes = read_tsv(table_by_name(artifacts, "biotype_counts.tsv") or Path(""))
    if biotypes:
        totals: defaultdict[str, int] = defaultdict(int)
        for record in biotypes:
            biotype = record.get("biotype") or record.get("feature_type") or record.get("gene_type") or "unclassified"
            count_text = record.get("count") or record.get("reads") or record.get("read_count") or "0"
            try:
                totals[biotype] += int(float(count_text))
            except ValueError:
                pass
        ordered = sorted(totals.items(), key=lambda item: (-item[1], item[0]))
        add_context(rows, "smallrna_residual", "biotype_count_rows", len(biotypes), "biotype_counts.tsv")
        add_context(
            rows,
            "smallrna_residual",
            "top_biotypes",
            ",".join(f"{name}:{count}" for name, count in ordered[:10]),
            "biotype_counts.tsv",
        )
    features = read_tsv(table_by_name(artifacts, "feature_counts.tsv") or Path(""))
    if features:
        add_context(rows, "smallrna_residual", "feature_count_rows", len(features), "feature_counts.tsv")


def summarize_strandedness(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    report = read_tsv(table_by_name(artifacts, "strandedness_report.tsv") or Path(""))
    if not report:
        return
    add_context(rows, "strandedness", "sample_count", len(report), "strandedness_report.tsv")
    if "inferred_strandness" in report[0]:
        add_context(rows, "strandedness", "inferred_counts", dict_count(report, "inferred_strandness"), "strandedness_report.tsv")
    warnings = [record for record in report if record.get("warning")]
    add_context(rows, "strandedness", "warning_count", len(warnings), "strandedness_report.tsv")
    if warnings:
        add_context(rows, "strandedness", "warnings", "; ".join(record.get("warning", "") for record in warnings[:5]), "strandedness_report.tsv")


def summarize_biotypes(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    counts = read_tsv(table_by_name(artifacts, "count_biotype_summary.tsv") or Path(""))
    if counts:
        add_context(rows, "biotypes", "count_summary_rows", len(counts), "count_biotype_summary.tsv")
        for level in sorted({record.get("level", "") for record in counts if record.get("level")}):
            selected = [record for record in counts if record.get("level") == level]
            ordered = sorted(
                selected,
                key=lambda record: -safe_float(record.get("total_count", "")),
            )
            add_context(
                rows,
                "biotypes",
                f"{level}_top_total_counts",
                ",".join(f"{record.get('biotype', 'unclassified')}:{record.get('total_count', '0')}" for record in ordered[:10]),
                "count_biotype_summary.tsv",
            )
    differential = read_tsv(table_by_name(artifacts, "differential_biotype_summary.tsv") or Path(""))
    if differential:
        add_context(rows, "biotypes", "differential_summary_rows", len(differential), "differential_biotype_summary.tsv")
        add_context(rows, "biotypes", "differential_levels", sorted({record.get("level", "") for record in differential}), "differential_biotype_summary.tsv")


def summarize_mirna_mrna(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    manifest = read_tsv(table_by_name(artifacts, "mirna_mrna_manifest.tsv") or Path(""))
    if not manifest:
        return
    add_context(rows, "mirna_mrna", "contrast_count", len(manifest), "mirna_mrna_manifest.tsv")
    if "status" in manifest[0]:
        add_context(rows, "mirna_mrna", "status_counts", dict_count(manifest, "status"), "mirna_mrna_manifest.tsv")
    for key in ("n_pairs", "n_inverse_pairs", "n_anticorrelated_pairs"):
        if key in manifest[0]:
            add_context(rows, "mirna_mrna", key, sum(int(safe_float(record.get(key, ""))) for record in manifest), "mirna_mrna_manifest.tsv")


def summarize_dtu(rows: list[dict[str, str]], artifacts: list[Path]) -> None:
    plan = read_tsv(table_by_name(artifacts, "dtu_plan.tsv") or Path(""))
    if not plan:
        return
    add_context(rows, "dtu", "plan_rows", len(plan), "dtu_plan.tsv")
    for key in ("method", "status", "reason", "candidate_methods"):
        if key in plan[0]:
            add_context(rows, "dtu", key, plan[0].get(key, ""), "dtu_plan.tsv")


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def summarize_config(rows: list[dict[str, str]], config: dict, assay: str) -> None:
    design = config.get("design", {}) or {}
    add_context(rows, "config", "condition_col", design.get("condition_col", ""), "config")
    add_context(rows, "config", "control_label", design.get("control_label", ""), "config")
    add_context(rows, "config", "covariates", design.get("covariates", []), "config")
    if assay == "rnaseq":
        alignment = config.get("rnaseq_alignment", {}) or {}
        quant = config.get("rnaseq_quantification", {}) or {}
        diff = config.get("rnaseq_differential", {}) or {}
        add_context(rows, "config", "rnaseq_aligner", alignment.get("aligner", ""), "config")
        add_context(rows, "config", "rnaseq_reference_fasta", alignment.get("reference_fasta", ""), "config")
        add_context(rows, "config", "rnaseq_annotation_gtf", alignment.get("annotation_gtf", ""), "config")
        add_context(rows, "config", "infer_strandedness", alignment.get("infer_strandedness", ""), "config")
        add_context(rows, "config", "rnaseq_design_formula", diff.get("design_formula", config.get("design", {}).get("model_formula", "")), "config")
        add_context(rows, "config", "transcriptome_mode", quant.get("transcriptome_mode", ""), "config")
        add_context(rows, "config", "biotype_summary", quant.get("biotype_summary", ""), "config")
        add_context(rows, "config", "differential_levels", diff.get("levels", []), "config")
        add_context(rows, "config", "contrast_by", diff.get("contrast_by", []), "config")
        add_context(rows, "config", "report_feature_sets", diff.get("report_feature_sets", ""), "config")
        add_context(rows, "config", "report_feature_set_tables", diff.get("report_feature_set_tables", ""), "config")
        dtu = config.get("rnaseq_dtu", {}) or {}
        add_context(rows, "config", "dtu_run", dtu.get("run", ""), "config")
        add_context(rows, "config", "dtu_candidate_methods", dtu.get("candidate_methods", ""), "config")
        add_file_context(rows, "reference", "rnaseq_reference_fasta", alignment.get("reference_fasta", ""), "config")
        add_file_context(rows, "reference", "rnaseq_annotation_gtf", alignment.get("annotation_gtf", ""), "config")
    if assay == "smallrna":
        smallrna = config.get("smallrna", {}) or {}
        add_context(rows, "config", "adapter", smallrna.get("adapter", ""), "config")
        add_context(rows, "config", "length_window", f"{smallrna.get('min_length', '')}-{smallrna.get('max_length', '')}", "config")
        add_context(rows, "config", "depletion_run", smallrna.get("depletion_run", ""), "config")
        add_context(rows, "config", "contaminant_fasta", smallrna.get("contaminant_fasta", ""), "config")
        add_context(rows, "config", "contaminant_mismatches", smallrna.get("contaminant_mismatches", ""), "config")
        add_context(rows, "config", "alignment_mismatches", smallrna.get("alignment_mismatches", ""), "config")
        add_context(rows, "config", "mirbase_species_prefix", smallrna.get("mirbase_species_prefix", ""), "config")
        add_context(rows, "config", "mirbase_fasta", smallrna.get("mirbase_fasta", ""), "config")
        add_context(rows, "config", "residual_run", smallrna.get("residual_run", ""), "config")
        add_context(rows, "config", "residual_genome_fasta", smallrna.get("residual_genome_fasta", ""), "config")
        add_context(rows, "config", "residual_annotation_gtf", smallrna.get("residual_annotation_gtf", ""), "config")
        add_context(rows, "config", "residual_mismatches", smallrna.get("residual_mismatches", ""), "config")
        add_context(rows, "config", "target_enrichment_mode", smallrna.get("target_enrichment_mode", ""), "config")
        add_context(rows, "config", "target_table", smallrna.get("target_table", ""), "config")
        add_context(rows, "config", "target_tables", smallrna.get("target_tables", ""), "config")
        add_context(rows, "config", "target_feature_sets", smallrna.get("target_feature_sets", ""), "config")
        add_context(rows, "config", "target_feature_set_tables", smallrna.get("target_feature_set_tables", ""), "config")
        integration = config.get("mirna_mrna_integration", {}) or {}
        add_context(rows, "config", "mirna_mrna_integration_run", integration.get("run", ""), "config")
        add_context(rows, "config", "mirna_mrna_match_columns", integration.get("match_columns", ""), "config")
        add_file_context(rows, "reference", "mirbase_fasta", smallrna.get("mirbase_fasta", ""), "config")
        add_file_context(rows, "reference", "contaminant_fasta", smallrna.get("contaminant_fasta", ""), "config")
        add_file_context(rows, "reference", "residual_genome_fasta", smallrna.get("residual_genome_fasta", ""), "config")
        add_file_context(rows, "reference", "residual_annotation_gtf", smallrna.get("residual_annotation_gtf", ""), "config")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    config_snapshot = copy_snapshot(args.configfile, Path(args.config_snapshot))
    intake_snapshot = copy_snapshot(args.intake, Path(args.intake_snapshot))

    artifact_paths = [Path(path) for path in args.artifacts if str(path).strip()]
    if args.preflight_report:
        artifact_paths.append(Path(args.preflight_report))
    if config_snapshot:
        artifact_paths.append(config_snapshot)
    if intake_snapshot:
        artifact_paths.append(intake_snapshot)

    seen = set()
    unique_artifacts = []
    for path in artifact_paths:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        unique_artifacts.append(path)

    manifest_rows = [
        manifest_row(path.name or "artifact", path)
        for path in unique_artifacts
    ]

    context_rows: list[dict[str, str]] = []
    add_context(context_rows, "run", "generated_at_utc", datetime.now(timezone.utc).isoformat(), "provenance")
    add_context(context_rows, "run", "assay", args.assay, "provenance")
    add_context(context_rows, "run", "project", args.project, "provenance")
    add_context(context_rows, "run", "artifact_count", len(unique_artifacts), "provenance")
    add_context(context_rows, "run", "configfile", args.configfile, "provenance")
    add_context(context_rows, "run", "intake", args.intake, "provenance")
    if args.preflight_report:
        add_context(context_rows, "run", "preflight_report", args.preflight_report, "provenance")

    summarize_config(context_rows, load_yaml(config_snapshot), args.assay)
    summarize_design(context_rows, unique_artifacts)
    summarize_plan(context_rows, unique_artifacts, args.project, args.assay)
    summarize_counts(context_rows, unique_artifacts)
    summarize_sample_qc(context_rows, unique_artifacts)
    summarize_differential(context_rows, unique_artifacts)
    summarize_strandedness(context_rows, unique_artifacts)
    summarize_biotypes(context_rows, unique_artifacts)
    summarize_mirna_mrna(context_rows, unique_artifacts)
    summarize_dtu(context_rows, unique_artifacts)
    if args.assay == "smallrna":
        summarize_smallrna_residual(context_rows, unique_artifacts)

    write_tsv(Path(args.manifest), manifest_rows, MANIFEST_COLUMNS)
    write_tsv(Path(args.summary), context_rows, CONTEXT_COLUMNS)
    write_tsv(
        Path(args.done),
        [
            {
                "status": "ok",
                "artifact_count": str(len(unique_artifacts)),
                "context_rows": str(len(context_rows)),
            }
        ],
        ["status", "artifact_count", "context_rows"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
