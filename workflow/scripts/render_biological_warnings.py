#!/usr/bin/env python3
"""Aggregate biological and design warnings into one branch-level report."""

from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


WARNING_COLUMNS = ["assay", "project", "severity", "category", "item", "message", "source"]
MANIFEST_COLUMNS = ["resource", "status", "path", "rows", "detail"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assay", required=True, choices=("rnaseq", "smallrna"))
    parser.add_argument("--project", required=True)
    parser.add_argument("--design", default="")
    parser.add_argument("--sample-qc-metrics", default="")
    parser.add_argument("--sample-correlations", default="")
    parser.add_argument("--strandedness-report", default="")
    parser.add_argument("--biotype-count-summary", default="")
    parser.add_argument("--biotype-differential-summary", default="")
    parser.add_argument("--transcript-discovery-summary", default="")
    parser.add_argument("--transcript-discovery-differential-summary", default="")
    parser.add_argument(
        "--deseq2-manifest",
        action="append",
        default=[],
        help="Optional DESeq2 contrast manifest TSV; may be repeated.",
    )
    parser.add_argument("--residual-manifest", default="")
    parser.add_argument("--residual-biotype-counts", default="")
    parser.add_argument("--length-stage-summary", default="")
    parser.add_argument("--arm-summary", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--warnings", required=True)
    parser.add_argument("--summary-html", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--min-detected-features", type=int, default=10)
    parser.add_argument("--min-library-size", type=int, default=100)
    parser.add_argument("--min-sample-correlation", type=float, default=0.6)
    parser.add_argument("--max-unclassified-biotype-fraction", type=float, default=0.5)
    parser.add_argument("--max-true-novel-transcript-fraction", type=float, default=0.2)
    parser.add_argument("--warn-high-true-novel-transcript-fraction", action="store_true")
    parser.add_argument("--max-residual-genome-fraction", type=float, default=0.5)
    parser.add_argument("--min-deseq2-replicates", type=int, default=2)
    parser.add_argument("--min-deseq2-tested-features", type=int, default=10)
    return parser.parse_args()


def read_table(path_text: str, required: set[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path_text:
        return [], []
    path = Path(path_text)
    if not path.exists():
        return [], []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return [], []
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return list(reader.fieldnames), [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_table(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def parse_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return 0.0


def parse_int(value: str) -> int:
    return int(parse_float(value))


def add_warning(rows: list[dict[str, str]], args: argparse.Namespace, severity: str, category: str, item: str, message: str, source: str) -> None:
    rows.append(
        {
            "assay": args.assay,
            "project": args.project,
            "severity": severity,
            "category": category,
            "item": item,
            "message": message,
            "source": source,
        }
    )


def split_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def simple_design_terms(formula: str) -> list[str]:
    text = formula.strip()
    if text.startswith("~"):
        text = text[1:]
    terms: list[str] = []
    for raw_term in text.split("+"):
        term = raw_term.strip()
        if not term or term in {"0", "1"}:
            continue
        if any(operator in term for operator in [":", "*", "/", "(", ")"]):
            continue
        terms.append(term)
    return terms


def count_matrix_qc(path_text: str, samples: list[str]) -> dict[str, dict[str, int]]:
    if not path_text or not samples:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    qc = {sample: {"library_size": 0, "detected_features": 0} for sample in samples}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return qc
        available = [sample for sample in samples if sample in reader.fieldnames]
        for row in reader:
            for sample in available:
                value = parse_float(row.get(sample, ""))
                qc[sample]["library_size"] += int(value)
                if value > 0:
                    qc[sample]["detected_features"] += 1
    return qc


def coldata_by_sample(path_text: str, samples: list[str]) -> dict[str, dict[str, str]]:
    if not path_text:
        return {}
    _columns, rows = read_table(path_text)
    if not rows:
        return {}
    id_column = "library_id" if "library_id" in rows[0] else ""
    if not id_column:
        for candidate in ["sample", "sample_id", "Sample", "SampleID"]:
            if candidate in rows[0]:
                id_column = candidate
                break
    if not id_column:
        return {}
    wanted = set(samples)
    return {row[id_column]: row for row in rows if row.get(id_column, "") in wanted}


def design_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    _columns, rows = read_table(args.design)
    for row in rows:
        if row.get("differential_status") == "blocked":
            add_warning(
                warnings,
                args,
                "error",
                "design",
                row.get("condition", "design"),
                row.get("reason", "") or "differential design is blocked",
                args.design,
            )
        if row.get("model_warning", ""):
            add_warning(warnings, args, "warning", "design", row.get("condition", "design"), row["model_warning"], args.design)


def deseq2_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    required = {"contrast_id", "status"}
    for manifest in args.deseq2_manifest:
        _columns, rows = read_table(manifest, required)
        for row in rows:
            contrast_id = row["contrast_id"]
            status = row.get("status", "")
            if status and status != "ok":
                message = row.get("reason", "") or f"DESeq2 contrast status is {status}"
                add_warning(warnings, args, "error", "deseq2_status", contrast_id, message, manifest)

            n_control = parse_int(row.get("n_control", ""))
            n_test = parse_int(row.get("n_test", ""))
            if (
                row.get("n_control", "").strip()
                and row.get("n_test", "").strip()
                and (n_control < args.min_deseq2_replicates or n_test < args.min_deseq2_replicates)
            ):
                add_warning(
                    warnings,
                    args,
                    "warning",
                    "deseq2_replicates",
                    contrast_id,
                    (
                        f"contrast has {n_control} control and {n_test} test samples; "
                        f"minimum expected is {args.min_deseq2_replicates} per group"
                    ),
                    manifest,
                )

            n_tested = parse_int(row.get("n_features_tested", ""))
            if row.get("n_features_tested", "").strip() and n_tested < args.min_deseq2_tested_features:
                add_warning(
                    warnings,
                    args,
                    "warning",
                    "deseq2_features",
                    contrast_id,
                    f"only {n_tested} features passed DESeq2 filtering; expected at least {args.min_deseq2_tested_features}",
                    manifest,
                )

            n_significant = parse_int(row.get("n_significant", ""))
            if n_tested and row.get("n_significant", "").strip() and n_significant == 0:
                padj = row.get("padj_threshold", "")
                log2fc = row.get("log2fc_threshold", "")
                add_warning(
                    warnings,
                    args,
                    "info",
                    "deseq2_signal",
                    contrast_id,
                    f"no significant features at padj <= {padj} and |log2FC| >= {log2fc}",
                    manifest,
                )

            samples = split_csv(row.get("samples", ""))
            for sample, metrics in count_matrix_qc(row.get("counts", ""), samples).items():
                library_size = metrics["library_size"]
                detected = metrics["detected_features"]
                if library_size < args.min_library_size:
                    add_warning(
                        warnings,
                        args,
                        "warning",
                        "deseq2_sample_qc",
                        f"{contrast_id}:{sample}",
                        f"contrast count matrix library size {library_size} is below {args.min_library_size}",
                        row.get("counts", ""),
                    )
                if detected < args.min_detected_features:
                    add_warning(
                        warnings,
                        args,
                        "warning",
                        "deseq2_sample_qc",
                        f"{contrast_id}:{sample}",
                        f"contrast count matrix detected features {detected} is below {args.min_detected_features}",
                        row.get("counts", ""),
                    )

            coldata = coldata_by_sample(row.get("coldata", ""), samples)
            condition_col = row.get("condition_col", "condition")
            formula = row.get("effective_design_formula", "") or row.get("design_formula", "")
            for term in simple_design_terms(formula):
                if term == condition_col:
                    continue
                values_by_condition: dict[str, set[str]] = {}
                condition_values: dict[str, set[str]] = {}
                term_values: set[str] = set()
                missing = False
                for sample in samples:
                    sample_row = coldata.get(sample, {})
                    if term not in sample_row:
                        missing = True
                        continue
                    condition = sample_row.get(condition_col, "")
                    value = sample_row.get(term, "")
                    if value:
                        term_values.add(value)
                    if condition and value:
                        values_by_condition.setdefault(value, set()).add(condition)
                        condition_values.setdefault(condition, set()).add(value)
                if missing:
                    add_warning(
                        warnings,
                        args,
                        "warning",
                        "deseq2_design",
                        f"{contrast_id}:{term}",
                        f"design covariate {term} is absent from one or more DESeq2 coldata rows",
                        row.get("coldata", ""),
                    )
                elif len(term_values) <= 1:
                    add_warning(
                        warnings,
                        args,
                        "info",
                        "deseq2_design",
                        f"{contrast_id}:{term}",
                        f"design covariate {term} is constant in this contrast",
                        row.get("coldata", ""),
                    )
                elif (
                    len(condition_values) > 1
                    and all(len(values) == 1 for values in condition_values.values())
                    and all(len(conditions) == 1 for conditions in values_by_condition.values())
                ):
                    add_warning(
                        warnings,
                        args,
                        "warning",
                        "deseq2_design",
                        f"{contrast_id}:{term}",
                        f"design covariate {term} is confounded with {condition_col} in this contrast",
                        row.get("coldata", ""),
                    )


def sample_qc_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    _columns, rows = read_table(args.sample_qc_metrics, {"library_id", "library_size", "detected_features"})
    for row in rows:
        library_id = row["library_id"]
        library_size = int(parse_float(row.get("library_size", "")))
        detected = int(parse_float(row.get("detected_features", "")))
        if library_size < args.min_library_size:
            add_warning(warnings, args, "warning", "sample_qc", library_id, f"library size {library_size} is below {args.min_library_size}", args.sample_qc_metrics)
        if detected < args.min_detected_features:
            add_warning(warnings, args, "warning", "sample_qc", library_id, f"detected features {detected} is below {args.min_detected_features}", args.sample_qc_metrics)
    _corr_columns, corr_rows = read_table(args.sample_correlations, {"sample_a", "sample_b", "pearson_log_cpm"})
    for row in corr_rows:
        if row["sample_a"] >= row["sample_b"]:
            continue
        correlation = parse_float(row.get("pearson_log_cpm", ""))
        if correlation < args.min_sample_correlation:
            add_warning(
                warnings,
                args,
                "warning",
                "sample_correlation",
                f"{row['sample_a']}:{row['sample_b']}",
                f"log-CPM correlation {correlation:.3g} is below {args.min_sample_correlation}",
                args.sample_correlations,
            )


def strandedness_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    _columns, rows = read_table(args.strandedness_report, {"library_id", "warning"})
    for row in rows:
        if row.get("warning", ""):
            add_warning(warnings, args, "warning", "strandedness", row["library_id"], row["warning"], args.strandedness_report)


def biotype_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    _columns, rows = read_table(args.biotype_count_summary, {"level", "biotype", "detected_features"})
    by_level: dict[str, dict[str, int]] = {}
    for row in rows:
        level = row["level"]
        by_level.setdefault(level, {})
        by_level[level][row["biotype"]] = by_level[level].get(row["biotype"], 0) + int(parse_float(row.get("detected_features", "")))
    for level, counts in by_level.items():
        total = sum(counts.values())
        unclassified = counts.get("unclassified", 0)
        if total and unclassified / total > args.max_unclassified_biotype_fraction:
            add_warning(
                warnings,
                args,
                "warning",
                "biotype",
                level,
                f"{unclassified}/{total} detected features are unclassified",
                args.biotype_count_summary,
            )
    _diff_columns, diff_rows = read_table(args.biotype_differential_summary, {"level", "contrast_id", "biotype", "tested"})
    if args.assay == "rnaseq" and not rows and not diff_rows:
        add_warning(warnings, args, "info", "biotype", "rnaseq", "biotype summary is not available", args.biotype_count_summary)


def transcript_discovery_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    if not args.warn_high_true_novel_transcript_fraction:
        return
    _columns, rows = read_table(
        args.transcript_discovery_summary,
        {"level", "transcript_discovery_class", "true_novel_candidate", "detected_features"},
    )
    if args.assay != "rnaseq" or not rows:
        return
    total = 0
    true_novel = 0
    by_class: dict[str, int] = {}
    for row in rows:
        detected = int(parse_float(row.get("detected_features", "")))
        total += detected
        if row.get("true_novel_candidate", "").lower() == "yes":
            true_novel += detected
            discovery_class = row.get("transcript_discovery_class", "unclassified")
            by_class[discovery_class] = by_class.get(discovery_class, 0) + detected
    fraction = true_novel / total if total else 0.0
    if total and fraction > args.max_true_novel_transcript_fraction:
        class_detail = ", ".join(f"{name}={count}" for name, count in sorted(by_class.items())) or "none"
        add_warning(
            warnings,
            args,
            "warning",
            "transcript_discovery",
            "true_novel_fraction",
            (
                f"{true_novel}/{total} detected transcripts ({fraction:.1%}) are true-novel candidates; "
                f"expected upper reference is {args.max_true_novel_transcript_fraction:.1%}. "
                f"Inspect annotation version, StringTie/gffcompare settings, and read support. Classes: {class_detail}"
            ),
            args.transcript_discovery_summary,
        )


def residual_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    _columns, rows = read_table(args.residual_manifest, {"library_id", "input_reads", "genome_aligned_reads"})
    for row in rows:
        input_reads = parse_float(row.get("input_reads", ""))
        aligned = parse_float(row.get("genome_aligned_reads", ""))
        fraction = aligned / input_reads if input_reads else 0.0
        if fraction > args.max_residual_genome_fraction:
            add_warning(
                warnings,
                args,
                "warning",
                "smallrna_residual",
                row["library_id"],
                f"{fraction:.1%} of miRBase-unmapped reads aligned to the residual genome",
                args.residual_manifest,
            )
    _bio_columns, biotype_rows = read_table(args.residual_biotype_counts, {"biotype"})
    if args.assay == "smallrna" and args.residual_manifest and not biotype_rows:
        add_warning(warnings, args, "info", "smallrna_residual", "biotypes", "residual biotype counts are not available", args.residual_biotype_counts)


def length_warnings(args: argparse.Namespace, warnings: list[dict[str, str]]) -> None:
    _columns, rows = read_table(args.length_stage_summary, {"stage", "library_id", "modal_length"})
    for row in rows:
        if row.get("stage") != "trimmed":
            continue
        modal = int(parse_float(row.get("modal_length", "")))
        if modal and (modal < 18 or modal > 30):
            add_warning(
                warnings,
                args,
                "warning",
                "smallrna_length",
                row["library_id"],
                f"trimmed modal read length {modal} is outside the usual smallRNA range 18-30 nt",
                args.length_stage_summary,
            )
    _arm_columns, arm_rows = read_table(args.arm_summary, {"arm", "fraction"})
    if arm_rows:
        arm_fractions = {row["arm"]: parse_float(row.get("fraction", "")) for row in arm_rows}
        if arm_fractions.get("unannotated_arm", 0.0) > 0.5:
            add_warning(warnings, args, "info", "smallrna_arm", "unannotated_arm", "more than half of miRNA counts lack a 5p/3p arm suffix", args.arm_summary)


def render_html(path: Path, rows: list[dict[str, str]], args: argparse.Namespace) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    severity_counts = {severity: sum(1 for row in rows if row["severity"] == severity) for severity in ["error", "warning", "info"]}
    body = []
    for row in rows:
        body.append(
            "<tr>"
            + "".join(f"<td>{html.escape(row.get(column, ''))}</td>" for column in WARNING_COLUMNS)
            + "</tr>"
        )
    if not body:
        body.append(f'<tr><td colspan="{len(WARNING_COLUMNS)}">No biological warnings were raised.</td></tr>')
    header = "".join(f"<th>{html.escape(column)}</th>" for column in WARNING_COLUMNS)
    metrics = " ".join(f"{key}: {value}" for key, value in severity_counts.items())
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(args.project)} {html.escape(args.assay)} biological warnings</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; vertical-align: top; }}
    th {{ background: #f2f2f2; }}
  </style>
</head>
<body>
  <h1>{html.escape(args.project)} {html.escape(args.assay)} biological warnings</h1>
  <p>{html.escape(metrics)}</p>
  <table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>
</body>
</html>
"""
    path.write_text(content, encoding="utf-8")


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    errors = sum(1 for row in rows if row["severity"] == "error")
    warnings = sum(1 for row in rows if row["severity"] == "warning")
    infos = sum(1 for row in rows if row["severity"] == "info")
    status = "error" if errors else "warning" if warnings else "ok"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\terrors\twarnings\tinfos\ttotal\n")
        handle.write(f"{status}\t{errors}\t{warnings}\t{infos}\t{len(rows)}\n")


def main() -> int:
    args = parse_args()
    warnings: list[dict[str, str]] = []
    design_warnings(args, warnings)
    deseq2_warnings(args, warnings)
    sample_qc_warnings(args, warnings)
    strandedness_warnings(args, warnings)
    biotype_warnings(args, warnings)
    transcript_discovery_warnings(args, warnings)
    residual_warnings(args, warnings)
    length_warnings(args, warnings)
    write_table(Path(args.warnings), WARNING_COLUMNS, warnings)
    render_html(Path(args.summary_html), warnings, args)
    write_table(
        Path(args.manifest),
        MANIFEST_COLUMNS,
        [
            {"resource": "warnings", "status": "ok", "path": args.warnings, "rows": str(len(warnings)), "detail": "aggregated biological warning table"},
            {"resource": "summary_html", "status": "ok", "path": args.summary_html, "rows": str(len(warnings)), "detail": "aggregated biological warning HTML"},
        ],
    )
    write_done(Path(args.done), warnings)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
