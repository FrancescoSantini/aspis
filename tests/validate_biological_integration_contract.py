#!/usr/bin/env python3
"""Exercise biological integration contracts that do not require real data."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("results/biological_integration_contract")
INPUT = BASE / "input"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)


def setup_common_inputs() -> dict[str, Path]:
    if BASE.exists():
        shutil.rmtree(BASE)
    INPUT.mkdir(parents=True)
    paths = {
        "small_plan": INPUT / "smallrna_plan.tsv",
        "small_manifest": INPUT / "smallrna_deseq2_manifest.tsv",
        "small_results": INPUT / "smallrna_results.tsv",
        "small_filtered": INPUT / "smallrna_filtered.tsv",
        "small_norm": INPUT / "smallrna_normalized.tsv",
        "small_samples": INPUT / "smallrna_samples.tsv",
        "rnaseq_manifest": INPUT / "rnaseq_gene_deseq2_manifest.tsv",
        "rnaseq_results": INPUT / "rnaseq_gene_results.tsv",
        "rnaseq_filtered": INPUT / "rnaseq_gene_filtered.tsv",
        "rnaseq_norm": INPUT / "rnaseq_gene_normalized.tsv",
        "rnaseq_samples": INPUT / "rnaseq_samples.tsv",
        "targets_validated": INPUT / "targets_validated.tsv",
        "targets_predicted": INPUT / "targets_predicted.tsv",
    }
    contrast = "treated_vs_control__time_h_24"
    write_tsv(paths["small_plan"], ["stage", "status", "reason"], [{"stage": "mirna_target_enrichment", "status": "ready", "reason": ""}])
    write_tsv(
        paths["small_results"],
        ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"],
        [
            {"Geneid": "hsa-miR-1-3p", "baseMean": "100", "log2FoldChange": "2.0", "pvalue": "0.001", "padj": "0.01"},
            {"Geneid": "hsa-miR-2-3p", "baseMean": "90", "log2FoldChange": "-2.0", "pvalue": "0.002", "padj": "0.02"},
        ],
    )
    write_tsv(paths["small_filtered"], ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"], read_tsv(paths["small_results"]))
    write_tsv(
        paths["small_norm"],
        ["Geneid", "sm1", "sm2"],
        [
            {"Geneid": "hsa-miR-1-3p", "sm1": "10", "sm2": "100"},
            {"Geneid": "hsa-miR-2-3p", "sm1": "90", "sm2": "20"},
        ],
    )
    write_tsv(
        paths["small_manifest"],
        ["contrast_id", "status", "reason", "results", "filtered", "normalized_counts"],
        [{"contrast_id": contrast, "status": "ok", "reason": "", "results": str(paths["small_results"]), "filtered": str(paths["small_filtered"]), "normalized_counts": str(paths["small_norm"])}],
    )
    write_tsv(
        paths["small_samples"],
        ["library_id", "condition", "replicate", "time_h", "biospecimen_id"],
        [
            {"library_id": "sm1", "condition": "control", "replicate": "1", "time_h": "24", "biospecimen_id": "bio1"},
            {"library_id": "sm2", "condition": "treated", "replicate": "1", "time_h": "24", "biospecimen_id": "bio2"},
        ],
    )
    write_tsv(
        paths["rnaseq_results"],
        ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"],
        [
            {"Geneid": "GENE1", "baseMean": "100", "log2FoldChange": "-2.0", "pvalue": "0.001", "padj": "0.01"},
            {"Geneid": "GENE2", "baseMean": "90", "log2FoldChange": "2.0", "pvalue": "0.002", "padj": "0.02"},
            {"Geneid": "GENE3", "baseMean": "80", "log2FoldChange": "-1.2", "pvalue": "0.01", "padj": "0.05"},
        ],
    )
    write_tsv(paths["rnaseq_filtered"], ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"], read_tsv(paths["rnaseq_results"]))
    write_tsv(
        paths["rnaseq_norm"],
        ["Geneid", "rna1", "rna2"],
        [
            {"Geneid": "GENE1", "rna1": "100", "rna2": "10"},
            {"Geneid": "GENE2", "rna1": "20", "rna2": "80"},
            {"Geneid": "GENE3", "rna1": "80", "rna2": "20"},
        ],
    )
    write_tsv(
        paths["rnaseq_manifest"],
        ["contrast_id", "status", "reason", "results", "filtered", "normalized_counts"],
        [{"contrast_id": contrast, "status": "ok", "reason": "", "results": str(paths["rnaseq_results"]), "filtered": str(paths["rnaseq_filtered"]), "normalized_counts": str(paths["rnaseq_norm"])}],
    )
    write_tsv(
        paths["rnaseq_samples"],
        ["library_id", "condition", "replicate", "time_h", "biospecimen_id"],
        [
            {"library_id": "rna1", "condition": "control", "replicate": "1", "time_h": "24", "biospecimen_id": "bio1"},
            {"library_id": "rna2", "condition": "treated", "replicate": "1", "time_h": "24", "biospecimen_id": "bio2"},
        ],
    )
    write_tsv(
        paths["targets_validated"],
        ["mirna_id", "target_id", "target_symbol", "source", "source_type", "evidence"],
        [
            {"mirna_id": "hsa-miR-1-3p", "target_id": "GENE1", "target_symbol": "Gene one", "source": "miRTarBase", "source_type": "validated", "evidence": "strong"},
            {"mirna_id": "hsa-miR-2-3p", "target_id": "GENE2", "target_symbol": "Gene two", "source": "miRTarBase", "source_type": "validated", "evidence": "strong"},
        ],
    )
    write_tsv(
        paths["targets_predicted"],
        ["mirna_id", "target_id", "target_symbol", "source", "source_type", "evidence"],
        [{"mirna_id": "hsa-miR-1-3p", "target_id": "GENE3", "target_symbol": "Gene three", "source": "TargetScan", "source_type": "predicted", "evidence": "context"}],
    )
    return paths


def exercise_target_and_integration(paths: dict[str, Path]) -> None:
    target_manifest = BASE / "target_enrichment" / "target_manifest.tsv"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_smallrna_target_enrichment.py",
            "--smallrna-plan",
            str(paths["small_plan"]),
            "--deseq2-manifest",
            str(paths["small_manifest"]),
            "--target-tables",
            f"{paths['targets_validated']},{paths['targets_predicted']}",
            "--outdir",
            str(BASE / "target_enrichment"),
            "--manifest",
            str(target_manifest),
            "--done",
            str(BASE / "target_enrichment" / "target_enrichment.done"),
            "--min-overlap",
            "1",
            "--top-n",
            "10",
        ]
    )
    target_rows = read_tsv(target_manifest, {"status", "mirna_targets", "target_source_summary"})
    if target_rows[0]["status"] != "ok":
        raise ValueError(f"target enrichment was not ok: {target_rows[0]}")
    source_summary = read_tsv(Path(target_rows[0]["target_source_summary"]), {"target_source", "target_source_type"})
    if not {"validated", "predicted"} <= {row["target_source_type"] for row in source_summary}:
        raise ValueError(f"target sources were not propagated: {source_summary}")

    integration_manifest = BASE / "mirna_mrna" / "mirna_mrna_manifest.tsv"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_mirna_mrna_integration.py",
            "--smallrna-samples",
            str(paths["small_samples"]),
            "--rnaseq-samples",
            str(paths["rnaseq_samples"]),
            "--smallrna-deseq2-manifest",
            str(paths["small_manifest"]),
            "--rnaseq-gene-manifest",
            str(paths["rnaseq_manifest"]),
            "--target-manifest",
            str(target_manifest),
            "--outdir",
            str(BASE / "mirna_mrna"),
            "--manifest",
            str(integration_manifest),
            "--done",
            str(BASE / "mirna_mrna" / "mirna_mrna.done"),
            "--match-columns",
            "biospecimen_id",
            "--min-pairs",
            "2",
        ]
    )
    integration_rows = read_tsv(integration_manifest, {"status", "n_pairs", "n_inverse_pairs", "n_anticorrelated_pairs"})
    row = integration_rows[0]
    if row["status"] != "ok" or int(row["n_pairs"]) < 2 or int(row["n_inverse_pairs"]) < 2:
        raise ValueError(f"unexpected miRNA-mRNA integration row: {row}")


def exercise_biotype_and_dtu(paths: dict[str, Path]) -> None:
    gtf = INPUT / "annotation.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\ttoy\tgene\t1\t100\t.\t+\t.\tgene_id "GENE1"; gene_biotype "protein_coding";',
                'chr1\ttoy\ttranscript\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; transcript_biotype "protein_coding";',
                'chr1\ttoy\tgene\t200\t300\t.\t-\t.\tgene_id "GENE2"; gene_biotype "lncRNA";',
                'chr1\ttoy\ttranscript\t200\t300\t.\t-\t.\tgene_id "GENE2"; transcript_id "TX2"; transcript_biotype "lncRNA";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    gene_counts = INPUT / "gene_counts.tsv"
    gene_metadata = INPUT / "gene_metadata.tsv"
    transcript_counts = INPUT / "transcript_counts.tsv"
    transcript_metadata = INPUT / "transcript_metadata.tsv"
    aligned_samples = INPUT / "aligned_samples.tsv"
    write_tsv(gene_counts, ["Geneid", "s1", "s2"], [{"Geneid": "GENE1", "s1": "10", "s2": "20"}, {"Geneid": "GENE2", "s1": "3", "s2": "0"}])
    write_tsv(gene_metadata, ["Geneid", "feature_type"], [{"Geneid": "GENE1", "feature_type": "protein_coding"}, {"Geneid": "GENE2", "feature_type": "lncRNA"}])
    write_tsv(transcript_counts, ["transcript_id", "s1", "s2"], [{"transcript_id": "TX1", "s1": "10", "s2": "20"}, {"transcript_id": "TX2", "s1": "3", "s2": "0"}])
    write_tsv(transcript_metadata, ["transcript_id", "gene_id"], [{"transcript_id": "TX1", "gene_id": "GENE1"}, {"transcript_id": "TX2", "gene_id": "GENE2"}])
    write_tsv(aligned_samples, ["library_id", "bam"], [{"library_id": "s1", "bam": "s1.bam"}, {"library_id": "s2", "bam": "s2.bam"}])
    biotype_dir = BASE / "biotypes"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_biotype_summary.py",
            "--annotation-gtf",
            str(gtf),
            "--gene-counts",
            str(gene_counts),
            "--gene-metadata",
            str(gene_metadata),
            "--transcript-counts",
            str(transcript_counts),
            "--transcript-metadata",
            str(transcript_metadata),
            "--gene-deseq2-manifest",
            str(paths["rnaseq_manifest"]),
            "--outdir",
            str(biotype_dir),
            "--manifest",
            str(biotype_dir / "biotype_manifest.tsv"),
            "--count-summary",
            str(biotype_dir / "count_biotype_summary.tsv"),
            "--differential-summary",
            str(biotype_dir / "differential_biotype_summary.tsv"),
            "--transcript-discovery-summary",
            str(biotype_dir / "transcript_discovery_summary.tsv"),
            "--transcript-discovery-differential-summary",
            str(biotype_dir / "transcript_discovery_differential_summary.tsv"),
            "--html",
            str(biotype_dir / "biotype_summary.html"),
            "--done",
            str(biotype_dir / "biotype_summary.done"),
        ]
    )
    count_rows = read_tsv(biotype_dir / "count_biotype_summary.tsv", {"level", "biotype", "detected_features"})
    if not {"protein_coding", "lncRNA"} <= {row["biotype"] for row in count_rows}:
        raise ValueError(f"unexpected biotype count rows: {count_rows}")

    dtu_dir = BASE / "dtu"
    run_command(
        [
            sys.executable,
            "workflow/scripts/plan_rnaseq_dtu.py",
            "--samples",
            str(paths["rnaseq_samples"]),
            "--transcript-counts",
            str(transcript_counts),
            "--transcript-metadata",
            str(transcript_metadata),
            "--annotation-gtf",
            str(gtf),
            "--output",
            str(dtu_dir / "dtu_plan.tsv"),
            "--done",
            str(dtu_dir / "dtu.done"),
            "--project",
            "ASPIS_CONTRACT",
        ]
    )
    dtu_rows = read_tsv(dtu_dir / "dtu_plan.tsv", {"status", "candidate_methods"})
    if dtu_rows[0]["status"] != "planned" or "DRIMSeq" not in dtu_rows[0]["candidate_methods"]:
        raise ValueError(f"unexpected DTU plan row: {dtu_rows[0]}")
    run_command(
        [
            sys.executable,
            "workflow/scripts/run_rnaseq_dtu_methods.py",
            "--plan",
            str(dtu_dir / "dtu_plan.tsv"),
            "--samples",
            str(paths["rnaseq_samples"]),
            "--aligned-samples",
            str(aligned_samples),
            "--transcript-counts",
            str(transcript_counts),
            "--transcript-metadata",
            str(transcript_metadata),
            "--annotation-gtf",
            str(gtf),
            "--outdir",
            str(dtu_dir / "methods"),
            "--manifest",
            str(dtu_dir / "dtu_method_manifest.tsv"),
            "--done",
            str(dtu_dir / "dtu_methods.done"),
            "--project",
            "ASPIS_CONTRACT",
        ]
    )
    method_rows = read_tsv(dtu_dir / "dtu_method_manifest.tsv", {"method", "status", "reason"})
    methods = {row["method"] for row in method_rows}
    if not {"DRIMSeq", "DEXSeq", "SUPPA2", "rMATS"} <= methods:
        raise ValueError(f"unexpected DTU method rows: {method_rows}")
    if {row["status"] for row in method_rows} != {"planned"}:
        raise ValueError(f"unconfigured DTU methods should be planned: {method_rows}")
    helper = INPUT / "write_fake_rmats.py"
    helper.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "outdir = Path(sys.argv[1])",
                "outdir.mkdir(parents=True, exist_ok=True)",
                "path = outdir / 'SE.MATS.JC.txt'",
                "path.write_text(",
                "    'ID\\tGeneID\\tgeneSymbol\\tPValue\\tFDR\\tIncLevelDifference\\n'",
                "    'EVENT1\\tGENE1\\tGene One\\t0.001\\t0.01\\t0.25\\n',",
                "    encoding='utf-8',",
                ")",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    executed_dir = BASE / "dtu_executed"
    run_command(
        [
            sys.executable,
            "workflow/scripts/run_rnaseq_dtu_methods.py",
            "--plan",
            str(dtu_dir / "dtu_plan.tsv"),
            "--samples",
            str(paths["rnaseq_samples"]),
            "--aligned-samples",
            str(aligned_samples),
            "--transcript-counts",
            str(transcript_counts),
            "--transcript-metadata",
            str(transcript_metadata),
            "--annotation-gtf",
            str(gtf),
            "--outdir",
            str(executed_dir / "methods"),
            "--manifest",
            str(executed_dir / "dtu_method_manifest.tsv"),
            "--done",
            str(executed_dir / "dtu_methods.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "rMATS",
            "--rmats-command",
            f"{sys.executable} {helper} {{outdir}}",
        ]
    )
    executed_rows = read_tsv(executed_dir / "dtu_method_manifest.tsv", {"method", "status", "standardized_results", "standardized_result_count", "standardized_status"})
    if executed_rows[0]["status"] != "completed" or executed_rows[0]["standardized_status"] != "ok":
        raise ValueError(f"executed DTU method was not standardized: {executed_rows}")
    standardized = read_tsv(Path(executed_rows[0]["standardized_results"]), {"method", "feature_id", "gene_id", "gene_name", "event_type", "pvalue", "padj", "delta_psi", "direction"})
    expected = standardized[0]
    if expected["feature_id"] != "EVENT1" or expected["gene_id"] != "GENE1" or expected["event_type"] != "SE":
        raise ValueError(f"standardized rMATS row lost identifiers: {standardized}")
    if expected["padj"] != "0.01" or expected["direction"] != "increased_usage":
        raise ValueError(f"standardized rMATS row lost statistics: {standardized}")
    exercise_rnaseq_report_dtu_assets(dtu_dir / "dtu_plan.tsv", executed_dir / "dtu_method_manifest.tsv")


def exercise_rnaseq_report_dtu_assets(dtu_plan: Path, dtu_manifest: Path) -> None:
    report_dir = BASE / "rnaseq_report_index"
    contrast = "treated_vs_control__time_h_24"
    base_row = {
        "project": "ASPIS_CONTRACT",
        "level": "transcript",
        "contrast_id": contrast,
        "status": "ok",
        "reason": "",
        "results": str(INPUT / "rnaseq_gene_results.tsv"),
        "filtered": str(INPUT / "rnaseq_gene_filtered.tsv"),
        "summary_html": str(report_dir / "summary.html"),
        "volcano_pdf": str(report_dir / "volcano.pdf"),
        "ma_pdf": str(report_dir / "ma.pdf"),
        "pca_pdf": str(report_dir / "pca.pdf"),
        "pca_metrics_tsv": str(report_dir / "pca_metrics.tsv"),
        "heatmap_pdf": str(report_dir / "heatmap.pdf"),
        "heatmap_panel_tsv": str(report_dir / "heatmap_panel.tsv"),
        "vst_tsv": str(report_dir / "vst.tsv"),
        "enrichment_manifest": str(report_dir / "enrichment_manifest.tsv"),
        "ranked_features": str(report_dir / "ranked_features.tsv"),
        "significant_features": str(report_dir / "significant_features.tsv"),
        "up_features": str(report_dir / "up_features.tsv"),
        "down_features": str(report_dir / "down_features.tsv"),
        "feature_set_universe": str(report_dir / "feature_set_universe.tsv"),
        "feature_set_results": str(report_dir / "feature_set_results.tsv"),
        "feature_set_plot": str(report_dir / "feature_set_plot.pdf"),
        "ranked_feature_set_results": str(report_dir / "ranked_feature_set_results.tsv"),
        "ranked_feature_set_plot": str(report_dir / "ranked_feature_set_plot.pdf"),
        "n_features": "2",
        "n_significant": "1",
        "n_up": "1",
        "n_down": "0",
        "n_ranked": "2",
        "n_feature_sets": "1",
        "n_feature_set_resources": "1",
        "n_feature_set_terms": "1",
        "n_ranked_feature_set_terms": "1",
    }
    plan_columns = [
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "results",
        "filtered",
        "summary_html",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "pca_metrics_tsv",
        "heatmap_pdf",
        "heatmap_panel_tsv",
        "vst_tsv",
        "enrichment_manifest",
    ]
    plots_columns = [
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "pca_metrics_tsv",
        "heatmap_pdf",
        "heatmap_panel_tsv",
        "vst_tsv",
        "n_features",
        "n_significant",
    ]
    enrichment_columns = [
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "enrichment_manifest",
        "ranked_features",
        "significant_features",
        "up_features",
        "down_features",
        "feature_set_universe",
        "feature_set_results",
        "feature_set_plot",
        "ranked_feature_set_results",
        "ranked_feature_set_plot",
        "n_ranked",
        "n_significant",
        "n_up",
        "n_down",
        "n_feature_sets",
        "n_feature_set_resources",
        "n_feature_set_terms",
        "n_ranked_feature_set_terms",
    ]
    summary_columns = [
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "summary_html",
        "results",
        "filtered",
        "ma_pdf",
        "pca_metrics_tsv",
        "n_features",
        "n_significant",
        "n_up",
        "n_down",
    ]
    write_tsv(report_dir / "report_plan.tsv", plan_columns, [base_row])
    write_tsv(report_dir / "plots_manifest.tsv", plots_columns, [base_row])
    write_tsv(report_dir / "enrichment_manifest.tsv", enrichment_columns, [base_row])
    write_tsv(report_dir / "summary_manifest.tsv", summary_columns, [base_row])
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_differential_report_index.py",
            "--plan",
            str(report_dir / "report_plan.tsv"),
            "--plots-manifest",
            str(report_dir / "plots_manifest.tsv"),
            "--enrichment-manifest",
            str(report_dir / "enrichment_manifest.tsv"),
            "--summary-manifest",
            str(report_dir / "summary_manifest.tsv"),
            "--dtu-plan",
            str(dtu_plan),
            "--dtu-method-manifest",
            str(dtu_manifest),
            "--asset-manifest",
            str(report_dir / "asset_manifest.tsv"),
            "--output",
            str(report_dir / "index.html"),
            "--done",
            str(report_dir / "report_index.done"),
        ]
    )
    assets = read_tsv(report_dir / "asset_manifest.tsv", {"asset_group", "asset_label", "path"})
    dtu_assets = [row for row in assets if row["asset_group"] == "dtu"]
    if not any(row["asset_label"] == "dtu_method_manifest" for row in dtu_assets):
        raise ValueError(f"DTU manifest was not exposed as a report asset: {dtu_assets}")
    if not any(row["asset_label"] == "rMATS_standardized_results" for row in dtu_assets):
        raise ValueError(f"standardized DTU results were not exposed as report assets: {dtu_assets}")
    html_text = (report_dir / "index.html").read_text(encoding="utf-8")
    if "DTU / splicing methods" not in html_text or "standardized rows: 1" not in html_text:
        raise ValueError("DTU summary was not rendered in the RNA-seq report index")


def main() -> int:
    paths = setup_common_inputs()
    exercise_target_and_integration(paths)
    exercise_biotype_and_dtu(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
