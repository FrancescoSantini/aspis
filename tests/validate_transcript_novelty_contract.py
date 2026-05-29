#!/usr/bin/env python3
"""Validate transcript discovery classes from StringTie/gffcompare outputs."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("results/transcript_novelty_contract")
INPUT = BASE / "input"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
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


def write_quant_gtf(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                'chr1\tStringTie\ttranscript\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_KNOWN"; gene_name "Gene1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding"; cov "10";',
                'chr1\tStringTie\texon\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_KNOWN";',
                'chr1\tStringTie\ttranscript\t1\t160\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_J"; gene_name "Gene1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding"; cov "8";',
                'chr1\tStringTie\texon\t1\t70\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_J";',
                'chr1\tStringTie\texon\t100\t160\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_J";',
                'chr1\tStringTie\ttranscript\t1\t90\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_C"; gene_name "Gene1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding"; cov "7";',
                'chr1\tStringTie\texon\t1\t90\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_C";',
                'chr1\tStringTie\ttranscript\t300\t380\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "TX_U"; gene_biotype "lncRNA"; transcript_biotype "lncRNA"; cov "5";',
                'chr1\tStringTie\texon\t300\t380\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "TX_U";',
                'chr1\tStringTie\ttranscript\t500\t580\t.\t+\t.\tgene_id "GENE2"; transcript_id "TX_I"; gene_name "Gene2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA"; cov "4";',
                'chr1\tStringTie\texon\t500\t580\t.\t+\t.\tgene_id "GENE2"; transcript_id "TX_I";',
                'chr1\tStringTie\ttranscript\t700\t780\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX_P"; gene_name "Gene3"; gene_biotype "processed_pseudogene"; transcript_biotype "processed_pseudogene"; cov "3";',
                'chr1\tStringTie\texon\t700\t780\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX_P";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def exercise_transcript_matrix() -> tuple[Path, Path]:
    if BASE.exists():
        shutil.rmtree(BASE)
    INPUT.mkdir(parents=True)
    quant_gtf = INPUT / "sample1.gtf"
    write_quant_gtf(quant_gtf)
    manifest = INPUT / "quant_manifest.tsv"
    plan = INPUT / "quantification_plan.tsv"
    tmap = INPUT / "merged.tmap"
    counts = BASE / "counts" / "transcript_counts.tsv"
    metadata = BASE / "counts" / "transcript_metadata.tsv"
    write_tsv(manifest, ["library_id", "quant_gtf", "status"], [{"library_id": "sample1", "quant_gtf": str(quant_gtf), "status": "ok"}])
    write_tsv(plan, ["status", "read_length"], [{"status": "ready", "read_length": "50"}])
    write_tsv(
        tmap,
        ["q_id", "class_code"],
        [
            {"q_id": "TX_KNOWN", "class_code": "="},
            {"q_id": "TX_J", "class_code": "j"},
            {"q_id": "TX_C", "class_code": "c"},
            {"q_id": "TX_U", "class_code": "u"},
            {"q_id": "TX_I", "class_code": "i"},
            {"q_id": "TX_P", "class_code": "p"},
        ],
    )
    run_command(
        [
            sys.executable,
            "workflow/scripts/build_stringtie_transcript_matrix.py",
            "--quant-manifest",
            str(manifest),
            "--plan",
            str(plan),
            "--tmap",
            str(tmap),
            "--counts",
            str(counts),
            "--metadata",
            str(metadata),
            "--done",
            str(BASE / "counts" / "transcript_counts.done"),
            "--known-codes-strict",
            "=",
            "--known-codes-lenient",
            "=,c,k,m,n,y",
        ]
    )
    rows = {
        row["transcript_id"]: row
        for row in read_tsv(
            metadata,
            {
                "transcript_id",
                "gene_biotype",
                "transcript_biotype",
                "class_code",
                "transcript_discovery_class",
                "transcript_novelty",
                "true_novel_candidate",
                "transcript_plot_group",
                "transcript_plot_label",
            },
        )
    }
    expected_discovery = {
        "TX_KNOWN": ("known_transcript", "known", "no", "known_compatible"),
        "TX_C": ("reference_contained_or_containing", "reference_overlap", "no", "known_compatible"),
        "TX_J": ("novel_isoform_known_gene", "novel_isoform", "yes", "novel_isoform"),
        "TX_U": ("intergenic_novel_locus", "novel_locus", "yes", "novel_locus"),
        "TX_I": ("intronic_novel_candidate", "ambiguous_overlap", "no", "ambiguous"),
        "TX_P": ("likely_artifact_or_repeat", "low_confidence", "no", "artifact"),
    }
    for tx_id, expected in expected_discovery.items():
        observed = rows[tx_id]
        actual = (
            observed["transcript_discovery_class"],
            observed["transcript_novelty"],
            observed["true_novel_candidate"],
            observed["transcript_plot_group"],
        )
        if actual != expected:
            raise ValueError(f"{tx_id} discovery classification mismatch: expected {expected}, observed {observed}")
    if rows["TX_KNOWN"]["transcript_biotype"] != "protein_coding" or rows["TX_U"]["transcript_biotype"] != "lncRNA":
        raise ValueError(f"transcript biotypes were not propagated from StringTie GTF: {rows}")
    if rows["TX_KNOWN"]["transcript_novelty"] != "known" or rows["TX_KNOWN"]["gene_type_strict"] != "Known":
        raise ValueError(f"exact match was not treated as known: {rows['TX_KNOWN']}")
    if rows["TX_KNOWN"]["transcript_plot_group"] != "known_compatible":
        raise ValueError(f"exact match was not assigned to known-compatible plots: {rows['TX_KNOWN']}")
    if rows["TX_J"]["transcript_novelty"] != "novel_isoform" or rows["TX_J"]["true_novel_candidate"] != "yes":
        raise ValueError(f"class j was not treated as a true novel isoform: {rows['TX_J']}")
    if rows["TX_J"]["transcript_plot_group"] != "novel_isoform":
        raise ValueError(f"class j was not assigned to novel-isoform plots: {rows['TX_J']}")
    if rows["TX_J"]["gene_type_strict"] != "Novel":
        raise ValueError(f"class j was still treated as strict Known: {rows['TX_J']}")
    if rows["TX_U"]["transcript_novelty"] != "novel_locus" or rows["TX_U"]["true_novel_candidate"] != "yes":
        raise ValueError(f"class u was not treated as a true novel locus: {rows['TX_U']}")
    if rows["TX_U"]["transcript_plot_group"] != "novel_locus":
        raise ValueError(f"class u was not assigned to novel-locus plots: {rows['TX_U']}")
    return counts, metadata


def exercise_grouped_plots() -> None:
    plot_dir = BASE / "plots"
    results = plot_dir / "transcript_results.tsv"
    filtered = plot_dir / "transcript_filtered.tsv"
    normalized = plot_dir / "normalized_counts.tsv"
    coldata = plot_dir / "coldata.tsv"
    heatmap_features = plot_dir / "heatmap_features.tsv"
    plan = plot_dir / "report_plan.tsv"
    rows = [
        {
            "transcript_id": "TX_KNOWN",
            "baseMean": "100",
            "log2FoldChange": "2.0",
            "pvalue": "0.001",
            "padj": "0.01",
            "transcript_plot_group": "known_compatible",
            "transcript_plot_label": "Known/reference-compatible",
            "transcript_biotype": "protein_coding",
        },
        {
            "transcript_id": "TX_J",
            "baseMean": "90",
            "log2FoldChange": "-2.0",
            "pvalue": "0.002",
            "padj": "0.02",
            "transcript_plot_group": "novel_isoform",
            "transcript_plot_label": "Novel isoform",
            "transcript_biotype": "protein_coding",
        },
        {
            "transcript_id": "TX_U",
            "baseMean": "80",
            "log2FoldChange": "1.5",
            "pvalue": "0.01",
            "padj": "0.05",
            "transcript_plot_group": "novel_locus",
            "transcript_plot_label": "Novel locus",
            "transcript_biotype": "lncRNA",
        },
        {
            "transcript_id": "TX_I",
            "baseMean": "70",
            "log2FoldChange": "-1.2",
            "pvalue": "0.03",
            "padj": "0.08",
            "transcript_plot_group": "ambiguous",
            "transcript_plot_label": "Ambiguous overlap",
            "transcript_biotype": "lncRNA",
        },
        {
            "transcript_id": "TX_P",
            "baseMean": "50",
            "log2FoldChange": "1.1",
            "pvalue": "0.04",
            "padj": "0.09",
            "transcript_plot_group": "artifact",
            "transcript_plot_label": "Artifact/repeat",
            "transcript_biotype": "processed_pseudogene",
        },
    ]
    write_tsv(results, list(rows[0]), rows)
    write_tsv(filtered, list(rows[0]), rows)
    write_tsv(
        normalized,
        ["transcript_id", "s1", "s2", "s3"],
        [
            {"transcript_id": "TX_KNOWN", "s1": "10", "s2": "20", "s3": "30"},
            {"transcript_id": "TX_J", "s1": "30", "s2": "20", "s3": "10"},
            {"transcript_id": "TX_U", "s1": "5", "s2": "15", "s3": "30"},
            {"transcript_id": "TX_I", "s1": "20", "s2": "10", "s3": "8"},
            {"transcript_id": "TX_P", "s1": "6", "s2": "10", "s3": "12"},
        ],
    )
    write_tsv(
        coldata,
        ["library_id", "condition"],
        [
            {"library_id": "s1", "condition": "control"},
            {"library_id": "s2", "condition": "treated"},
            {"library_id": "s3", "condition": "treated"},
        ],
    )
    write_tsv(
        heatmap_features,
        ["feature_id", "feature_list", "level", "plot_group"],
        [
            {"feature_id": "TX_KNOWN", "feature_list": "curated_switch_candidates", "level": "transcript", "plot_group": "all"},
            {"feature_id": "TX_J", "feature_list": "curated_switch_candidates", "level": "transcript", "plot_group": "all"},
        ],
    )
    write_tsv(
        plan,
        [
            "project",
            "level",
            "contrast_id",
            "status",
            "reason",
            "results",
            "filtered",
            "normalized_counts",
            "coldata",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "heatmap_panel_tsv",
            "vst_tsv",
        ],
        [
            {
                "project": "ASPIS_TRANSCRIPT_NOVELTY",
                "level": "transcript",
                "contrast_id": "treated_vs_control",
                "status": "ready",
                "reason": "",
                "results": str(results),
                "filtered": str(filtered),
                "normalized_counts": str(normalized),
                "coldata": str(coldata),
                "volcano_pdf": str(plot_dir / "volcano.pdf"),
                "ma_pdf": str(plot_dir / "ma.pdf"),
                "pca_pdf": str(plot_dir / "pca.pdf"),
                "heatmap_pdf": str(plot_dir / "heatmap.pdf"),
                "heatmap_panel_tsv": str(plot_dir / "heatmap_panels.tsv"),
                "vst_tsv": str(plot_dir / "vst.tsv"),
            }
        ],
    )
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript is required for grouped transcript plot validation")
    run_command(
        [
            rscript,
            "workflow/scripts/render_rnaseq_differential_plots.R",
            "--plan",
            str(plan),
            "--manifest",
            str(plot_dir / "plots_manifest.tsv"),
            "--done",
            str(plot_dir / "plots.done"),
            "--top-n",
            "3",
            "--heatmap-modes",
            "significant,variable,feature_list",
            "--heatmap-feature-lists",
            str(heatmap_features),
            "--padj",
            "0.1",
            "--log2fc",
            "1.0",
            "--transcript-plot-groups",
            "all,known_compatible,novel_isoform,novel_locus,ambiguous,artifact",
            "--transcript-biotype-plot-groups",
            "protein_coding,lncRNA,processed_pseudogene",
        ]
    )
    manifest = read_tsv(
        plot_dir / "plots_manifest.tsv",
        {"status", "volcano_pdf", "heatmap_pdf", "heatmap_panel_tsv", "plot_group_tsv"},
    )
    if manifest[0]["status"] != "ok":
        raise ValueError(f"grouped plot rendering failed: {manifest[0]}")
    for key in ["volcano_pdf", "heatmap_pdf"]:
        output = Path(manifest[0][key])
        if not output.exists() or output.stat().st_size == 0:
            raise ValueError(f"grouped plot output was not written: {output}")
    panel_rows = read_tsv(
        Path(manifest[0]["heatmap_panel_tsv"]),
        {"plot_group", "heatmap_mode", "status", "n_plotted_features", "features"},
    )
    if not any(row["heatmap_mode"] == "feature_list" and row["status"] == "ok" for row in panel_rows):
        raise ValueError(f"configured heatmap feature-list panel was not rendered: {panel_rows}")
    if not any(row["heatmap_mode"] == "variable" and row["status"] == "ok" for row in panel_rows):
        raise ValueError(f"variable-feature heatmap panel was not rendered: {panel_rows}")
    group_rows = read_tsv(
        Path(manifest[0]["plot_group_tsv"]),
        {"plot_group", "plot_group_type", "plot_label", "n_features"},
    )
    groups = {row["plot_group"] for row in group_rows}
    expected_groups = {
        "transcript_biotype__protein_coding",
        "transcript_biotype__lncrna",
        "transcript_biotype__processed_pseudogene",
        "ambiguous",
        "artifact",
        "transcript_novelty_biotype__known_compatible__protein_coding",
        "transcript_novelty_biotype__novel_locus__lncrna",
        "transcript_novelty_biotype__ambiguous__lncrna",
        "transcript_novelty_biotype__artifact__processed_pseudogene",
    }
    if not expected_groups <= groups:
        raise ValueError(f"missing transcript biotype/novelty plot groups: {group_rows}")


def exercise_transcript_novelty_report_summary() -> None:
    report_dir = BASE / "report_summary"
    results = report_dir / "transcript_results.tsv"
    filtered = report_dir / "transcript_filtered.tsv"
    deseq2_summary = report_dir / "deseq2_summary.tsv"
    pca_metrics = report_dir / "pca_metrics.tsv"
    enrichment_resources = report_dir / "enrichment_resources.tsv"
    report_plan = report_dir / "report_plan.tsv"
    summary_manifest = report_dir / "summary_manifest.tsv"
    novelty_summary = report_dir / "novelty_summary.tsv"
    rows = [
        {
            "transcript_id": "TX_KNOWN",
            "baseMean": "100",
            "log2FoldChange": "2.0",
            "pvalue": "0.001",
            "padj": "0.01",
            "class_code": "=",
            "transcript_novelty": "known",
            "transcript_plot_group": "known_compatible",
            "transcript_plot_label": "Known/reference-compatible",
            "true_novel_candidate": "no",
        },
        {
            "transcript_id": "TX_J",
            "baseMean": "90",
            "log2FoldChange": "-2.0",
            "pvalue": "0.002",
            "padj": "0.02",
            "class_code": "j",
            "transcript_novelty": "novel_isoform",
            "transcript_plot_group": "novel_isoform",
            "transcript_plot_label": "Novel isoform",
            "true_novel_candidate": "yes",
        },
        {
            "transcript_id": "TX_U",
            "baseMean": "80",
            "log2FoldChange": "0.5",
            "pvalue": "0.1",
            "padj": "0.2",
            "class_code": "u",
            "transcript_novelty": "novel_locus",
            "transcript_plot_group": "novel_locus",
            "transcript_plot_label": "Novel locus",
            "true_novel_candidate": "yes",
        },
    ]
    write_tsv(results, list(rows[0]), rows)
    write_tsv(filtered, list(rows[0]), rows[:2])
    write_tsv(
        deseq2_summary,
        ["status", "padj_threshold", "log2fc_threshold"],
        [{"status": "ok", "padj_threshold": "0.1", "log2fc_threshold": "1.0"}],
    )
    write_tsv(
        pca_metrics,
        ["status", "pc1_variance_percent", "pc2_variance_percent"],
        [{"status": "ok", "pc1_variance_percent": "40.0", "pc2_variance_percent": "20.0"}],
    )
    write_tsv(
        enrichment_resources,
        ["resource", "status", "path"],
        [{"resource": "feature_set_results", "status": "ok", "path": str(report_dir / "feature_sets.tsv")}],
    )
    write_tsv(
        report_plan,
        [
            "project",
            "level",
            "contrast_id",
            "status",
            "reason",
            "results",
            "filtered",
            "deseq2_summary",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "pca_metrics_tsv",
            "sample_distance_pdf",
            "heatmap_pdf",
            "heatmap_panel_tsv",
            "plot_group_tsv",
            "novelty_summary_tsv",
            "vst_tsv",
            "enrichment_manifest",
            "summary_html",
        ],
        [
            {
                "project": "ASPIS_TRANSCRIPT_NOVELTY",
                "level": "transcript",
                "contrast_id": "treated_vs_control",
                "status": "ready",
                "reason": "",
                "results": str(results),
                "filtered": str(filtered),
                "deseq2_summary": str(deseq2_summary),
                "volcano_pdf": str(report_dir / "volcano.pdf"),
                "ma_pdf": str(report_dir / "ma.pdf"),
                "pca_pdf": str(report_dir / "pca.pdf"),
                "pca_metrics_tsv": str(pca_metrics),
                "sample_distance_pdf": str(report_dir / "sample_distance.pdf"),
                "heatmap_pdf": str(report_dir / "heatmap.pdf"),
                "heatmap_panel_tsv": str(report_dir / "heatmap_panels.tsv"),
                "plot_group_tsv": str(report_dir / "plot_groups.tsv"),
                "novelty_summary_tsv": str(novelty_summary),
                "vst_tsv": str(report_dir / "vst.tsv"),
                "enrichment_manifest": str(enrichment_resources),
                "summary_html": str(report_dir / "summary.html"),
            }
        ],
    )
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_differential_summary.py",
            "--plan",
            str(report_plan),
            "--manifest",
            str(summary_manifest),
            "--done",
            str(report_dir / "summary.done"),
            "--top-n",
            "3",
        ]
    )
    summary_rows = read_tsv(
        novelty_summary,
        {
            "transcript_novelty",
            "transcript_plot_group",
            "n_tested",
            "fraction_tested",
            "n_significant",
            "fraction_significant",
            "n_up",
            "n_down",
            "n_true_novel_candidates",
            "n_significant_true_novel_candidates",
        },
    )
    by_novelty = {row["transcript_novelty"]: row for row in summary_rows}
    if by_novelty["known"]["n_tested"] != "1" or by_novelty["known"]["fraction_tested"] != "0.333333":
        raise ValueError(f"known transcript novelty counts were not summarized: {summary_rows}")
    if by_novelty["novel_isoform"]["n_down"] != "1" or by_novelty["novel_isoform"]["n_significant_true_novel_candidates"] != "1":
        raise ValueError(f"novel-isoform transcript novelty counts were not summarized: {summary_rows}")
    if by_novelty["novel_locus"]["n_significant"] != "0" or by_novelty["novel_locus"]["n_true_novel_candidates"] != "1":
        raise ValueError(f"novel-locus transcript novelty counts were not summarized: {summary_rows}")
    html_text = (report_dir / "summary.html").read_text(encoding="utf-8")
    if "Transcript Novelty Summary" not in html_text or "Known/reference-compatible" not in html_text:
        raise ValueError("transcript novelty summary was not rendered in the contrast HTML")

    plots_manifest = report_dir / "plots_manifest.tsv"
    index_enrichment_manifest = report_dir / "index_enrichment_manifest.tsv"
    report_index = report_dir / "report_index.html"
    asset_manifest = report_dir / "asset_manifest.tsv"
    write_tsv(
        plots_manifest,
        [
            "project",
            "level",
            "contrast_id",
            "status",
            "reason",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "pca_metrics_tsv",
            "sample_distance_pdf",
            "heatmap_pdf",
            "heatmap_panel_tsv",
            "plot_group_tsv",
            "vst_tsv",
            "n_features",
            "n_significant",
        ],
        [
            {
                "project": "ASPIS_TRANSCRIPT_NOVELTY",
                "level": "transcript",
                "contrast_id": "treated_vs_control",
                "status": "ok",
                "reason": "",
                "volcano_pdf": str(report_dir / "volcano.pdf"),
                "ma_pdf": str(report_dir / "ma.pdf"),
                "pca_pdf": str(report_dir / "pca.pdf"),
                "pca_metrics_tsv": str(pca_metrics),
                "sample_distance_pdf": str(report_dir / "sample_distance.pdf"),
                "heatmap_pdf": str(report_dir / "heatmap.pdf"),
                "heatmap_panel_tsv": str(report_dir / "heatmap_panels.tsv"),
                "plot_group_tsv": str(report_dir / "plot_groups.tsv"),
                "vst_tsv": str(report_dir / "vst.tsv"),
                "n_features": "3",
                "n_significant": "2",
            }
        ],
    )
    write_tsv(
        index_enrichment_manifest,
        [
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
        ],
        [
            {
                "project": "ASPIS_TRANSCRIPT_NOVELTY",
                "level": "transcript",
                "contrast_id": "treated_vs_control",
                "status": "ok",
                "reason": "",
                "enrichment_manifest": str(enrichment_resources),
                "ranked_features": str(report_dir / "ranked.tsv"),
                "significant_features": str(report_dir / "significant.tsv"),
                "up_features": str(report_dir / "up.tsv"),
                "down_features": str(report_dir / "down.tsv"),
                "feature_set_universe": "",
                "feature_set_results": "",
                "feature_set_plot": "",
                "ranked_feature_set_results": "",
                "ranked_feature_set_plot": "",
                "n_ranked": "3",
                "n_significant": "2",
                "n_up": "1",
                "n_down": "1",
                "n_feature_sets": "0",
                "n_feature_set_resources": "0",
                "n_feature_set_terms": "0",
                "n_ranked_feature_set_terms": "0",
            }
        ],
    )
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_differential_report_index.py",
            "--plan",
            str(report_plan),
            "--plots-manifest",
            str(plots_manifest),
            "--enrichment-manifest",
            str(index_enrichment_manifest),
            "--summary-manifest",
            str(summary_manifest),
            "--asset-manifest",
            str(asset_manifest),
            "--output",
            str(report_index),
            "--done",
            str(report_dir / "report_index.done"),
        ]
    )
    index_html = report_index.read_text(encoding="utf-8")
    if "novelty" not in index_html:
        raise ValueError("report index did not link the transcript novelty summary")
    asset_rows = read_tsv(asset_manifest, {"asset_label", "path"})
    if not any(row["asset_label"] == "novelty_summary_tsv" and row["path"] == str(novelty_summary) for row in asset_rows):
        raise ValueError(f"report asset manifest did not include the transcript novelty summary: {asset_rows}")


def exercise_gene_biotype_plots() -> None:
    plot_dir = BASE / "gene_plots"
    results = plot_dir / "gene_results.tsv"
    filtered = plot_dir / "gene_filtered.tsv"
    normalized = plot_dir / "normalized_counts.tsv"
    coldata = plot_dir / "coldata.tsv"
    plan = plot_dir / "report_plan.tsv"
    rows = [
        {"Geneid": "GENE1", "baseMean": "100", "log2FoldChange": "2.0", "pvalue": "0.001", "padj": "0.01", "gene_biotype": "protein_coding"},
        {"Geneid": "GENE2", "baseMean": "80", "log2FoldChange": "-1.5", "pvalue": "0.01", "padj": "0.04", "gene_biotype": "lncRNA"},
        {"Geneid": "GENE3", "baseMean": "60", "log2FoldChange": "0.5", "pvalue": "0.2", "padj": "0.5", "gene_biotype": "processed_pseudogene"},
    ]
    write_tsv(results, list(rows[0]), rows)
    write_tsv(filtered, list(rows[0]), rows[:2])
    write_tsv(
        normalized,
        ["Geneid", "s1", "s2", "s3"],
        [
            {"Geneid": "GENE1", "s1": "10", "s2": "20", "s3": "30"},
            {"Geneid": "GENE2", "s1": "30", "s2": "20", "s3": "10"},
            {"Geneid": "GENE3", "s1": "5", "s2": "15", "s3": "30"},
        ],
    )
    write_tsv(
        coldata,
        ["library_id", "condition"],
        [
            {"library_id": "s1", "condition": "control"},
            {"library_id": "s2", "condition": "treated"},
            {"library_id": "s3", "condition": "treated"},
        ],
    )
    write_tsv(
        plan,
        [
            "project",
            "level",
            "contrast_id",
            "status",
            "reason",
            "results",
            "filtered",
            "normalized_counts",
            "coldata",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "heatmap_panel_tsv",
            "vst_tsv",
        ],
        [
            {
                "project": "ASPIS_GENE_BIOTYPE",
                "level": "gene",
                "contrast_id": "treated_vs_control",
                "status": "ready",
                "reason": "",
                "results": str(results),
                "filtered": str(filtered),
                "normalized_counts": str(normalized),
                "coldata": str(coldata),
                "volcano_pdf": str(plot_dir / "volcano.pdf"),
                "ma_pdf": str(plot_dir / "ma.pdf"),
                "pca_pdf": str(plot_dir / "pca.pdf"),
                "heatmap_pdf": str(plot_dir / "heatmap.pdf"),
                "heatmap_panel_tsv": str(plot_dir / "heatmap_panels.tsv"),
                "vst_tsv": str(plot_dir / "vst.tsv"),
            }
        ],
    )
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript is required for gene biotype plot validation")
    run_command(
        [
            rscript,
            "workflow/scripts/render_rnaseq_differential_plots.R",
            "--plan",
            str(plan),
            "--manifest",
            str(plot_dir / "plots_manifest.tsv"),
            "--done",
            str(plot_dir / "plots.done"),
            "--top-n",
            "3",
            "--padj",
            "0.1",
            "--log2fc",
            "1.0",
            "--gene-biotype-plot-groups",
            "protein_coding,lncRNA,pseudogene",
        ]
    )
    manifest = read_tsv(plot_dir / "plots_manifest.tsv", {"status", "plot_group_tsv", "heatmap_panel_tsv"})[0]
    if manifest["status"] != "ok":
        raise ValueError(f"gene biotype plot rendering failed: {manifest}")
    panel_rows = read_tsv(Path(manifest["heatmap_panel_tsv"]), {"heatmap_mode", "status"})
    if not any(row["heatmap_mode"] == "significant" for row in panel_rows):
        raise ValueError(f"gene heatmap significant panel was not planned: {panel_rows}")
    group_rows = read_tsv(Path(manifest["plot_group_tsv"]), {"plot_group", "plot_group_type"})
    groups = {row["plot_group"] for row in group_rows}
    expected_groups = {"gene_biotype__protein_coding", "gene_biotype__lncrna", "gene_biotype__processed_pseudogene"}
    if not expected_groups <= groups:
        raise ValueError(f"missing gene biotype plot groups: {group_rows}")


def exercise_discovery_reports(counts: Path, metadata: Path) -> None:
    outdir = BASE / "biotypes"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_biotype_summary.py",
            "--transcript-counts",
            str(counts),
            "--transcript-metadata",
            str(metadata),
            "--outdir",
            str(outdir),
            "--manifest",
            str(outdir / "biotype_manifest.tsv"),
            "--count-summary",
            str(outdir / "count_biotype_summary.tsv"),
            "--differential-summary",
            str(outdir / "differential_biotype_summary.tsv"),
            "--transcript-discovery-summary",
            str(outdir / "transcript_discovery_summary.tsv"),
            "--transcript-discovery-differential-summary",
            str(outdir / "transcript_discovery_differential_summary.tsv"),
            "--html",
            str(outdir / "biotype_summary.html"),
            "--done",
            str(outdir / "biotype_summary.done"),
        ]
    )
    discovery_rows = read_tsv(
        outdir / "transcript_discovery_summary.tsv",
        {
            "transcript_discovery_class",
            "transcript_novelty",
            "true_novel_candidate",
            "transcript_plot_group",
            "detected_features",
            "detected_feature_fraction",
            "true_novel_reference_fraction",
        },
    )
    classes = {row["transcript_discovery_class"] for row in discovery_rows}
    expected_classes = {
        "known_transcript",
        "reference_contained_or_containing",
        "novel_isoform_known_gene",
        "intergenic_novel_locus",
        "intronic_novel_candidate",
        "likely_artifact_or_repeat",
    }
    if not expected_classes <= classes:
        raise ValueError(f"missing transcript discovery classes: {discovery_rows}")
    groups = {row["transcript_plot_group"] for row in discovery_rows}
    if not {"known_compatible", "novel_isoform", "novel_locus", "ambiguous", "artifact"} <= groups:
        raise ValueError(f"missing transcript plot groups: {discovery_rows}")

    warnings_dir = BASE / "warnings"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_biological_warnings.py",
            "--assay",
            "rnaseq",
            "--project",
            "ASPIS_TRANSCRIPT_NOVELTY",
            "--transcript-discovery-summary",
            str(outdir / "transcript_discovery_summary.tsv"),
            "--outdir",
            str(warnings_dir),
            "--warnings",
            str(warnings_dir / "warnings.tsv"),
            "--summary-html",
            str(warnings_dir / "warnings.html"),
            "--manifest",
            str(warnings_dir / "warnings_manifest.tsv"),
            "--done",
            str(warnings_dir / "warnings.done"),
            "--max-true-novel-transcript-fraction",
            "0.2",
        ]
    )
    warning_rows = read_tsv(warnings_dir / "warnings.tsv", {"category", "item"})
    if any(row["category"] == "transcript_discovery" and row["item"] == "true_novel_fraction" for row in warning_rows):
        raise ValueError(f"true-novel fraction warning should be opt-in: {warning_rows}")

    opt_in_dir = BASE / "warnings_opt_in"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_biological_warnings.py",
            "--assay",
            "rnaseq",
            "--project",
            "ASPIS_TRANSCRIPT_NOVELTY",
            "--transcript-discovery-summary",
            str(outdir / "transcript_discovery_summary.tsv"),
            "--outdir",
            str(opt_in_dir),
            "--warnings",
            str(opt_in_dir / "warnings.tsv"),
            "--summary-html",
            str(opt_in_dir / "warnings.html"),
            "--manifest",
            str(opt_in_dir / "warnings_manifest.tsv"),
            "--done",
            str(opt_in_dir / "warnings.done"),
            "--warn-high-true-novel-transcript-fraction",
            "--max-true-novel-transcript-fraction",
            "0.2",
        ]
    )
    opt_in_warnings = read_tsv(opt_in_dir / "warnings.tsv", {"category", "item"})
    if not any(row["category"] == "transcript_discovery" and row["item"] == "true_novel_fraction" for row in opt_in_warnings):
        raise ValueError(f"opt-in true-novel fraction warning was not emitted: {opt_in_warnings}")


def main() -> int:
    counts, metadata = exercise_transcript_matrix()
    exercise_grouped_plots()
    exercise_transcript_novelty_report_summary()
    exercise_gene_biotype_plots()
    exercise_discovery_reports(counts, metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
