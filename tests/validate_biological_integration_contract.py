#!/usr/bin/env python3
"""Exercise biological integration contracts that do not require real data."""

from __future__ import annotations

import csv
import re
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
        "targets_multimir_cache": INPUT / "targets_multimir_cache.tsv",
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
    write_tsv(
        paths["targets_multimir_cache"],
        ["mature_mirna_id", "target_symbol", "target_entrez", "database", "type", "support_type", "pubmed_id"],
        [
            {
                "mature_mirna_id": "hsa-miR-2-3p",
                "target_symbol": "GENE3",
                "target_entrez": "333",
                "database": "miRTarBase",
                "type": "validated",
                "support_type": "functional MTI",
                "pubmed_id": "123456",
            }
        ],
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
            "--target-cache",
            str(paths["targets_multimir_cache"]),
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
    if "miRTarBase" not in {row["target_source"] for row in source_summary}:
        raise ValueError(f"cached multiMiR-style source was not propagated: {source_summary}")
    if "validated" not in {row["target_evidence_type"] for row in source_summary}:
        raise ValueError(f"cached multiMiR-style evidence was not propagated: {source_summary}")

    match_table = BASE / "mirna_mrna_match_table.tsv"
    write_tsv(
        match_table,
        ["pair_id", "smallrna_library_id", "rnaseq_library_id"],
        [
            {"pair_id": "pair_1", "smallrna_library_id": "sm1", "rnaseq_library_id": "rna1"},
            {"pair_id": "pair_2", "smallrna_library_id": "sm2", "rnaseq_library_id": "rna2"},
        ],
    )
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
            "--match-table",
            str(match_table),
            "--min-pairs",
            "2",
        ]
    )
    integration_rows = read_tsv(integration_manifest, {"status", "sample_pairing", "n_sample_pairs", "n_pairs", "n_inverse_pairs", "n_anticorrelated_pairs"})
    row = integration_rows[0]
    if row["status"] != "ok" or int(row["n_sample_pairs"]) != 2 or int(row["n_pairs"]) < 2 or int(row["n_inverse_pairs"]) < 2:
        raise ValueError(f"unexpected miRNA-mRNA integration row: {row}")
    pairing_rows = read_tsv(Path(row["sample_pairing"]), {"pair_id", "smallrna_library_id", "rnaseq_library_id", "match_source"})
    if {item["pair_id"] for item in pairing_rows} != {"pair_1", "pair_2"}:
        raise ValueError(f"explicit pairing provenance was not preserved: {pairing_rows}")


def exercise_biotype_and_dtu(paths: dict[str, Path]) -> None:
    gtf = INPUT / "annotation.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\ttoy\tgene\t1\t100\t.\t+\t.\tgene_id "GENE1"; gene_biotype "protein_coding";',
                'chr1\ttoy\ttranscript\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; transcript_biotype "protein_coding";',
                'chr1\ttoy\ttranscript\t20\t120\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1B"; transcript_biotype "protein_coding";',
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
    write_tsv(gene_counts, ["Geneid", "rna1", "rna2"], [{"Geneid": "GENE1", "rna1": "20", "rna2": "35"}, {"Geneid": "GENE2", "rna1": "3", "rna2": "0"}])
    write_tsv(gene_metadata, ["Geneid", "feature_type"], [{"Geneid": "GENE1", "feature_type": "protein_coding"}, {"Geneid": "GENE2", "feature_type": "lncRNA"}])
    write_tsv(transcript_counts, ["transcript_id", "rna1", "rna2"], [{"transcript_id": "TX1", "rna1": "12", "rna2": "30"}, {"transcript_id": "TX1B", "rna1": "8", "rna2": "5"}, {"transcript_id": "TX2", "rna1": "3", "rna2": "0"}])
    write_tsv(transcript_metadata, ["transcript_id", "gene_id"], [{"transcript_id": "TX1", "gene_id": "GENE1"}, {"transcript_id": "TX1B", "gene_id": "GENE1"}, {"transcript_id": "TX2", "gene_id": "GENE2"}])
    write_tsv(aligned_samples, ["library_id", "bam"], [{"library_id": "rna1", "bam": "rna1.bam"}, {"library_id": "rna2", "bam": "rna2.bam"}])
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
            "--method",
            "DRIMSeq",
            "--contrast-by",
            "time_h",
            "--min-replicates",
            "1",
        ]
    )
    dtu_rows = read_tsv(dtu_dir / "dtu_plan.tsv", {"status", "candidate_methods", "contrast_id", "method", "n_multi_isoform_genes"})
    if dtu_rows[0]["status"] != "ready" or dtu_rows[0]["method"] != "DRIMSeq":
        raise ValueError(f"unexpected DTU plan row: {dtu_rows[0]}")
    if dtu_rows[0]["contrast_id"] != "treated_vs_control__time_h_24" or dtu_rows[0]["n_multi_isoform_genes"] != "1":
        raise ValueError(f"unexpected DTU contrast planning: {dtu_rows[0]}")

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
            str(dtu_dir / "suppa2_plan.tsv"),
            "--done",
            str(dtu_dir / "suppa2_plan.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "SUPPA2",
            "--contrast-by",
            "time_h",
            "--min-replicates",
            "1",
        ]
    )
    suppa2_plan_rows = read_tsv(dtu_dir / "suppa2_plan.tsv", {"status", "method"})
    if suppa2_plan_rows[0]["status"] != "ready" or suppa2_plan_rows[0]["method"] != "SUPPA2":
        raise ValueError(f"SUPPA2 should be planned as a native ready method: {suppa2_plan_rows}")

    fake_drimseq = INPUT / "fake_drimseq_runner.py"
    fake_drimseq.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "def value(flag):",
                "    return Path(args[args.index(flag) + 1])",
                "gene = value('--gene-results')",
                "tx = value('--transcript-results')",
                "summary = value('--summary')",
                "gene.parent.mkdir(parents=True, exist_ok=True)",
                "gene.write_text('gene_id\\tpvalue\\tpadj\\tstatus\\nGENE1\\t0.002\\t0.01\\tok\\n', encoding='utf-8')",
                "tx.write_text('gene_id\\tfeature_id\\tmean_usage_control\\tmean_usage_test\\tdelta_usage\\tstatus\\nGENE1\\tTX1\\t0.2\\t0.8\\t0.6\\tok\\n', encoding='utf-8')",
                "summary.write_text('status\\treason\\tn_tested_genes\\tn_usage_transcripts\\nok\\t\\t1\\t1\\n', encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
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
            "--method",
            "DRIMSeq",
            "--contrast-id",
            "treated_vs_control__time_h_24",
            "--rscript",
            sys.executable,
            "--drimseq-script",
            str(fake_drimseq),
            "--dtu-min-samples",
            "1",
        ]
    )
    drimseq_rows = read_tsv(dtu_dir / "dtu_method_manifest.tsv", {"method", "contrast_id", "status", "standardized_results", "standardized_result_count", "standardized_status"})
    if drimseq_rows[0]["status"] != "completed" or drimseq_rows[0]["standardized_status"] != "ok":
        raise ValueError(f"native DRIMSeq DTU was not standardized: {drimseq_rows}")
    drimseq_standardized = read_tsv(Path(drimseq_rows[0]["standardized_results"]), {"method", "contrast_id", "feature_id", "gene_id", "event_type", "pvalue", "padj"})
    if drimseq_standardized[0]["method"] != "DRIMSeq" or drimseq_standardized[0]["event_type"] != "transcript_usage":
        raise ValueError(f"standardized DRIMSeq row lost method/event type: {drimseq_standardized}")

    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_dtu_plots.py",
            "--method-manifest",
            str(dtu_dir / "dtu_method_manifest.tsv"),
            "--outdir",
            str(dtu_dir / "plots"),
            "--manifest",
            str(dtu_dir / "plots" / "dtu_plot_manifest.tsv"),
            "--done",
            str(dtu_dir / "plots" / "dtu_plots.done"),
            "--top-n",
            "5",
            "--max-points",
            "100",
        ]
    )
    dtu_plot_rows = read_tsv(
        dtu_dir / "plots" / "dtu_plot_manifest.tsv",
        {"method", "status", "reason", "transcript_metadata", "annotation_gtf", "overview_plot", "usage_plot", "feature_plot", "plot_qa_status", "plot_file_count"},
    )
    if dtu_plot_rows[0]["status"] != "ok":
        raise ValueError(f"DTU plot rendering was not ok: {dtu_plot_rows}")
    if dtu_plot_rows[0]["plot_qa_status"] != "ok" or int(dtu_plot_rows[0]["plot_file_count"]) < 2:
        raise ValueError(f"DTU plot QA did not confirm rendered SVGs: {dtu_plot_rows}")
    if dtu_plot_rows[0]["transcript_metadata"] != drimseq_rows[0]["transcript_metadata"]:
        raise ValueError(f"DTU plot manifest lost transcript metadata provenance: {dtu_plot_rows}")
    if dtu_plot_rows[0]["annotation_gtf"] != drimseq_rows[0]["annotation_gtf"]:
        raise ValueError(f"DTU plot manifest lost annotation GTF provenance: {dtu_plot_rows}")
    if (
        not dtu_plot_rows[0]["overview_plot"]
        or not Path(dtu_plot_rows[0]["overview_plot"]).exists()
        or not dtu_plot_rows[0]["usage_plot"]
        or not Path(dtu_plot_rows[0]["usage_plot"]).exists()
    ):
        raise ValueError(f"DTU plot files were not written: {dtu_plot_rows}")
    usage_svg = Path(dtu_plot_rows[0]["usage_plot"]).read_text(encoding="utf-8")
    if "Top transcript-usage genes" not in usage_svg:
        raise ValueError(f"DTU transcript-usage detail plot was not written: {dtu_plot_rows}")
    if "DRIMSeq reports gene-level significance" not in dtu_plot_rows[0].get("reason", ""):
        raise ValueError(f"DRIMSeq missing ranked-candidate reason was not recorded: {dtu_plot_rows}")
    if dtu_plot_rows[0]["feature_plot"]:
        feature_svg = Path(dtu_plot_rows[0]["feature_plot"]).read_text(encoding="utf-8")
        if "Ranked transcript-usage candidates" not in feature_svg:
            raise ValueError(f"DTU transcript-usage candidate plot was not written: {dtu_plot_rows}")

    run_command(
        [
            sys.executable,
            "workflow/scripts/merge_rnaseq_dtu_consensus.py",
            "--method-manifest",
            str(dtu_dir / "dtu_method_manifest.tsv"),
            "--gene-summary",
            str(dtu_dir / "consensus" / "dtu_consensus_gene_summary.tsv"),
            "--method-detail",
            str(dtu_dir / "consensus" / "dtu_consensus_method_detail.tsv"),
            "--done",
            str(dtu_dir / "consensus" / "dtu_consensus.done"),
            "--padj",
            "0.05",
        ]
    )
    consensus_rows = read_tsv(
        dtu_dir / "consensus" / "dtu_consensus_gene_summary.tsv",
        {"gene_id", "methods_detected", "methods_significant", "support_class"},
    )
    if consensus_rows[0]["gene_id"] != "GENE1" or consensus_rows[0]["methods_significant"] != "DRIMSeq":
        raise ValueError(f"DTU consensus gene summary was not populated from DRIMSeq: {consensus_rows}")

    fake_dexseq = INPUT / "fake_dexseq_runner.py"
    fake_dexseq.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "def value(flag):",
                "    return Path(args[args.index(flag) + 1])",
                "gene = value('--gene-results')",
                "tx = value('--feature-results')",
                "summary = value('--summary')",
                "gene.parent.mkdir(parents=True, exist_ok=True)",
                "gene.write_text('gene_id\\tn_features\\tmin_pvalue\\tmin_padj\\tstatus\\nGENE1\\t2\\t0.004\\t0.02\\tok\\n', encoding='utf-8')",
                "tx.write_text('gene_id\\tfeature_id\\tstatistic\\tlog2_fold_change\\tpvalue\\tpadj\\tevent_type\\tmean_usage_control\\tmean_usage_test\\tdelta_usage\\tstatus\\nGENE1\\tTX1\\t5.2\\t1.1\\t0.004\\t0.02\\ttranscript_feature_usage\\t0.2\\t0.8\\t0.6\\tok\\n', encoding='utf-8')",
                "summary.write_text('status\\treason\\tn_tested_genes\\tn_usage_transcripts\\nok\\ttranscript-feature fake DEXSeq\\t1\\t1\\n', encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    dexseq_dir = BASE / "dtu_dexseq"
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
            str(dexseq_dir / "methods"),
            "--manifest",
            str(dexseq_dir / "dtu_method_manifest.tsv"),
            "--done",
            str(dexseq_dir / "dtu_methods.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "DEXSeq",
            "--contrast-id",
            "treated_vs_control__time_h_24",
            "--rscript",
            sys.executable,
            "--dexseq-script",
            str(fake_dexseq),
            "--dtu-min-samples",
            "1",
        ]
    )
    dexseq_rows = read_tsv(dexseq_dir / "dtu_method_manifest.tsv", {"method", "status", "standardized_results", "standardized_status"})
    if dexseq_rows[0]["method"] != "DEXSeq" or dexseq_rows[0]["status"] != "completed" or dexseq_rows[0]["standardized_status"] != "ok":
        raise ValueError(f"native DEXSeq DTU was not standardized: {dexseq_rows}")
    dexseq_standardized = read_tsv(Path(dexseq_rows[0]["standardized_results"]), {"method", "event_type", "feature_id", "gene_id", "pvalue", "padj"})
    if dexseq_standardized[0]["method"] != "DEXSeq" or dexseq_standardized[0]["event_type"] != "transcript_feature_usage":
        raise ValueError(f"standardized DEXSeq row lost method/event type: {dexseq_standardized}")

    fake_suppa = INPUT / "fake_suppa.py"
    fake_suppa.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "cmd = args[0]",
                "def value(flag):",
                "    return Path(args[args.index(flag) + 1])",
                "if cmd == 'generateEvents':",
                "    out = value('-o')",
                "    out.parent.mkdir(parents=True, exist_ok=True)",
                "    Path(str(out) + '.ioi').write_text('gene_id\\tevent_id\\nGENE1\\tGENE1;TX1\\n', encoding='utf-8')",
                "elif cmd == 'psiPerIsoform':",
                "    out = value('-o')",
                "    out.parent.mkdir(parents=True, exist_ok=True)",
                "    Path(str(out) + '.psi').write_text('event\\ts1\\nGENE1;TX1\\t0.5\\n', encoding='utf-8')",
                "elif cmd == 'diffSplice':",
                "    out = value('-o')",
                "    out.parent.mkdir(parents=True, exist_ok=True)",
                "    Path(str(out) + '.dpsi.temp.0').write_text('control_treated_dPSI\\tcontrol_treated_p-val\\nGENE1;TX1\\t0.45\\t0\\nGENE1;TX1\\t0.45\\t0\\nGENE1;TX2\\t0.35\\t0\\nGENE1;TX3\\t0.25\\t0\\n', encoding='utf-8')",
                "    Path(str(out) + '.psivec').write_text('s1\\nGENE1;TX1\\t0.5\\n', encoding='utf-8')",
                "else:",
                "    raise SystemExit(2)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    suppa2_dir = BASE / "dtu_suppa2"
    python_executable = Path(sys.executable).as_posix()
    fake_suppa_path = fake_suppa.as_posix()
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
            str(suppa2_dir / "methods"),
            "--manifest",
            str(suppa2_dir / "dtu_method_manifest.tsv"),
            "--done",
            str(suppa2_dir / "dtu_methods.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "SUPPA2",
            "--contrast-id",
            "treated_vs_control__time_h_24",
            "--suppa2-executable",
            f"{python_executable} {fake_suppa_path}",
        ]
    )
    suppa2_rows = read_tsv(suppa2_dir / "dtu_method_manifest.tsv", {"method", "status", "standardized_results", "standardized_status", "transcript_results", "summary"})
    if suppa2_rows[0]["method"] != "SUPPA2" or suppa2_rows[0]["status"] != "completed" or suppa2_rows[0]["standardized_status"] != "ok":
        raise ValueError(f"native SUPPA2 DTU was not standardized: {suppa2_rows}")
    suppa2_standardized = read_tsv(Path(suppa2_rows[0]["standardized_results"]), {"method", "event_type", "feature_id", "gene_id", "delta_psi", "pvalue", "padj"})
    if len(suppa2_standardized) != 3:
        raise ValueError(f"SUPPA2 duplicate events were not deduplicated: {suppa2_standardized}")
    if suppa2_standardized[0]["method"] != "SUPPA2" or suppa2_standardized[0]["event_type"] != "transcript_event":
        raise ValueError(f"standardized SUPPA2 row lost method/event type: {suppa2_standardized}")
    if suppa2_standardized[0]["gene_id"] != "GENE1" or suppa2_standardized[0]["delta_psi"] != "0.45":
        raise ValueError(f"standardized SUPPA2 row lost identifiers/statistics: {suppa2_standardized}")
    suppa2_events = read_tsv(Path(suppa2_rows[0]["transcript_results"]), {"event_id", "event_type", "delta_psi", "pvalue"})
    if len(suppa2_events) != 3 or suppa2_events[0]["event_id"] != "GENE1;TX1" or suppa2_events[0]["pvalue"] != "0":
        raise ValueError(f"SUPPA2 event result table was not written: {suppa2_events}")
    suppa2_summary = read_tsv(Path(suppa2_rows[0]["summary"]), {"n_tested_genes", "n_usage_transcripts", "n_events"})
    if suppa2_summary[0]["n_tested_genes"] != "1" or suppa2_summary[0]["n_usage_transcripts"] != "3" or suppa2_summary[0]["n_events"] != "3":
        raise ValueError(f"SUPPA2 summary did not expose event counts: {suppa2_summary}")
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_dtu_plots.py",
            "--method-manifest",
            str(suppa2_dir / "dtu_method_manifest.tsv"),
            "--outdir",
            str(suppa2_dir / "plots"),
            "--manifest",
            str(suppa2_dir / "plots" / "dtu_plot_manifest.tsv"),
            "--done",
            str(suppa2_dir / "plots" / "dtu_plots.done"),
        ]
    )
    suppa2_plot_rows = read_tsv(
        suppa2_dir / "plots" / "dtu_plot_manifest.tsv",
        {"method", "status", "usage_plot", "plot_qa_status", "plot_file_count"},
    )
    if suppa2_plot_rows[0]["plot_qa_status"] != "ok" or int(suppa2_plot_rows[0]["plot_file_count"]) < 2:
        raise ValueError(f"SUPPA2 plot QA did not confirm rendered SVGs: {suppa2_plot_rows}")
    suppa2_usage_plot = Path(suppa2_plot_rows[0]["usage_plot"])
    suppa2_overview_svg = Path(suppa2_plot_rows[0]["overview_plot"]).read_text(encoding="utf-8")
    if "Exact zero SUPPA2 p-values are displayed at a finite floor" not in suppa2_overview_svg:
        raise ValueError("SUPPA2 exact-zero overview did not use its method-specific display floor")
    if ">300</text>" in suppa2_overview_svg:
        raise ValueError("SUPPA2 exact-zero overview retained the shared 1e-300 display scale")
    significant_x = re.findall(r'<circle class="sig" cx="([^"]+)"', suppa2_overview_svg)
    if len(significant_x) != 3 or len(set(significant_x)) != 3:
        raise ValueError(f"SUPPA2 exact-zero points were not visibly separated: {significant_x}")
    suppa2_usage_svg = suppa2_usage_plot.read_text(encoding="utf-8")
    if "Top SUPPA2 genes: event detail" not in suppa2_usage_svg or "delta PSI" not in suppa2_usage_svg:
        raise ValueError(f"SUPPA2 delta-PSI usage plot was not rendered correctly: {suppa2_usage_plot}")

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
    exercise_rnaseq_report_dtu_assets(
        dtu_dir / "dtu_plan.tsv",
        dtu_dir / "dtu_method_manifest.tsv",
        dtu_dir / "consensus" / "dtu_consensus_gene_summary.tsv",
        dtu_dir / "consensus" / "dtu_consensus_method_detail.tsv",
        dtu_dir / "consensus" / "dtu_consensus.done",
        dtu_dir / "plots" / "dtu_plot_manifest.tsv",
    )

def exercise_rnaseq_report_dtu_assets(
    dtu_plan: Path,
    dtu_manifest: Path,
    dtu_consensus_gene_summary: Path,
    dtu_consensus_method_detail: Path,
    dtu_consensus_done: Path,
    dtu_plot_manifest: Path,
) -> None:
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
    isoform_dtu_evidence = report_dir / "isoform_dtu_evidence.tsv"
    isoform_dtu_summary = report_dir / "isoform_dtu_evidence_summary.tsv"
    isoform_interpretation = report_dir / "isoform_interpretation_consensus.tsv"
    isoform_interpretation_summary = report_dir / "isoform_interpretation_consensus_summary.tsv"
    write_tsv(
        isoform_dtu_evidence,
        [
            "event_id",
            "contrast_id",
            "gene_id",
            "isoform_id",
            "dtu_evidence_status",
            "dtu_methods_detected",
            "dtu_methods_significant",
        ],
        [
            {
                "event_id": "switch1",
                "contrast_id": contrast,
                "gene_id": "GENE1",
                "isoform_id": "TX1",
                "dtu_evidence_status": "supported_significant",
                "dtu_methods_detected": "DRIMSeq,SUPPA2",
                "dtu_methods_significant": "DRIMSeq,SUPPA2",
            }
        ],
    )
    write_tsv(
        isoform_dtu_summary,
        [
            "status",
            "isoform_candidates",
            "switch_events",
            "genes_with_dtu_support",
            "genes_with_significant_dtu_support",
            "candidate_rows_with_dtu_support",
            "candidate_rows_with_significant_dtu_support",
            "dtu_methods_seen",
            "dtu_methods_significant",
            "reason",
        ],
        [
            {
                "status": "ok",
                "isoform_candidates": "1",
                "switch_events": "1",
                "genes_with_dtu_support": "1",
                "genes_with_significant_dtu_support": "1",
                "candidate_rows_with_dtu_support": "1",
                "candidate_rows_with_significant_dtu_support": "1",
                "dtu_methods_seen": "DRIMSeq,SUPPA2",
                "dtu_methods_significant": "DRIMSeq,SUPPA2",
                "reason": "contract fixture",
            }
        ],
    )
    write_tsv(
        isoform_interpretation,
        [
            "event_id",
            "contrast_id",
            "gene_id",
            "gene_name",
            "isoform_id",
            "switch_role",
            "interpretation_status",
            "interpretation_priority",
            "interpretation_label",
            "review_reason",
            "dtu_support_class",
        ],
        [
            {
                "event_id": "switch1",
                "contrast_id": contrast,
                "gene_id": "GENE1",
                "gene_name": "Gene One",
                "isoform_id": "TX1",
                "switch_role": "switch_in",
                "interpretation_status": "supported",
                "interpretation_priority": "high",
                "interpretation_label": "isoform switch with multi-method DTU/splicing support",
                "review_reason": "contract fixture",
                "dtu_support_class": "multi_method_significant",
            }
        ],
    )
    write_tsv(
        isoform_interpretation_summary,
        [
            "status",
            "interpretation_rows",
            "high_priority_rows",
            "medium_priority_rows",
            "low_priority_rows",
            "multi_method_supported_rows",
            "single_method_supported_rows",
            "no_dtu_support_rows",
            "reason",
        ],
        [
            {
                "status": "ok",
                "interpretation_rows": "1",
                "high_priority_rows": "1",
                "medium_priority_rows": "0",
                "low_priority_rows": "0",
                "multi_method_supported_rows": "1",
                "single_method_supported_rows": "0",
                "no_dtu_support_rows": "0",
                "reason": "contract fixture",
            }
        ],
    )
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
            "--dtu-consensus-gene-summary",
            str(dtu_consensus_gene_summary),
            "--dtu-consensus-method-detail",
            str(dtu_consensus_method_detail),
            "--dtu-consensus-done",
            str(dtu_consensus_done),
            "--dtu-plot-manifest",
            str(dtu_plot_manifest),
            "--isoform-dtu-evidence",
            str(isoform_dtu_evidence),
            "--isoform-dtu-evidence-summary",
            str(isoform_dtu_summary),
            "--isoform-interpretation-consensus",
            str(isoform_interpretation),
            "--isoform-interpretation-consensus-summary",
            str(isoform_interpretation_summary),
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
    if not any(row["asset_label"] == "DRIMSeq_standardized_results" for row in dtu_assets):
        raise ValueError(f"standardized DRIMSeq results were not exposed as report assets: {dtu_assets}")
    if not any(row["asset_label"] == "dtu_plot_manifest" for row in dtu_assets):
        raise ValueError(f"DTU plot manifest was not exposed as a report asset: {dtu_assets}")
    if not any(row["asset_label"] == "dtu_consensus_gene_summary" for row in dtu_assets):
        raise ValueError(f"DTU consensus gene summary was not exposed as a report asset: {dtu_assets}")
    if not any(row["asset_label"] == "dtu_consensus_method_detail" for row in dtu_assets):
        raise ValueError(f"DTU consensus method detail was not exposed as a report asset: {dtu_assets}")
    if not any(row["asset_label"].endswith("_overview_plot") for row in dtu_assets):
        raise ValueError(f"DTU overview plots were not exposed as report assets: {dtu_assets}")
    isoform_assets = [row for row in assets if row["asset_group"] == "isoform_switch"]
    if not any(row["asset_label"] == "isoform_dtu_evidence" for row in isoform_assets):
        raise ValueError(f"isoform/DTU evidence was not exposed as a report asset: {isoform_assets}")
    if not any(row["asset_label"] == "isoform_interpretation_consensus" for row in isoform_assets):
        raise ValueError(f"isoform interpretation consensus was not exposed as a report asset: {isoform_assets}")
    html_text = (report_dir / "index.html").read_text(encoding="utf-8")
    if "DTU / splicing methods" not in html_text or "standardized rows: 1" not in html_text:
        raise ValueError("DTU summary was not rendered in the RNA-seq report index")
    if "padj&lt;0.05 rows: 1" not in html_text or "usage table" not in html_text or "overview plot" not in html_text:
        raise ValueError("DTU contrast table was not rendered in the RNA-seq report index")
    if "DRIMSeq gene-level differential transcript usage" not in html_text or "link-grid" not in html_text:
        raise ValueError("DTU method table did not render compact links and fallback method reasons")
    if "Cross-method DTU consensus" not in html_text or "single_method_significant" not in html_text:
        raise ValueError("DTU consensus summary was not rendered in the RNA-seq report index")
    if "Isoform-switch / DTU interpretation" not in html_text or "candidate rows with significant DTU support: 1" not in html_text:
        raise ValueError("isoform/DTU evidence summary was not rendered in the RNA-seq report index")
    if "high priority: 1" not in html_text or "interpretation consensus" not in html_text:
        raise ValueError("isoform interpretation consensus summary was not rendered in the RNA-seq report index")


def main() -> int:
    paths = setup_common_inputs()
    exercise_target_and_integration(paths)
    exercise_biotype_and_dtu(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
